# Native in-memory static-archive build for unifrac-binaries.
#
# Activated by the top-level `inmem_static` target. Produces
# libssu_inmem.a — the same in-memory subset shipped as libssu_wasm.a,
# but compiled with the host toolchain so OpenMP is enabled. Intended
# for downstream native projects that want to embed the in-memory
# UniFrac surface without dragging in HDF5, lz4, or the dlopen-based
# CPU dispatch.
#
# Excluded vs. the standard libssu.so:
#   biom.o, tsv.o, cmd.o — same exclusions as the WASM build (see
#   wasm/emscripten_build.mk for rationale).
#
# Skbb is treated as the embedder's responsibility: this archive
# contains only unifrac translation units, and the staged skbb header
# tree exists purely to satisfy `<scikit-bio-binaries/...>` includes
# at compile time.
#
# This Makefile fragment is included from src/Makefile and must be
# run with src/ as the working directory.

# --------------------------------------------------------------------------
# skbb dependency
# --------------------------------------------------------------------------
# SKBB_DIR overlaps with the WASM build's variable on purpose: a
# single sibling checkout serves both. Override on the command line
# for non-sibling layouts:
#   make inmem_static SKBB_DIR=/some/other/path
SKBB_DIR ?= $(abspath $(CURDIR)/../../scikit-bio-binaries)
INMEM_SKBB_EXTERN := $(SKBB_DIR)/src/extern

INMEM_SKBB_INC_STAGE := .inmem-skbb-include
INMEM_SKBB_STAGED_HS := $(INMEM_SKBB_INC_STAGE)/scikit-bio-binaries/util.h \
                       $(INMEM_SKBB_INC_STAGE)/scikit-bio-binaries/distance.h \
                       $(INMEM_SKBB_INC_STAGE)/scikit-bio-binaries/ordination.h

$(INMEM_SKBB_INC_STAGE)/scikit-bio-binaries/%.h: $(INMEM_SKBB_EXTERN)/%.h
	@mkdir -p $(INMEM_SKBB_INC_STAGE)/scikit-bio-binaries
	cp $< $@

# --------------------------------------------------------------------------
# Compiler vars
# --------------------------------------------------------------------------
INMEM_CXX      ?= $(CXX)
INMEM_AR       ?= ar
INMEM_MPFLAG   ?= -fopenmp
INMEM_CXXFLAGS := -std=c++17 -O3 -Wall -fPIC -I. -I$(INMEM_SKBB_INC_STAGE) \
                  -DUNIFRAC_WASM=1 \
                  -DSKIP_MMAP=1 \
                  -DNOGPU=1 \
                  $(INMEM_MPFLAG) \
                  -Wno-unknown-pragmas

# --------------------------------------------------------------------------
# Object list
# --------------------------------------------------------------------------
# Same translation units as libssu_wasm.a — see wasm/emscripten_build.mk
# for the rationale on which TUs are excluded. The .inmem.o suffix keeps
# these objects distinct from the regular .o (full libssu.so) and .wasm.o
# objects in the same directory.
INMEM_OBJS := \
    tree.inmem.o \
    biom_inmem.inmem.o \
    biom_subsampled.inmem.o \
    unifrac.inmem.o \
    unifrac_internal.inmem.o \
    unifrac_accapi_cpu.inmem.o \
    unifrac_task_cpu.inmem.o \
    unifrac_cmp_cpu.inmem.o \
    skbio_alt.inmem.o \
    api.inmem.o

# Generated cpp sources need to exist before their .inmem.o rules fire.
# Reuse the native Makefile's generators by listing them as prerequisites.
unifrac_accapi_cpu.inmem.o: unifrac_accapi_cpu.cpp unifrac_accapi.hpp unifrac_accapi_impl.hpp
	$(INMEM_CXX) $(INMEM_CXXFLAGS) -DSUCMP_NM=su_cpu -c $< -o $@

unifrac_task_cpu.inmem.o: unifrac_task_noclass_cpu.cpp unifrac_task_noclass.hpp unifrac_task_impl.hpp
	$(INMEM_CXX) $(INMEM_CXXFLAGS) -DSUCMP_NM=su_cpu -c $< -o $@

unifrac_cmp_cpu.inmem.o: unifrac_cmp.cpp unifrac_cmp.hpp unifrac_internal.hpp unifrac.hpp unifrac_task.hpp unifrac_task_noclass.hpp biom_interface.hpp tree.hpp
	$(INMEM_CXX) $(INMEM_CXXFLAGS) -DSUCMP_NM=su_cpu -c $< -o $@

# Plain-cpp rules. Each TU gets a tracked-prereq list so that header
# changes trigger rebuilds. The skbio_alt.inmem.o rule depends on the
# staged skbb headers so the `<scikit-bio-binaries/...>` include resolves.
tree.inmem.o: tree.cpp tree.hpp
	$(INMEM_CXX) $(INMEM_CXXFLAGS) -c $< -o $@

biom_inmem.inmem.o: biom_inmem.cpp biom_inmem.hpp biom_interface.hpp
	$(INMEM_CXX) $(INMEM_CXXFLAGS) -c $< -o $@

biom_subsampled.inmem.o: biom_subsampled.cpp biom_subsampled.hpp biom_inmem.hpp omp_stub.h
	$(INMEM_CXX) $(INMEM_CXXFLAGS) -c $< -o $@

unifrac.inmem.o: unifrac.cpp unifrac.hpp unifrac_internal.hpp unifrac_task.hpp tree.hpp
	$(INMEM_CXX) $(INMEM_CXXFLAGS) -c $< -o $@

unifrac_internal.inmem.o: unifrac_internal.cpp unifrac_internal.hpp tree.hpp biom_interface.hpp
	$(INMEM_CXX) $(INMEM_CXXFLAGS) -c $< -o $@

skbio_alt.inmem.o: skbio_alt.cpp skbio_alt.hpp $(INMEM_SKBB_STAGED_HS)
	$(INMEM_CXX) $(INMEM_CXXFLAGS) -c $< -o $@

api.inmem.o: api.cpp api.hpp api_compat.hpp unifrac.hpp skbio_alt.hpp biom_inmem.hpp biom_subsampled.hpp tree.hpp
	$(INMEM_CXX) $(INMEM_CXXFLAGS) -c $< -o $@

# --------------------------------------------------------------------------
# Archive
# --------------------------------------------------------------------------
libssu_inmem.a: $(INMEM_OBJS)
	rm -f $@
	$(INMEM_AR) rcs $@ $(INMEM_OBJS)

inmem_static: libssu_inmem.a

# --------------------------------------------------------------------------
# Install (archive + public headers under a stable prefix layout).
# Embedders pick up libssu_inmem.a + the unifrac/ header tree via
# vcpkg / CMake / pkg-config conventions.
# --------------------------------------------------------------------------
install_inmem: libssu_inmem.a
	mkdir -p ${PREFIX}/lib ${PREFIX}/include/unifrac
	rm -f ${PREFIX}/lib/libssu_inmem.a; cp libssu_inmem.a ${PREFIX}/lib/
	rm -f ${PREFIX}/include/unifrac/api.hpp; cp api.hpp ${PREFIX}/include/unifrac/
	rm -f ${PREFIX}/include/unifrac/task_parameters.hpp; cp task_parameters.hpp ${PREFIX}/include/unifrac/
	rm -f ${PREFIX}/include/unifrac/status_enum.hpp; cp status_enum.hpp ${PREFIX}/include/unifrac/

# --------------------------------------------------------------------------
# Cleanup
# --------------------------------------------------------------------------
inmem_clean:
	rm -f libssu_inmem.a *.inmem.o
	rm -rf $(INMEM_SKBB_INC_STAGE)

.PHONY: inmem_static install_inmem inmem_clean
