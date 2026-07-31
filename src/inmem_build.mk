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
# Where the public skbb headers are staged from. Defaults to a source
# checkout; override to consume an installed skbb instead, e.g. a conda
# package, so headers and library come from the same version:
#   make inmem_static INMEM_SKBB_EXTERN=$CONDA_PREFIX/include/scikit-bio-binaries
INMEM_SKBB_EXTERN ?= $(SKBB_DIR)/src/extern

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
# Test
# --------------------------------------------------------------------------
# Concurrency coverage for the archive as embedders link it. The suite in
# test_su.cpp exercises the same entry points, but only as built for
# libssu.so; this build defines UNIFRAC_WASM (no signal handler, CPU_SETSIZE
# fallback) while still being multi-threaded, so it is a distinct
# configuration.
#
# Linking is the embedder's problem in general -- the archive carries no skbb
# -- but the test has to resolve those symbols somehow. It links whatever skbb
# is installed under PREFIX, which is also where INMEM_SKBB_EXTERN should point
# so the headers match the library.
#
# The rpath is what lets the test run straight out of the build directory:
# conda does not put its lib dir on LD_LIBRARY_PATH, so without it the binary
# links fine and then fails to start.
INMEM_SKBB_LIB     ?= -lskbb
INMEM_TEST_LDFLAGS ?= -L$(PREFIX)/lib -Wl,-rpath,$(PREFIX)/lib

test_concurrency_inmem: tests/inmem/test_concurrency_inmem.cpp libssu_inmem.a \
                        tests/wasm/fixtures.hpp tests/wasm/check_macros.hpp \
                        api.hpp $(INMEM_SKBB_STAGED_HS)
	$(INMEM_CXX) $(INMEM_CXXFLAGS) $< -o $@ libssu_inmem.a \
	    $(INMEM_TEST_LDFLAGS) $(INMEM_SKBB_LIB) -lpthread

inmem_test: test_concurrency_inmem
	./test_concurrency_inmem

# ASan variant. A plain run only catches a fault that happens to land; the
# n_substeps cases in particular corrupted the heap silently before they were
# clamped, and only ASan called it. Note that the report_status use-after-free
# that motivated the concurrency work cannot fire in *this* configuration --
# UNIFRAC_WASM means no SIGUSR1 handler is installed, so no flag is ever set --
# so what this gate covers is the memory-safety class generally, not that
# specific bug. Notes:
#   - The runtime has to be preloaded even though the binary links it:
#     libskbb.so gets initialized ahead of it and ASan then refuses to start.
#     -static-libasan would sidestep the preload, but conda-forge's
#     libsanitizer package ships no libasan.a, so it will not link there.
#   - Leak detection is off. The target here is memory safety in unifrac's own
#     code under concurrency, not allocation hygiene in whatever skbb build
#     happens to be installed.
#   - The archive keeps its shipped -O3; the n_substeps overflow this covers
#     was confirmed to report at that level.
#   - Objects are rebuilt instrumented, so this cleans on the way in and back
#     out again; otherwise an instrumented libssu_inmem.a would sit there
#     looking up to date and get shipped.
inmem_test_asan:
	$(MAKE) inmem_clean
	$(MAKE) test_concurrency_inmem INMEM_MPFLAG="-fopenmp -fsanitize=address -g"
	@asan_rt=`$(INMEM_CXX) -print-file-name=libasan.so`; \
	    test -f "$$asan_rt" || { echo "ERROR: no ASan runtime from '$(INMEM_CXX) -print-file-name=libasan.so' (got '$$asan_rt')"; exit 1; }; \
	    echo "LD_PRELOAD=$$asan_rt ASAN_OPTIONS=detect_leaks=0 ./test_concurrency_inmem"; \
	    LD_PRELOAD=$$asan_rt ASAN_OPTIONS=detect_leaks=0 ./test_concurrency_inmem; \
	    rc=$$?; $(MAKE) inmem_clean; exit $$rc

# --------------------------------------------------------------------------
# Install (archive + public headers under a stable prefix layout).
# Embedders pick up libssu_inmem.a + the unifrac/ header tree via
# vcpkg / CMake / pkg-config conventions.
#
# Guards against the empty-PREFIX footgun: the top-level Makefile
# falls back to CONDA_PREFIX, but if both are unset `mkdir -p /lib`
# would silently target the root filesystem.
# --------------------------------------------------------------------------
install_inmem: libssu_inmem.a
	@test -n "$(PREFIX)" || { echo "ERROR: PREFIX is unset (and CONDA_PREFIX is unset). Pass PREFIX=/path or activate a conda env."; exit 1; }
	mkdir -p ${PREFIX}/lib ${PREFIX}/include/unifrac
	rm -f ${PREFIX}/lib/libssu_inmem.a; cp libssu_inmem.a ${PREFIX}/lib/
	rm -f ${PREFIX}/include/unifrac/api.hpp; cp api.hpp ${PREFIX}/include/unifrac/
	rm -f ${PREFIX}/include/unifrac/task_parameters.hpp; cp task_parameters.hpp ${PREFIX}/include/unifrac/
	rm -f ${PREFIX}/include/unifrac/status_enum.hpp; cp status_enum.hpp ${PREFIX}/include/unifrac/

# --------------------------------------------------------------------------
# Cleanup
# --------------------------------------------------------------------------
inmem_clean:
	rm -f libssu_inmem.a *.inmem.o test_concurrency_inmem
	rm -rf $(INMEM_SKBB_INC_STAGE)

.PHONY: inmem_static inmem_test inmem_test_asan install_inmem inmem_clean
