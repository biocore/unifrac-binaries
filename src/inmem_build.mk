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
# Objects and archives
# --------------------------------------------------------------------------
# Same translation units as libssu_wasm.a -- see wasm/emscripten_build.mk for
# the rationale on which TUs are excluded. The .inmem.o suffix keeps these
# objects distinct from the regular .o (full libssu.so) and .wasm.o objects in
# the same directory; .inmem-asan.o does the same for the instrumented build, so
# the two do not invalidate each other and an instrumented archive can never be
# mistaken for a shippable one.
#
# One call emits both rules for a TU: $(1) object stem, $(2) source,
# $(3) extra flags, $(4) extra prereqs. Same shape as skbb's skbb_cpu_tu.
define inmem_tu
$(1).inmem.o: $(2) $(4)
	$$(INMEM_CXX) $$(INMEM_CXXFLAGS) $(3) -c $$< -o $$@
$(1).inmem-asan.o: $(2) $(4)
	$$(INMEM_CXX) $$(INMEM_CXXFLAGS) $$(INMEM_ASAN_FLAGS) $(3) -c $$< -o $$@
endef

# unifrac_accapi_cpu.cpp and unifrac_task_noclass_cpu.cpp are generated; the
# native Makefile's rules for them fire first because they are listed as
# prerequisites here. skbio_alt needs the staged skbb headers so its
# <scikit-bio-binaries/...> includes resolve.
$(eval $(call inmem_tu,tree,tree.cpp,,tree.hpp))
$(eval $(call inmem_tu,biom_inmem,biom_inmem.cpp,,biom_inmem.hpp biom_interface.hpp))
$(eval $(call inmem_tu,biom_subsampled,biom_subsampled.cpp,,biom_subsampled.hpp biom_inmem.hpp omp_stub.h))
$(eval $(call inmem_tu,unifrac,unifrac.cpp,,unifrac.hpp unifrac_internal.hpp unifrac_task.hpp tree.hpp))
$(eval $(call inmem_tu,unifrac_internal,unifrac_internal.cpp,,unifrac_internal.hpp tree.hpp biom_interface.hpp))
$(eval $(call inmem_tu,skbio_alt,skbio_alt.cpp,,skbio_alt.hpp $(INMEM_SKBB_STAGED_HS)))
$(eval $(call inmem_tu,api,api.cpp,,api.hpp api_compat.hpp unifrac.hpp skbio_alt.hpp biom_inmem.hpp biom_subsampled.hpp tree.hpp))
$(eval $(call inmem_tu,unifrac_accapi_cpu,unifrac_accapi_cpu.cpp,-DSUCMP_NM=su_cpu,unifrac_accapi.hpp unifrac_accapi_impl.hpp))
$(eval $(call inmem_tu,unifrac_task_cpu,unifrac_task_noclass_cpu.cpp,-DSUCMP_NM=su_cpu,unifrac_task_noclass.hpp unifrac_task_impl.hpp))
$(eval $(call inmem_tu,unifrac_cmp_cpu,unifrac_cmp.cpp,-DSUCMP_NM=su_cpu,unifrac_cmp.hpp unifrac_internal.hpp unifrac.hpp unifrac_task.hpp unifrac_task_noclass.hpp biom_interface.hpp tree.hpp))

INMEM_STEMS := tree biom_inmem biom_subsampled unifrac unifrac_internal \
               unifrac_accapi_cpu unifrac_task_cpu unifrac_cmp_cpu skbio_alt api
INMEM_OBJS      := $(INMEM_STEMS:=.inmem.o)
INMEM_ASAN_OBJS := $(INMEM_STEMS:=.inmem-asan.o)

libssu_inmem.a: $(INMEM_OBJS)
	rm -f $@
	$(INMEM_AR) rcs $@ $(INMEM_OBJS)

libssu_inmem_asan.a: $(INMEM_ASAN_OBJS)
	rm -f $@
	$(INMEM_AR) rcs $@ $(INMEM_ASAN_OBJS)

inmem_static: libssu_inmem.a

# --------------------------------------------------------------------------
# Test
# --------------------------------------------------------------------------
# Why this suite exists at all, and what it can and cannot cover, is documented
# in the header of tests/inmem/test_concurrency_inmem.cpp.
#
# Linking is the embedder's problem in general -- the archive carries no skbb --
# but the test has to resolve those symbols somehow. It links whatever skbb is
# installed under PREFIX, which is also where INMEM_SKBB_EXTERN should point so
# the headers match the library. The rpath is what lets the test run straight
# out of the build directory: conda does not put its lib dir on
# LD_LIBRARY_PATH, so without it the binary links fine and then fails to start.
INMEM_SKBB_LIB     ?= -lskbb
INMEM_TEST_LDFLAGS ?= -L$(PREFIX)/lib -Wl,-rpath,$(PREFIX)/lib
INMEM_TEST_DEPS    := tests/inmem/test_concurrency_inmem.cpp \
                      tests/wasm/fixtures.hpp tests/wasm/check_macros.hpp \
                      api.hpp $(INMEM_SKBB_STAGED_HS)

test_concurrency_inmem: $(INMEM_TEST_DEPS) libssu_inmem.a
	$(INMEM_CXX) $(INMEM_CXXFLAGS) $< -o $@ libssu_inmem.a \
	    $(INMEM_TEST_LDFLAGS) $(INMEM_SKBB_LIB) -lpthread

inmem_test: test_concurrency_inmem
	./test_concurrency_inmem

# ASan variant. A plain run only catches a fault that happens to land, so heap
# corruption under concurrency can pass silently without this. Three notes:
#   - The runtime has to be preloaded even though the binary links it:
#     libskbb.so gets initialized ahead of it and ASan then refuses to start.
#     -static-libasan would sidestep the preload, but conda-forge's
#     libsanitizer package ships no libasan.a, so it will not link there.
#   - Leak detection is off. The target is memory safety in unifrac's own code
#     under concurrency, not allocation hygiene in whatever skbb is installed.
#   - gcc and clang only. Both -fsanitize=address and the -print-file-name
#     lookup below are gcc/clang spellings; the NVIDIA HPC SDK and the AMD
#     offload compilers do not accept them. That costs nothing in practice --
#     this archive is CPU-only (-DUNIFRAC_WASM=1, no accelerator TUs) and
#     INMEM_CXX defaults to $(CXX) -- but it must fail clearly rather than emit
#     a binary that was never instrumented, so the target checks first.
#     Override INMEM_ASAN_FLAGS if a compiler spells it differently.
# The archive keeps its shipped -O3, so what is instrumented is what ships.
INMEM_ASAN_FLAGS ?= -fsanitize=address -g

# gcc reports "gcc"/"g++", clang and its derivatives report "clang" in --version
INMEM_ASAN_CXX_OK := $(shell $(INMEM_CXX) --version 2>/dev/null | head -1 | grep -ciE 'gcc|g\+\+|clang')

inmem_asan_supported:
	@test "$(INMEM_ASAN_CXX_OK)" != "0" || { \
	    echo "ERROR: the ASan targets need gcc or clang; INMEM_CXX='$(INMEM_CXX)' reports:"; \
	    $(INMEM_CXX) --version 2>&1 | head -1; \
	    echo "       Build the plain 'inmem_test' target, or set INMEM_CXX to gcc/clang."; \
	    exit 1; }

test_concurrency_inmem_asan: $(INMEM_TEST_DEPS) libssu_inmem_asan.a | inmem_asan_supported
	$(INMEM_CXX) $(INMEM_CXXFLAGS) $(INMEM_ASAN_FLAGS) $< -o $@ libssu_inmem_asan.a \
	    $(INMEM_TEST_LDFLAGS) $(INMEM_SKBB_LIB) -lpthread

inmem_test_asan: test_concurrency_inmem_asan
	@asan_rt=`$(INMEM_CXX) -print-file-name=libasan.so`; \
	    test -f "$$asan_rt" || { echo "ERROR: no ASan runtime from '$(INMEM_CXX) -print-file-name=libasan.so' (got '$$asan_rt')"; exit 1; }; \
	    echo "LD_PRELOAD=$$asan_rt ASAN_OPTIONS=detect_leaks=0 ./test_concurrency_inmem_asan"; \
	    LD_PRELOAD=$$asan_rt ASAN_OPTIONS=detect_leaks=0 ./test_concurrency_inmem_asan

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
	rm -f libssu_inmem.a libssu_inmem_asan.a *.inmem.o *.inmem-asan.o \
	      test_concurrency_inmem test_concurrency_inmem_asan
	rm -rf $(INMEM_SKBB_INC_STAGE)

.PHONY: inmem_static inmem_test inmem_test_asan inmem_asan_supported install_inmem inmem_clean
