# WebAssembly build rules for unifrac-binaries.
#
# Activated by the top-level `wasm` target. Produces libssu_wasm.a, a
# single-threaded static archive intended to be linked into downstream
# emscripten projects (e.g. duckdb-miint).
#
# Backend: scikit-bio-binaries WASM build (libskbb_wasm.a from skbb's
# `make wasm` target). No HDF5, no lz4, no GPU, no OpenMP, no pthread,
# no CPU-arch dispatch.
#
# Expected toolchain (activated emsdk on PATH):
#   emcc, em++, emar
#
# This Makefile fragment is included from src/Makefile and must be run
# with src/ as the working directory.

# --------------------------------------------------------------------------
# skbb dependency
# --------------------------------------------------------------------------
# SKBB_DIR is the path to a scikit-bio-binaries checkout that has been
# built with `make wasm` (i.e. libskbb_wasm.a exists at its repo root).
# Override on the make command line for non-sibling layouts:
#   make wasm SKBB_DIR=/some/other/path
SKBB_DIR ?= $(abspath $(CURDIR)/../../scikit-bio-binaries)
WASM_SKBB_LIB := $(SKBB_DIR)/src/libskbb_wasm.a
WASM_SKBB_EXTERN := $(SKBB_DIR)/src/extern

# skbio_alt.cpp #include's <scikit-bio-binaries/util.h>. Stage skbb's
# public headers under .wasm-skbb-include/scikit-bio-binaries/ so that
# prefix resolves identically to a `make install` layout.
WASM_SKBB_INC_STAGE := .wasm-skbb-include
WASM_SKBB_STAGED_HS := $(WASM_SKBB_INC_STAGE)/scikit-bio-binaries/util.h \
                      $(WASM_SKBB_INC_STAGE)/scikit-bio-binaries/distance.h \
                      $(WASM_SKBB_INC_STAGE)/scikit-bio-binaries/ordination.h

$(WASM_SKBB_INC_STAGE)/scikit-bio-binaries/%.h: $(WASM_SKBB_EXTERN)/%.h
	@mkdir -p $(WASM_SKBB_INC_STAGE)/scikit-bio-binaries
	cp $< $@

# --------------------------------------------------------------------------
# Compiler vars
# --------------------------------------------------------------------------
WASM_CXX      := em++
WASM_AR       := emar
WASM_CXXFLAGS := -std=c++17 -O3 -Wall -I. -I$(WASM_SKBB_INC_STAGE) \
                 -DUNIFRAC_WASM=1 \
                 -DSKIP_MMAP=1 \
                 -DNOGPU=1 \
                 -fno-exceptions \
                 -Wno-unknown-pragmas

# --------------------------------------------------------------------------
# WASM object list
# --------------------------------------------------------------------------
# Mirrors the in-memory subset of the native build: same sources, just
# compiled with em++ to .wasm.o suffix so the suffixes don't clash with
# the native object files in the same directory.
#
# Excluded vs. native libssu.so:
#   biom.o      — HDF5-only BIOM v2 reader. WASM uses biom_inmem only.
#   tsv.o       — TSV grouping-file parser used only by the file-based
#                 PERMANOVA path (gated under UNIFRAC_WASM in api.cpp).
#   cmd.o       — driver for the ssu/faithpd CLI binaries. Library only.
WASM_OBJS := \
    tree.wasm.o \
    biom_inmem.wasm.o \
    biom_subsampled.wasm.o \
    unifrac.wasm.o \
    unifrac_internal.wasm.o \
    unifrac_accapi_cpu.wasm.o \
    unifrac_task_cpu.wasm.o \
    unifrac_cmp_cpu.wasm.o \
    skbio_alt.wasm.o \
    api.wasm.o

# Generated cpp sources need to exist before their .wasm.o rules fire.
# Reuse the native Makefile's generators by listing them as prerequisites.
unifrac_accapi_cpu.wasm.o: unifrac_accapi_cpu.cpp unifrac_accapi.hpp unifrac_accapi_impl.hpp
	$(WASM_CXX) $(WASM_CXXFLAGS) -DSUCMP_NM=su_cpu -c $< -o $@

unifrac_task_cpu.wasm.o: unifrac_task_noclass_cpu.cpp unifrac_task_noclass.hpp unifrac_task_impl.hpp
	$(WASM_CXX) $(WASM_CXXFLAGS) -DSUCMP_NM=su_cpu -c $< -o $@

unifrac_cmp_cpu.wasm.o: unifrac_cmp.cpp unifrac_cmp.hpp unifrac_internal.hpp unifrac.hpp unifrac_task.hpp unifrac_task_noclass.hpp biom_interface.hpp tree.hpp
	$(WASM_CXX) $(WASM_CXXFLAGS) -DSUCMP_NM=su_cpu -c $< -o $@

# Plain-cpp rules. Each TU gets a tracked-prereq list so that header
# changes trigger rebuilds. The skbio_alt.wasm.o rule depends on the
# staged skbb headers so the `<scikit-bio-binaries/...>` include resolves.
tree.wasm.o: tree.cpp tree.hpp
	$(WASM_CXX) $(WASM_CXXFLAGS) -c $< -o $@

biom_inmem.wasm.o: biom_inmem.cpp biom_inmem.hpp biom_interface.hpp
	$(WASM_CXX) $(WASM_CXXFLAGS) -c $< -o $@

biom_subsampled.wasm.o: biom_subsampled.cpp biom_subsampled.hpp biom_inmem.hpp omp_stub.h
	$(WASM_CXX) $(WASM_CXXFLAGS) -c $< -o $@

unifrac.wasm.o: unifrac.cpp unifrac.hpp unifrac_internal.hpp unifrac_task.hpp tree.hpp
	$(WASM_CXX) $(WASM_CXXFLAGS) -c $< -o $@

unifrac_internal.wasm.o: unifrac_internal.cpp unifrac_internal.hpp tree.hpp biom_interface.hpp
	$(WASM_CXX) $(WASM_CXXFLAGS) -c $< -o $@

skbio_alt.wasm.o: skbio_alt.cpp skbio_alt.hpp $(WASM_SKBB_STAGED_HS)
	$(WASM_CXX) $(WASM_CXXFLAGS) -c $< -o $@

api.wasm.o: api.cpp api.hpp api_compat.hpp unifrac.hpp skbio_alt.hpp biom_inmem.hpp biom_subsampled.hpp tree.hpp
	$(WASM_CXX) $(WASM_CXXFLAGS) -c $< -o $@

# --------------------------------------------------------------------------
# Archive
# --------------------------------------------------------------------------
libssu_wasm.a: $(WASM_OBJS)
	rm -f $@
	$(WASM_AR) rcs $@ $(WASM_OBJS)

wasm: libssu_wasm.a

# --------------------------------------------------------------------------
# WASM test infrastructure
# --------------------------------------------------------------------------
# Emscripten emits a sibling .js + .wasm pair; node runs the .js, which
# auto-loads the .wasm. NODERAWFS=0 is the default — WASM tests do not
# touch the filesystem.
WASM_TEST_LDFLAGS := -sEXIT_RUNTIME=1 \
                     -sALLOW_MEMORY_GROWTH=1 \
                     -sENVIRONMENT=node \
                     -sNODERAWFS=0

test_smoke_wasm.js: libssu_wasm.a $(WASM_SKBB_LIB) tests/wasm/test_smoke_wasm.cpp \
                   tests/wasm/check_macros.hpp $(WASM_SKBB_STAGED_HS)
	$(WASM_CXX) $(WASM_CXXFLAGS) tests/wasm/test_smoke_wasm.cpp \
	    libssu_wasm.a $(WASM_SKBB_LIB) \
	    $(WASM_TEST_LDFLAGS) -o $@

wasm_test: test_smoke_wasm.js
	@echo "--- smoke ---"
	node test_smoke_wasm.js

# --------------------------------------------------------------------------
# Cleanup
# --------------------------------------------------------------------------
wasm_clean:
	rm -f libssu_wasm.a *.wasm.o *_wasm.js *_wasm.wasm
	rm -rf $(WASM_SKBB_INC_STAGE)

.PHONY: wasm wasm_test wasm_clean
