/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2026, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

/*
 * Shared test fixtures + tiny helpers for the WASM test suite and its
 * native ground-truth generators. Single source of truth; consumed by
 * test_*_wasm.cpp and generate_*_expected.cpp.
 *
 * Synthetic 5-OTU / 6-sample CSR table identical to the data in
 * test/capi_inmem_test.c:21-27, whose unweighted-UniFrac matrix is
 * already validated by the existing native CI pass.
 *
 * Per-sample observed OTU sets (decoded from FIXTURE_INDICES /
 * FIXTURE_INDPTR):
 *   Sample 0: {GG_OTU_2, GG_OTU_4}
 *   Sample 1: {GG_OTU_2, GG_OTU_4, GG_OTU_5}
 *   Sample 2: {GG_OTU_1, GG_OTU_3, GG_OTU_4, GG_OTU_5}
 *   Sample 3: {GG_OTU_2, GG_OTU_3}
 *   Sample 4: {GG_OTU_2}
 *   Sample 5: {GG_OTU_2, GG_OTU_3, GG_OTU_4}
 *
 * Tree topology: multifurcating root with three children — GG_OTU_1,
 * an inner clade {GG_OTU_2, GG_OTU_3}, and an inner clade
 * {GG_OTU_5, GG_OTU_4}. Branch lengths set to 1.0 so Faith's PD and
 * the UniFrac variants are non-trivial.
 *
 * The 6x6 distance matrix FIXTURE_UNWEIGHTED_DIST is the unweighted
 * UniFrac result on this fixture, hardcoded at
 * test/capi_inmem_test.c:99-104. Used as input to the PCoA and
 * PERMANOVA tests.
 */

#ifndef UNIFRAC_WASM_FIXTURES_HPP
#define UNIFRAC_WASM_FIXTURES_HPP

#include <cstdint>
#include <cstdio>

// ---- table (CSR) --------------------------------------------------------
static const int        FIXTURE_N_OBS  = 5;
static const int        FIXTURE_N_SAMP = 6;
static const char* const FIXTURE_OBS_IDS[]  = {"GG_OTU_1", "GG_OTU_2", "GG_OTU_3", "GG_OTU_4", "GG_OTU_5"};
static const char* const FIXTURE_SAMP_IDS[] = {"Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6"};
static const uint32_t    FIXTURE_INDICES[]  = {2, 0, 1, 3, 4, 5, 2, 3, 5, 0, 1, 2, 5, 1, 2};
static const uint32_t    FIXTURE_INDPTR[]   = {0, 1, 6, 9, 13, 15};
static const double      FIXTURE_DATA[]     = {1., 5., 1., 2., 3., 1., 1., 4., 2., 2., 1., 1., 1., 1., 1.};

// ---- tree (balanced parens, all branch lengths 1.0) ---------------------
static const unsigned int FIXTURE_NPARENS = 16;
static const bool   FIXTURE_STRUCTURE[] = { true, true, false, true,
                                            true, false, true, false,
                                            false, true, true, false,
                                            true, false, false, false };
static const double FIXTURE_LENGTHS[]   = { 0.0, 1.0, 0.0, 1.0,
                                            1.0, 0.0, 1.0, 0.0,
                                            0.0, 1.0, 1.0, 0.0,
                                            1.0, 0.0, 0.0, 0.0 };
static const char* const FIXTURE_NAMES[] = {"", "GG_OTU_1", "", "",
                                            "GG_OTU_2", "", "GG_OTU_3", "",
                                            "", "", "GG_OTU_5", "",
                                            "GG_OTU_4", "", "", ""};

// ---- pre-computed unweighted-UniFrac distance matrix --------------------
// 6x6 row-major. From test/capi_inmem_test.c:99-104.
static const double FIXTURE_UNWEIGHTED_DIST[36] = {
    0.0,         0.2,         0.57142857,  0.6,         0.5,         0.2,
    0.2,         0.0,         0.42857143,  0.66666667,  0.6,         0.33333333,
    0.57142857,  0.42857143,  0.0,         0.71428571,  0.85714286,  0.42857143,
    0.6,         0.66666667,  0.71428571,  0.0,         0.33333333,  0.4,
    0.5,         0.6,         0.85714286,  0.33333333,  0.0,         0.6,
    0.2,         0.33333333,  0.42857143,  0.4,         0.6,         0.0
};

// 2-group split for PERMANOVA: arbitrary but fixed.
static const uint32_t FIXTURE_GROUPING[6] = {0, 0, 1, 1, 1, 0};

// ---- generator-only helpers ---------------------------------------------
// Both helpers are header-only and `static inline`; only the generator
// binaries reference them. Tests that don't need these are unaffected.

#ifdef UNIFRAC_WASM_FIXTURES_USE_SKBB_SEED
#include <random>
#include <scikit-bio-binaries/util.h>

// Replicates su::set_random_seed: seed mt19937 with `s`, draw one
// uint32, pass that to skbb_set_random_seed. Generators must mirror
// this so the seed seen by skbb matches what ssu_set_random_seed
// produces under WASM.
static inline void su_compatible_seed(uint32_t s) {
    std::mt19937 g(s);
    uint32_t skbb_seed = g();
    skbb_set_random_seed(skbb_seed);
}
#endif

// Emit a C-array initializer of hex-encoded doubles.
static inline void emit_array_double(const char* name, const double* arr, unsigned int n) {
    std::printf("static const double %s[%u] = {\n", name, n);
    for (unsigned int i = 0; i < n; i++) {
        std::printf("    %a%s\n", arr[i], (i + 1 == n) ? "" : ",");
    }
    std::printf("};\n\n");
}

#endif /* UNIFRAC_WASM_FIXTURES_HPP */
