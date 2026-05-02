/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2026, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

/*
 * Stage 5 WASM correctness test for PCoA via the existing pcoa() C
 * wrapper in api.hpp (which calls su::pcoa -> skbb_pcoa_fsvd_fp64).
 *
 * Input: 6x6 unweighted-UniFrac distance matrix from
 *   test/capi_inmem_test.c:99-104 (native-CI-validated).
 * Seed: 42, propagated through ssu_set_random_seed -> su::set_random_seed
 *   -> skbb_set_random_seed.
 *
 * Ground truth produced natively by
 *   src/tests/wasm/generate_pcoa_expected.cpp linked against skbb-build
 *   conda env's libskbb.so. Stored as hex-encoded doubles in
 *   src/tests/wasm/expected/pcoa_expected.h.
 *
 * Tolerances (mirrors skbb's WASM PCoA test):
 *   - eigenvalues:          1e-6 absolute
 *   - proportion_explained: 1e-6 absolute
 *   - samples (coords):     1e-3 absolute, with sign-per-axis flip to
 *                           handle the eigenvector-sign ambiguity.
 */

#include "tests/wasm/check_macros.hpp"
#include "tests/wasm/expected/pcoa_expected.h"
#include "api.hpp"

#include <cstdint>
#include <cstdlib>

// Same 6x6 unweighted UniFrac distance matrix the generator hardcodes.
static const double DIST[36] = {
    0.0,         0.2,         0.57142857,  0.6,         0.5,         0.2,
    0.2,         0.0,         0.42857143,  0.66666667,  0.6,         0.33333333,
    0.57142857,  0.42857143,  0.0,         0.71428571,  0.85714286,  0.42857143,
    0.6,         0.66666667,  0.71428571,  0.0,         0.33333333,  0.4,
    0.5,         0.6,         0.85714286,  0.33333333,  0.0,         0.6,
    0.2,         0.33333333,  0.42857143,  0.4,         0.6,         0.0
};

// Sample-coord layout in skbb's output is row-major n_samples × n_eighs
// (each row = one sample's K coordinates; each column = one eigenvector
// axis). See skbb principal_coordinate_analysis.cpp:574-578. Compare
// each axis (column) under both sign conventions.
static bool samples_match_with_sign_flip(const double* obs, const double* exp,
                                          unsigned int K, unsigned int N,
                                          double tol) {
    for (unsigned int col = 0; col < K; col++) {
        double max_abs_diff_pos = 0.0;
        double max_abs_diff_neg = 0.0;
        for (unsigned int row = 0; row < N; row++) {
            double g = obs[row * K + col];
            double w = exp[row * K + col];
            double dpos = g - w;  if (dpos < 0) dpos = -dpos;
            double dneg = g + w;  if (dneg < 0) dneg = -dneg;
            if (dpos > max_abs_diff_pos) max_abs_diff_pos = dpos;
            if (dneg > max_abs_diff_neg) max_abs_diff_neg = dneg;
        }
        // Accept the sign that minimizes per-axis L_inf error.
        double ax_err = (max_abs_diff_pos < max_abs_diff_neg) ? max_abs_diff_pos
                                                              : max_abs_diff_neg;
        if (ax_err > tol) {
            std::fprintf(stderr,
                "FAIL pcoa samples axis %u: max diff (pos sign) = %.3e, "
                "(neg sign) = %.3e; tol = %.3e\n",
                col, max_abs_diff_pos, max_abs_diff_neg, tol);
            return false;
        }
    }
    return true;
}

int main(void) {
    const unsigned int N = PCOA_EXPECTED_N;       // 6
    const unsigned int K = PCOA_EXPECTED_K;       // 3
    CHECK_EQ(N, 6u);
    CHECK_EQ(K, 3u);

    ssu_set_random_seed(PCOA_EXPECTED_SEED);

    double *eigenvalues = nullptr;
    double *samples     = nullptr;
    double *prop_exp    = nullptr;

    // pcoa() in api.hpp:1074 wraps su::pcoa, which seeds skbb (-1 -> use
    // skbb's mt19937 which we just set) and mallocs all three output
    // buffers.
    pcoa(DIST, N, K, &eigenvalues, &samples, &prop_exp);

    CHECK(eigenvalues != nullptr);
    CHECK(samples != nullptr);
    CHECK(prop_exp != nullptr);

    // (1) Eigenvalues to 1e-6.
    for (unsigned int k = 0; k < K; k++) {
        if (!almost_equal<double>(eigenvalues[k], PCOA_EXPECTED_EIGENVALUES[k], 1e-6)) {
            std::fprintf(stderr,
                "FAIL eigenvalues[%u] = %.12g; expected %.12g\n",
                k, eigenvalues[k], PCOA_EXPECTED_EIGENVALUES[k]);
            std::exit(1);
        }
    }

    // (2) Proportion explained to 1e-6.
    // Note: for non-Euclidean distances like UniFrac, the centered matrix
    // can have negative eigenvalues, and skbb normalizes by the sum of
    // positive eigenvalues only, so the top-K sum is not bounded by 1.0.
    // The bit-equality check against the native-generated expected is
    // the authoritative assertion.
    for (unsigned int k = 0; k < K; k++) {
        if (!almost_equal<double>(prop_exp[k], PCOA_EXPECTED_PROPORTION_EXPLAINED[k], 1e-6)) {
            std::fprintf(stderr,
                "FAIL proportion_explained[%u] = %.12g; expected %.12g\n",
                k, prop_exp[k], PCOA_EXPECTED_PROPORTION_EXPLAINED[k]);
            std::exit(1);
        }
    }

    // (3) Samples (coords) with axis-wise sign flip, tol 1e-3.
    if (!samples_match_with_sign_flip(samples, PCOA_EXPECTED_SAMPLES, K, N, 1e-3)) {
        std::exit(1);
    }

    // (4) Cleanup. su::pcoa allocates with malloc; free() is correct.
    std::free(eigenvalues);
    std::free(samples);
    std::free(prop_exp);

    std::fprintf(stdout, "OK test_pcoa_wasm\n");
    return 0;
}
