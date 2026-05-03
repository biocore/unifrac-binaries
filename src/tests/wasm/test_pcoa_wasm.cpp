/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2026, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

/*
 * WASM correctness test for PCoA via the existing pcoa() C wrapper
 * (api.hpp:1074, calls su::pcoa -> skbb_pcoa_fsvd_fp64).
 *
 * Tolerances (mirrors skbb's WASM PCoA test):
 *   - eigenvalues:          1e-6 absolute
 *   - proportion_explained: 1e-6 absolute
 *   - samples (coords):     1e-3 absolute, with sign-per-axis flip to
 *                           handle the eigenvector-sign ambiguity.
 *
 * Layout: skbb writes samples row-major n_samples × n_eighs (each row
 * is one sample, each column an axis). See skbb
 * principal_coordinate_analysis.cpp:574-578.
 */

#include "tests/wasm/check_macros.hpp"
#include "tests/wasm/fixtures.hpp"
#include "tests/wasm/expected/pcoa_expected.h"
#include "api.hpp"

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
    pcoa(FIXTURE_UNWEIGHTED_DIST, N, K, &eigenvalues, &samples, &prop_exp);

    CHECK(eigenvalues != nullptr);
    CHECK(samples != nullptr);
    CHECK(prop_exp != nullptr);

    for (unsigned int k = 0; k < K; k++) {
        if (!almost_equal<double>(eigenvalues[k], PCOA_EXPECTED_EIGENVALUES[k], 1e-6)) {
            std::fprintf(stderr,
                "FAIL eigenvalues[%u] = %.12g; expected %.12g\n",
                k, eigenvalues[k], PCOA_EXPECTED_EIGENVALUES[k]);
            std::exit(1);
        }
    }

    // For non-Euclidean distances like UniFrac, skbb normalizes
    // proportion_explained by the sum of positive eigenvalues only,
    // so the top-K sum is not bounded by 1.0. Bit-equality against
    // the native-generated expected is the authoritative assertion.
    for (unsigned int k = 0; k < K; k++) {
        if (!almost_equal<double>(prop_exp[k], PCOA_EXPECTED_PROPORTION_EXPLAINED[k], 1e-6)) {
            std::fprintf(stderr,
                "FAIL proportion_explained[%u] = %.12g; expected %.12g\n",
                k, prop_exp[k], PCOA_EXPECTED_PROPORTION_EXPLAINED[k]);
            std::exit(1);
        }
    }

    if (!samples_match_with_sign_flip(samples, PCOA_EXPECTED_SAMPLES, K, N, 1e-3)) {
        std::exit(1);
    }

    std::free(eigenvalues);
    std::free(samples);
    std::free(prop_exp);

    std::fprintf(stdout, "OK test_pcoa_wasm\n");
    return 0;
}
