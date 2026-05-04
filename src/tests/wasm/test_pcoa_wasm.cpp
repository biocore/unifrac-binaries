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
 * (api.hpp, calls su::pcoa -> skbb_pcoa_fsvd_fp64).
 *
 * Asserts structural invariants of the decomposition rather than pinning
 * to a native oracle:
 *   - eigenvalues are returned in non-increasing order;
 *   - the leading eigenvalue is positive;
 *   - prop_exp values are non-negative and share a single denominator with
 *     the eigenvalues (so prop_exp[i]/prop_exp[j] == eigval[i]/eigval[j]);
 *   - samples are deterministic across re-seeded re-runs (bit equal).
 *
 * Layout: skbb writes samples row-major n_samples × n_eighs.
 */

#include "tests/wasm/check_macros.hpp"
#include "tests/wasm/fixtures.hpp"
#include "api.hpp"

int main(void) {
    const unsigned int N = 6;
    const unsigned int K = 3;
    const unsigned int seed = 42;

    ssu_set_random_seed(seed);

    double *eigenvalues = nullptr;
    double *samples     = nullptr;
    double *prop_exp    = nullptr;
    pcoa(FIXTURE_UNWEIGHTED_DIST, N, K, &eigenvalues, &samples, &prop_exp);

    CHECK(eigenvalues != nullptr);
    CHECK(samples != nullptr);
    CHECK(prop_exp != nullptr);

    CHECK(eigenvalues[0] > 0.0);
    for (unsigned int k = 1; k < K; k++) {
        CHECK(eigenvalues[k] <= eigenvalues[k - 1]);
    }

    // skbb computes prop_exp[k] = eigval[k] / trace(centered_gram), which
    // shares a single positive scalar denominator across all axes. Verify
    // that ratio invariant without needing access to the trace.
    for (unsigned int k = 0; k < K; k++) {
        CHECK(prop_exp[k] >= 0.0);
    }
    CHECK(prop_exp[0] > 0.0);
    for (unsigned int k = 1; k < K; k++) {
        // eigval[k] / eigval[0] should equal prop_exp[k] / prop_exp[0].
        double lhs = eigenvalues[k] / eigenvalues[0];
        double rhs = prop_exp[k]    / prop_exp[0];
        if (!almost_equal<double>(lhs, rhs, 1e-9)) {
            std::fprintf(stderr,
                "FAIL prop_exp ratio mismatch at k=%u: eig %.12g vs prop %.12g\n",
                k, lhs, rhs);
            std::exit(1);
        }
    }

    // Determinism: same seed → identical eigenvalues and samples.
    ssu_set_random_seed(seed);
    double *eig2 = nullptr, *samp2 = nullptr, *pe2 = nullptr;
    pcoa(FIXTURE_UNWEIGHTED_DIST, N, K, &eig2, &samp2, &pe2);
    for (unsigned int k = 0; k < K; k++) {
        CHECK(eig2[k] == eigenvalues[k]);
    }
    for (unsigned int i = 0; i < N * K; i++) {
        CHECK(samp2[i] == samples[i]);
    }

    std::free(eigenvalues);
    std::free(samples);
    std::free(prop_exp);
    std::free(eig2);
    std::free(samp2);
    std::free(pe2);

    std::fprintf(stdout, "OK test_pcoa_wasm\n");
    return 0;
}
