/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2026, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

/*
 * WASM correctness test for compute_permanova_inmem_fp64.
 *
 * Asserts structural invariants of the result rather than pinning to a
 * native oracle: with the FIXTURE_UNWEIGHTED_DIST + FIXTURE_GROUPING two-
 * group split and 999 permutations under a fixed seed, the computed
 * F-statistic must be strictly positive and the p-value must lie in
 * (0, 1]. Within-binary determinism is asserted exact: re-seeding and
 * re-running reproduces fstat and pvalue bitwise.
 *
 * The native-side capi_inmem_test.c covers the same call with the same
 * invariants and also exercises the dispatcher routing path.
 */

#include "tests/wasm/check_macros.hpp"
#include "tests/wasm/fixtures.hpp"
#include "api.hpp"

int main(void) {
    const unsigned int n_perm = 999;
    const unsigned int seed   = 42;

    ssu_set_random_seed(seed);

    double fstat  = 0.0;
    double pvalue = 0.0;
    ComputeStatus status = compute_permanova_inmem_fp64(
        FIXTURE_UNWEIGHTED_DIST, 6, FIXTURE_GROUPING, n_perm, &fstat, &pvalue);

    CHECK(status == okay);
    CHECK(fstat > 0.0);
    CHECK(pvalue > 0.0);
    CHECK(pvalue <= 1.0);

    // Determinism: re-seed and re-run; same fstat / pvalue.
    ssu_set_random_seed(seed);
    double fstat2 = 0.0, pvalue2 = 0.0;
    compute_permanova_inmem_fp64(FIXTURE_UNWEIGHTED_DIST, 6, FIXTURE_GROUPING,
                                 n_perm, &fstat2, &pvalue2);
    CHECK(fstat2 == fstat);
    CHECK(pvalue2 == pvalue);

    // NULL-input rejection.
    {
        double f = 0.0, p = 0.0;
        CHECK(compute_permanova_inmem_fp64(NULL, 6, FIXTURE_GROUPING, 999, &f, &p) != okay);
        CHECK(compute_permanova_inmem_fp64(FIXTURE_UNWEIGHTED_DIST, 6, NULL, 999, &f, &p) != okay);
    }

    std::fprintf(stdout, "OK test_permanova_wasm fstat=%.6f pvalue=%.6f\n", fstat, pvalue);
    return 0;
}
