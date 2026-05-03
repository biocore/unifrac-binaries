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
 * Tolerance:
 *   - fstat / pvalue: 1e-6 abs vs native-generated expected. The
 *     unpermuted F is computed via sum-over-distance-matrix
 *     reductions; native LAPACK and WASM Eigen builds can diverge by
 *     ~8 ULPs due to floating-point associativity in parallel
 *     reductions even at OMP_NUM_THREADS=1, well below any meaningful
 *     PERMANOVA resolution.
 *   - within-binary determinism (re-seed + re-run): exact-equality.
 */

#include "tests/wasm/check_macros.hpp"
#include "tests/wasm/fixtures.hpp"
#include "tests/wasm/expected/permanova_expected.h"
#include "api.hpp"

int main(void) {
    CHECK_EQ(PERMANOVA_EXPECTED_N, 6u);

    ssu_set_random_seed(PERMANOVA_EXPECTED_SEED);

    double fstat  = 0.0;
    double pvalue = 0.0;
    ComputeStatus status = compute_permanova_inmem_fp64(
        FIXTURE_UNWEIGHTED_DIST, PERMANOVA_EXPECTED_N, FIXTURE_GROUPING,
        PERMANOVA_EXPECTED_N_PERM, &fstat, &pvalue);

    CHECK(status == okay);

    if (!almost_equal<double>(fstat, PERMANOVA_EXPECTED_FSTAT, 1e-6)) {
        std::fprintf(stderr,
            "FAIL fstat = %a (%.12g); expected %a (%.12g)\n",
            fstat, fstat, PERMANOVA_EXPECTED_FSTAT, PERMANOVA_EXPECTED_FSTAT);
        std::exit(1);
    }
    if (!almost_equal<double>(pvalue, PERMANOVA_EXPECTED_PVALUE, 1e-6)) {
        std::fprintf(stderr,
            "FAIL pvalue = %a (%.12g); expected %a (%.12g)\n",
            pvalue, pvalue, PERMANOVA_EXPECTED_PVALUE, PERMANOVA_EXPECTED_PVALUE);
        std::exit(1);
    }

    // Determinism: re-seed and re-run; same fstat / pvalue.
    ssu_set_random_seed(PERMANOVA_EXPECTED_SEED);
    double fstat2 = 0.0, pvalue2 = 0.0;
    compute_permanova_inmem_fp64(FIXTURE_UNWEIGHTED_DIST, PERMANOVA_EXPECTED_N,
                                 FIXTURE_GROUPING, PERMANOVA_EXPECTED_N_PERM,
                                 &fstat2, &pvalue2);
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
