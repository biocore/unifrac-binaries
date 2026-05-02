/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2026, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

/*
 * Stage 6 WASM correctness test for compute_permanova_inmem_fp64.
 *
 * Tolerance:
 *   - fstat:  1e-6 abs. The unpermuted F is computed from
 *     sum_upper_square reductions over the n×n distance matrix; native
 *     LAPACK and WASM Eigen reductions can diverge by ~8 ULPs due to
 *     associativity differences even at OMP_NUM_THREADS=1, which is
 *     numerically irrelevant for PERMANOVA.
 *   - pvalue: 1e-6 abs. skbb's portable Fisher-Yates (commit 14b5306)
 *     keeps the permutation sequence identical, so cross-toolchain
 *     drift in pvalue is bounded by the rounding in s_W per permutation
 *     — well below any meaningful pvalue resolution at n_perm=999.
 *   - determinism (re-seed + re-run): exact-equality (==) within a
 *     single binary, since identical RNG state -> identical FP ops.
 */

#include "tests/wasm/check_macros.hpp"
#include "tests/wasm/expected/permanova_expected.h"
#include "api.hpp"

#include <cstdint>
#include <cstdlib>

static const double DIST[36] = {
    0.0,         0.2,         0.57142857,  0.6,         0.5,         0.2,
    0.2,         0.0,         0.42857143,  0.66666667,  0.6,         0.33333333,
    0.57142857,  0.42857143,  0.0,         0.71428571,  0.85714286,  0.42857143,
    0.6,         0.66666667,  0.71428571,  0.0,         0.33333333,  0.4,
    0.5,         0.6,         0.85714286,  0.33333333,  0.0,         0.6,
    0.2,         0.33333333,  0.42857143,  0.4,         0.6,         0.0
};
static const uint32_t GROUPING[6] = {0, 0, 1, 1, 1, 0};

int main(void) {
    CHECK_EQ(PERMANOVA_EXPECTED_N, 6u);

    ssu_set_random_seed(PERMANOVA_EXPECTED_SEED);

    double fstat  = 0.0;
    double pvalue = 0.0;
    ComputeStatus status = compute_permanova_inmem_fp64(
        DIST, PERMANOVA_EXPECTED_N, GROUPING,
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

    // Determinism check: re-seed and re-run; same fstat / pvalue.
    ssu_set_random_seed(PERMANOVA_EXPECTED_SEED);
    double fstat2 = 0.0, pvalue2 = 0.0;
    compute_permanova_inmem_fp64(DIST, PERMANOVA_EXPECTED_N, GROUPING,
                                 PERMANOVA_EXPECTED_N_PERM, &fstat2, &pvalue2);
    CHECK(fstat2 == fstat);
    CHECK(pvalue2 == pvalue);

    // NULL-input rejection.
    {
        double f = 0.0, p = 0.0;
        CHECK(compute_permanova_inmem_fp64(NULL, 6, GROUPING, 999, &f, &p) != okay);
        CHECK(compute_permanova_inmem_fp64(DIST, 6, NULL, 999, &f, &p) != okay);
        CHECK(compute_permanova_inmem_fp64(DIST, 0, GROUPING, 999, &f, &p) != okay);
        CHECK(compute_permanova_inmem_fp64(DIST, 6, GROUPING, 0, &f, &p) != okay);
    }

    std::fprintf(stdout, "OK test_permanova_wasm fstat=%.6f pvalue=%.6f\n", fstat, pvalue);
    return 0;
}
