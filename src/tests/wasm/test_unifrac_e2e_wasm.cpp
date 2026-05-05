/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2026, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

/*
 * End-to-end WASM correctness test.
 *
 * Runs one_off_matrix_inmem_v3 under WASM for the unweighted method and
 * verifies the resulting 6x6 distance matrix matches FIXTURE_UNWEIGHTED_DIST,
 * which is the same oracle hardcoded in test/capi_inmem_test.c:99-104 and
 * thus exercised by the existing native CI.
 *
 * For weighted_normalized, weighted_unnormalized, and generalized α=1,
 * structural invariants are checked: zero diagonal, symmetry, finite values.
 * Per-method bit-precise correctness for these variants is covered by the
 * existing native test_su.cpp suite.
 *
 * Tolerance for unweighted matches src/testdata/validation_tests.sh (1e-5).
 */

#include "tests/wasm/check_macros.hpp"
#include "tests/wasm/fixtures.hpp"
#include "api.hpp"

#include <cmath>

static void check_structural(const char* label, const double* m, unsigned int n) {
    for (unsigned int i = 0; i < n; i++) {
        if (m[i * n + i] != 0.0) {
            std::fprintf(stderr, "FAIL %s diag[%u] = %.10g (expected 0)\n",
                label, i, m[i * n + i]);
            std::exit(1);
        }
        for (unsigned int j = i + 1; j < n; j++) {
            double a = m[i * n + j];
            double b = m[j * n + i];
            if (!std::isfinite(a) || !std::isfinite(b)) {
                std::fprintf(stderr, "FAIL %s [%u,%u] non-finite\n", label, i, j);
                std::exit(1);
            }
            if (!almost_equal<double>(a, b, 1e-12)) {
                std::fprintf(stderr,
                    "FAIL %s asymmetric [%u,%u]=%.10g vs [%u,%u]=%.10g\n",
                    label, i, j, a, j, i, b);
                std::exit(1);
            }
            if (a < 0.0) {
                std::fprintf(stderr, "FAIL %s [%u,%u] = %.10g (negative)\n",
                    label, i, j, a);
                std::exit(1);
            }
        }
    }
}

static void run_method(const char* label, const char* method_arg, double alpha,
                       const double* oracle /* may be nullptr */) {
    const support_biom_t   table = {(char**) FIXTURE_OBS_IDS, (char**) FIXTURE_SAMP_IDS,
                                    (uint32_t*) FIXTURE_INDICES, (uint32_t*) FIXTURE_INDPTR,
                                    (double*) FIXTURE_DATA, FIXTURE_N_OBS, FIXTURE_N_SAMP, 0};
    const support_bptree_t tree  = {(bool*) FIXTURE_STRUCTURE, (double*) FIXTURE_LENGTHS,
                                    (char**) FIXTURE_NAMES, (int) FIXTURE_NPARENS};

    mat_full_fp64_t* result = nullptr;
    ComputeStatus status = one_off_matrix_inmem_v3(
        &table, &tree, method_arg,
        /*variance_adjust*/ false, alpha,
        /*bypass_tips*/ false, /*normalize_sample_counts*/ true,
        /*n_substeps*/ 1,
        /*subsample_depth*/ 0, /*subsample_with_replacement*/ false,
        /*mmap_dir*/ NULL,
        &result);

    CHECK(status == okay);
    CHECK(result != nullptr);
    CHECK_EQ(result->n_samples, (unsigned int) FIXTURE_N_SAMP);

    check_structural(label, result->matrix, FIXTURE_N_SAMP);

    if (oracle != nullptr) {
        const unsigned int total = FIXTURE_N_SAMP * FIXTURE_N_SAMP;
        for (unsigned int i = 0; i < total; i++) {
            if (!almost_equal<double>(result->matrix[i], oracle[i], 1e-5)) {
                std::fprintf(stderr,
                    "FAIL %s entry [%u/%u] = %.10g; expected %.10g\n",
                    label, i / FIXTURE_N_SAMP, i % FIXTURE_N_SAMP,
                    result->matrix[i], oracle[i]);
                std::exit(1);
            }
        }
    }

    destroy_mat_full_fp64(&result);
    std::fprintf(stdout, "  %-22s OK\n", label);
}

int main(void) {
    run_method("unweighted",            "unweighted_fp64",            1.0, FIXTURE_UNWEIGHTED_DIST);
    run_method("weighted_normalized",   "weighted_normalized_fp64",   1.0, nullptr);
    run_method("weighted_unnormalized", "weighted_unnormalized_fp64", 1.0, nullptr);
    run_method("generalized",           "generalized_fp64",           1.0, nullptr);

    std::fprintf(stdout, "OK test_unifrac_e2e_wasm\n");
    return 0;
}
