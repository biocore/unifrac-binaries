/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2026, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

/*
 * WASM correctness test for faith_pd_inmem.
 *
 * Hand-derived from the multifurcating tree in fixtures.hpp (root has three
 * children — GG_OTU_1, an inner clade {GG_OTU_2, GG_OTU_3}, and an inner
 * clade {GG_OTU_5, GG_OTU_4}; all branches = 1.0) and the per-sample OTU
 * sets implied by FIXTURE_INDICES / FIXTURE_INDPTR.
 *
 * Tolerance 1e-6 (matches src/test_su.cpp:1641).
 */

#include "tests/wasm/check_macros.hpp"
#include "tests/wasm/fixtures.hpp"
#include "api.hpp"

int main(void) {
    const support_biom_t   table = {
        (char**) FIXTURE_OBS_IDS, (char**) FIXTURE_SAMP_IDS,
        (uint32_t*) FIXTURE_INDICES, (uint32_t*) FIXTURE_INDPTR,
        (double*) FIXTURE_DATA,
        FIXTURE_N_OBS, FIXTURE_N_SAMP, 0
    };
    const support_bptree_t tree  = {
        (bool*) FIXTURE_STRUCTURE, (double*) FIXTURE_LENGTHS,
        (char**) FIXTURE_NAMES, (int) FIXTURE_NPARENS
    };

    const double expected[6] = {4., 5., 6., 3., 2., 5.};

    r_vec* result = nullptr;
    ComputeStatus status = faith_pd_inmem(&table, &tree, &result);

    CHECK(status == okay);
    CHECK(result != nullptr);
    CHECK_EQ(result->n_samples, (unsigned int) FIXTURE_N_SAMP);

    for (unsigned int i = 0; i < (unsigned int) FIXTURE_N_SAMP; i++) {
        if (!almost_equal<double>(result->values[i], expected[i], 1e-6)) {
            std::fprintf(stderr,
                "FAIL faith_pd[%u] = %.9f; expected %.9f\n",
                i, result->values[i], expected[i]);
            std::exit(1);
        }
    }

    destroy_results_vec(&result);

    std::fprintf(stdout, "OK test_faith_pd_wasm\n");
    return 0;
}
