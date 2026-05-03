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
 * Tolerance 1e-6 (matches src/test_su.cpp:1641).
 */

#include "tests/wasm/check_macros.hpp"
#include "tests/wasm/fixtures.hpp"
#include "tests/wasm/expected/faith_pd_expected.h"
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

    r_vec* result = nullptr;
    ComputeStatus status = faith_pd_inmem(&table, &tree, &result);

    CHECK(status == okay);
    CHECK(result != nullptr);
    CHECK_EQ(result->n_samples, FAITH_PD_EXPECTED_N);

    for (unsigned int i = 0; i < FAITH_PD_EXPECTED_N; i++) {
        if (!almost_equal<double>(result->values[i], FAITH_PD_EXPECTED[i], 1e-6)) {
            std::fprintf(stderr,
                "FAIL faith_pd[%u] = %.9f; expected %.9f\n",
                i, result->values[i], FAITH_PD_EXPECTED[i]);
            std::exit(1);
        }
    }

    destroy_results_vec(&result);

    std::fprintf(stdout, "OK test_faith_pd_wasm\n");
    return 0;
}
