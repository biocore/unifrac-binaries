/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2026, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

/*
 * Stage 3 WASM correctness test for faith_pd_inmem.
 *
 * Fixture: the 5-OTU / 6-sample synthetic table from
 * test/capi_inmem_test.c (whose unweighted-UniFrac correctness is
 * guaranteed by the existing native CI pass), paired with the same
 * tree topology but with all branch lengths set to 1.0 so Faith's PD
 * is non-trivial.
 *
 * Expected values committed at src/tests/wasm/expected/faith_pd_expected.h
 * with full provenance documentation. Tolerance 1e-6 (matches the native
 * test_su.cpp:1641 tolerance for the same family of computations).
 */

#include "tests/wasm/check_macros.hpp"
#include "tests/wasm/expected/faith_pd_expected.h"
#include "api.hpp"

#include <cstdint>

// ---- biom (CSR format) ---------------------------------------------------
// Identical to test/capi_inmem_test.c:21-27. Decoded per-sample OTU sets:
//   Sample 0: GG_OTU_2, GG_OTU_4
//   Sample 1: GG_OTU_2, GG_OTU_4, GG_OTU_5
//   Sample 2: GG_OTU_1, GG_OTU_3, GG_OTU_4, GG_OTU_5
//   Sample 3: GG_OTU_2, GG_OTU_3
//   Sample 4: GG_OTU_2
//   Sample 5: GG_OTU_2, GG_OTU_3, GG_OTU_4
static const int n_obs  = 5;
static const int n_samp = 6;
static const char* const obs_ids[]   = {"GG_OTU_1", "GG_OTU_2", "GG_OTU_3", "GG_OTU_4", "GG_OTU_5"};
static const char* const samp_ids[]  = {"Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6"};
static const uint32_t    indices[]   = {2, 0, 1, 3, 4, 5, 2, 3, 5, 0, 1, 2, 5, 1, 2};
static const uint32_t    indptr[]    = {0, 1, 6, 9, 13, 15};
static const double      data[]      = {1., 5., 1., 2., 3., 1., 1., 4., 2., 2., 1., 1., 1., 1., 1.};

// ---- tree ---------------------------------------------------------------
// Topology identical to test/capi_inmem_test.c:30-42 (multifurcating root
// with three children: GG_OTU_1, inner clade {GG_OTU_2, GG_OTU_3}, inner
// clade {GG_OTU_5, GG_OTU_4}). Branch lengths overridden to 1.0 so Faith's
// PD is non-trivial.
static const unsigned int nparens = 16;
static const bool structure[] = { true, true, false, true,
                                  true, false, true, false,
                                  false, true, true, false,
                                  true, false, false, false };
// One length per paren slot; only the open-paren slots carry a meaningful
// edge length. The native BPTree constructor uses the value at the
// position of the corresponding open-paren node.
static const double lengths[] = { 0.0, 1.0, 0.0, 1.0,
                                  1.0, 0.0, 1.0, 0.0,
                                  0.0, 1.0, 1.0, 0.0,
                                  1.0, 0.0, 0.0, 0.0 };
static const char* const names[] = {"", "GG_OTU_1", "", "",
                                    "GG_OTU_2", "", "GG_OTU_3", "",
                                    "", "", "GG_OTU_5", "",
                                    "GG_OTU_4", "", "", ""};

int main(void) {
    const support_biom_t   table = {
        (char**) obs_ids, (char**) samp_ids,
        (uint32_t*) indices, (uint32_t*) indptr,
        (double*) data,
        n_obs, n_samp, 0
    };
    const support_bptree_t tree  = {
        (bool*) structure, (double*) lengths, (char**) names, (int) nparens
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
