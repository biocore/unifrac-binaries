/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2026, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

/*
 * Stage 4 WASM test for subsample_table_inmem.
 *
 * No external ground-truth oracle is needed — assertions are structural
 * (per-sample rarefied count == depth, conservation, seed determinism).
 *
 * Fixture: same 5-OTU / 6-sample CSR table from test/capi_inmem_test.c
 * used by Stage 3. Per-sample input totals are
 *   {7, 3, 7, 4, 3, 4}; subsampling depth 2 keeps all 6 samples.
 *
 * Determinism is governed by the global skbb-side mt19937; we re-seed
 * before each subsampling call via ssu_set_random_seed and assert that
 * two calls with the same seed produce bit-identical per-sample dense
 * vectors.
 */

#include "tests/wasm/check_macros.hpp"
#include "api.hpp"

#include <cstdint>
#include <cstring>

static const int n_obs  = 5;
static const int n_samp = 6;
static const char* const obs_ids[]   = {"GG_OTU_1", "GG_OTU_2", "GG_OTU_3", "GG_OTU_4", "GG_OTU_5"};
static const char* const samp_ids[]  = {"Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6"};
static const uint32_t    indices[]   = {2, 0, 1, 3, 4, 5, 2, 3, 5, 0, 1, 2, 5, 1, 2};
static const uint32_t    indptr[]    = {0, 1, 6, 9, 13, 15};
static const double      data[]      = {1., 5., 1., 2., 3., 1., 1., 4., 2., 2., 1., 1., 1., 1., 1.};

static void run_subsample(opaque_biom_inmem_t **out, unsigned int seed) {
    ssu_set_random_seed(seed);
    const support_biom_t table = {
        (char**) obs_ids, (char**) samp_ids,
        (uint32_t*) indices, (uint32_t*) indptr,
        (double*) data,
        n_obs, n_samp, 0
    };
    ComputeStatus status = subsample_table_inmem(&table, /*depth=*/2u,
                                                 /*with_replacement=*/false,
                                                 out);
    CHECK(status == okay);
    CHECK(*out != nullptr);
}

// Sum each surviving sample's counts across all OTUs into per_sample_sum.
// Returns the number of surviving samples (== subsampled_n_samples).
static unsigned int collect_per_sample_sums(const opaque_biom_inmem_t *sub,
                                            double *per_sample_sum) {
    unsigned int n_s = subsampled_n_samples(sub);
    unsigned int n_o = subsampled_n_obs(sub);
    for (unsigned int j = 0; j < n_s; j++) per_sample_sum[j] = 0.0;

    double row[64];                                // n_samp <= 64 fixture-bound
    CHECK(n_s <= 64);
    for (unsigned int i = 0; i < n_o; i++) {
        const char *obs_id = subsampled_get_obs_id(sub, i);
        CHECK(obs_id != nullptr);
        bool ok = subsampled_get_obs_data(sub, obs_id, row);
        CHECK(ok);
        for (unsigned int j = 0; j < n_s; j++) {
            per_sample_sum[j] += row[j];
        }
    }
    return n_s;
}

// Pull the full dense subsampled table into a flat buffer (row-major,
// n_obs rows × n_samples cols) for byte-exact comparison across runs.
static void collect_dense(const opaque_biom_inmem_t *sub,
                          double *flat, unsigned int max_cells) {
    unsigned int n_s = subsampled_n_samples(sub);
    unsigned int n_o = subsampled_n_obs(sub);
    CHECK(n_s * n_o <= max_cells);
    double row[64];
    for (unsigned int i = 0; i < n_o; i++) {
        const char *obs_id = subsampled_get_obs_id(sub, i);
        CHECK(obs_id != nullptr);
        CHECK(subsampled_get_obs_data(sub, obs_id, row));
        for (unsigned int j = 0; j < n_s; j++) {
            flat[i * n_s + j] = row[j];
        }
    }
}

int main(void) {
    constexpr unsigned int DEPTH    = 2;
    constexpr unsigned int SEED     = 42;
    constexpr unsigned int FLAT_CAP = 1024;

    // (1) Run subsampling once. All 6 input samples have totals >= depth,
    //     so all 6 survive.
    opaque_biom_inmem_t *sub_a = nullptr;
    run_subsample(&sub_a, SEED);
    CHECK_EQ(subsampled_n_samples(sub_a), 6);
    CHECK(subsampled_n_obs(sub_a) >= 1);
    CHECK(subsampled_n_obs(sub_a) <= 5);

    // (2) Each surviving sample's rarefied count must equal DEPTH.
    {
        double per_sample[6] = {0};
        unsigned int n_s = collect_per_sample_sums(sub_a, per_sample);
        CHECK_EQ(n_s, 6);
        for (unsigned int j = 0; j < n_s; j++) {
            if (per_sample[j] != (double) DEPTH) {
                std::fprintf(stderr,
                    "FAIL sample %u count = %.1f; expected %u\n",
                    j, per_sample[j], DEPTH);
                std::exit(1);
            }
        }
    }

    // (3) Re-seed and run again; output must be bit-identical.
    opaque_biom_inmem_t *sub_b = nullptr;
    run_subsample(&sub_b, SEED);
    CHECK_EQ(subsampled_n_samples(sub_b), subsampled_n_samples(sub_a));
    CHECK_EQ(subsampled_n_obs(sub_b),     subsampled_n_obs(sub_a));

    {
        double flat_a[FLAT_CAP], flat_b[FLAT_CAP];
        collect_dense(sub_a, flat_a, FLAT_CAP);
        collect_dense(sub_b, flat_b, FLAT_CAP);
        unsigned int cells = subsampled_n_samples(sub_a) * subsampled_n_obs(sub_a);
        for (unsigned int k = 0; k < cells; k++) {
            if (flat_a[k] != flat_b[k]) {
                std::fprintf(stderr,
                    "FAIL determinism break at cell %u: %.6f vs %.6f\n",
                    k, flat_a[k], flat_b[k]);
                std::exit(1);
            }
        }
    }

    // (4) A different seed must (almost certainly) produce different
    //     output. We don't strictly need this — the bit-equality check
    //     above already proves seed plumbing — but it guards against a
    //     no-op stub silently passing both branches.
    opaque_biom_inmem_t *sub_c = nullptr;
    run_subsample(&sub_c, SEED + 1);
    {
        double flat_a[FLAT_CAP], flat_c[FLAT_CAP];
        collect_dense(sub_a, flat_a, FLAT_CAP);
        collect_dense(sub_c, flat_c, FLAT_CAP);
        unsigned int cells = subsampled_n_samples(sub_a) * subsampled_n_obs(sub_a);
        unsigned int diffs = 0;
        for (unsigned int k = 0; k < cells; k++) {
            if (flat_a[k] != flat_c[k]) diffs++;
        }
        // Our fixture is small; "different" means at least 1 differing
        // cell across the whole 5x6 table.
        CHECK(diffs >= 1);
    }

    destroy_subsampled_inmem(&sub_a);
    destroy_subsampled_inmem(&sub_b);
    destroy_subsampled_inmem(&sub_c);
    CHECK(sub_a == nullptr);
    CHECK(sub_b == nullptr);
    CHECK(sub_c == nullptr);

    std::fprintf(stdout, "OK test_subsample_wasm\n");
    return 0;
}
