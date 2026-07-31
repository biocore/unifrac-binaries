/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2026, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

/*
 * Concurrency test for libssu_inmem.a -- the native in-memory static archive
 * that embedders link (see inmem_build.mk). The equivalent test in
 * src/test_su.cpp covers the same entry points, but only as built for
 * libssu.so; this configuration differs in ways that matter here:
 *
 *   - UNIFRAC_WASM is defined, so no SIGUSR1 handler is installed and
 *     CPU_SETSIZE falls back to 32 (unifrac_internal.cpp) -- in a build that,
 *     unlike WASM, is genuinely multi-threaded. One consequence worth being
 *     explicit about: the report_status use-after-free that motivated this
 *     work is inert here, because no flag is ever set. What is being pinned
 *     is that concurrent computes agree with the serial answer, not that
 *     particular bug.
 *   - OpenMP is on, so each caller thread spawns its own team.
 *   - Only the in-memory API surface is declared; a file-based entry point
 *     slipping into this test would not compile.
 *
 * Fixtures are shared with the WASM suite (tests/wasm/fixtures.hpp): the same
 * 5-OTU / 6-sample CSR table whose unweighted matrix the native CI already
 * validates.
 *
 * Worker threads never call CHECK -- it exits the process, which would race the
 * other workers and lose the diagnostic. Each records into its own slot and the
 * main thread checks every slot after joining.
 */

#include "tests/wasm/check_macros.hpp"
#include "tests/wasm/fixtures.hpp"
#include "api.hpp"

#include <thread>
#include <vector>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

static const unsigned int N_THREADS = 4;
static const unsigned int N_ITERS   = 25;

/* Subsampling depth 2 keeps every sample: per-sample input totals are
 * {7, 3, 7, 4, 3, 4}. Same depth the WASM subsample test uses.
 */
static const unsigned int SUBSAMPLE_DEPTH = 2;
static const int          SUBSAMPLE_SEED  = 42;

struct outcome {
    unsigned int n_ok       = 0;  // computes that returned okay
    unsigned int n_status   = 0;  // computes that returned something else
    unsigned int n_mismatch = 0;  // computes that disagreed with the reference
};

static support_biom_t make_table() {
    support_biom_t table = {
        (char**) FIXTURE_OBS_IDS, (char**) FIXTURE_SAMP_IDS,
        (uint32_t*) FIXTURE_INDICES, (uint32_t*) FIXTURE_INDPTR,
        (double*) FIXTURE_DATA,
        FIXTURE_N_OBS, FIXTURE_N_SAMP, 0
    };
    return table;
}

static support_bptree_t make_tree() {
    support_bptree_t tree = {
        (bool*) FIXTURE_STRUCTURE, (double*) FIXTURE_LENGTHS,
        (char**) FIXTURE_NAMES, (int) FIXTURE_NPARENS
    };
    return tree;
}

/* One unweighted compute, collected into a vector. seed < 0 selects the
 * unsubsampled v3 path; seed >= 0 the subsampled v4 path.
 */
static ComputeStatus run_matrix(std::vector<float> &out,
                                unsigned int n_substeps,
                                int seed) {
    const support_biom_t   table = make_table();
    const support_bptree_t tree  = make_tree();

    mat_full_fp32_t* mat = nullptr;
    ComputeStatus rc = (seed < 0)
        ? one_off_matrix_inmem_fp32_v3(&table, &tree, "unweighted_fp32",
                                       false, 1.0, false, true, n_substeps,
                                       0, false,                       // no subsampling
                                       nullptr, &mat)
        : one_off_matrix_inmem_fp32_v4(&table, &tree, "unweighted_fp32",
                                       false, 1.0, false, true, n_substeps,
                                       SUBSAMPLE_DEPTH, false, seed,
                                       nullptr, &mat);
    if (rc != okay) return rc;

    const size_t n_els = size_t(mat->n_samples) * size_t(mat->n_samples);
    out.assign(mat->matrix, mat->matrix + n_els);
    destroy_mat_full_fp32(&mat);
    return okay;
}

/* The unweighted compute is deterministic, so a concurrent result must be
 * bit-identical to the serial one -- no tolerance needed.
 */
static void matrix_worker(int seed, const std::vector<float>* reference, outcome* out) {
    for (unsigned int i = 0; i < N_ITERS; i++) {
        std::vector<float> got;
        if (run_matrix(got, 1, seed) != okay) {
            out->n_status++;
            continue;
        }
        if (got != *reference) out->n_mismatch++;
        out->n_ok++;
    }
}

static void faith_pd_worker(outcome* out) {
    const support_biom_t   table = make_table();
    const support_bptree_t tree  = make_tree();
    const double expected[6] = {4., 5., 6., 3., 2., 5.};

    for (unsigned int i = 0; i < N_ITERS; i++) {
        r_vec* res = nullptr;
        if (faith_pd_inmem(&table, &tree, &res) != okay) {
            out->n_status++;
            continue;
        }
        if (res->n_samples != FIXTURE_N_SAMP) {
            out->n_mismatch++;
        } else {
            for (int j = 0; j < FIXTURE_N_SAMP; j++) {
                if (std::fabs(res->values[j] - expected[j]) > 1e-6) {
                    out->n_mismatch++;
                    break;
                }
            }
        }
        destroy_results_vec(&res);
        out->n_ok++;
    }
}

static void check_all(const std::vector<outcome> &results, const char* what) {
    for (unsigned int t = 0; t < results.size(); t++) {
        if (results[t].n_status != 0 || results[t].n_mismatch != 0 ||
            results[t].n_ok != N_ITERS) {
            std::fprintf(stderr,
                "FAIL %s: thread %u ok=%u status=%u mismatch=%u (expected ok=%u)\n",
                what, t, results[t].n_ok, results[t].n_status,
                results[t].n_mismatch, N_ITERS);
            std::exit(1);
        }
    }
}

template<typename F>
static void run_workers(F worker, const char* what) {
    std::vector<outcome> results(N_THREADS);
    std::vector<std::thread> workers;
    for (unsigned int t = 0; t < N_THREADS; t++)
        workers.emplace_back(worker, &results[t]);
    for (unsigned int t = 0; t < N_THREADS; t++)
        workers[t].join();
    check_all(results, what);
}

int main(void) {
#ifdef _OPENMP
    std::fprintf(stdout, "omp_get_max_threads=%d\n", omp_get_max_threads());
#else
#error "libssu_inmem.a is built with OpenMP; this test must be too"
#endif

    // ---- unsubsampled, concurrent -------------------------------------
    std::vector<float> reference;
    CHECK(run_matrix(reference, 1, -1) == okay);
    CHECK(reference.size() == size_t(FIXTURE_N_SAMP) * size_t(FIXTURE_N_SAMP));

    // the serial answer is the one the native CI already validates
    for (int i = 0; i < FIXTURE_N_SAMP * FIXTURE_N_SAMP; i++)
        CHECK(std::fabs(reference[i] - FIXTURE_UNWEIGHTED_DIST[i]) < 1e-6);

    run_workers([&reference](outcome* o) { matrix_worker(-1, &reference, o); },
                "concurrent one_off_matrix_inmem_fp32_v3");

    // ---- subsampled with a per-call seed, concurrent -------------------
    /* Bit-exactness here relies on every call seeing the same OpenMP width,
     * since the subsample draw is distributed across the team. That holds
     * within one process -- nthreads-var is a per-thread ICV and no one
     * changes it -- but is not a claim of reproducibility across widths.
     */
    std::vector<float> seeded_reference;
    CHECK(run_matrix(seeded_reference, 1, SUBSAMPLE_SEED) == okay);
    CHECK(seeded_reference.size() == size_t(FIXTURE_N_SAMP) * size_t(FIXTURE_N_SAMP));

    run_workers([&seeded_reference](outcome* o) {
                    matrix_worker(SUBSAMPLE_SEED, &seeded_reference, o);
                },
                "concurrent one_off_matrix_inmem_fp32_v4 (seeded)");

    // ---- faith_pd, concurrent ------------------------------------------
    run_workers([](outcome* o) { faith_pd_worker(o); },
                "concurrent faith_pd_inmem");

    // ---- n_substeps outside the stripe range ---------------------------
    /* 6 samples -> 3 stripes. 0 used to divide by zero and 4 used to run off
     * the end of the stripe array; both are clamped now, and neither changes
     * the answer.
     */
    for (unsigned int n_substeps : {0u, 1u, 3u, 4u}) {
        std::vector<float> got;
        CHECK(run_matrix(got, n_substeps, -1) == okay);
        CHECK(got == reference);
    }

    std::fprintf(stdout, "OK test_concurrency_inmem\n");
    return 0;
}
