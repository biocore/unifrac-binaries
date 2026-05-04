#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "api.hpp"

#ifndef bool
#define bool char
#define true 1
#define false 0
#endif

void err(bool condition, const char* msg) {
    if(condition) {
        fprintf(stderr, "%s\n", msg);
        exit(1);
    }
} 

// biom
const unsigned int n_obs = 5;
const unsigned int n_samp = 6;
const char* const obs_ids[] = {"GG_OTU_1", "GG_OTU_2", "GG_OTU_3", "GG_OTU_4", "GG_OTU_5"};
const char* const samp_ids[] = {"Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6"};
const uint32_t    indices[] = {2, 0, 1, 3, 4, 5, 2, 3, 5, 0, 1, 2, 5, 1, 2};
const uint32_t    indptr[] = {0,  1,  6,  9, 13, 15};
const double      data[]= {1., 5., 1., 2., 3., 1., 1., 4., 2., 2., 1., 1., 1., 1., 1.};

// tree
const unsigned int nparens = 16;
const double lengths[] = { 0, 0, 0, 0,
	                   0, 0, 0, 0,
	                   0, 0, 0, 0,
	                   0, 0, 0, 0 };
const bool structure[] = { true, true, false, true,
	                   true, false, true, false,
			   false, true, true, false,
			   true, false, false, false };
const char * const names[] = {"", "GG_OTU_1", "", "",
                              "GG_OTU_2", "", "GG_OTU_3", "",
			      "", "", "GG_OTU_5", "",
			      "GG_OTU_4", "", "", ""};

void test_su_dense(int num_cores){
    double result = 0.0;
    const support_bptree_t tree = {(bool*) structure, (double*) lengths, (char**) names, nparens};
    const char* table_oids[] = { "GG_OTU_1", "GG_OTU_2", "GG_OTU_3", "GG_OTU_4", "GG_OTU_5" };
    const double sample1[] = { 0, 5, 0, 2, 0 };
    const double sample2[] = { 0, 1, 0, 1, 1 };
    const double sample3[] = { 1, 0, 1, 1, 1 };
    const char* method = "unweighted";
    float exp13 = 0.57142857;
    float exp21 = 0.2;
    
    opaque_bptree_t *tree_data;
    convert_bptree_opaque(&tree, &tree_data);

    ComputeStatus status;
    status = one_dense_pair_v2t(5, table_oids, sample1, sample3,
	                        tree_data, method,
                                false, 1.0, false,
				&result);

    err(status != okay, "Compute failed");

    err(fabs(exp13 - result) > 0.00001, "Result is wrong");

    status = one_dense_pair_v2t(5, table_oids, sample2, sample1,
	                        tree_data, method,
                                false, 1.0, false,
				&result);

    err(status != okay, "Compute failed");

    err(fabs(exp21 - result) > 0.00001, "Result is wrong");

    destroy_bptree_opaque(&tree_data);

    // check direct, expanded tree structure, too
    status = one_dense_pair_v2(5, table_oids, sample1, sample3,
	                       &tree, method,
                               false, 1.0, false,
			       &result);

    err(status != okay, "Compute failed");

    err(fabs(exp13 - result) > 0.00001, "Result is wrong");

}

void test_su_matrix(int num_cores){
    mat_full_fp64_t* result = NULL;
    mat_full_fp32_t* result_fp32 = NULL;
    mat_full_fp64_t* result2 = NULL;
    mat_full_fp32_t* result2_fp32 = NULL;
    const support_biom_t table = {(char**) obs_ids, (char**) samp_ids, (uint32_t*) indices, (uint32_t*) indptr, (double*) data, n_obs, n_samp, 0};
    const support_bptree_t tree = {(bool*) structure, (double*) lengths, (char**) names, nparens};

    float exp[] = { 0.0,        0.2,        0.57142857, 0.6,        0.5,        0.2, 
                    0.2,        0.0,        0.42857143, 0.66666667, 0.6,        0.33333333,
		    0.57142857, 0.42857143, 0.0,        0.71428571, 0.85714286, 0.42857143, 
		    0.6,        0.66666667, 0.71428571, 0.0,        0.33333333, 0.4,
		    0.5,        0.6,        0.85714286, 0.33333333, 0.0,        0.6,
                    0.2,        0.33333333, 0.42857143, 0.4,        0.6,        0.0        };
    
    ComputeStatus status;

    status = one_off_matrix_inmem_v2(&table, &tree, "unweighted_fp64",
                                    false, 1.0, false, num_cores,
                                    0, true, NULL,
				    &result);

    err(status != okay, "Compute failed");
    err(result == NULL, "Empty result");
    err(result->n_samples != 6, "Wrong number of samples");

    for(unsigned int i = 0; i < (result->n_samples*result->n_samples); i++) {
        err(fabs(exp[i] - result->matrix[i]) > 0.00001, "Result is wrong");
    }


    destroy_mat_full_fp64(&result);

    status = one_off_matrix_inmem_fp32_v2(&table, &tree, "unweighted_fp32",
                                    false, 1.0, false, num_cores,
                                    0, true, NULL,
				    &result_fp32);

    err(status != okay, "Compute failed");
    err(result == NULL, "Empty result");
    err(result_fp32->n_samples != 6, "Wrong number of samples");

    for(unsigned int i = 0; i < (result_fp32->n_samples*result_fp32->n_samples); i++) {
        err(fabs(exp[i] - result_fp32->matrix[i]) > 0.00001, "Result is wrong");
    }


    destroy_mat_full_fp32(&result_fp32);

    // exercise the old interface, too
    status = one_off_inmem(&table, &tree, "unweighted_fp64",
                                  false, 1.0, false, num_cores,
				  &result2);

    err(status != okay, "Compute failed");
    err(result2 == NULL, "Empty result");
    err(result2->n_samples != 6, "Wrong number of samples");

    for(unsigned int i = 0; i < (result2->n_samples*result2->n_samples); i++) {
        err(fabs(exp[i] - result2->matrix[i]) > 0.00001, "Result is wrong");
    }


    destroy_mat_full_fp64(&result2);

    status = one_off_inmem_fp32(&table, &tree, "unweighted_fp32",
                                    false, 1.0, false, num_cores,
				    &result2_fp32);

    err(status != okay, "Compute failed");
    err(result2 == NULL, "Empty result");
    err(result2_fp32->n_samples != 6, "Wrong number of samples");

    for(unsigned int i = 0; i < (result2_fp32->n_samples*result2_fp32->n_samples); i++) {
        err(fabs(exp[i] - result2_fp32->matrix[i]) > 0.00001, "Result is wrong");
    }


    destroy_mat_full_fp32(&result2_fp32);

}


// Branch lengths used by the in-memory Faith PD / subsample / permanova tests.
// The shared `lengths[]` above is all-zero (which produces all-zero Faith PD);
// override to unit lengths so the tree carries meaningful path information.
const double lengths_unit[] = { 0., 1., 0., 1.,
                                1., 0., 1., 0.,
                                0., 1., 1., 0.,
                                1., 0., 0., 0. };

void test_faith_pd_inmem_capi(void) {
    const support_biom_t   table = {(char**) obs_ids, (char**) samp_ids,
                                    (uint32_t*) indices, (uint32_t*) indptr,
                                    (double*) data, n_obs, n_samp, 0};
    const support_bptree_t tree  = {(bool*) structure, (double*) lengths_unit,
                                    (char**) names, nparens};

    // Multifurcating tree with unit branches; expected per-sample PD computed
    // by hand from the OTU sets implied by the CSR table.
    const double exp_pd[] = {4., 5., 6., 3., 2., 5.};

    r_vec* result = NULL;
    ComputeStatus status = faith_pd_inmem(&table, &tree, &result);
    err(status != okay, "faith_pd_inmem failed");
    err(result == NULL, "faith_pd_inmem returned NULL");
    err(result->n_samples != n_samp, "faith_pd_inmem n_samples mismatch");
    for (unsigned int i = 0; i < n_samp; i++) {
        err(fabs(result->values[i] - exp_pd[i]) > 1e-6, "faith_pd_inmem value mismatch");
    }
    destroy_results_vec(&result);

    // Error path: NULL tree.
    r_vec* tmp = NULL;
    err(faith_pd_inmem(&table, NULL, &tmp) != tree_missing,
        "faith_pd_inmem should reject NULL tree");
    // Error path: NULL table.
    err(faith_pd_inmem(NULL, &tree, &tmp) != table_missing,
        "faith_pd_inmem should reject NULL table");
}

void test_subsample_inmem_capi(void) {
    const support_biom_t table = {(char**) obs_ids, (char**) samp_ids,
                                  (uint32_t*) indices, (uint32_t*) indptr,
                                  (double*) data, n_obs, n_samp, 0};

    // All six samples have counts >= 3 ({7,3,4,6,3,3}), so depth=3 retains
    // every sample. Each retained sample sums to exactly depth.
    const unsigned int depth = 3;

    ssu_set_random_seed(42);
    opaque_biom_inmem_t* sub = NULL;
    err(subsample_table_inmem(&table, depth, false, &sub) != okay,
        "subsample_table_inmem failed");
    err(sub == NULL, "subsample_table_inmem returned NULL handle");

    err(subsampled_n_samples(sub) != n_samp, "subsampled_n_samples mismatch");

    // Sum each sample column across the dense observation rows; must equal depth.
    unsigned int n_sub_obs = subsampled_n_obs(sub);
    unsigned int n_sub_samp = subsampled_n_samples(sub);
    double col_sums[6] = {0., 0., 0., 0., 0., 0.};
    double row_buf[6];
    for (unsigned int i = 0; i < n_sub_obs; i++) {
        const char* oid = subsampled_get_obs_id(sub, i);
        err(oid == NULL, "subsampled_get_obs_id returned NULL");
        for (unsigned int j = 0; j < n_sub_samp; j++) row_buf[j] = 0.0;
        err(!subsampled_get_obs_data(sub, oid, row_buf),
            "subsampled_get_obs_data failed for valid obs id");
        for (unsigned int j = 0; j < n_sub_samp; j++) col_sums[j] += row_buf[j];
    }
    for (unsigned int j = 0; j < n_sub_samp; j++) {
        err(fabs(col_sums[j] - (double) depth) > 1e-9,
            "subsample column sum != depth");
    }

    // Sample IDs round-trip.
    for (unsigned int j = 0; j < n_sub_samp; j++) {
        err(subsampled_get_sample_id(sub, j) == NULL,
            "subsampled_get_sample_id returned NULL for in-range idx");
    }
    err(subsampled_get_sample_id(sub, n_sub_samp) != NULL,
        "subsampled_get_sample_id should return NULL for out-of-range idx");
    err(subsampled_get_obs_data(sub, "NOT_A_REAL_OTU", row_buf),
        "subsampled_get_obs_data should reject unknown obs id");

    // Determinism: identical seed reproduces the same data byte-for-byte.
    opaque_biom_inmem_t* sub2 = NULL;
    ssu_set_random_seed(42);
    err(subsample_table_inmem(&table, depth, false, &sub2) != okay,
        "subsample_table_inmem (second call) failed");
    err(subsampled_n_obs(sub2) != n_sub_obs, "determinism: n_obs differs");
    for (unsigned int i = 0; i < n_sub_obs; i++) {
        const char* oid = subsampled_get_obs_id(sub, i);
        double r1[6] = {0.}, r2[6] = {0.};
        subsampled_get_obs_data(sub, oid, r1);
        subsampled_get_obs_data(sub2, oid, r2);
        for (unsigned int j = 0; j < n_sub_samp; j++) {
            err(r1[j] != r2[j], "determinism: subsampled cell differs");
        }
    }

    destroy_subsampled_inmem(&sub);
    err(sub != NULL, "destroy_subsampled_inmem did not NULL the handle");
    destroy_subsampled_inmem(&sub2);
}

void test_permanova_inmem_capi(int num_cores) {
    // Build a real (non-degenerate) unweighted UniFrac matrix to feed into
    // PERMANOVA. The shared `lengths[]` is all-zero, which would yield an
    // all-zero distance matrix and a degenerate F-statistic — use unit
    // branch lengths so the test exercises a meaningful pairwise structure.
    const support_biom_t   table = {(char**) obs_ids, (char**) samp_ids,
                                    (uint32_t*) indices, (uint32_t*) indptr,
                                    (double*) data, n_obs, n_samp, 0};
    const support_bptree_t tree  = {(bool*) structure, (double*) lengths_unit,
                                    (char**) names, nparens};

    mat_full_fp64_t* dm = NULL;
    err(one_off_matrix_inmem_v2(&table, &tree, "unweighted_fp64",
                                false, 1.0, false, num_cores,
                                0, true, NULL, &dm) != okay,
        "one_off_matrix_inmem_v2 failed");

    const uint32_t grouping[] = {0, 0, 1, 1, 1, 0};

    ssu_set_random_seed(42);
    double fstat = 0.0, pvalue = 0.0;
    err(compute_permanova_inmem_fp64(dm->matrix, dm->n_samples, grouping,
                                     999, &fstat, &pvalue) != okay,
        "compute_permanova_inmem_fp64 failed");
    err(!(fstat > 0.0), "permanova fstat must be positive");
    err(!(pvalue > 0.0 && pvalue <= 1.0), "permanova pvalue out of range");

    // Re-seed and rerun. Under multi-threaded OMP, parallel reductions can
    // differ by ULPs across runs, so don't require bit-exactness — just
    // require the result reproduces within a tight numeric tolerance.
    ssu_set_random_seed(42);
    double fstat2 = 0.0, pvalue2 = 0.0;
    err(compute_permanova_inmem_fp64(dm->matrix, dm->n_samples, grouping,
                                     999, &fstat2, &pvalue2) != okay,
        "compute_permanova_inmem_fp64 (rerun) failed");
    err(fabs(fstat2 - fstat) > 1e-6, "permanova fstat differs across reruns");
    err(fabs(pvalue2 - pvalue) > 1e-2, "permanova pvalue differs across reruns");

    // Error paths.
    err(compute_permanova_inmem_fp64(NULL, dm->n_samples, grouping, 9, &fstat, &pvalue) == okay,
        "permanova should reject NULL matrix");
    err(compute_permanova_inmem_fp64(dm->matrix, dm->n_samples, NULL, 9, &fstat, &pvalue) == okay,
        "permanova should reject NULL grouping");

    // fp32 variant on a fp32 matrix.
    mat_full_fp32_t* dm32 = NULL;
    err(one_off_matrix_inmem_fp32_v2(&table, &tree, "unweighted_fp32",
                                     false, 1.0, false, num_cores,
                                     0, true, NULL, &dm32) != okay,
        "one_off_matrix_inmem_fp32_v2 failed");

    ssu_set_random_seed(42);
    float fstat32 = 0.0f, pvalue32 = 0.0f;
    err(compute_permanova_inmem_fp32(dm32->matrix, dm32->n_samples, grouping,
                                     999, &fstat32, &pvalue32) != okay,
        "compute_permanova_inmem_fp32 failed");
    err(!(fstat32 > 0.0f), "permanova fp32 fstat must be positive");

    destroy_mat_full_fp64(&dm);
    destroy_mat_full_fp32(&dm32);
}

int main(int argc, char** argv) {
    int num_cores = strtol(argv[1], NULL, 10);

    printf("Testing Striped UniFrac one_dense_pair...\n");
    test_su_dense(num_cores);
    printf("Testing Striped UniFrac one_off matrix...\n");
    test_su_matrix(num_cores);
    printf("Testing faith_pd_inmem...\n");
    test_faith_pd_inmem_capi();
    printf("Testing subsample_table_inmem + accessors...\n");
    test_subsample_inmem_capi();
    printf("Testing compute_permanova_inmem_fp64/fp32...\n");
    test_permanova_inmem_capi(num_cores);
    printf("Tests passed.\n");
    return 0;
}

