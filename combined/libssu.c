/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2023-2025, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

/*
 * Implement wrappers around all the EXTERN functions
 * defined in api.hpp.
 *
 */

/*********************************************************************/

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* Pick the right libssu implementation */
/* Must be before the include of ssu_ld */
#ifndef BASIC_ONLY
static const char *ssu_get_lib_name() {
#if defined(__x86_64__)
   __builtin_cpu_init ();
   bool has_v2  = __builtin_cpu_supports ("x86-64-v2");
   bool has_v3  = __builtin_cpu_supports ("x86-64-v3");
   bool has_v4  = __builtin_cpu_supports ("x86-64-v4");
#else
   bool has_v2  = false;
   bool has_v3  = false;
   bool has_v4  = false;
#endif
   const char* env_max_cpu = getenv("UNIFRAC_MAX_CPU");

   if (env_max_cpu!=NULL) {
    if (strcmp(env_max_cpu,"basic")==0) {
      has_v2 = false;
      has_v3 = false;
      has_v4 = false;
    } else if ((strcmp(env_max_cpu,"x86-64-v2")==0) || (strcmp(env_max_cpu,"avx")==0)) {
      has_v3 = false;
      has_v4 = false;
    } else if ((strcmp(env_max_cpu,"x86-64-v3")==0) || (strcmp(env_max_cpu,"avx2")==0)) {
      has_v4 = false;
    }
   }

   const char *ssu = "unknown"; // just to have a default
   if (has_v4) {
     ssu = "libssu_cpu_x86_v4.so";
   } else if (has_v3) {
     ssu = "libssu_cpu_x86_v3.so";
   } else if (has_v2) {
     ssu = "libssu_cpu_x86_v2.so";
   } else {
     ssu = "libssu_cpu_basic.so";
   }

   return ssu;
}

#else
static const char *ssu_get_lib_name() {
   const char *ssu = "libssu_cpu_basic.so";

   const char* env_gpu_info = getenv("UNIFRAC_GPU_INFO");
   if ((env_gpu_info!=NULL) && (env_gpu_info[0]=='Y')) {
         printf("INFO (unifrac): No GPU support in this version\n");
   }

   return ssu;
}


#endif

/* Import the actual implementation of ld handling */
#include "../src/ssu_ld.c"

/* And we will need the headers for implemenation of the wrappers */
#include "../src/api.hpp"

/*********************************************************************/
/* All the functons below are wrappers
 * and each has its own function pointer
 * that is initialized on first use */
/*********************************************************************/

static void (*dl_ssu_set_random_seed)(unsigned int) = NULL;
void ssu_set_random_seed(unsigned int new_seed) {
   cond_ssu_load("ssu_set_random_seed", (void **) &dl_ssu_set_random_seed);

   (*dl_ssu_set_random_seed)(new_seed);
}

/*********************************************************************/

static void (*dl_destroy_mat)(mat_t**) = NULL;
static void (*dl_destroy_mat_full_fp64)(mat_full_fp64_t**) = NULL;
static void (*dl_destroy_mat_full_fp32)(mat_full_fp32_t**) = NULL;
static void (*dl_destroy_partial_mat)(partial_mat_t**) = NULL;
static void (*dl_destroy_partial_dyn_mat)(partial_dyn_mat_t**) = NULL;
static void (*dl_destroy_results_vec)(r_vec**) = NULL;
static void (*dl_destroy_bptree_opaque)(opaque_bptree_t**) = NULL;

static IOStatus (*dl_read_bptree_opaque)(const char*, opaque_bptree_t**) = NULL;
static void (*dl_load_bptree_opaque)(const char*, opaque_bptree_t**) = NULL;
static void (*dl_convert_bptree_opaque)(const support_bptree_t*, opaque_bptree_t**) = NULL;
static int (*dl_get_bptree_opaque_els)(opaque_bptree_t*) = NULL;

void destroy_mat(mat_t** result) {
   cond_ssu_load("destroy_mat", (void **) &dl_destroy_mat);

   (*dl_destroy_mat)(result);
}

void destroy_mat_full_fp64(mat_full_fp64_t** result) {
   cond_ssu_load("destroy_mat_full_fp64", (void **) &dl_destroy_mat_full_fp64);

   (*dl_destroy_mat_full_fp64)(result);
}

void destroy_mat_full_fp32(mat_full_fp32_t** result) {
   cond_ssu_load("destroy_mat_full_fp32", (void **) &dl_destroy_mat_full_fp32);

   (*dl_destroy_mat_full_fp32)(result);
}

void destroy_partial_mat(partial_mat_t** result) {
   cond_ssu_load("destroy_partial_mat", (void **) &dl_destroy_partial_mat);

   (*dl_destroy_partial_mat)(result);
}

void destroy_partial_dyn_mat(partial_dyn_mat_t** result) {
   cond_ssu_load("destroy_partial_dyn_mat", (void **) &dl_destroy_partial_dyn_mat);

   (*dl_destroy_partial_dyn_mat)(result);
}

void destroy_results_vec(r_vec** result) {
   cond_ssu_load("destroy_results_vec", (void **) &dl_destroy_results_vec);

   (*dl_destroy_results_vec)(result);
}

void destroy_bptree_opaque(opaque_bptree_t** tree_data) {
   cond_ssu_load("destroy_bptree_opaque", (void **) &dl_destroy_bptree_opaque);

   (*dl_destroy_bptree_opaque)(tree_data);
}

IOStatus read_bptree_opaque(const char* tree_filename, opaque_bptree_t** tree_data) {
   cond_ssu_load("read_bptree_opaque", (void **) &dl_read_bptree_opaque);

   return (*dl_read_bptree_opaque)(tree_filename,tree_data);
}

void load_bptree_opaque(const char* newick, opaque_bptree_t** tree_data) {
   cond_ssu_load("load_bptree_opaque", (void **) &dl_load_bptree_opaque);

   (*dl_load_bptree_opaque)(newick,tree_data);
}

void convert_bptree_opaque(const support_bptree_t* in_tree, opaque_bptree_t** tree_data) {
   cond_ssu_load("convert_bptree_opaque", (void **) &dl_convert_bptree_opaque);

   (*dl_convert_bptree_opaque)(in_tree,tree_data);
}

int get_bptree_opaque_els(opaque_bptree_t* tree_data) {
   cond_ssu_load("get_bptree_opaque_els", (void **) &dl_get_bptree_opaque_els);

   return (*dl_get_bptree_opaque_els)(tree_data);
}

/*********************************************************************/

static ComputeStatus (*dl_one_off_v3)(const char*, const char*, const char*, bool, double, bool, bool, unsigned int, mat_t**) = NULL;
static ComputeStatus (*dl_one_off_wtree_v3)(const char*, const opaque_bptree_t*, const char*, bool, double, bool, bool, unsigned int, mat_t**) = NULL;

ComputeStatus one_off_v3(const char* biom_filename, const char* tree_filename,
                             const char* unifrac_method, bool variance_adjust, double alpha,
                             bool bypass_tips, bool normalize_sample_counts, unsigned int n_substeps, mat_t** result) {
   cond_ssu_load("one_off_v3", (void **) &dl_one_off_v3);

   return (*dl_one_off_v3)(biom_filename, tree_filename, unifrac_method, variance_adjust, alpha, bypass_tips, normalize_sample_counts, n_substeps, result);
}

ComputeStatus one_off_wtree_v3(const char* biom_filename, const opaque_bptree_t* tree_data,
                                   const char* unifrac_method, bool variance_adjust, double alpha,
                                   bool bypass_tips, bool normalize_sample_counts, unsigned int n_substeps, mat_t** result) {
   cond_ssu_load("one_off_wtree_v3", (void **) &dl_one_off_wtree_v3);

   return (*dl_one_off_wtree_v3)(biom_filename, tree_data, unifrac_method, variance_adjust, alpha, bypass_tips, normalize_sample_counts, n_substeps, result);
}

/*********************************************************************/

static ComputeStatus (*dl_one_off_matrix_inmem_v3)(const support_biom_t *, const support_bptree_t *, const char*, bool, double, 
                                                   bool, bool, unsigned int, unsigned int, bool, const char *, mat_full_fp64_t**) = NULL;
static ComputeStatus (*dl_one_off_matrix_inmem_fp32_v3)(const support_biom_t *, const support_bptree_t *, const char*, bool, double,
                                                        bool, bool, unsigned int, unsigned int, bool, const char *, mat_full_fp32_t**) = NULL;

ComputeStatus one_off_matrix_inmem_v3(const support_biom_t *table_data, const support_bptree_t *tree_data,
                                             const char* unifrac_method, bool variance_adjust, double alpha,
                                             bool bypass_tips, bool normalize_sample_counts, unsigned int n_substeps,
                                             unsigned int subsample_depth, bool subsample_with_replacement, const char *mmap_dir,
                                             mat_full_fp64_t** result) {
   cond_ssu_load("one_off_matrix_inmem_v3", (void **) &dl_one_off_matrix_inmem_v3);

   return (*dl_one_off_matrix_inmem_v3)(table_data, tree_data, unifrac_method, variance_adjust, alpha,
                                 bypass_tips, normalize_sample_counts, n_substeps, subsample_depth, subsample_with_replacement, mmap_dir, result);
}

ComputeStatus one_off_matrix_inmem_fp32_v3(const support_biom_t *table_data, const support_bptree_t *tree_data,
                                                  const char* unifrac_method, bool variance_adjust, double alpha,
                                                  bool bypass_tips, bool normalize_sample_counts, unsigned int n_substeps,
                                                  unsigned int subsample_depth, bool subsample_with_replacement, const char *mmap_dir,
                                                  mat_full_fp32_t** result) {
   cond_ssu_load("one_off_matrix_inmem_fp32_v3", (void **) &dl_one_off_matrix_inmem_fp32_v3);

   return (*dl_one_off_matrix_inmem_fp32_v3)(table_data, tree_data, unifrac_method, variance_adjust, alpha,
                                      bypass_tips, normalize_sample_counts, n_substeps, subsample_depth, subsample_with_replacement, mmap_dir, result);
}

/*********************************************************************/

static ComputeStatus (*dl_one_off_matrix_v3)(const char*, const char*, const char*, bool, double,
                                             bool, bool, unsigned int, unsigned int, bool, const char *, mat_full_fp64_t**) = NULL;
static ComputeStatus (*dl_one_off_matrix_v3t)(const char*, const opaque_bptree_t*, const char*, bool, double,
                                             bool, bool, unsigned int, unsigned int, bool, const char *, mat_full_fp64_t**) = NULL;
static ComputeStatus (*dl_one_off_matrix_fp32_v3)(const char*, const char*, const char*, bool, double,
                                                  bool, bool, unsigned int, unsigned int, bool, const char *, mat_full_fp32_t**) = NULL;
static ComputeStatus (*dl_one_off_matrix_fp32_v3t)(const char*, const opaque_bptree_t*, const char*, bool, double,
                                                  bool, bool, unsigned int, unsigned int, bool, const char *, mat_full_fp32_t**) = NULL;

ComputeStatus one_off_matrix_v3(const char* biom_filename, const char* tree_filename,
                                       const char* unifrac_method, bool variance_adjust, double alpha,
                                       bool bypass_tips, bool normalize_sample_counts, unsigned int n_substeps,
                                       unsigned int subsample_depth, bool subsample_with_replacement, const char *mmap_dir,
                                       mat_full_fp64_t** result) {
   cond_ssu_load("one_off_matrix_v3", (void **) &dl_one_off_matrix_v3);

   return (*dl_one_off_matrix_v3)(biom_filename, tree_filename, unifrac_method, variance_adjust, alpha,
                           bypass_tips, normalize_sample_counts, n_substeps, subsample_depth, subsample_with_replacement, mmap_dir, result);
}

ComputeStatus one_off_matrix_v3t(const char* biom_filename, const opaque_bptree_t* tree_data,
                                        const char* unifrac_method, bool variance_adjust, double alpha,
                                        bool bypass_tips, bool normalize_sample_counts, unsigned int n_substeps,
                                        unsigned int subsample_depth, bool subsample_with_replacement, const char *mmap_dir,
                                        mat_full_fp64_t** result) {
   cond_ssu_load("one_off_matrix_v3t", (void **) &dl_one_off_matrix_v3t);

   return (*dl_one_off_matrix_v3t)(biom_filename, tree_data, unifrac_method, variance_adjust, alpha,
                           bypass_tips, normalize_sample_counts, n_substeps, subsample_depth, subsample_with_replacement, mmap_dir, result);
}

ComputeStatus one_off_matrix_fp32_v3(const char* biom_filename, const char* tree_filename,
                                            const char* unifrac_method, bool variance_adjust, double alpha,
                                            bool bypass_tips, bool normalize_sample_counts, unsigned int n_substeps,
                                            unsigned int subsample_depth, bool subsample_with_replacement, const char *mmap_dir,
                                            mat_full_fp32_t** result) {
   cond_ssu_load("one_off_matrix_fp32_v3", (void **) &dl_one_off_matrix_fp32_v3);

   return (*dl_one_off_matrix_fp32_v3)(biom_filename, tree_filename, unifrac_method, variance_adjust, alpha,
                                bypass_tips, normalize_sample_counts, n_substeps, subsample_depth, subsample_with_replacement, mmap_dir, result);
}

ComputeStatus one_off_matrix_fp32_v3t(const char* biom_filename, const opaque_bptree_t* tree_data,
                                             const char* unifrac_method, bool variance_adjust, double alpha,
                                             bool bypass_tips, bool normalize_sample_counts, unsigned int n_substeps,
                                             unsigned int subsample_depth, bool subsample_with_replacement, const char *mmap_dir,
                                             mat_full_fp32_t** result) {
   cond_ssu_load("one_off_matrix_fp32_v3t", (void **) &dl_one_off_matrix_fp32_v3t);

   return (*dl_one_off_matrix_fp32_v3t)(biom_filename, tree_data, unifrac_method, variance_adjust, alpha,
                                bypass_tips, normalize_sample_counts, n_substeps, subsample_depth, subsample_with_replacement, mmap_dir, result);
}

/*********************************************************************/

static ComputeStatus (*dl_faith_pd_one_off)(const char*, const char*, r_vec**) = NULL;
ComputeStatus faith_pd_one_off(const char* biom_filename, const char* tree_filename,
                                      r_vec** result) {
   cond_ssu_load("faith_pd_one_off", (void **) &dl_faith_pd_one_off);

   return (*dl_faith_pd_one_off)(biom_filename, tree_filename, result);
}

/*********************************************************************/

static ComputeStatus (*dl_unifrac_to_txt_file_v3)(const char*, const char*, const char*,
                                              const char*, bool, double,
                                              bool, bool, unsigned int,
                                              const char *) = NULL;

ComputeStatus unifrac_to_txt_file_v3(const char* biom_filename, const char* tree_filename, const char* out_filename,
                                        const char* unifrac_method, bool variance_adjust, double alpha,
                                        bool bypass_tips, bool normalize_sample_counts, unsigned int n_substeps,
                                        const char *mmap_dir){
   cond_ssu_load("unifrac_to_txt_file_v3", (void **) &dl_unifrac_to_txt_file_v3);

   return (*dl_unifrac_to_txt_file_v3)(biom_filename, tree_filename, out_filename, unifrac_method, variance_adjust, alpha,
                            bypass_tips, normalize_sample_counts, n_substeps,
                            mmap_dir);
}

/*********************************************************************/

static ComputeStatus (*dl_unifrac_to_file_v3)(const char*, const char*, const char*, const char*, bool, double,
                                              bool, bool, unsigned int, const char*, unsigned int, bool, 
                                              unsigned int, unsigned int, const char *, const char *, const char *) = NULL;
static ComputeStatus (*dl_unifrac_multi_to_file_v3)(const char*, const char*, const char*, const char*, bool, double,
                                              bool, bool, unsigned int, const char*, unsigned int, unsigned int, bool, 
                                              unsigned int, unsigned int, const char *, const char *, const char *) = NULL;

ComputeStatus unifrac_to_file_v3(const char* biom_filename, const char* tree_filename, const char* out_filename,
                                        const char* unifrac_method, bool variance_adjust, double alpha,
                                        bool bypass_tips, bool normalize_sample_counts, unsigned int n_substeps, const char* format,
                                        unsigned int subsample_depth, bool subsample_with_replacement, 
                                        unsigned int pcoa_dims,
                                        unsigned int permanova_perms, const char *grouping_filename, const char *grouping_columns,
                                        const char *mmap_dir){
   cond_ssu_load("unifrac_to_file_v3", (void **) &dl_unifrac_to_file_v3);

   return (*dl_unifrac_to_file_v3)(biom_filename, tree_filename, out_filename, unifrac_method, variance_adjust, alpha,
                            bypass_tips, normalize_sample_counts, n_substeps, format, subsample_depth, subsample_with_replacement,
                            pcoa_dims, permanova_perms, grouping_filename, grouping_columns, mmap_dir);
}

ComputeStatus unifrac_multi_to_file_v3(const char* biom_filename, const char* tree_filename, const char* out_filename,
                                              const char* unifrac_method, bool variance_adjust, double alpha,
                                              bool bypass_tips, bool normalize_sample_counts, unsigned int n_substeps, const char* format,
                                              unsigned int n_subsamples, unsigned int subsample_depth, bool subsample_with_replacement, 
                                              unsigned int pcoa_dims,
                                              unsigned int permanova_perms, const char *grouping_filename, const char *grouping_columns,
                                              const char *mmap_dir) {
   cond_ssu_load("unifrac_multi_to_file_v3", (void **) &dl_unifrac_multi_to_file_v3);

   return (*dl_unifrac_multi_to_file_v3)(biom_filename, tree_filename, out_filename, unifrac_method, variance_adjust, alpha,
                                  bypass_tips, normalize_sample_counts, n_substeps, format, n_subsamples, subsample_depth, subsample_with_replacement,
                                  pcoa_dims, permanova_perms, grouping_filename, grouping_columns, mmap_dir);
}


/*********************************************************************/

static ComputeStatus (*dl_compute_permanova_fp64)(const char *, unsigned int, const char**, mat_full_fp64_t *, unsigned int, double *, double *) = NULL;
static ComputeStatus (*dl_compute_permanova_fp32)(const char *, unsigned int, const char**, mat_full_fp32_t *, unsigned int, float *, float *) = NULL;

ComputeStatus compute_permanova_fp64(const char *grouping_filename, unsigned int n_columns, const char* *columns,
                                            mat_full_fp64_t * result, unsigned int permanova_perms,
                                            double *fstats, double *pvalues) {
   cond_ssu_load("compute_permanova_fp64", (void **) &dl_compute_permanova_fp64);

   return (*dl_compute_permanova_fp64)(grouping_filename, n_columns, columns, result, permanova_perms, fstats, pvalues);
}

ComputeStatus compute_permanova_fp32(const char *grouping_filename, unsigned int n_columns, const char* * columns,
                                            mat_full_fp32_t * result, unsigned int permanova_perms,
                                            float *fstats, float *pvalues) {
   cond_ssu_load("compute_permanova_fp32", (void **) &dl_compute_permanova_fp32);

   return (*dl_compute_permanova_fp32)(grouping_filename, n_columns, columns, result, permanova_perms, fstats, pvalues);
}

/*********************************************************************/

static IOStatus (*dl_write_mat)(const char*, mat_t*) = NULL;
static IOStatus (*dl_write_mat_from_matrix)(const char*, mat_full_fp64_t*) = NULL;
static IOStatus (*dl_write_mat_from_matrix_fp32)(const char*, mat_full_fp32_t*) = NULL;
static IOStatus (*dl_write_vec)(const char*, r_vec*) = NULL;

IOStatus write_mat(const char* filename, mat_t* result) {
   cond_ssu_load("write_mat", (void **) &dl_write_mat);

   return (*dl_write_mat)(filename, result);
}

IOStatus write_mat_from_matrix(const char* filename, mat_full_fp64_t* result) {
   cond_ssu_load("write_mat_from_matrix", (void **) &dl_write_mat_from_matrix);

   return (*dl_write_mat_from_matrix)(filename, result);
}

IOStatus write_mat_from_matrix_fp32(const char* filename, mat_full_fp32_t* result) {
   cond_ssu_load("write_mat_from_matrix_fp32", (void **) &dl_write_mat_from_matrix_fp32);

   return (*dl_write_mat_from_matrix_fp32)(filename, result);
}

IOStatus write_vec(const char* filename, r_vec* result) {
   cond_ssu_load("write_vec", (void **) &dl_write_vec);

   return (*dl_write_vec)(filename, result);
}

/*********************************************************************/

static IOStatus (*dl_write_mat_from_matrix_hdf5_fp64_v2)(const char*, mat_full_fp64_t*, unsigned int, int, unsigned int,
                                                         const char* *, const char**, const double *, const double *, const unsigned int *,
                                                         const char**, const unsigned int *) = NULL;
static IOStatus (*dl_write_mat_from_matrix_hdf5_fp32_v2)(const char*, mat_full_fp32_t*, unsigned int, int, unsigned int,
                                                         const char**, const char**, const float *, const float *, const unsigned int *,
                                                         const char**, const unsigned int *) = NULL;

IOStatus write_mat_from_matrix_hdf5_fp64_v2(const char* output_filename, mat_full_fp64_t* result,
                                                   unsigned int pcoa_dims, int save_dist,
                                                   unsigned int stat_n_vals,
                                                   const char*  *stat_method_arr,     const char*        *stat_name_arr,
                                                   const double *stat_val_arr,        const double       *stat_pval_arr, const unsigned int *stat_perm_count_arr,
                                                   const char*  *stat_group_name_arr, const unsigned int *stat_group_count_arr) {
   cond_ssu_load("write_mat_from_matrix_hdf5_fp64_v2", (void **) &dl_write_mat_from_matrix_hdf5_fp64_v2);

   return (*dl_write_mat_from_matrix_hdf5_fp64_v2)(output_filename, result, pcoa_dims, save_dist, stat_n_vals,
                                            stat_method_arr, stat_name_arr, stat_val_arr, stat_pval_arr, stat_perm_count_arr,
                                            stat_group_name_arr, stat_group_count_arr);
}

IOStatus write_mat_from_matrix_hdf5_fp32_v2(const char* output_filename, mat_full_fp32_t* result,
                                                   unsigned int pcoa_dims, int save_dist,
                                                   unsigned int stat_n_vals,
                                                   const char*  *stat_method_arr,     const char*        *stat_name_arr,
                                                   const float  *stat_val_arr,        const float        *stat_pval_arr, const unsigned int *stat_perm_count_arr,
                                                   const char*  *stat_group_name_arr, const unsigned int *stat_group_count_arr) {
   cond_ssu_load("write_mat_from_matrix_hdf5_fp32_v2", (void **) &dl_write_mat_from_matrix_hdf5_fp32_v2);

   return (*dl_write_mat_from_matrix_hdf5_fp32_v2)(output_filename, result, pcoa_dims, save_dist, stat_n_vals,
                                            stat_method_arr, stat_name_arr, stat_val_arr, stat_pval_arr, stat_perm_count_arr,
                                            stat_group_name_arr, stat_group_count_arr);
}

/*********************************************************************/

static ComputeStatus (*dl_one_dense_pair_v3t)(unsigned int, const char **, const double*,const double*,const opaque_bptree_t*,const char*, bool, double, bool, bool, double*) = NULL;
static ComputeStatus (*dl_one_dense_pair_v3)(unsigned int, const char **, const double*,const double*,const support_bptree_t*,const char*, bool, double, bool, bool, double*) = NULL;

ComputeStatus one_dense_pair_v3t(unsigned int n_obs, const char ** obs_ids, const double* sample1, const double* sample2,
		                        const opaque_bptree_t* tree_data,
                                        const char* unifrac_method, bool variance_adjust, double alpha,
                                        bool bypass_tips, bool normalize_sample_counts, double* result) {
   cond_ssu_load("one_dense_pair_v3t", (void **) &dl_one_dense_pair_v3t);

   return (*dl_one_dense_pair_v3t)(n_obs,obs_ids,sample1,sample2,tree_data,unifrac_method,variance_adjust,alpha,bypass_tips,normalize_sample_counts,result);
}

ComputeStatus one_dense_pair_v3(unsigned int n_obs, const char ** obs_ids, const double* sample1, const double* sample2,
		                       const support_bptree_t* tree_data,
                                       const char* unifrac_method, bool variance_adjust, double alpha,
                                       bool bypass_tips, bool normalize_sample_counts, double* result) {
   cond_ssu_load("one_dense_pair_v3", (void **) &dl_one_dense_pair_v3);

   return (*dl_one_dense_pair_v3)(n_obs,obs_ids,sample1,sample2,tree_data,unifrac_method,variance_adjust,alpha,bypass_tips,normalize_sample_counts,result);
}

/*********************************************************************/

static ComputeStatus (*dl_partial_v3)(const char*, const char*, const char*, bool, double, bool, bool, unsigned int, unsigned int, unsigned int, partial_mat_t**) = NULL;
static MergeStatus (*dl_merge_partial_to_mmap_matrix)(partial_dyn_mat_t**, int, const char *, mat_full_fp64_t**) = NULL;
static MergeStatus (*dl_merge_partial_to_mmap_matrix_fp32)(partial_dyn_mat_t**, int, const char *, mat_full_fp32_t**) = NULL;
static MergeStatus (*dl_validate_partial)(const partial_dyn_mat_t* const *, int);
static IOStatus (*dl_read_partial)(const char*, partial_mat_t**);
static IOStatus (*dl_read_partial_header)(const char*, partial_dyn_mat_t**);
static IOStatus (*dl_read_partial_one_stripe)(partial_dyn_mat_t*, uint32_t);
static IOStatus (*dl_write_partial)(const char*, const partial_mat_t*);


ComputeStatus partial_v3(const char* biom_filename, const char* tree_filename,
                             const char* unifrac_method, bool variance_adjust, double alpha,
                             bool bypass_tips, bool normalize_sample_counts, unsigned int n_substeps, unsigned int stripe_start,
                             unsigned int stripe_stop, partial_mat_t** result) {
   cond_ssu_load("partial_v3", (void **) &dl_partial_v3);

   return (*dl_partial_v3)(biom_filename,tree_filename,unifrac_method,variance_adjust,alpha,
		   bypass_tips,normalize_sample_counts,n_substeps,stripe_start,stripe_stop,result);
}

MergeStatus merge_partial_to_mmap_matrix(partial_dyn_mat_t* * partial_mats, int n_partials, const char *mmap_dir, mat_full_fp64_t** result) {
   cond_ssu_load("merge_partial_to_mmap_matrix", (void **) &dl_merge_partial_to_mmap_matrix);

   return (*dl_merge_partial_to_mmap_matrix)(partial_mats,n_partials,mmap_dir,result);
}

MergeStatus merge_partial_to_mmap_matrix_fp32(partial_dyn_mat_t* * partial_mats, int n_partials, const char *mmap_dir, mat_full_fp32_t** result) {
   cond_ssu_load("merge_partial_to_mmap_matrix_fp32", (void **) &dl_merge_partial_to_mmap_matrix_fp32);

   return (*dl_merge_partial_to_mmap_matrix_fp32)(partial_mats,n_partials,mmap_dir,result);
}

MergeStatus validate_partial(const partial_dyn_mat_t* const * partial_mats, int n_partials) {
   cond_ssu_load("validate_partial", (void **) &dl_validate_partial);

   return (*dl_validate_partial)(partial_mats,n_partials);
}

IOStatus read_partial(const char* input_filename, partial_mat_t** result_out) {
   cond_ssu_load("read_partial", (void **) &dl_read_partial);

   return (*dl_read_partial)(input_filename,result_out);
}

IOStatus read_partial_header(const char* input_filename, partial_dyn_mat_t** result_out) {
   cond_ssu_load("read_partial_header", (void **) &dl_read_partial_header);

   return (*dl_read_partial_header)(input_filename,result_out);
}

IOStatus read_partial_one_stripe(partial_dyn_mat_t* result, uint32_t stripe_idx) {
   cond_ssu_load("read_partial_one_stripe", (void **) &dl_read_partial_one_stripe);

   return (*dl_read_partial_one_stripe)(result,stripe_idx);
}

IOStatus write_partial(const char* filename, const partial_mat_t* result) {
   cond_ssu_load("write_partial", (void **) &dl_write_partial);

   return (*dl_write_partial)(filename,result);
}

// compat versions

#include "../src/api_compat.hpp"


