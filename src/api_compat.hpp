/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2025, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

/*
 * This file contains deprecated functions that we still want to support for backwards compatibility
 * but are really just wrappers around more modern versions of the same.
 *
 * Meant to be included alongide the modern implementation source code.
 *
 * What belongs here versus in api.cpp: a forwarder whose old name is
 * *superseded* by a newer version of the same entry point goes here. A
 * non-seeded peer that remains first-class -- subsample_table_inmem,
 * compute_permanova_inmem_fp64/fp32, pcoa* -- stays in api.cpp. The difference
 * is not cosmetic: ../combined/libssu.c includes this header, so anything
 * defined here needs no dlsym stub there, while anything in api.cpp needs one.
 *
 */

#ifndef UNIFRAC_WASM
ComputeStatus partial(const char* biom_filename, const char* tree_filename,
                       const char* unifrac_method, bool variance_adjust, double alpha, bool bypass_tips,
                       unsigned int n_substeps, unsigned int stripe_start, unsigned int stripe_stop,
                       partial_mat_t** result) {
  const bool normalize_sample_counts = true;
  return partial_v3(biom_filename, tree_filename, unifrac_method,
		    variance_adjust, alpha, bypass_tips, normalize_sample_counts, n_substeps,
		    stripe_start, stripe_stop, result);
}

ComputeStatus one_off(const char* biom_filename, const char* tree_filename,
                       const char* unifrac_method, bool variance_adjust, double alpha,
                       bool bypass_tips, unsigned int n_substeps, mat_t** result) {
    const bool normalize_sample_counts = true;
    return one_off_v3(biom_filename,tree_filename,unifrac_method,variance_adjust,alpha,bypass_tips,normalize_sample_counts,n_substeps,result);
}

ComputeStatus one_off_wtree(const char* biom_filename, const opaque_bptree_t* tree_data,
                             const char* unifrac_method, bool variance_adjust, double alpha,
                             bool bypass_tips, unsigned int n_substeps, mat_t** result) {
    const bool normalize_sample_counts = true;
    return one_off_wtree_v3(biom_filename,tree_data,unifrac_method,variance_adjust,alpha,bypass_tips,normalize_sample_counts,n_substeps,result);
}

ComputeStatus one_off_matrix_v2(const char* biom_filename, const char* tree_filename,
                                 const char* unifrac_method, bool variance_adjust, double alpha,
                                 bool bypass_tips, unsigned int n_substeps,
                                 unsigned int subsample_depth, bool subsample_with_replacement, const char *mmap_dir,
                                 mat_full_fp64_t** result) {
    bool normalize_sample_counts = true;
    return one_off_matrix_v3(biom_filename,tree_filename,unifrac_method,variance_adjust,alpha,bypass_tips,normalize_sample_counts,n_substeps,subsample_depth,subsample_with_replacement,mmap_dir,result);
}

ComputeStatus one_off_matrix_fp32_v2(const char* biom_filename, const char* tree_filename,
                                      const char* unifrac_method, bool variance_adjust, double alpha,
                                      bool bypass_tips, unsigned int n_substeps,
                                      unsigned int subsample_depth, bool subsample_with_replacement, const char *mmap_dir,
                                      mat_full_fp32_t** result) {
    bool normalize_sample_counts = true;
    return one_off_matrix_fp32_v3(biom_filename,tree_filename,unifrac_method,variance_adjust,alpha,bypass_tips,normalize_sample_counts,n_substeps,subsample_depth,subsample_with_replacement,mmap_dir,result);
}

ComputeStatus one_off_matrix_v2t(const char* biom_filename, const opaque_bptree_t* tree_data,
                                 const char* unifrac_method, bool variance_adjust, double alpha,
                                 bool bypass_tips, unsigned int n_substeps,
                                 unsigned int subsample_depth, bool subsample_with_replacement, const char *mmap_dir,
                                 mat_full_fp64_t** result) {
    bool normalize_sample_counts = true;
    return one_off_matrix_v3t(biom_filename,tree_data,unifrac_method,variance_adjust,alpha,bypass_tips,normalize_sample_counts,n_substeps,subsample_depth,subsample_with_replacement,mmap_dir,result);
}

ComputeStatus one_off_matrix_fp32_v2t(const char* biom_filename, const opaque_bptree_t* tree_data,
                                      const char* unifrac_method, bool variance_adjust, double alpha,
                                      bool bypass_tips, unsigned int n_substeps,
                                      unsigned int subsample_depth, bool subsample_with_replacement, const char *mmap_dir,
                                      mat_full_fp32_t** result) {
    bool normalize_sample_counts = true;
    return one_off_matrix_fp32_v3t(biom_filename,tree_data,unifrac_method,variance_adjust,alpha,bypass_tips,normalize_sample_counts,n_substeps,subsample_depth,subsample_with_replacement,mmap_dir,result);
}

// Old interface
ComputeStatus one_off_matrix(const char* biom_filename, const char* tree_filename,
                              const char* unifrac_method, bool variance_adjust, double alpha,
                              bool bypass_tips, unsigned int nthreads,
                              const char *mmap_dir,
                              mat_full_fp64_t** result) {
    return one_off_matrix_v2(biom_filename,tree_filename,unifrac_method,variance_adjust,alpha,bypass_tips,nthreads,0,true,mmap_dir,result);
}

// Old interface
ComputeStatus one_off_matrix_fp32(const char* biom_filename, const char* tree_filename,
                                   const char* unifrac_method, bool variance_adjust, double alpha,
                                   bool bypass_tips, unsigned int nthreads,
                                   const char *mmap_dir,
                                   mat_full_fp32_t** result) {
    return one_off_matrix_fp32_v2(biom_filename,tree_filename,unifrac_method,variance_adjust,alpha,bypass_tips,nthreads,0,true,mmap_dir,result);
}
#endif // UNIFRAC_WASM (file-based v2 compat wrappers)

ComputeStatus one_off_matrix_inmem_v3(const support_biom_t *table_data, const support_bptree_t *tree_data,
                                       const char* unifrac_method, bool variance_adjust, double alpha,
                                       bool bypass_tips, bool normalize_sample_counts, unsigned int n_substeps,
                                       unsigned int subsample_depth, bool subsample_with_replacement, const char *mmap_dir,
                                       mat_full_fp64_t** result) {
    return one_off_matrix_inmem_v4(table_data,tree_data,unifrac_method,variance_adjust,alpha,bypass_tips,normalize_sample_counts,n_substeps,subsample_depth,subsample_with_replacement,/*seed*/ -1,mmap_dir,result);
}

ComputeStatus one_off_matrix_inmem_v2(const support_biom_t *table_data, const support_bptree_t *tree_data,
                                       const char* unifrac_method, bool variance_adjust, double alpha,
                                       bool bypass_tips, unsigned int n_substeps,
                                       unsigned int subsample_depth, bool subsample_with_replacement, const char *mmap_dir,
                                       mat_full_fp64_t** result) {
    bool normalize_sample_counts = true;
    return one_off_matrix_inmem_v3(table_data,tree_data,unifrac_method,variance_adjust,alpha,bypass_tips,normalize_sample_counts,n_substeps,subsample_depth,subsample_with_replacement,mmap_dir,result);
}

// Old interface
ComputeStatus one_off_inmem(const support_biom_t *table_data, const support_bptree_t *tree_data,
                             const char* unifrac_method, bool variance_adjust, double alpha,
                             bool bypass_tips, unsigned int nthreads, mat_full_fp64_t** result) {
    return one_off_matrix_inmem_v2(table_data, tree_data, unifrac_method, variance_adjust, alpha, bypass_tips, nthreads,
                                   0, true,  NULL,
                                   result);
}

ComputeStatus one_off_matrix_inmem_fp32_v3(const support_biom_t *table_data, const support_bptree_t *tree_data,
                                            const char* unifrac_method, bool variance_adjust, double alpha,
                                            bool bypass_tips, bool normalize_sample_counts, unsigned int n_substeps,
                                            unsigned int subsample_depth, bool subsample_with_replacement, const char *mmap_dir,
                                            mat_full_fp32_t** result) {
    return one_off_matrix_inmem_fp32_v4(table_data,tree_data,unifrac_method,variance_adjust,alpha,bypass_tips,normalize_sample_counts,n_substeps,subsample_depth,subsample_with_replacement,/*seed*/ -1,mmap_dir,result);
}

ComputeStatus one_off_matrix_inmem_fp32_v2(const support_biom_t *table_data, const support_bptree_t *tree_data,
                                            const char* unifrac_method, bool variance_adjust, double alpha,
                                            bool bypass_tips, unsigned int n_substeps,
                                            unsigned int subsample_depth, bool subsample_with_replacement, const char *mmap_dir,
                                            mat_full_fp32_t** result) {
    bool normalize_sample_counts = true;
    return one_off_matrix_inmem_fp32_v3(table_data,tree_data,unifrac_method,variance_adjust,alpha,bypass_tips,normalize_sample_counts,n_substeps,subsample_depth,subsample_with_replacement,mmap_dir,result);
}

// Old interface
ComputeStatus one_off_inmem_fp32(const support_biom_t *table_data, const support_bptree_t *tree_data,
                                  const char* unifrac_method, bool variance_adjust, double alpha,
                                  bool bypass_tips, unsigned int nthreads, mat_full_fp32_t** result) {
    return one_off_matrix_inmem_fp32_v2(table_data, tree_data, unifrac_method, variance_adjust, alpha, bypass_tips, nthreads,
                                        0, true,  NULL,
                                        result);
}
ComputeStatus one_dense_pair_v2t(unsigned int n_obs, const char ** obs_ids, const double* sample1, const double* sample2,
		                  const opaque_bptree_t* tree_data,
                                  const char* unifrac_method, bool variance_adjust, double alpha,
                                  bool bypass_tips, double* result) {
    bool normalize_sample_counts = true;
    return one_dense_pair_v3t(n_obs,obs_ids,sample1,sample2,tree_data,unifrac_method,variance_adjust,alpha,bypass_tips,normalize_sample_counts,result);
}

ComputeStatus one_dense_pair_v2(unsigned int n_obs, const char ** obs_ids, const double* sample1, const double* sample2,
		                 const support_bptree_t* tree_data,
                                 const char* unifrac_method, bool variance_adjust, double alpha,
                                 bool bypass_tips, double* result) {
    bool normalize_sample_counts = true;
    return one_dense_pair_v3(n_obs,obs_ids,sample1,sample2,tree_data,unifrac_method,variance_adjust,alpha,bypass_tips,normalize_sample_counts,result);
}

#ifndef UNIFRAC_WASM
ComputeStatus unifrac_to_file_v2(const char* biom_filename, const char* tree_filename, const char* out_filename,
                                  const char* unifrac_method, bool variance_adjust, double alpha,
                                  bool bypass_tips, unsigned int n_substeps, const char* format,
                                  unsigned int subsample_depth, bool subsample_with_replacement,
                                  unsigned int pcoa_dims,
                                  unsigned int permanova_perms, const char *grouping_filename, const char *grouping_columns,
                                  const char *mmap_dir) {
  bool normalize_sample_counts = true;
  return unifrac_to_file_v3(biom_filename,tree_filename,out_filename,unifrac_method,variance_adjust,alpha,
		            bypass_tips,normalize_sample_counts,n_substeps,format,subsample_depth,subsample_with_replacement,
			    pcoa_dims,pcoa_dims,grouping_filename,grouping_columns,mmap_dir);
}

ComputeStatus unifrac_to_file(const char* biom_filename, const char* tree_filename, const char* out_filename,
                               const char* unifrac_method, bool variance_adjust, double alpha,
                               bool bypass_tips, unsigned int threads, const char* format,
                               unsigned int pcoa_dims, const char *mmap_dir) {
  return unifrac_to_file_v2(biom_filename,tree_filename,out_filename,unifrac_method,variance_adjust,alpha,bypass_tips,
                            threads,format,0,true,pcoa_dims,0,NULL,NULL,mmap_dir);
}

ComputeStatus unifrac_multi_to_file_v2(const char* biom_filename, const char* tree_filename, const char* out_filename,
                                        const char* unifrac_method, bool variance_adjust, double alpha,
                                        bool bypass_tips, unsigned int nsubsteps, const char* format,
                                        unsigned int n_subsamples, unsigned int subsample_depth, bool subsample_with_replacement,
                                        unsigned int pcoa_dims,
                                        unsigned int permanova_perms, const char *grouping_filename, const char *grouping_columns,
                                        const char *mmap_dir) {
    bool normalize_sample_counts = true;
    return unifrac_multi_to_file_v3(biom_filename,tree_filename,out_filename,unifrac_method,variance_adjust,alpha,bypass_tips,normalize_sample_counts,nsubsteps,format,
		                    n_subsamples,subsample_depth,subsample_with_replacement,pcoa_dims,permanova_perms,grouping_filename,grouping_columns,mmap_dir);
}

IOStatus write_mat_from_matrix_hdf5_fp64(const char* output_filename, mat_full_fp64_t* result, unsigned int pcoa_dims, int save_dist) {
  return write_mat_from_matrix_hdf5_fp64_v2(output_filename,result,pcoa_dims,save_dist,
                             0,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
}

IOStatus write_mat_from_matrix_hdf5_fp32(const char* output_filename, mat_full_fp32_t* result, unsigned int pcoa_dims, int save_dist) {
  return write_mat_from_matrix_hdf5_fp32_v2(output_filename,result,pcoa_dims,save_dist,
                             0,NULL,NULL,NULL,NULL,NULL,NULL,NULL);
}
#endif // UNIFRAC_WASM (unifrac_to_file_v2 + write_mat_from_matrix_hdf5 compat)

