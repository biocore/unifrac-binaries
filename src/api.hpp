#include "task_parameters.hpp"
#include "status_enum.hpp"

#ifdef __cplusplus
#include <vector>
#define EXTERN extern "C"


#else
#include <stdbool.h>
#define EXTERN
#endif

/*
 *
 * Note: Each function declared EXTERN must both have 
 *       an implementation in api.cpp, AND
 *       a wrapper in ../combined/libssu.c
 *
 */

#define PARTIAL_MAGIC "SSU-PARTIAL-01"
#define PARTIAL_MAGIC_V2 0x088ABA02

/*
 * Set random seed used by this library.
 */
EXTERN void ssu_set_random_seed(unsigned int new_seed);

/* a result matrix
 *
 * n_samples <uint> the number of samples.
 * cf_size <uint> the size of the condensed form.
 * is_upper_triangle <bool> if true, indicates condensed_form represents a square
 *      matrix, and only the upper triangle is contained. if false,
 *      condensed_form represents the lower triangle of a matrix.
 * condensed_form <double*> the matrix values of length cf_size.
 * sample_ids <char**> the sample IDs of length n_samples.
 */
typedef struct mat {
    unsigned int n_samples;
    unsigned int cf_size;
    bool is_upper_triangle;
    double* condensed_form;
    char** sample_ids;
} mat_t;

/* a result matrix, full, fp64
 *
 * n_samples <uint> the number of samples.
 * matrix <double*> the matrix values, n_sample**2 size
 * sample_ids <char**> the sample IDs of length n_samples.
 */
typedef struct mat_full_fp64 {
    uint32_t n_samples;
    uint32_t flags; //opaque, 0 for default behavior
    double* matrix;
    char** sample_ids;
} mat_full_fp64_t;

/* a result matrix, full, fp32
 *
 * n_samples <uint> the number of samples.
 * matrix <float*> the matrix values, n_sample**2 size
 * sample_ids <char**> the sample IDs of length n_samples.
 */
typedef struct mat_full_fp32 {
    uint32_t n_samples;
    uint32_t flags; //opaque, 0 for default behavior
    float* matrix;
    char** sample_ids;
} mat_full_fp32_t;



/* a result vector
 *
 * n_samples <uint> the number of samples.
 * values <double*> the score values of length n_samples.
 * sample_ids <char**> the sample IDs of length n_samples.
 */
typedef struct results_vec{
    unsigned int n_samples;
    double* values;
    char** sample_ids;
} r_vec;

/* a partial result containing stripe data
 *
 * n_samples <uint> the number of samples.
 * sample_ids <char**> the sample IDs of length n_samples.
 * stripes <double**> the stripe data of dimension (stripe_stop - stripe_start, n_samples)
 * stripe_start <uint> the logical starting stripe in the final matrix.
 * stripe_stop <uint> the logical stopping stripe in the final matrix.
 * stripe_total <uint> the total number of stripes present in the final matrix.
 * is_upper_triangle <bool> whether the stripes correspond to the upper triangle of the resulting matrix.
 *      This is useful for asymmetric unifrac metrics.
 */
typedef struct partial_mat {
    uint32_t n_samples;
    char** sample_ids;
    double** stripes;
    uint32_t stripe_start;
    uint32_t stripe_stop;
    uint32_t stripe_total;
    bool is_upper_triangle;
} partial_mat_t;

/* a partial resuly, can be populated dynamically
 *
 * n_samples <uint> the number of samples.
 * sample_ids <char**> the sample IDs of length n_samples.
 * offsets <uint64_t*> offsets to the stripes in the file; 0 means unknown
 * stripes <double**> the stripe data of dimension (stripe_stop - stripe_start, n_samples)
 * stripe_start <uint> the logical starting stripe in the final matrix.
 * stripe_stop <uint> the logical stopping stripe in the final matrix.
 * stripe_total <uint> the total number of stripes present in the final matrix.
 * is_upper_triangle <bool> whether the stripes correspond to the upper triangle of the resulting matrix.
 *      This is useful for asymmetric unifrac metrics.
 * filename <char*> Name of the file from which to read
 */
typedef struct partial_dyn_mat {
    uint32_t n_samples;
    char** sample_ids;
    uint64_t* offsets;
    double** stripes;
    uint32_t stripe_start;
    uint32_t stripe_stop;
    uint32_t stripe_total;
    bool is_upper_triangle;
    char* filename;
} partial_dyn_mat_t;

/* support structure to carry in biom table information 
 *
 * obs_ids <char**> the observation IDs
 * sample_ids <char**> the sample IDs
 * indices <int32_t*> the indices of the data values
 * indptr <int32_t*> the row offset of the data values
 * data <double*> the actual matrix values
 * n_obs <int> the number of observations, corresponding to length of obs_ids
 * n_samples <int> the number of samples, corresponding to the length of sample_ids
 * nnz <int> the number of nonzero values, corresponding to the length of data and indices (DEPRECATED)
 */
typedef struct support_biom {
    char** obs_ids;
    char** sample_ids;
    uint32_t* indices;
    uint32_t* indptr;
    double* data;
    int n_obs;
    int n_samples;
    int nnz;
} support_biom_t;

/* support structure to carry in bptree information
 *
 * structure <bool*> the topology of the tree
 * lengths <double*> the branch lengths of the tree
 * names <char**> the names of the tips and internal nodes of hte tree
 * n_parens <int> the length of the structure array
 */
typedef struct support_bptree {
    bool* structure;
    double* lengths;
    char** names;
    int n_parens;
} support_bptree_t;

/* Opaque BTTree structure for externalizing the full BPTree object
 * Do not assume anything about the internals of the pointer
 */
typedef struct opaque_bptree {
    int dummy;
} opaque_bptree_t;

/* Read tree from file and fill tree_data
 *
 * tree_filename <const char*> the filename to the correspodning tree.
 * tree_data <opaque_bptree_t**> the resulting tree data object, this is initialized within the method so using **
 *
 * read_bptree_opaque returns the following error codes:
 *
 * read_okay      : no problems encountered
 * open_error     : the filename for the tree does not exist
 * unexpected_end : any other error.
 */
EXTERN IOStatus read_bptree_opaque(const char* tree_filename, opaque_bptree_t** tree_data);

/* Load tree from newick string and fill tree_data */
EXTERN void load_bptree_opaque(const char* newick, opaque_bptree_t** tree_data);

/* Convert nexplicit tree structure into the opaque form */
EXTERN void convert_bptree_opaque(const support_bptree_t* in_tree, opaque_bptree_t** tree_data);

EXTERN void destroy_mat(mat_t** result);
EXTERN void destroy_mat_full_fp64(mat_full_fp64_t** result);
EXTERN void destroy_mat_full_fp32(mat_full_fp32_t** result);
EXTERN void destroy_partial_mat(partial_mat_t** result);
EXTERN void destroy_partial_dyn_mat(partial_dyn_mat_t** result);
EXTERN void destroy_results_vec(r_vec** result);

EXTERN void destroy_bptree_opaque(opaque_bptree_t** tree_data);

/* Return number of elements in BPTree, equvalent to n_parens */
EXTERN int get_bptree_opaque_els(opaque_bptree_t* tree_data);

/* Opaque subsampled-table handle for externalizing su::biom_subsampled.
 * Do not assume anything about the internals of the pointer.
 */
typedef struct opaque_biom_inmem {
    int dummy;
} opaque_biom_inmem_t;

/* Subsample a feature table in memory.
 *
 * Wraps su::skbio_biom_subsampled. Determinism is governed by the global
 * RNG: call ssu_set_random_seed(seed) before invoking this function to
 * obtain reproducible subsampling.
 *
 * Samples whose total count is less than `depth` are dropped from the
 * output. OTUs that end up with zero counts in every surviving sample
 * are likewise dropped.
 *
 * table_data       <support_biom_t*> CSR-encoded input table.
 * depth            <unsigned int>    target per-sample read depth.
 * with_replacement <bool>            true = multinomial; false = permute.
 * out              <opaque_biom_inmem_t**> opaque handle to the
 *                  subsampled table; query via subsampled_* accessors;
 *                  release via destroy_subsampled_inmem.
 *
 * Returns table_empty if the input has no samples or no observations,
 * okay otherwise.
 */
EXTERN ComputeStatus subsample_table_inmem(const support_biom_t *table_data,
                                           unsigned int depth,
                                           bool with_replacement,
                                           opaque_biom_inmem_t **out);

EXTERN unsigned int subsampled_n_samples(const opaque_biom_inmem_t *t);
EXTERN unsigned int subsampled_n_obs(const opaque_biom_inmem_t *t);

/* Fetch a dense per-OTU vector of counts. Returns false if obs_id is not
 * present (out is left untouched in that case). */
EXTERN bool subsampled_get_obs_data(const opaque_biom_inmem_t *t,
                                    const char *obs_id,
                                    double *out);

/* Returns NULL if idx is out of range. The returned pointer remains
 * valid until destroy_subsampled_inmem is called. */
EXTERN const char* subsampled_get_sample_id(const opaque_biom_inmem_t *t, unsigned int idx);
EXTERN const char* subsampled_get_obs_id(const opaque_biom_inmem_t *t, unsigned int idx);

EXTERN void destroy_subsampled_inmem(opaque_biom_inmem_t **t);

/* Compute UniFrac - condensed form
 *
 * biom_filename <const char*> the filename to the biom table.
 * tree_filename <const char*> the filename to the correspodning tree.
 * unifrac_method <const char*> the requested unifrac method.
 * variance_adjust <bool> whether to apply variance adjustment.
 * alpha <double> GUniFrac alpha, only relevant if method == generalized.
 * bypass_tips <bool> disregard tips, reduces compute by about 50%
 * normalize_sample_counts <bool> normalize sample counts, use false for absolute quants mode
 * n_substeps <uint> the number of substeps to use.
 * result <mat_t**> the resulting distance matrix in condensed form, this is initialized within the method so using **
 *
 * one_off returns the following error codes:
 *
 * okay           : no problems encountered
 * table_missing  : the filename for the table does not exist
 * tree_missing   : the filename for the tree does not exist
 * unknown_method : the requested method is unknown.
 * table_empty    : the table does not have any entries
 */
EXTERN ComputeStatus one_off_v3(const char* biom_filename, const char* tree_filename,
                                const char* unifrac_method, bool variance_adjust, double alpha,
                                bool bypass_tips, bool normalize_sample_counts, unsigned int n_substeps, mat_t** result);

/* As above, but from a pre-loaded tree object */
EXTERN ComputeStatus one_off_wtree_v3(const char* biom_filename, const opaque_bptree_t* tree_data,
                                      const char* unifrac_method, bool variance_adjust, double alpha,
                                      bool bypass_tips, bool normalize_sample_counts, unsigned int n_substeps, mat_t** result);

/* Older version, will be deprecated in the future */
EXTERN ComputeStatus one_off(const char* biom_filename, const char* tree_filename,
                             const char* unifrac_method, bool variance_adjust, double alpha,
                             bool bypass_tips, unsigned int n_substeps, mat_t** result);

EXTERN ComputeStatus one_off_wtree(const char* biom_filename, const opaque_bptree_t* tree_data,
                                   const char* unifrac_method, bool variance_adjust, double alpha,
                                   bool bypass_tips, unsigned int n_substeps, mat_t** result);


/* Compute UniFrac - against in-memory objects returning full form matrix
 *
 * table <biom> a constructed BIOM object
 * tree <BPTree> a constructed BPTree object
 * unifrac_method <const char*> the requested unifrac method.
 * variance_adjust <bool> whether to apply variance adjustment.
 * alpha <double> GUniFrac alpha, only relevant if method == generalized.
 * bypass_tips <bool> disregard tips, reduces compute by about 50%
 * normalize_sample_counts <bool> normalize sample counts, use false for absolute quants mode
 * n_substeps <uint> the number of substeps to use.
 * subsample_depth <uint> Depth of subsampling, if >0
 * subsample_with_replacement <bool> Use subsampling with replacement? (only True supported)
 * mmap_dir <const char*> If not NULL, area to use for temp memory storage
 * result <mat_full_fp64_t**> the resulting distance matrix in full form, this is initialized within the method so using **
 *
 * one_off_inmem returns the following error codes:
 *
 * okay           : no problems encountered
 * unknown_method : the requested method is unknown.
 * table_empty    : the table does not have any entries
 */
EXTERN ComputeStatus one_off_matrix_inmem_v3(const support_biom_t *table_data, const support_bptree_t *tree_data,
                                             const char* unifrac_method, bool variance_adjust, double alpha,
                                             bool bypass_tips, bool normalize_sample_counts, unsigned int n_substeps,
                                             unsigned int subsample_depth, bool subsample_with_replacement, const char *mmap_dir,
                                             mat_full_fp64_t** result);

/* Older version, will be deprecated in the future */
EXTERN ComputeStatus one_off_matrix_inmem_v2(const support_biom_t *table_data, const support_bptree_t *tree_data,
                                             const char* unifrac_method, bool variance_adjust, double alpha,
                                             bool bypass_tips, unsigned int n_substeps,
                                             unsigned int subsample_depth, bool subsample_with_replacement, const char *mmap_dir,
                                             mat_full_fp64_t** result);

/* Older version, will be deprecated in the future */
EXTERN ComputeStatus one_off_inmem(const support_biom_t *table_data, const support_bptree_t *tree_data,
                                   const char* unifrac_method, bool variance_adjust, double alpha,
                                   bool bypass_tips, unsigned int n_substeps, mat_full_fp64_t** result);

/* Compute UniFrac - against in-memory objects returning full form matrix, fp32
 *
 * table <biom> a constructed BIOM object
 * tree <BPTree> a constructed BPTree object
 * unifrac_method <const char*> the requested unifrac method.
 * variance_adjust <bool> whether to apply variance adjustment.
 * alpha <double> GUniFrac alpha, only relevant if method == generalized.
 * bypass_tips <bool> disregard tips, reduces compute by about 50%
 * normalize_sample_counts <bool> normalize sample counts, use false for absolute quants mode
 * n_substeps <uint> the number of substeps to use.
 * subsample_depth <uint> Depth of subsampling, if >0
 * subsample_with_replacement <bool> Use subsampling with replacement? (only True supported)
 * mmap_dir <const char*> If not NULL, area to use for temp memory storage
 * result <mat_full_fp32_t**> the resulting distance matrix in full form, this is initialized within the method so using **
 *
 * one_off_inmem returns the following error codes:
 *
 * okay           : no problems encountered
 * unknown_method : the requested method is unknown.
 * table_empty    : the table does not have any entries
 */
EXTERN ComputeStatus one_off_matrix_inmem_fp32_v3(const support_biom_t *table_data, const support_bptree_t *tree_data,
                                                  const char* unifrac_method, bool variance_adjust, double alpha,
                                                  bool bypass_tips, bool normalize_sample_counts, unsigned int n_substeps,
                                                  unsigned int subsample_depth, bool subsample_with_replacement, const char *mmap_dir,
                                                  mat_full_fp32_t** result);

/* Older version, will be deprecated in the future */
EXTERN ComputeStatus one_off_matrix_inmem_fp32_v2(const support_biom_t *table_data, const support_bptree_t *tree_data,
                                                  const char* unifrac_method, bool variance_adjust, double alpha,
                                                  bool bypass_tips, unsigned int n_substeps,
                                                  unsigned int subsample_depth, bool subsample_with_replacement, const char *mmap_dir,
                                                  mat_full_fp32_t** result);

/* Older version, will be deprecated in the future */
EXTERN ComputeStatus one_off_inmem_fp32(const support_biom_t *table_data, const support_bptree_t *tree_data,
                                        const char* unifrac_method, bool variance_adjust, double alpha,
                                        bool bypass_tips, unsigned int n_substeps, mat_full_fp32_t** result);

/* Compute UniFrac - matrix form
 *
 * biom_filename <const char*> the filename to the biom table.
 * tree_filename <const char*> the filename to the correspodning tree.
 * unifrac_method <const char*> the requested unifrac method.
 * variance_adjust <bool> whether to apply variance adjustment.
 * alpha <double> GUniFrac alpha, only relevant if method == generalized.
 * bypass_tips <bool> disregard tips, reduces compute by about 50%
 * normalize_sample_counts <bool> normalize sample counts, use false for absolute quants mode
 * n_substeps <uint> the number of substeps/blocks to use.
 * subsample_depth <uint> Depth of subsampling, if >0
 * subsample_with_replacement <bool> Use subsampling with replacement? (only True supported)
 * mmap_dir <const char*> If not NULL, area to use for temp memory storage
 * result <mat_full_fp64_t**> the resulting distance matrix in matrix form, this is initialized within the method so using **
 *
 * one_off_matrix returns the following error codes:
 *
 * okay           : no problems encountered
 * table_missing  : the filename for the table does not exist
 * tree_missing   : the filename for the tree does not exist
 * unknown_method : the requested method is unknown.
 * table_empty    : the table does not have any entries
 */
EXTERN ComputeStatus one_off_matrix_v3(const char* biom_filename, const char* tree_filename,
                                       const char* unifrac_method, bool variance_adjust, double alpha,
                                       bool bypass_tips, bool normalize_sample_counts, unsigned int n_substeps,
                                       unsigned int subsample_depth, bool subsample_with_replacement, const char *mmap_dir,
                                       mat_full_fp64_t** result);

/* As above, but from a pre-loaded tree object */
EXTERN ComputeStatus one_off_matrix_v3t(const char* biom_filename, const opaque_bptree_t* tree_data,
                                        const char* unifrac_method, bool variance_adjust, double alpha,
                                        bool bypass_tips, bool normalize_sample_counts, unsigned int n_substeps,
                                        unsigned int subsample_depth, bool subsample_with_replacement, const char *mmap_dir,
                                        mat_full_fp64_t** result);

/* Older version, will be deprecated in the future */
EXTERN ComputeStatus one_off_matrix_v2(const char* biom_filename, const char* tree_filename,
                                       const char* unifrac_method, bool variance_adjust, double alpha,
                                       bool bypass_tips, unsigned int n_substeps,
                                       unsigned int subsample_depth, bool subsample_with_replacement, const char *mmap_dir,
                                       mat_full_fp64_t** result);

/* As above, but from a pre-loaded tree object */
EXTERN ComputeStatus one_off_matrix_v2t(const char* biom_filename, const opaque_bptree_t* tree_data,
                                        const char* unifrac_method, bool variance_adjust, double alpha,
                                        bool bypass_tips, unsigned int n_substeps,
                                        unsigned int subsample_depth, bool subsample_with_replacement, const char *mmap_dir,
                                        mat_full_fp64_t** result);

EXTERN ComputeStatus one_off_matrix(const char* biom_filename, const char* tree_filename,
                                    const char* unifrac_method, bool variance_adjust, double alpha,
                                    bool bypass_tips, unsigned int n_substeps,
                                    const char *mmap_dir,
                                    mat_full_fp64_t** result);

/* Compute UniFrac - matrix form, fp32 variant
 *
 * biom_filename <const char*> the filename to the biom table.
 * tree_filename <const char*> the filename to the correspodning tree.
 * unifrac_method <const char*> the requested unifrac method.
 * variance_adjust <bool> whether to apply variance adjustment.
 * alpha <double> GUniFrac alpha, only relevant if method == generalized.
 * bypass_tips <bool> disregard tips, reduces compute by about 50%
 * normalize_sample_counts <bool> normalize sample counts, use false for absolute quants mode
 * n_substeps <uint> the number of substeps/blocks to use.
 * subsample_depth <uint> Depth of subsampling, if >0
 * subsample_with_replacement <bool> Use subsampling with replacement? (only True supported)
 * mmap_dir <const char*> If not NULL, area to use for temp memory storage
 * result <mat_full_fp32_t**> the resulting distance matrix in matrix form, this is initialized within the method so using **
 *
 * one_off_matrix_fp32 returns the following error codes:
 *
 * okay           : no problems encountered
 * table_missing  : the filename for the table does not exist
 * tree_missing   : the filename for the tree does not exist
 * unknown_method : the requested method is unknown.
 * table_empty    : the table does not have any entries
 */
EXTERN ComputeStatus one_off_matrix_fp32_v3(const char* biom_filename, const char* tree_filename,
                                            const char* unifrac_method, bool variance_adjust, double alpha,
                                            bool bypass_tips, bool normalize_sample_counts, unsigned int n_substeps,
                                            unsigned int subsample_depth, bool subsample_with_replacement, const char *mmap_dir,
                                            mat_full_fp32_t** result);

/* As above, but from a pre-loaded tree object */
EXTERN ComputeStatus one_off_matrix_fp32_v3t(const char* biom_filename, const opaque_bptree_t* tree_data,
                                             const char* unifrac_method, bool variance_adjust, double alpha,
                                             bool bypass_tips, bool normalize_sample_counts, unsigned int n_substeps,
                                             unsigned int subsample_depth, bool subsample_with_replacement, const char *mmap_dir,
                                             mat_full_fp32_t** result);

/* Older version, will be deprecated in the future */
EXTERN ComputeStatus one_off_matrix_fp32_v2(const char* biom_filename, const char* tree_filename,
                                            const char* unifrac_method, bool variance_adjust, double alpha,
                                            bool bypass_tips, unsigned int n_substeps,
                                            unsigned int subsample_depth, bool subsample_with_replacement, const char *mmap_dir,
                                            mat_full_fp32_t** result);

/* As above, but from a pre-loaded tree object */
EXTERN ComputeStatus one_off_matrix_fp32_v2t(const char* biom_filename, const opaque_bptree_t* tree_data,
                                             const char* unifrac_method, bool variance_adjust, double alpha,
                                             bool bypass_tips, unsigned int n_substeps,
                                             unsigned int subsample_depth, bool subsample_with_replacement, const char *mmap_dir,
                                             mat_full_fp32_t** result);

/* Older version, will be deprecated in the future */
EXTERN ComputeStatus one_off_matrix_fp32(const char* biom_filename, const char* tree_filename,
                                         const char* unifrac_method, bool variance_adjust, double alpha,
                                         bool bypass_tips, unsigned int n_substeps,
                                         const char *mmap_dir,
                                         mat_full_fp32_t** result);

/* Compute UniFrac from a pair of dense vectors 
 *
 * n_obs <unsigned int> the number of observations, corresponding to length of obs_ids, sample1 and sample2
 * obs_ids <const char**> the observation IDs
 * sample1 <const double*> sample values
 * sample2 <const double*> sample values
 * unifrac_method <const char*> the requested unifrac method.
 * variance_adjust <bool> whether to apply variance adjustment.
 * alpha <double> GUniFrac alpha, only relevant if method == generalized.
 * bypass_tips <bool> disregard tips, reduces compute by about 50%
 * normalize_sample_counts <bool> normalize sample counts, use false for absolute quants mode
 * result <double*> the resulting distance
 *
 * one_dense_pair_v3t returns the following error codes:
 *
 * okay           : no problems encountered
 * table_empty    : the table does not have any entries
 */
EXTERN ComputeStatus one_dense_pair_v3t(unsigned int n_obs, const char ** obs_ids, const double* sample1, const double* sample2,
		                        const opaque_bptree_t* tree_data,
                                        const char* unifrac_method, bool variance_adjust, double alpha,
                                        bool bypass_tips, bool normalize_sample_counts, double* result);

/* Same as above, but usign the explicit treee structure... slightly less efficient */
EXTERN ComputeStatus one_dense_pair_v3(unsigned int n_obs, const char ** obs_ids, const double* sample1, const double* sample2,
		                       const support_bptree_t* tree_data,
                                       const char* unifrac_method, bool variance_adjust, double alpha,
                                       bool bypass_tips, bool normalize_sample_counts, double* result);

/* Old variants, deprecated */

EXTERN ComputeStatus one_dense_pair_v2t(unsigned int n_obs, const char ** obs_ids, const double* sample1, const double* sample2,
		                        const opaque_bptree_t* tree_data,
                                        const char* unifrac_method, bool variance_adjust, double alpha,
                                        bool bypass_tips, double* result);

/* Same as above, but usign the explicit treee structure... slightly less efficient */
EXTERN ComputeStatus one_dense_pair_v2(unsigned int n_obs, const char ** obs_ids, const double* sample1, const double* sample2,
		                       const support_bptree_t* tree_data,
                                       const char* unifrac_method, bool variance_adjust, double alpha,
                                       bool bypass_tips, double* result);

/* compute Faith PD
 * biom_filename <const char*> the filename to the biom table.
 * tree_filename <const char*> the filename to the correspodning tree.
 * result <r_vec**> the resulting vector of computed Faith PD values
 *
 * faith_pd_one_off returns the following error codes:
 *
 * okay           : no problems encountered
 * table_missing  : the filename for the table does not exist
 * tree_missing   : the filename for the tree does not exist
 * table_empty    : the table does not have any entries
 */
EXTERN ComputeStatus faith_pd_one_off(const char* biom_filename, const char* tree_filename,
                                      r_vec** result);

/* compute Faith PD from in-memory inputs.
 *
 * In-memory analogue of faith_pd_one_off. Accepts a sparse table and a
 * pre-built tree directly, with no filesystem access. Required for the
 * WASM build and useful for any embedding context where BIOM v2 (HDF5)
 * is unavailable or undesired.
 *
 * table_data <support_biom_t*> CSR-encoded feature table.
 * tree_data  <support_bptree_t*> balanced-parens tree structure.
 * result     <r_vec**> the resulting vector of computed Faith PD values.
 *
 * Returns the same error codes as faith_pd_one_off, plus tree_missing
 * if tree_data is NULL.
 */
EXTERN ComputeStatus faith_pd_inmem(const support_biom_t *table_data,
                                    const support_bptree_t *tree_data,
                                    r_vec** result);

/* Compute UniFrac and save to file
 *
 * biom_filename <const char*> the filename to the biom table.
 * tree_filename <const char*> the filename to the correspodning tree.
 * out_filename <const char*> the filename of the output file.
 * unifrac_method <const char*> the requested unifrac method.
 * variance_adjust <bool> whether to apply variance adjustment.
 * alpha <double> GUniFrac alpha, only relevant if method == generalized.
 * bypass_tips <bool> disregard tips, reduces compute by about 50%
 * normalize_sample_counts <bool> normalize sample counts, use false for absolute quants mode
 * n_substeps <uint> the number of substeps to use.
 * format <const char*> output format to use.
 * subsample_depth <uint> Depth of subsampling, if >0
 * subsample_with_replacement <bool> Use subsampling with replacement? (only True supported)
 * pcoa_dims <uint> if not 0, number of dimensions to use or PCoA
 * permanova_perms <uint> If not 0, compute PERMANOVA using that many permutations
 * grouping_filename <const char*> the TSV filename containing grouping information
 * grouping_columns <const char *> the columns to use for grouping
 * mmap_dir <const char*> if not empty, temp dir to use for disk-based memory 
 *
 * unifrac_to_file_v3 returns the following error codes:
 *
 * okay           : no problems encountered
 * table_missing  : the filename for the table does not exist
 * tree_missing   : the filename for the tree does not exist
 * unknown_method : the requested method is unknown.
 * table_empty    : the table does not have any entries
 * output_error   : failed to properly write the output file
 */
EXTERN ComputeStatus unifrac_to_file_v3(const char* biom_filename, const char* tree_filename, const char* out_filename,
                                        const char* unifrac_method, bool variance_adjust, double alpha,
                                        bool bypass_tips, bool normalize_sample_counts, unsigned int n_substeps, const char* format,
                                        unsigned int subsample_depth, bool subsample_with_replacement, 
                                        unsigned int pcoa_dims,
                                        unsigned int permanova_perms, const char *grouping_filename, const char *grouping_columns,
                                        const char *mmap_dir);

/* As above, but write using the text format */
EXTERN ComputeStatus unifrac_to_txt_file_v3(const char* biom_filename, const char* tree_filename, const char* out_filename,
                                            const char* unifrac_method, bool variance_adjust, double alpha,
                                            bool bypass_tips, bool normalize_sample_counts, unsigned int n_substeps,
                                            const char *mmap_dir);

/* Older version, will be deprecated in the future */
EXTERN ComputeStatus unifrac_to_file_v2(const char* biom_filename, const char* tree_filename, const char* out_filename,
                                        const char* unifrac_method, bool variance_adjust, double alpha,
                                        bool bypass_tips, unsigned int n_substeps, const char* format,
                                        unsigned int subsample_depth, bool subsample_with_replacement, 
                                        unsigned int pcoa_dims,
                                        unsigned int permanova_perms, const char *grouping_filename, const char *grouping_columns,
                                        const char *mmap_dir);

// backwards compatible version
EXTERN ComputeStatus unifrac_to_file(const char* biom_filename, const char* tree_filename, const char* out_filename,
                                     const char* unifrac_method, bool variance_adjust, double alpha,
                                     bool bypass_tips, unsigned int n_substeps, const char* format,
                                     unsigned int pcoa_dims, const char *mmap_dir);

/* Compute UniFrac with subsampling several times and save to file
 *
 * biom_filename <const char*> the filename to the biom table.
 * tree_filename <const char*> the filename to the corresponding tree.
 * out_filename <const char*> the filename of the output file.
 * unifrac_method <const char*> the requested unifrac method.
 * variance_adjust <bool> whether to apply variance adjustment.
 * alpha <double> GUniFrac alpha, only relevant if method == generalized.
 * bypass_tips <bool> disregard tips, reduces compute by about 50%
 * normalize_sample_counts <bool> normalize sample counts, use false for absolute quants mode
 * n_substeps <uint> the number of substeps to use.
 * format <const char*> output format to use.
 * n_subsamples <uint> Number of subsamples to compute.
 * subsample_depth <uint> Depth of subsampling.
 * subsample_with_replacement <bool> Use subsampling with replacement? (only True supported)
 * pcoa_dims <uint> if not 0, number of dimensions to use or PCoA
 * permanova_perms <uint> If not 0, compute PERMANOVA using that many permutations
 * grouping_filename <const char*> the TSV filename containing grouping information
 * grouping_columns <const char *> the columns to use for grouping
 * mmap_dir <const char*> if not empty, temp dir to use for disk-based memory 
 *
 * unifrac_multi_to_file_v3 returns the following error codes:
 *
 * okay           : no problems encountered
 * table_missing  : the filename for the table does not exist
 * tree_missing   : the filename for the tree does not exist
 * unknown_method : the requested method is unknown.
 * table_empty    : the table does not have any entries
 * output_error   : failed to properly write the output file
 */
EXTERN ComputeStatus unifrac_multi_to_file_v3(const char* biom_filename, const char* tree_filename, const char* out_filename,
                                              const char* unifrac_method, bool variance_adjust, double alpha,
                                              bool bypass_tips, bool normalize_sample_counts, unsigned int n_substeps, const char* format,
                                              unsigned int n_subsamples, unsigned int subsample_depth, bool subsample_with_replacement, 
                                              unsigned int pcoa_dims,
                                              unsigned int permanova_perms, const char *grouping_filename, const char *grouping_columns,
                                              const char *mmap_dir);

/* Older version, will be deprecated in the future */
EXTERN ComputeStatus unifrac_multi_to_file_v2(const char* biom_filename, const char* tree_filename, const char* out_filename,
                                              const char* unifrac_method, bool variance_adjust, double alpha,
                                              bool bypass_tips, unsigned int n_substeps, const char* format,
                                              unsigned int n_subsamples, unsigned int subsample_depth, bool subsample_with_replacement, 
                                              unsigned int pcoa_dims,
                                              unsigned int permanova_perms, const char *grouping_filename, const char *grouping_columns,
                                              const char *mmap_dir);

/* Compute PERMANOVA - fp64 variant
 *
 * grouping_filename <const char*> the filename to the grouping TSV file
 * n_columns <uint> the number of grouping columns to process
 * columns <const char**> the grouping columns, array of size n_columns
 * mmap_dir <const char*> If not NULL, area to use for temp memory storage
 * result <mat_full_fp64_t*> the distance matrix in matrix form
 * permanova_perms <uint> the number of permutations, should be >0
 * fstats <double *> the resulting PERMANOVA fstats, array of size n_columns
 * pvalues <double *> the resulting PERMANOVA pvalues, array of size n_columns
 *
 * compute_permanova_fp64 returns the following error codes:
 *
 * okay               : no problems encountered
 * grouping_missing  : the filename for the grouping does not exist or is not valid
 */
EXTERN ComputeStatus compute_permanova_fp64(const char *grouping_filename, unsigned int n_columns, const char* *columns,
                                            mat_full_fp64_t * result, unsigned int permanova_perms,
                                            double *fstats, double *pvalues);

/* Compute PERMANOVA - fp32 variant
 *
 * grouping_filename <const char*> the filename to the grouping TSV file
 * n_columns <uint> the number of grouping columns to process
 * columns <const char**> the grouping columns, array of size n_columns
 * mmap_dir <const char*> If not NULL, area to use for temp memory storage
 * result <mat_full_fp32_t*> the distance matrix in matrix form
 * permanova_perms <uint> the number of permutations, should be >0
 * fstats <float *> the resulting PERMANOVA fstats, array of size n_columns
 * pvalues <float *> the resulting PERMANOVA pvalues, array of size n_columns
 *
 * compute_permanova_fp32 returns the following error codes:
 *
 * okay               : no problems encountered
 * grouping_missing  : the filename for the grouping does not exist or is not valid
 */
EXTERN ComputeStatus compute_permanova_fp32(const char *grouping_filename, unsigned int n_columns, const char* * columns,
                                            mat_full_fp32_t * result, unsigned int permanova_perms,
                                            float *fstats, float *pvalues);

/* Compute PERMANOVA from in-memory inputs.
 *
 * In-memory analogue of compute_permanova_fp64 that bypasses the TSV
 * grouping parser. Required for the WASM build (no filesystem access)
 * and useful in any embedding context where labels are already in
 * memory. Determinism is governed by the global RNG: call
 * ssu_set_random_seed(seed) before invoking this function for
 * reproducible p-values; skbb's portable Fisher-Yates means results
 * are bit-exact across native and WASM builds at the same seed.
 *
 * mat              <const double*> distance matrix, n_dims x n_dims, row-major.
 * n_dims           <unsigned int>  size of the matrix.
 * grouping         <const uint32_t*> length-n_dims integer group labels;
 *                  caller-built. Group ids do not need to be contiguous.
 * permanova_perms  <unsigned int>  number of permutations, > 0.
 * fstat            <double*>       out, computed F statistic.
 * pvalue           <double*>       out, computed p-value.
 *
 * Returns okay on success; permanova_failed on bad input.
 */
EXTERN ComputeStatus compute_permanova_inmem_fp64(const double *mat, unsigned int n_dims,
                                                  const uint32_t *grouping,
                                                  unsigned int permanova_perms,
                                                  double *fstat, double *pvalue);

EXTERN ComputeStatus compute_permanova_inmem_fp32(const float *mat, unsigned int n_dims,
                                                  const uint32_t *grouping,
                                                  unsigned int permanova_perms,
                                                  float *fstat, float *pvalue);

/* Write a matrix object using the text format
 *
 * filename <const char*> the file to write into
 * result <mat_t*> the results object
 *
 * The following error codes are returned:
 *
 * write_okay : no problems
 */
EXTERN IOStatus write_mat(const char* filename, mat_t* result);

/* Write a matrix object using hdf5 format, using fp64 precision
 *
 * filename <const char*> the file to write into
 * result <mat_t*> the results object
 * pcoa_dims <uint> PCoAdimensions to compute, if >0
 * save_dist <bool> If false, do not same the distance matrix data
 *
 * The following error codes are returned:
 *
 * write_okay : no problems
 */
// backwards compatible version, deprecated
IOStatus write_mat_hdf5_fp64(const char* filename, mat_t* result, unsigned int pcoa_dims, int save_dist);

/* Write a matrix object using hdf5 format, using fp32 precision
 *
 * filename <const char*> the file to write into
 * result <mat_t*> the results object
 * pcoa_dims <uint> PCoAdimensions to compute, if >0
 * save_dist <bool> If false, do not same the distance matrix data
 *
 * The following error codes are returned:
 *
 * write_okay : no problems
 */
// backwards compatible version, deprecated
IOStatus write_mat_hdf5_fp32(const char* filename, mat_t* result, unsigned int pcoa_dims, int save_dist);

/* Write a matrix object using the text format
 *
 * filename <const char*> the file to write into
 * result <mat_full_t*> the results object
 *
 * The following error codes are returned:
 *
 * write_okay : no problems
 */
EXTERN IOStatus write_mat_from_matrix(const char* filename, mat_full_fp64_t* result);

/* as above but fp32 */
EXTERN IOStatus write_mat_from_matrix_fp32(const char* filename, mat_full_fp32_t* result);


/* Write a matrix object from buffer using hdf5 format, using fp64 precision
 *
 * filename <const char*> the file to write into
 * result <mat_full_t*> the results object
 * pcoa_dims <uint> PCoA dimensions to compute, if >0
 * save_dist <bool> If false, do not same the distance matrix data
 * stat_n_vals <uint> Number of statistics to save
 * stat_method_arr <const char**> Array of stat. methods names
 * stat_name_arr   <const char**> Array of stat. test names
 * stat_val_arr   <float*>        Array of stat. values
 * stat_pval_arr  <float*>        Array of stat. p-values
 * stat_perm_count_arr <uint*>    Array of the number of permutations used
 * stat_group_name_arr <const char**> Array group names used for the grouping
 * stat_group_count_arr <uint*>       Array number of groups used for the statistics
 *
 * The following error codes are returned:
 *
 * write_okay : no problems
 * write_error : something went wrong
 */
EXTERN IOStatus write_mat_from_matrix_hdf5_fp64_v2(const char* output_filename, mat_full_fp64_t* result,
                                                   unsigned int pcoa_dims, int save_dist,
                                                   unsigned int stat_n_vals,
                                                   const char*  *stat_method_arr,     const char*        *stat_name_arr,
                                                   const double *stat_val_arr,        const double       *stat_pval_arr, const unsigned int *stat_perm_count_arr,
                                                   const char*  *stat_group_name_arr, const unsigned int *stat_group_count_arr);

// backwards compatible version, deprecated
EXTERN IOStatus write_mat_from_matrix_hdf5_fp64(const char* filename, mat_full_fp64_t* result, unsigned int pcoa_dims, int save_dist);

/* Write a matrix object from buffer using hdf5 format, using fp32 precision
 *
 * filename <const char*> the file to write into
 * result <mat_full_fp32_t*> the results object
 * pcoa_dims <uint> PCoAdimensions to compute, if >0
 * save_dist <bool> If false, do not same the distance matrix data
 * stat_n_vals <uint> Number of statistics to save
 * stat_method_arr <const char**> Array of stat. methods names
 * stat_name_arr   <const char**> Array of stat. test names
 * stat_val_arr   <float*>        Array of stat. values
 * stat_pval_arr  <float*>        Array of stat. p-values
 * stat_perm_count_arr <uint*>    Array of the number of permutations used
 * stat_group_name_arr <const char**> Array group names used for the grouping
 * stat_group_count_arr <uint*>       Array number of groups used for the statistics
 *
 * The following error codes are returned:
 *
 * write_okay : no problems
 * write_error : something went wrong
 */
EXTERN IOStatus write_mat_from_matrix_hdf5_fp32_v2(const char* output_filename, mat_full_fp32_t* result,
                                                   unsigned int pcoa_dims, int save_dist,
                                                   unsigned int stat_n_vals,
                                                   const char*  *stat_method_arr,     const char*        *stat_name_arr,
                                                   const float  *stat_val_arr,        const float        *stat_pval_arr, const unsigned int *stat_perm_count_arr,
                                                   const char*  *stat_group_name_arr, const unsigned int *stat_group_count_arr);

// backwards compatible version, deprecated
EXTERN IOStatus write_mat_from_matrix_hdf5_fp32(const char* filename, mat_full_fp32_t* result, unsigned int pcoa_dims, int save_dist);

/* Write a series
 *
 * filename <const char*> the file to write into
 * result <r_vec*> the results object
 *
 * The following error codes are returned:
 *
 * write_okay : no problems
 */
EXTERN IOStatus write_vec(const char* filename, r_vec* result);

/* Read a matrix object
 *
 * filename <const char*> the file to write into
 * result <mat_t**> the results object
 *
 * The following error codes are returned:
 *
 * read_okay : no problems
 * open_error : could not open the file
 * magic_incompatible : format magic not found or incompatible
 * unexpected_end     : format end not found in expected location
 */
//EXTERN IOStatus read_mat(const char* filename, mat_t** result);

/* Compute a subset of a UniFrac distance matrix
 *
 * biom_filename <const char*> the filename to the biom table.
 * tree_filename <const char*> the filename to the correspodning tree.
 * unifrac_method <const char*> the requested unifrac method.
 * variance_adjust <bool> whether to apply variance adjustment.
 * alpha <double> GUniFrac alpha, only relevant if method == generalized.
 * bypass_tips <bool> disregard tips, reduces compute by about 50%
 * normalize_sample_counts <bool> normalize sample counts, use false for absolute quants mode
 * n_substeps <uint> the number of substeps to use.
 * stripe_start <uint> the starting stripe to compute
 * stripe_stop <uint> the last stripe to compute
 * dm_stripes <vector of double*> the unique branch length stripes. This is expected to be
 *      uninitialized, and is an output parameter.
 * dm_stripes_total <vector of double*> the total branch length stripes. This is expected to be
 *      uninitialized, and is an output parameter.
 * result <partial_mat_t**> the resulting distance matrix in condensed form, this is initialized within the method so using **
 *
 * partial returns the following error codes:
 *
 * okay           : no problems encountered
 * table_missing  : the filename for the table does not exist
 * tree_missing   : the filename for the tree does not exist
 * unknown_method : the requested method is unknown.
 */

EXTERN ComputeStatus partial_v3(const char* biom_filename, const char* tree_filename,
                                const char* unifrac_method, bool variance_adjust, double alpha,
                                bool bypass_tips, bool normalize_sample_count, unsigned int n_substeps,
				unsigned int stripe_start, unsigned int stripe_stop, partial_mat_t** result);

/* Older version, will be deprecated in the future */
EXTERN ComputeStatus partial(const char* biom_filename, const char* tree_filename,
                             const char* unifrac_method, bool variance_adjust, double alpha,
                             bool bypass_tips, unsigned int n_substeps, unsigned int stripe_start,
                             unsigned int stripe_stop, partial_mat_t** result);

/* Write a partial matrix object
 *
 * filename <const char*> the file to write into
 * result <partial_mat_t*> the partial results object
 *
 * The following error codes are returned:
 *
 * write_okay : no problems
 * open_error : could not open the file
 *
 * The structure of the binary output file is as follows. Newlines added for clarity, but are not stored.
 * The file has logical blocks, but are not explicitly denoted in the format. These logical blocks are
 * just used to improve readability here, and are denoted by ### marks.
 *
 * ### HEADER ###
 * <MAGIC_LEN>          : uint16_t, the length of the magic
 * <MAGIC>              : char, e.g., SSU-PARTIAL-01
 * <N_SAMPLES>          : uint32_t, the number of samples
 * <N_STRIPES>          : uint32_t, the number of stripes represented in this file
 * <STRIPE_START>       : uint32_t, the starting stripe number
 * <STRIPES_TOTAL>      : uint32_t, the total number of stripes in the full matrix
 * <IS_UPPER_TRIANGLE>  : uint8_t, zero is false, nonzero is true
 *
 * ### SAMPLE IDS ###
 * <LEN>                : uint16_t, the length of the next sample ID
 * <SAMPLE_ID[0]>       : LEN bytes, char
 * ...                  : ... repeated <LEN><SAMPLE_ID[i]>
 * <LEN>                : uint16_t, the length of the next sample ID
 * <SAMPLE_ID[N]>       : LEN bytes, char
 *
 * ### STRIPE VALUES; SS -> STRIPE_START, NS -> N_STRIPES
 * <STRIPE[SS][0]>      : double, the first value in the 0th stripe
 * ...                  : ... repeated for N_SAMPLES values
 * <STRIPE[SS][N]>      : double, the last value in the 0th stripe
 * <STRIPE[SS + NS][0]> : double, the first value in the Kth stripe
 * ...                  : ... repeated for N_SAMPLES values
 * <STRIPE[SS + NS][N]> : double, the last value in the Kth stripe
 *
 * ### FOOTER ###
 * <MAGIC>              : char, e.g., SSU-PARTIAL-01, same as starting magic
 */
EXTERN IOStatus write_partial(const char* filename, const partial_mat_t* result);

/* Read a partial matrix object
 *
 * filename <const char*> the file to write into
 * result <partial_mat_t**> the partial results object, output parameter
 *
 * The following error codes are returned:
 *
 * read_okay          : no problems
 * open_error         : could not open the file
 * magic_incompatible : format magic not found or incompatible
 * bad_header         : header seems malformed
 * unexpected_end     : format end not found in expected location
 */
EXTERN IOStatus read_partial(const char* filename, partial_mat_t** result);

/* Read a partial matrix object header
 *
 * filename <const char*> the file to write into
 * result <partial_dyn_mat_t**> the partial results object, output parameter
 *
 * The following error codes are returned:
 *
 * read_okay          : no problems
 * open_error         : could not open the file
 * magic_incompatible : format magic not found or incompatible
 * bad_header         : header seems malformed
 * unexpected_end     : format end not found in expected location
 */
EXTERN IOStatus read_partial_header(const char* input_filename, partial_dyn_mat_t** result_out);

/* Read a stripe of a partial matrix
 *
 * filename <const char*> the file to write into
 * result <partial_dyn_mat_t*> the partial results object
 * stripe_idx <uint> relative stripe number 
 *
 * The following error codes are returned:
 *
 * read_okay          : no problems
 * open_error         : could not open the file
 * magic_incompatible : format magic not found or incompatible
 * bad_header         : header seems malformed
 * unexpected_end     : format end not found in expected location
 */
EXTERN IOStatus read_partial_one_stripe(partial_dyn_mat_t* result, uint32_t stripe_idx);

/*
 * Description TBD
 */
EXTERN MergeStatus validate_partial(const partial_dyn_mat_t* const * partial_mats, int n_partials);

#if 0
/* Merge partial results
 * DEPRECATED, not used by anyone anymore
 *
 * results <partial_mat_t**> an array of partial_mat_t*, the buffers will be destroyed in the process
 * n_partials <int> number of partial mats
 * merged <mat_t**> the full matrix, output parameters, this is initialized in the method so using **
 *
 * The following error codes are returned:
 *
 * merge_okay            : no problems
 * incomplete_stripe_set : not all stripes needed to create a full matrix were foun
 * sample_id_consistency : samples described by stripes are inconsistent
 * square_mismatch       : inconsistency on denotation of square matrix
 */
MergeStatus merge_partial(partial_mat_t** partial_mats, int n_partials, unsigned int dummy, mat_t** result);
#endif

/* Merge partial results
 * DEPRECATED, for backwards compatibility only, use the mmap version instead
 *
 * partial_mats <partial_dyn_mat_t**> an array of partial_dyn_mat_t*
 * n_partials <int> number of partial mats
 * result <mat_full_fp64_t**> the full matrix, output parameters, this is initialized in the method so using **
 *
 * The following error codes are returned:
 *
 * merge_okay            : no problems
 * incomplete_stripe_set : not all stripes needed to create a full matrix were foun
 * sample_id_consistency : samples described by stripes are inconsistent
 * square_mismatch       : inconsistency on denotation of square matrix
 */
MergeStatus merge_partial_to_matrix(partial_dyn_mat_t* * partial_mats, int n_partials, mat_full_fp64_t** result);

/* Merge partial results
 * DEPRECATED, for backwards compatibility only, use the mmap version instead
 *
 * partial_mats <partial_dyn_mat_t**> an array of partial_dyn_mat_t*
 * n_partials <int> number of partial mats
 * result <mat_full_fp32_t**> the full matrix, output parameters, this is initialized in the method so using **
 *
 * The following error codes are returned:
 *             
 * merge_okay            : no problems
 * incomplete_stripe_set : not all stripes needed to create a full matrix were foun
 * sample_id_consistency : samples described by stripes are inconsistent
 * square_mismatch       : inconsistency on denotation of square matrix
 */
MergeStatus merge_partial_to_matrix_fp32(partial_dyn_mat_t* * partial_mats, int n_partials, mat_full_fp32_t** result);


/* Merge partial results
 *
 * partial_mats <partial_dyn_mat_t**> an array of partial_dyn_mat_t*
 * n_partials <int> number of partial mats
 * mmap_dir <const char *> Where to host the mmap file, use NULL if no mmap
 * result <mat_full_fp64_t**> the full matrix, output parameters, this is initialized in the method so using **
 *
 * The following error codes are returned:
 *
 * merge_okay            : no problems
 * incomplete_stripe_set : not all stripes needed to create a full matrix were foun
 * sample_id_consistency : samples described by stripes are inconsistent
 * square_mismatch       : inconsistency on denotation of square matrix
 */
EXTERN MergeStatus merge_partial_to_mmap_matrix(partial_dyn_mat_t* * partial_mats, int n_partials, const char *mmap_dir, mat_full_fp64_t** result);

/* Merge partial results
 *
 * partial_mats <partial_dyn_mat_t**> an array of partial_dyn_mat_t*
 * n_partials <int> number of partial mats
 * mmap_dir <const char *> Where to host the mmap file, use NULL if no mmap
 * result <mat_full_fp32_t**> the full matrix, output parameters, this is initialized in the method so using **
 *
 * The following error codes are returned:
 *
 * merge_okay            : no problems
 * incomplete_stripe_set : not all stripes needed to create a full matrix were foun
 * sample_id_consistency : samples described by stripes are inconsistent
 * square_mismatch       : inconsistency on denotation of square matrix
 */
EXTERN MergeStatus merge_partial_to_mmap_matrix_fp32(partial_dyn_mat_t* * partial_mats, int n_partials, const char *mmap_dir, mat_full_fp32_t** result);


// Find eigen values and vectors
// Based on N. Halko, P.G. Martinsson, Y. Shkolnisky, and M. Tygert.
//     Original Paper: https://arxiv.org/abs/1007.5510
// centered == n x n, must be symmetric, Note: will be used in-place as temp buffer
void find_eigens_fast(const uint32_t n_samples, const uint32_t n_dims, double * centered, double **eigenvalues, double **eigenvectors);
void find_eigens_fast_p32(const uint32_t n_samples, const uint32_t n_dims,float * centered, float **eigenvalues, float **eigenvectors);

// Perform Principal Coordinate Analysis
// mat       - in, result of unifrac compute
// n_samples - in, size of the matrix (n x n)
// n_dims    - in, Dimensions to reduce the distance matrix to. This number determines how many eigenvectors and eigenvalues will be returned.
// eigenvalues - out, alocated buffer of size n_dims
// samples     - out, alocated buffer of size n_dims x n_samples
// proportion_explained - out, allocated buffer of size n_dims
void pcoa(const double * mat, const uint32_t n_samples, const uint32_t n_dims, double **eigenvalues, double **samples, double **proportion_explained);
void pcoa_fp32(const float * mat, const uint32_t n_samples, const uint32_t n_dims, float * *eigenvalues, float * *samples, float * *proportion_explained);
void pcoa_mixed(const double * mat, const uint32_t n_samples, const uint32_t n_dims, float * *eigenvalues, float * *samples, float * *proportion_explained);


#ifdef __cplusplus
// TODO: only needed for testing, should be encased in a macro
void set_tasks(std::vector<su::task_parameters> &tasks,
               double alpha,
               unsigned int n_samples,
               unsigned int stripe_start,
               unsigned int stripe_stop,
               bool bypass_tips,
               bool normalize_sample_counts,
               unsigned int n_tasks);

#endif
