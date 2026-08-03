#include "api.hpp"
#ifndef UNIFRAC_WASM
#include "biom.hpp"
#endif
#include "tree.hpp"
#include "tsv.hpp"
#include "unifrac.hpp"
#include "skbio_alt.hpp"
#include <fstream>
#include <iomanip>
#include <memory>
#include <thread>
#include <cstring>
#include <stdlib.h> 
#include <string.h> 
#include <string> 
#include <vector>
#include <stdexcept>
#include <charconv>

#ifndef UNIFRAC_WASM
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <lz4.h>
#endif
#include <time.h>
#if defined(_OPENMP)
#include <omp.h>
#endif

#define MMAP_FD_MASK 0x0fff
#define MMAP_FLAG    0x1000

// Note: Threading is now full controlled by OpenMP.
// Any threads variable is really referring to n_substeps.
// The old naming was retained to minimize code refactoring.

#define SETUP_TDBG(method) const char *tdbg_method=method; \
                          bool print_tdbg = false;\
                          if (const char* env_p = std::getenv("UNIFRAC_TIMING_INFO")) { \
                            print_tdbg = true; \
                            std::string env_s(env_p); \
                            if ((env_s=="NO") || (env_s=="N") || (env_s=="no") || (env_s=="n") || \
                               (env_s=="NEVER") || (env_s=="never")) print_tdbg = false; \
                          } \
                          time_t tgdb_t0; time(&tgdb_t0); \
                          if(print_tdbg) printf("INFO (unifrac): Starting %s\n",tdbg_method);

#define TDBG_STEP(sname) if(print_tdbg) {\
                           time_t tgdb_t1; time(&tgdb_t1); \
                           printf("INFO (unifrac): dt %4i : Completed %s.%s\n",(int)(tgdb_t1-tgdb_t0),tdbg_method,sname); \
                           tgdb_t0 = tgdb_t1; \
                         }

#define CHECK_FILE(filename, err) if(!is_file_exists(filename)) { \
                                      return err;                 \
                                  }

#define SET_METHOD(requested_method, err) Method method;                                                            \
                                          if(std::strcmp(requested_method, "unweighted") == 0)                      \
                                              method = unweighted_fp32;                                             \
                                          else if(std::strcmp(requested_method, "weighted_normalized") == 0)        \
                                              method = weighted_normalized_fp32;                                    \
                                          else if(std::strcmp(requested_method, "weighted_unnormalized") == 0)      \
                                              method = weighted_unnormalized_fp32;                                  \
                                          else if(std::strcmp(requested_method, "unweighted_unnormalized") == 0)    \
                                              method = unweighted_unnormalized_fp32;                                \
                                          else if(std::strcmp(requested_method, "generalized") == 0)                \
                                              method = generalized_fp32;                                            \
                                          else if(std::strcmp(requested_method, "unweighted_fp64") == 0)            \
                                              method = unweighted;                                                  \
                                          else if(std::strcmp(requested_method, "weighted_normalized_fp64") == 0)   \
                                              method = weighted_normalized;                                         \
                                          else if(std::strcmp(requested_method, "weighted_unnormalized_fp64") == 0) \
                                              method = weighted_unnormalized;                                       \
                                          else if(std::strcmp(requested_method, "unweighted_unnormalized_fp64") == 0) \
                                              method = unweighted_unnormalized;                                     \
                                          else if(std::strcmp(requested_method, "generalized_fp64") == 0)           \
                                              method = generalized;                                                 \
                                          else if(std::strcmp(requested_method, "unweighted_fp32") == 0)            \
                                              method = unweighted_fp32;                                             \
                                          else if(std::strcmp(requested_method, "weighted_normalized_fp32") == 0)   \
                                              method = weighted_normalized_fp32;                                    \
                                          else if(std::strcmp(requested_method, "weighted_unnormalized_fp32") == 0) \
                                              method = weighted_unnormalized_fp32;                                  \
                                          else if(std::strcmp(requested_method, "unweighted_unnormalized_fp32") == 0) \
                                              method = unweighted_unnormalized_fp32;                                \
                                          else if(std::strcmp(requested_method, "generalized_fp32") == 0)           \
                                              method = generalized_fp32;                                            \
                                          else {                                                               \
                                              return err;                                                      \
                                          }

#define SYNC_TREE_TABLE(tree, table) std::unordered_set<std::string> to_keep(table.get_obs_ids().begin(),           \
                                                                             table.get_obs_ids().end());            \
                                     su::BPTree tree_sheared = tree.shear(to_keep).collapse();

#define VALIDATE_TREE_TABLE(tree, table)                if( (table.n_samples <= 0) || (table.n_obs <= 0) ) {                     \
                                                            return table_empty;                                                  \
                                                        }                                                                        \
                                                        std::string bad_id = su::test_table_ids_are_subset_of_tree(table, tree); \
                                                        if(bad_id != "") {                                                       \
                                                            return table_and_tree_do_not_overlap;                                \
                                                        }   

#define PARSE_TREE_TABLE(tree_filename, table_filename) su::BPTree tree(get_tree_content(tree_filename)); \
                                                        su::biom table(biom_filename);                    \
							VALIDATE_TREE_TABLE(tree, table)

#define PARSE_SYNC_TREE_TABLE(tree_filename, table_filename) PARSE_TREE_TABLE(tree_filename, table_filename) \
                                                             SYNC_TREE_TABLE(tree, table)

using namespace su;
using namespace std;

void ssu_set_random_seed(unsigned int new_seed) {
  su::set_random_seed(new_seed);
}

#ifndef UNIFRAC_WASM
// https://stackoverflow.com/a/19841704/19741
bool is_file_exists(const char *fileName) {
    std::ifstream infile(fileName);
        return infile.good();
}

inline std::string get_tree_content(const char* tree_filename) {
    std::ifstream ifs(tree_filename);
    std::string content = std::string(std::istreambuf_iterator<char>(ifs),
                                      std::istreambuf_iterator<char>());
    return content;
}

/* Read tree from file and fill tree_data */
IOStatus read_bptree_opaque(const char* tree_filename, opaque_bptree_t** tree_data) {
    SETUP_TDBG("read_bptree_opaque")
    if(tree_data==NULL) return unexpected_end;
    CHECK_FILE(tree_filename, open_error)
    TDBG_STEP("load_files")
    *tree_data = (opaque_bptree_t*) new su::BPTree(get_tree_content(tree_filename));
    return read_okay;
}
#endif // UNIFRAC_WASM (file-reading helpers; WASM uses load_bptree_opaque from a Newick string)

void load_bptree_opaque(const char* newick, opaque_bptree_t** tree_data) {
    SETUP_TDBG("load_bptree_opaque")
    *tree_data = (opaque_bptree_t*) new su::BPTree(newick);
}

void convert_bptree_opaque(const support_bptree_t* in_tree, opaque_bptree_t** tree_data) {
    SETUP_TDBG("convert_bptree_opaque")
    *tree_data = (opaque_bptree_t*) new su::BPTree(
                    in_tree->structure,
                    in_tree->lengths,
                    in_tree->names,
                    in_tree->n_parens);
}

void destroy_bptree_opaque(opaque_bptree_t** tree_data) {
	if (tree_data!=NULL) {
		su::BPTree *tree = (su::BPTree *) (*tree_data);
		*tree_data = NULL;
		delete tree;
	}
}

/* Return number of elements in BPTree, equvalent to n_parens */
int get_bptree_opaque_els(opaque_bptree_t* tree_data) {
	if (tree_data!=NULL) {
		su::BPTree *tree = (su::BPTree *) tree_data;
		return tree->get_tip_names().size();
	}
	return 0; // just to have a reasonable default
}

void destroy_stripes(vector<double*> &dm_stripes, vector<double*> &dm_stripes_total, unsigned int n_samples,
                     unsigned int stripe_start, unsigned int stripe_stop) {
    unsigned int n_rotations = (n_samples + 1) / 2;

    if(stripe_stop == 0) {
        for(unsigned int i = 0; i < n_rotations; i++) {
            free(dm_stripes[i]);
            if(dm_stripes_total[i] != NULL)
                free(dm_stripes_total[i]);
        }
    } else {
        // if a stripe_stop is specified, and if we're in the stripe window, do not free
        // dm_stripes. this is done as the pointers in dm_stripes are assigned to the partial_mat_t
        // and subsequently freed in destroy_partial_mat. but, we do need to free dm_stripes_total
        // if appropriate
        for(unsigned int i = stripe_start; i < stripe_stop; i++) {
            if(dm_stripes_total[i] != NULL)
                free(dm_stripes_total[i]);
        }
    }
}


void initialize_mat(mat_t* &result, biom_interface &table, bool is_upper_triangle) {
    result = (mat_t*)malloc(sizeof(mat));
    result->n_samples = table.n_samples;

    result->cf_size = su::comb_2(table.n_samples);
    result->is_upper_triangle = is_upper_triangle;
    result->sample_ids = (char**)malloc(sizeof(char*) * result->n_samples);
    result->condensed_form = (double*)malloc(sizeof(double) * su::comb_2(table.n_samples));

    const std::vector<std::string> &table_sample_ids = table.get_sample_ids();
    for(unsigned int i = 0; i < result->n_samples; i++) {
        size_t len = table_sample_ids[i].length();
        result->sample_ids[i] = (char*)malloc(sizeof(char) * len + 1);
        table_sample_ids[i].copy(result->sample_ids[i], len);
        result->sample_ids[i][len] = '\0';
    }
}

void initialize_results_vec(r_vec* &result, biom_interface &table){
    // Stores results for Faith PD
    result = (r_vec*)malloc(sizeof(results_vec));
    result->n_samples = table.n_samples;
    result->values = (double*)malloc(sizeof(double) * result->n_samples);
    result->sample_ids = (char**)malloc(sizeof(char*) * result->n_samples);

    const std::vector<std::string> &table_sample_ids = table.get_sample_ids();
    for(unsigned int i = 0; i < result->n_samples; i++) {
        size_t len = table_sample_ids[i].length();
        result->sample_ids[i] = (char*)malloc(sizeof(char) * len + 1);
        table_sample_ids[i].copy(result->sample_ids[i], len);
        result->sample_ids[i][len] = '\0';
        result->values[i] = 0;
    }

}

void initialize_mat_no_biom(mat_t* &result, char** sample_ids, unsigned int n_samples, bool is_upper_triangle) {
    result = (mat_t*)malloc(sizeof(mat));
    result->n_samples = n_samples;

    result->cf_size = su::comb_2(n_samples);
    result->is_upper_triangle = is_upper_triangle;
    result->sample_ids = (char**)malloc(sizeof(char*) * result->n_samples);
    result->condensed_form = (double*)malloc(sizeof(double) * su::comb_2(n_samples));

    for(unsigned int i = 0; i < n_samples; i++) {
        result->sample_ids[i] = strdup(sample_ids[i]);
    }
}

inline compute_status is_fp64_method(const std::string &method_string, bool &fp64) {
    if ((method_string=="unweighted") || (method_string=="weighted_normalized") || (method_string=="weighted_unnormalized") || (method_string=="unweighted_unnormalized") || (method_string=="generalized")) {
        fp64 = false;
    } else if ((method_string=="unweighted_fp64") || (method_string=="weighted_normalized_fp64") || (method_string=="weighted_unnormalized_fp64") || (method_string=="unweighted_unnormalized_fp64") || (method_string=="generalized_fp64")) {
       fp64 = true;
    } else if ((method_string=="unweighted_fp32") || (method_string=="weighted_normalized_fp32") || (method_string=="weighted_unnormalized_fp32") || (method_string=="unweighted_unnormalized_fp32") || (method_string=="generalized_fp32")) {
       fp64 = false;
    } else {
        return unknown_method;
    }

    return okay;
}


inline compute_status is_fp64(const std::string &method_string, const std::string &format_string, bool &fp64, bool &save_dist) {
  if (format_string == "hdf5_fp32") {
    fp64 = false;
    save_dist = true;
  } else if (format_string == "hdf5_fp64") {
    fp64 = true;
    save_dist = true;
  } else if (format_string == "hdf5") {
    save_dist = true;
    return is_fp64_method(method_string, fp64);
  } else if (format_string == "hdf5_nodist") {
    save_dist = false;
    return is_fp64_method(method_string, fp64);
  } else {
    return unknown_method; 
  }

  return okay;
}


template<class TReal, class TMat>
void initialize_mat_full_no_biom_T(TMat* &result, const char* const * sample_ids, unsigned int n_samples, 
                                   const char *mmap_dir /* if NULL or "", use malloc */) {
    result = (TMat*)malloc(sizeof(mat));
    result->n_samples = n_samples;

    uint64_t n_samples_64 = result->n_samples; // force 64bit to avoit overflow problems

    result->sample_ids = (char**)malloc(sizeof(char*) * n_samples_64);
    result->flags=0;

    if (mmap_dir!=NULL) {
     if (mmap_dir[0]==0) mmap_dir = NULL; // easier to have a simple test going on
    }

    uint64_t msize = sizeof(TReal) * n_samples_64 * n_samples_64;
    if (mmap_dir==NULL) {
      result->matrix = (TReal*)malloc(msize);
    }
#ifndef UNIFRAC_WASM
    else {
      std::string mmap_template(mmap_dir);
      mmap_template+="/su_mmap_XXXXXX";
      // note: mkstemp/mkostemp will update mmap_template in place
#ifdef O_NOATIME
      int fd=mkostemp((char *) mmap_template.c_str(), O_NOATIME );
#else
      int fd=mkstemp((char *) mmap_template.c_str() );
#endif
      if (fd<0) {
         result->matrix = NULL;
         // leave error handling to the caller
      } else {
        // remove the file name, so it will be destroyed on close
        unlink(mmap_template.c_str());
        // make it big enough
        if (ftruncate(fd,msize)!=0) {
           // fail
           result->matrix = NULL;
           // leave error handling to the caller
	} else {
           // now can be used, just like a malloc-ed buffer
           result->matrix = (TReal*)mmap(NULL, msize,PROT_READ|PROT_WRITE, MAP_SHARED|MAP_NORESERVE, fd, 0);
           result->flags=(uint32_t(fd) & MMAP_FD_MASK) | MMAP_FLAG;
        }
      }
   }
#endif // UNIFRAC_WASM (mmap path; WASM always takes the malloc branch)

    for(unsigned int i = 0; i < n_samples; i++) {
        result->sample_ids[i] = strdup(sample_ids[i]);
    }
}

void initialize_partial_mat(partial_mat_t* &result, biom_interface &table, std::vector<double*> &dm_stripes,
                            unsigned int stripe_start, unsigned int stripe_stop, bool is_upper_triangle) {
    result = (partial_mat_t*)malloc(sizeof(partial_mat));
    result->n_samples = table.n_samples;

    result->sample_ids = (char**)malloc(sizeof(char*) * result->n_samples);
    const std::vector<std::string> &table_sample_ids = table.get_sample_ids();
    for(unsigned int i = 0; i < result->n_samples; i++) {
        size_t len = table_sample_ids[i].length();
        result->sample_ids[i] = (char*)malloc(sizeof(char) * len + 1);
        table_sample_ids[i].copy(result->sample_ids[i], len);
        result->sample_ids[i][len] = '\0';
    }

    result->stripes = (double**)malloc(sizeof(double*) * (stripe_stop - stripe_start));
    result->stripe_start = stripe_start;
    result->stripe_stop = stripe_stop;
    result->is_upper_triangle = is_upper_triangle;
    result->stripe_total = dm_stripes.size();

    for(unsigned int i = stripe_start; i < stripe_stop; i++) {
        result->stripes[i - stripe_start] = dm_stripes[i];
    }
}

void destroy_results_vec(r_vec** result) {
    // for Faith PD
    for(unsigned int i = 0; i < (*result)->n_samples; i++) {
        free((*result)->sample_ids[i]);
    };
    free((*result)->sample_ids);
    free((*result)->values);
    free(*result);
}

void destroy_mat(mat_t** result) {
    for(unsigned int i = 0; i < (*result)->n_samples; i++) {
        free((*result)->sample_ids[i]);
    };
    free((*result)->sample_ids);
    if (((*result)->condensed_form)!=NULL) {
      free((*result)->condensed_form);
    }
    free(*result);
}

template<class TMat, class TReal>
inline void destroy_mat_full_T(TMat** result) {
    for(uint32_t i = 0; i < (*result)->n_samples; i++) {
        free((*result)->sample_ids[i]);   
    };                                        
    free((*result)->sample_ids);          
    if (((*result)->matrix)!=NULL) {
      if (((*result)->flags & MMAP_FLAG) == 0)  {
         free((*result)->matrix);
      }
#ifndef UNIFRAC_WASM
      else {
         uint64_t n_samples = (*result)->n_samples;
         munmap((*result)->matrix, sizeof(TReal)*n_samples*n_samples);

         int fd = (*result)->flags & MMAP_FD_MASK;
         close(fd);
      }
#endif // UNIFRAC_WASM (no mmap-backed matrices in WASM build)
      (*result)->matrix=NULL;
    }
    free(*result);
}


void destroy_mat_full_fp64(mat_full_fp64_t** result) {
    destroy_mat_full_T<mat_full_fp64_t,double>(result);
}

void destroy_mat_full_fp32(mat_full_fp32_t** result) {
    destroy_mat_full_T<mat_full_fp32_t,float>(result);
}

void destroy_partial_mat(partial_mat_t** result) {
    for(unsigned int i = 0; i < (*result)->n_samples; i++) {
        if((*result)->sample_ids[i] != NULL)
            free((*result)->sample_ids[i]);
    };
    if((*result)->sample_ids != NULL)
        free((*result)->sample_ids);

    unsigned int n_stripes = (*result)->stripe_stop - (*result)->stripe_start;
    for(unsigned int i = 0; i < n_stripes; i++)
        if((*result)->stripes[i] != NULL)
            free((*result)->stripes[i]);
    if((*result)->stripes != NULL)
        free((*result)->stripes);

    free(*result);
}

void destroy_partial_dyn_mat(partial_dyn_mat_t** result) {
    for(unsigned int i = 0; i < (*result)->n_samples; i++) {
        if((*result)->sample_ids[i] != NULL)
            free((*result)->sample_ids[i]);
    };
    if((*result)->sample_ids != NULL)
        free((*result)->sample_ids);

    unsigned int n_stripes = (*result)->stripe_stop - (*result)->stripe_start;
    for(unsigned int i = 0; i < n_stripes; i++)
        if((*result)->stripes[i] != NULL)
            free((*result)->stripes[i]);
    if((*result)->stripes != NULL)
        free((*result)->stripes);
    if((*result)->offsets != NULL)
        free((*result)->offsets);
    if((*result)->filename != NULL)
        free((*result)->filename);

    free(*result);
}


/* A stripe_stop that does not make sense means "through the last stripe".
 * set_tasks applies that default internally, and callers need the same answer to
 * work out how many stripes their tasks are going to divide, so the rule lives
 * in exactly one place.
 */
static unsigned int effective_stripe_stop(unsigned int n_samples,
                                          unsigned int stripe_start,
                                          unsigned int stripe_stop) {
    return (stripe_stop <= stripe_start) ? ((n_samples + 1) / 2) : stripe_stop;
}

/* Force a caller-supplied n_substeps into the range set_tasks can actually
 * divide the stripes into. n_substeps only says how to split that range across
 * tasks, so out-of-range values are a request to be normalized rather than an
 * error:
 *
 *   0                        divides by zero in set_tasks below
 *   > n_stripes_in_range     leaves trailing tasks with an empty stripe range
 *                            whose start indexes one past the end of dm_stripes
 *                            (dereferenced unconditionally in
 *                            UnifracTaskVector, unifrac_task.hpp)
 *
 * n_stripes_in_range is stripe_stop - stripe_start, the number of stripes these
 * tasks will divide -- NOT the size of dm_stripes. The two differ for partial
 * computes, which allocate the full stripe vector but compute a sub-range.
 *
 * Must be called before sizing the tasks vector, since that is sized by the
 * returned value.
 */
static unsigned int clamp_substeps(unsigned int n_substeps, unsigned int n_stripes_in_range) {
    if (n_stripes_in_range < 1) n_stripes_in_range = 1;  // an empty table is rejected upstream
    if (n_substeps > n_stripes_in_range) {
        fprintf(stderr, "More substeps were requested than stripes. Using %u substeps.\n", n_stripes_in_range);
        return n_stripes_in_range;
    }
    if (n_substeps < 1) return 1;
    return n_substeps;
}

void set_tasks(std::vector<su::task_parameters> &tasks,
               double alpha,
               unsigned int n_samples,
               unsigned int stripe_start,
               unsigned int stripe_stop,
               bool bypass_tips,
	       bool normalize_sample_counts,
               unsigned int n_tasks) {

    // compute from start to the max possible stripe if stop doesn't make sense
    stripe_stop = effective_stripe_stop(n_samples, stripe_start, stripe_stop);

    /* chunking strategy is to balance as much as possible. eg if there are 15 stripes
     * and 4 threads, our goal is to assign 4 stripes to 3 threads, and 3 stripes to one thread.
     *
     * we use the remaining the chunksize for bins which cannot be full maximally
     */
    unsigned int fullchunk = ((stripe_stop - stripe_start) + n_tasks - 1) / n_tasks;  // this computes the ceiling
    unsigned int smallchunk = (stripe_stop - stripe_start) / n_tasks;

    unsigned int n_fullbins = (stripe_stop - stripe_start) % n_tasks;
    if(n_fullbins == 0)
        n_fullbins = n_tasks;

    unsigned int start = stripe_start;

    for(unsigned int tid = 0; tid < n_tasks; tid++) {
        tasks[tid].tid = tid;
        tasks[tid].start = start; // stripe start
        tasks[tid].bypass_tips = bypass_tips;
        tasks[tid].normalize_sample_counts = normalize_sample_counts;

        if(tid < n_fullbins) {
            tasks[tid].stop = start + fullchunk;  // stripe end
            start = start + fullchunk;
        } else {
            tasks[tid].stop = start + smallchunk;  // stripe end
            start = start + smallchunk;
        }

        tasks[tid].n_samples = n_samples;
        tasks[tid].g_unifrac_alpha = alpha;
    }
}

compute_status one_off_inmem_cpp(su::biom_interface &table, const su::BPTree &tree,
                             const char* unifrac_method, bool variance_adjust, double alpha,
                             bool bypass_tips, bool normalize_sample_counts, unsigned int n_substeps, mat_t** result) {
    SETUP_TDBG("one_off_inmem")
    SYNC_TREE_TABLE(tree, table)
    TDBG_STEP("sync_tree_table")
    SET_METHOD(unifrac_method, unknown_method)

    const unsigned int stripe_stop = (table.n_samples + 1) / 2;
    std::vector<double*> dm_stripes(stripe_stop);
    std::vector<double*> dm_stripes_total(stripe_stop);

    // whole range, so stripe_start is 0
    n_substeps = clamp_substeps(n_substeps, stripe_stop);

    std::vector<su::task_parameters> tasks(n_substeps);

    set_tasks(tasks, alpha, table.n_samples, 0, stripe_stop, bypass_tips, normalize_sample_counts, n_substeps);
    su::process_stripes(table, tree_sheared, method, variance_adjust, dm_stripes, dm_stripes_total, tasks);

    TDBG_STEP("process_stripes")
    initialize_mat(*result, table, true);  // true -> is_upper_triangle
    for(unsigned int tid = 0; tid < n_substeps; tid++) {
        su::stripes_to_condensed_form(dm_stripes,table.n_samples,(*result)->condensed_form,tasks[tid].start,tasks[tid].stop);
    }

    TDBG_STEP("stripes_to_condensed_form")
    destroy_stripes(dm_stripes, dm_stripes_total, table.n_samples, 0, 0);

    return okay;
}

#ifndef UNIFRAC_WASM
compute_status partial_v3(const char* biom_filename, const char* tree_filename,
                          const char* unifrac_method, bool variance_adjust, double alpha, bool bypass_tips, bool normalize_sample_counts,
                          unsigned int n_substeps, unsigned int stripe_start, unsigned int stripe_stop,
                          partial_mat_t** result) {

    SETUP_TDBG("partial")
    CHECK_FILE(biom_filename, table_missing)
    CHECK_FILE(tree_filename, tree_missing)
    SET_METHOD(unifrac_method, unknown_method)
    PARSE_SYNC_TREE_TABLE(tree_filename, table_filename)

    TDBG_STEP("load_files")
    // we resize to the largest number of possible stripes even if only computing
    // partial, however we do not allocate arrays for non-computed stripes so
    // there is a little memory waste here but should be on the order of
    // 8 bytes * N samples per vector.
    std::vector<double*> dm_stripes((table.n_samples + 1) / 2);
    std::vector<double*> dm_stripes_total((table.n_samples + 1) / 2);

    /* dm_stripes covers every stripe, but the tasks only divide the requested
     * sub-range, so that -- not dm_stripes.size() -- is what bounds n_substeps.
     */
    n_substeps = clamp_substeps(n_substeps,
                                effective_stripe_stop(table.n_samples, stripe_start, stripe_stop) - stripe_start);

    std::vector<su::task_parameters> tasks(n_substeps);

    if(((table.n_samples + 1) / 2) < stripe_stop) {
        fprintf(stderr, "Stopping stripe is out-of-bounds, max %d\n", (table.n_samples + 1) / 2);
        exit(EXIT_FAILURE);
    }

    set_tasks(tasks, alpha, table.n_samples, stripe_start, stripe_stop, bypass_tips, normalize_sample_counts, n_substeps);
    su::process_stripes(table, tree_sheared, method, variance_adjust, dm_stripes, dm_stripes_total, tasks);

    TDBG_STEP("process_stripes")
    initialize_partial_mat(*result, table, dm_stripes, stripe_start, stripe_stop, true);  // true -> is_upper_triangle
    TDBG_STEP("partial_mat")
    destroy_stripes(dm_stripes, dm_stripes_total, table.n_samples, stripe_start, stripe_stop);

    return okay;
}

compute_status faith_pd_one_off(const char* biom_filename, const char* tree_filename,
                                r_vec** result){
    SETUP_TDBG("faith_pd_one_off")
    CHECK_FILE(biom_filename, table_missing)
    CHECK_FILE(tree_filename, tree_missing)
    PARSE_TREE_TABLE(tree_filename, biom_filename)
    TDBG_STEP("load_files")

    // Filter out any elements with zero counts
    su::biom_inmem table_nz(table,1.0);
    if ((table_nz.n_samples==0) || (table_nz.n_obs==0)) {
      fprintf(stderr, "WARNING: All samples had zero counts. Forcing zero result.\n");
      SYNC_TREE_TABLE(tree, table)

      TDBG_STEP("sync_tree_table")
      initialize_results_vec(*result, table);
      // nothing else to do... results already initialized to 0
      TDBG_STEP("faith_pd")
    } else {
      if ((table_nz.n_samples!=table.n_samples) || (table_nz.n_obs!=table.n_obs)) {
        fprintf(stderr, "WARNING: Some samples had zero counts and were filtered out.\n");
      }
      SYNC_TREE_TABLE(tree, table_nz)

      TDBG_STEP("sync_tree_table")
    
      initialize_results_vec(*result, table_nz);

      // compute faithpd
      su::faith_pd(table_nz, tree_sheared, std::ref((*result)->values));
      TDBG_STEP("faith_pd")
    }

    return okay;
}

compute_status one_off_v3(const char* biom_filename, const char* tree_filename,
                          const char* unifrac_method, bool variance_adjust, double alpha,
                          bool bypass_tips, bool normalize_sample_counts, unsigned int n_substeps, mat_t** result) {
    SETUP_TDBG("one_off")
    CHECK_FILE(biom_filename, table_missing)
    CHECK_FILE(tree_filename, tree_missing)
    PARSE_TREE_TABLE(tree_filename, table_filename)

    TDBG_STEP("load_files")
    // condensed form
    return one_off_inmem_cpp(table, tree, unifrac_method, variance_adjust, alpha, bypass_tips, normalize_sample_counts, n_substeps, result);
}

compute_status one_off_wtree_v3(const char* biom_filename, const opaque_bptree_t* tree_data,
                                const char* unifrac_method, bool variance_adjust, double alpha,
                                bool bypass_tips, bool normalize_sample_counts, unsigned int n_substeps, mat_t** result) {
    SETUP_TDBG("one_off_wtree")
    if (tree_data==NULL) return tree_missing;
    CHECK_FILE(biom_filename, table_missing)
    const su::BPTree &tree = *( (const su::BPTree*) tree_data);
    su::biom table(biom_filename);
    VALIDATE_TREE_TABLE(tree, table)

    TDBG_STEP("load_files")
    // condensed form
    return one_off_inmem_cpp(table, tree, unifrac_method, variance_adjust, alpha, bypass_tips, normalize_sample_counts, n_substeps, result);
}
#endif // UNIFRAC_WASM (file-based v3 entries)

/*
 * ==============================   one_off_matrix
 */

// TMat mat_full_fp32_t
template<class TReal, class TMat>
compute_status one_off_matrix_T(su::biom_interface &table, const su::BPTree &tree,
                                const char* unifrac_method, bool variance_adjust, double alpha,
                                bool bypass_tips, bool normalize_sample_counts, unsigned int n_substeps,
                                const char *mmap_dir,
                                TMat** result) {
    SETUP_TDBG("one_off_matrix_inmem")
#ifdef UNIFRAC_WASM
    // mmap-backed matrices are not supported under WASM (no real fd
    // semantics in the browser sandbox). Force malloc path regardless of
    // what the caller passed.
    mmap_dir = NULL;
#endif
    if (mmap_dir!=NULL) {
     if (mmap_dir[0]==0) mmap_dir = NULL; // easier to have a simple test going on
    }

    SET_METHOD(unifrac_method, unknown_method)
    SYNC_TREE_TABLE(tree, table)

    TDBG_STEP("sync_tree_table")
    const unsigned int stripe_stop = (table.n_samples + 1) / 2;
    partial_mat_t *partial_mat = NULL;

    {
      std::vector<double*> dm_stripes(stripe_stop);
      std::vector<double*> dm_stripes_total(stripe_stop);

      // whole range, so stripe_start is 0
      n_substeps = clamp_substeps(n_substeps, stripe_stop);

      std::vector<su::task_parameters> tasks(n_substeps);

      set_tasks(tasks, alpha, table.n_samples, 0, stripe_stop, bypass_tips, normalize_sample_counts, n_substeps);
      su::process_stripes(table, tree_sheared, method, variance_adjust, dm_stripes, dm_stripes_total, tasks);

      TDBG_STEP("process_stripes")
      initialize_partial_mat(partial_mat, table, dm_stripes, 0, stripe_stop, true);  // true -> is_upper_triangle
      if ((partial_mat==NULL) || (partial_mat->stripes==NULL) || (partial_mat->sample_ids==NULL) ) {
          fprintf(stderr, "Memory allocation error! (initialize_partial_mat)\n");
          exit(EXIT_FAILURE);
      }
      destroy_stripes(dm_stripes, dm_stripes_total, table.n_samples, 0, stripe_stop);
    }

    // allow the caller to allocate the memory
    if((*result) == NULL) {
        initialize_mat_full_no_biom_T<TReal,TMat>(*result, partial_mat->sample_ids, partial_mat->n_samples,mmap_dir);
    }

    if (((*result)==NULL) || ((*result)->matrix==NULL) || ((*result)->sample_ids==NULL) ) {
        fprintf(stderr, "Memory allocation error! (initialize_mat)\n");
        exit(EXIT_FAILURE);
    }


    {
      MemoryStripes ps(partial_mat->stripes);
      const uint32_t tile_size = (mmap_dir==NULL) ? \
                                  (128/sizeof(TReal)) : /* keep it small for memory access, to fit in chip cache */ \
                                  (4096/sizeof(TReal)); /* make it larger for mmap, as the limiting factor is swapping */
      su::stripes_to_matrix_T<TReal>(ps, partial_mat->n_samples, partial_mat->stripe_total, (*result)->matrix, tile_size);
    }
    TDBG_STEP("stripes_to_matrix")
    destroy_partial_mat(&partial_mat);

    return okay;
}


template<class TReal, class TMat>
compute_status one_off_matrix_v4_T(su::biom_inmem &table, const su::BPTree &tree,
                                   const char* unifrac_method, bool variance_adjust, double alpha,
                                   bool bypass_tips, bool normalize_sample_counts, unsigned int n_substeps,
                                   unsigned int subsample_depth, bool subsample_with_replacement, int seed,
                                   const char *mmap_dir,
                                   TMat** result) {
    SETUP_TDBG("one_off_matrix_inmem_v4")
    if (subsample_depth>0) {
        /* Same seeding rule as subsample_table_inmem_seeded():
         *   seed >= 0  deterministic draw from the explicit seed, touching no
         *              shared state, so concurrent callers do not need a lock
         *   seed <  0  skbio_biom_subsampled draws from the process-global
         *              mt19937 that ssu_set_random_seed() sets, which is the
         *              legacy behaviour
         * The two branches differ in type, so hold the result by base pointer;
         * one_off_matrix_T takes a biom_interface&.
         */
        std::unique_ptr<su::biom_inmem> table_subsampled;
        if (seed < 0) {
            table_subsampled.reset(new su::skbio_biom_subsampled(table, subsample_with_replacement, subsample_depth));
        } else {
            table_subsampled.reset(new su::biom_subsampled(table, subsample_with_replacement, subsample_depth, (uint32_t) seed));
        }
        if ((table_subsampled->n_samples==0) || (table_subsampled->n_obs==0)) {
           return table_empty;
        }
        TDBG_STEP("subsample")
        return one_off_matrix_T<TReal,TMat>(*table_subsampled,tree,unifrac_method,variance_adjust,alpha,bypass_tips,normalize_sample_counts,n_substeps,mmap_dir,result);
    } else {
        return one_off_matrix_T<TReal,TMat>(table,tree,unifrac_method,variance_adjust,alpha,bypass_tips,normalize_sample_counts,n_substeps,mmap_dir,result);
    }
}

#ifndef UNIFRAC_WASM
compute_status one_off_matrix_v3(const char* biom_filename, const char* tree_filename,
                                 const char* unifrac_method, bool variance_adjust, double alpha,
                                 bool bypass_tips, bool normalize_sample_counts, unsigned int n_substeps,
                                 unsigned int subsample_depth, bool subsample_with_replacement, const char *mmap_dir,
                                 mat_full_fp64_t** result) {
    SETUP_TDBG("one_off_matrix")
    CHECK_FILE(biom_filename, table_missing)
    CHECK_FILE(tree_filename, tree_missing)
    PARSE_TREE_TABLE(tree_filename, biom_filename)
    TDBG_STEP("load_files")
    return one_off_matrix_v4_T<double,mat_full_fp64_t>(table,tree,unifrac_method,variance_adjust,alpha,bypass_tips,normalize_sample_counts,n_substeps,subsample_depth,subsample_with_replacement,/*seed*/ -1,mmap_dir,result);
}

compute_status one_off_matrix_fp32_v3(const char* biom_filename, const char* tree_filename,
                                      const char* unifrac_method, bool variance_adjust, double alpha,
                                      bool bypass_tips, bool normalize_sample_counts, unsigned int n_substeps,
                                      unsigned int subsample_depth, bool subsample_with_replacement, const char *mmap_dir,
                                      mat_full_fp32_t** result) {
    SETUP_TDBG("one_off_matrix_fp32")
    CHECK_FILE(biom_filename, table_missing)
    CHECK_FILE(tree_filename, tree_missing)
    PARSE_TREE_TABLE(tree_filename, biom_filename)
    TDBG_STEP("load_files")
    return one_off_matrix_v4_T<float,mat_full_fp32_t>(table,tree,unifrac_method,variance_adjust,alpha,bypass_tips,normalize_sample_counts,n_substeps,subsample_depth,subsample_with_replacement,/*seed*/ -1,mmap_dir,result);
}

/* As above, but from a pre-loaded tree object */
compute_status one_off_matrix_v3t(const char* biom_filename, const opaque_bptree_t* tree_data,
                                 const char* unifrac_method, bool variance_adjust, double alpha,
                                 bool bypass_tips, bool normalize_sample_counts, unsigned int n_substeps,
                                 unsigned int subsample_depth, bool subsample_with_replacement, const char *mmap_dir,
                                 mat_full_fp64_t** result) {
    SETUP_TDBG("one_off_matrix_wtree")
    if (tree_data==NULL) return tree_missing;
    CHECK_FILE(biom_filename, table_missing)
    const su::BPTree &tree = *( (const su::BPTree*) tree_data);
    su::biom table(biom_filename);
    VALIDATE_TREE_TABLE(tree, table)
    TDBG_STEP("load_files")
    return one_off_matrix_v4_T<double,mat_full_fp64_t>(table,tree,unifrac_method,variance_adjust,alpha,bypass_tips,normalize_sample_counts,n_substeps,subsample_depth,subsample_with_replacement,/*seed*/ -1,mmap_dir,result);
}

compute_status one_off_matrix_fp32_v3t(const char* biom_filename, const opaque_bptree_t* tree_data,
                                      const char* unifrac_method, bool variance_adjust, double alpha,
                                      bool bypass_tips, bool normalize_sample_counts, unsigned int n_substeps,
                                      unsigned int subsample_depth, bool subsample_with_replacement, const char *mmap_dir,
                                      mat_full_fp32_t** result) {
    SETUP_TDBG("one_off_matrix_fp32_wtree")
    if (tree_data==NULL) return tree_missing;
    CHECK_FILE(biom_filename, table_missing)
    const su::BPTree &tree = *( (const su::BPTree*) tree_data);
    su::biom table(biom_filename);
    VALIDATE_TREE_TABLE(tree, table)
    TDBG_STEP("load_files")
    return one_off_matrix_v4_T<float,mat_full_fp32_t>(table,tree,unifrac_method,variance_adjust,alpha,bypass_tips,normalize_sample_counts,n_substeps,subsample_depth,subsample_with_replacement,/*seed*/ -1,mmap_dir,result);
}
#endif // UNIFRAC_WASM (file-based one_off_matrix wrappers)

compute_status one_off_matrix_inmem_v4(const support_biom_t *table_data, const support_bptree_t *tree_data,
                                       const char* unifrac_method, bool variance_adjust, double alpha,
                                       bool bypass_tips, bool normalize_sample_counts, unsigned int n_substeps,
                                       unsigned int subsample_depth, bool subsample_with_replacement, int seed,
                                       int device_id, const char *mmap_dir,
                                       mat_full_fp64_t** result) {
    SETUP_TDBG("one_off_matrix_inmem")
    /* Device-resident input and output are not implemented. Rejected before any
     * work so a caller that asks cannot mistake a host-computed answer for a
     * device-computed one.
     */
    if (device_id >= 0) return unsupported_device;

    bool fp64;
    compute_status rc = is_fp64_method(unifrac_method, fp64);

    if (rc == okay) {
        if (!fp64) {
            return invalid_method;
        }
    } else {
        return rc;
    }

    if( (table_data->n_samples <= 0) || (table_data->n_obs <= 0) ) {
	return table_empty;
    }

    su::biom_inmem table(table_data->obs_ids, 
                         table_data->sample_ids, 
                         table_data->indices,
                         table_data->indptr,
                         table_data->data,
                         table_data->n_obs,
                         table_data->n_samples);

    su::BPTree tree(tree_data->structure,
                    tree_data->lengths,
                    tree_data->names,
                    tree_data->n_parens);

    VALIDATE_TREE_TABLE(tree,table)

    return one_off_matrix_v4_T<double,mat_full_fp64_t>(table,tree,unifrac_method,variance_adjust,alpha,bypass_tips,normalize_sample_counts,n_substeps,subsample_depth,subsample_with_replacement,seed,mmap_dir,result);
}

compute_status one_off_matrix_inmem_fp32_v4(const support_biom_t *table_data, const support_bptree_t *tree_data,
                                            const char* unifrac_method, bool variance_adjust, double alpha,
                                            bool bypass_tips, bool normalize_sample_counts, unsigned int n_substeps,
                                            unsigned int subsample_depth, bool subsample_with_replacement, int seed,
                                            int device_id, const char *mmap_dir,
                                            mat_full_fp32_t** result) {
    SETUP_TDBG("one_off_matrix_inmem_fp32")
    // see one_off_matrix_inmem_v4
    if (device_id >= 0) return unsupported_device;

    bool fp64;
    compute_status rc = is_fp64_method(unifrac_method, fp64);

    if (rc == okay) {
        if (fp64) {
            return invalid_method;
        }
    } else {
        return rc;
    }

    if( (table_data->n_samples <= 0) || (table_data->n_obs <= 0) ) {
	return table_empty;
    }

    su::biom_inmem table(table_data->obs_ids, 
                         table_data->sample_ids, 
                         table_data->indices,
                         table_data->indptr,
                         table_data->data,
                         table_data->n_obs,
                         table_data->n_samples);

    su::BPTree tree(tree_data->structure,
                    tree_data->lengths,
                    tree_data->names,
                    tree_data->n_parens);

    VALIDATE_TREE_TABLE(tree,table)

    return one_off_matrix_v4_T<float,mat_full_fp32_t>(table,tree,unifrac_method,variance_adjust,alpha,bypass_tips,normalize_sample_counts,n_substeps,subsample_depth,subsample_with_replacement,seed,mmap_dir,result);
}

compute_status faith_pd_inmem(const support_biom_t *table_data,
                              const support_bptree_t *tree_data,
                              r_vec** result) {
    SETUP_TDBG("faith_pd_inmem")
    if (tree_data == NULL) return tree_missing;
    if (table_data == NULL) return table_missing;

    su::biom_inmem table(table_data->obs_ids,
                         table_data->sample_ids,
                         table_data->indices,
                         table_data->indptr,
                         table_data->data,
                         table_data->n_obs,
                         table_data->n_samples);
    su::BPTree tree(tree_data->structure,
                    tree_data->lengths,
                    tree_data->names,
                    tree_data->n_parens);

    VALIDATE_TREE_TABLE(tree, table)
    TDBG_STEP("load_files")

    // Filter out any elements with zero counts (mirrors faith_pd_one_off).
    su::biom_inmem table_nz(table, 1.0);
    if ((table_nz.n_samples == 0) || (table_nz.n_obs == 0)) {
        fprintf(stderr, "WARNING: All samples had zero counts. Forcing zero result.\n");
        SYNC_TREE_TABLE(tree, table)
        TDBG_STEP("sync_tree_table")
        initialize_results_vec(*result, table);
    } else {
        if ((table_nz.n_samples != table.n_samples) || (table_nz.n_obs != table.n_obs)) {
            fprintf(stderr, "WARNING: Some samples had zero counts and were filtered out.\n");
        }
        SYNC_TREE_TABLE(tree, table_nz)
        TDBG_STEP("sync_tree_table")

        initialize_results_vec(*result, table_nz);
        su::faith_pd(table_nz, tree_sheared, std::ref((*result)->values));
        TDBG_STEP("faith_pd")
    }

    return okay;
}

compute_status subsample_table_inmem_seeded(const support_biom_t *table_data,
                                            unsigned int depth,
                                            bool with_replacement,
                                            int seed,
                                            opaque_biom_inmem_t **out) {
    SETUP_TDBG("subsample_table_inmem_seeded")
    if (table_data == NULL) return table_missing;
    if ((table_data->n_samples <= 0) || (table_data->n_obs <= 0)) {
        return table_empty;
    }

    su::biom_inmem table(table_data->obs_ids,
                         table_data->sample_ids,
                         table_data->indices,
                         table_data->indptr,
                         table_data->data,
                         table_data->n_obs,
                         table_data->n_samples);
    TDBG_STEP("load_table")

    // seed >= 0: deterministic draw from the explicit seed.
    // seed <  0: skbio_biom_subsampled draws from the global skbb-side
    //            mt19937 (re-seedable via ssu_set_random_seed) — matches
    //            the legacy non-seeded subsample_table_inmem behavior.
    // Heap-allocate either way so the object survives the function
    // return as an opaque handle. biom_inmem has a virtual destructor,
    // so destroy_subsampled_inmem's `delete` chains correctly through
    // either subclass.
    su::biom_inmem *sub;
    if (seed < 0) {
        sub = new su::skbio_biom_subsampled(table, with_replacement, depth);
    } else {
        sub = new su::biom_subsampled(table, with_replacement, depth, (uint32_t) seed);
    }
    *out = (opaque_biom_inmem_t*) sub;
    TDBG_STEP("subsample")
    return okay;
}

compute_status subsample_table_inmem(const support_biom_t *table_data,
                                     unsigned int depth,
                                     bool with_replacement,
                                     opaque_biom_inmem_t **out) {
    return subsample_table_inmem_seeded(table_data, depth, with_replacement, -1, out);
}

// The opaque handle returned by subsample_table_inmem* is a
// su::biom_inmem* — the concrete type is either skbio_biom_subsampled
// (global-RNG path) or biom_subsampled (explicit-seed path). All the
// accessors below only depend on the biom_inmem interface, so the
// casts target the common base. biom_inmem has a virtual destructor,
// so destroy_subsampled_inmem's `delete` chains correctly through
// either subclass.

unsigned int subsampled_n_samples(const opaque_biom_inmem_t *t) {
    if (t == NULL) return 0;
    const su::biom_inmem *sub = (const su::biom_inmem*) t;
    return sub->n_samples;
}

unsigned int subsampled_n_obs(const opaque_biom_inmem_t *t) {
    if (t == NULL) return 0;
    const su::biom_inmem *sub = (const su::biom_inmem*) t;
    return sub->n_obs;
}

bool subsampled_get_obs_data(const opaque_biom_inmem_t *t,
                             const char *obs_id,
                             double *out) {
    if (t == NULL || obs_id == NULL || out == NULL) return false;
    const su::biom_inmem *sub = (const su::biom_inmem*) t;
    std::string id_str(obs_id);
    if (!sub->has_obs_id(id_str)) return false;
    sub->get_obs_data(id_str, out);
    return true;
}

const char* subsampled_get_sample_id(const opaque_biom_inmem_t *t, unsigned int idx) {
    if (t == NULL) return NULL;
    const su::biom_inmem *sub = (const su::biom_inmem*) t;
    const std::vector<std::string> &ids = sub->get_sample_ids();
    if (idx >= ids.size()) return NULL;
    return ids[idx].c_str();
}

const char* subsampled_get_obs_id(const opaque_biom_inmem_t *t, unsigned int idx) {
    if (t == NULL) return NULL;
    const su::biom_inmem *sub = (const su::biom_inmem*) t;
    const std::vector<std::string> &ids = sub->get_obs_ids();
    if (idx >= ids.size()) return NULL;
    return ids[idx].c_str();
}

void destroy_subsampled_inmem(opaque_biom_inmem_t **t) {
    if (t == NULL || *t == NULL) return;
    su::biom_inmem *sub = (su::biom_inmem*) (*t);
    *t = NULL;
    delete sub;
}

compute_status compute_permanova_inmem_fp64_seeded(const double *mat, unsigned int n_dims,
                                                    const uint32_t *grouping,
                                                    unsigned int permanova_perms,
                                                    int seed,
                                                    double *fstat, double *pvalue) {
    if (mat == NULL || grouping == NULL) return grouping_missing;
    su::permanova(mat, n_dims, grouping, permanova_perms, *fstat, *pvalue, seed);
    return okay;
}

compute_status compute_permanova_inmem_fp32_seeded(const float *mat, unsigned int n_dims,
                                                    const uint32_t *grouping,
                                                    unsigned int permanova_perms,
                                                    int seed,
                                                    float *fstat, float *pvalue) {
    if (mat == NULL || grouping == NULL) return grouping_missing;
    su::permanova(mat, n_dims, grouping, permanova_perms, *fstat, *pvalue, seed);
    return okay;
}

compute_status compute_permanova_inmem_fp64(const double *mat, unsigned int n_dims,
                                             const uint32_t *grouping,
                                             unsigned int permanova_perms,
                                             double *fstat, double *pvalue) {
    return compute_permanova_inmem_fp64_seeded(mat, n_dims, grouping, permanova_perms,
                                               -1, fstat, pvalue);
}

compute_status compute_permanova_inmem_fp32(const float *mat, unsigned int n_dims,
                                             const uint32_t *grouping,
                                             unsigned int permanova_perms,
                                             float *fstat, float *pvalue) {
    return compute_permanova_inmem_fp32_seeded(mat, n_dims, grouping, permanova_perms,
                                               -1, fstat, pvalue);
}

/*
 * ==============================   one_dense_pair
 */

/* Compute UniFrac from a pair of dense vectors */
static inline compute_status one_dense_pair(unsigned int n_obs, const char ** obs_ids, const double* sample1, const double* sample2,
		                  const su::BPTree &tree,
                                  const char* unifrac_method, bool variance_adjust, double alpha,
                                  bool bypass_tips, bool normalize_sample_counts, double* result) {
    SETUP_TDBG("one_dense_pair")
    bool fp64;
    compute_status rc = is_fp64_method(unifrac_method, fp64);

    if (rc != okay) return rc;

    const char* sample_ids[2] = { "Sample1", "Sample2" };
    const double* dense_data[2] = { sample1, sample2 };
    su::biom_inmem table(obs_ids, 
                         sample_ids,
                         dense_data, 
                         n_obs,
                         2);

    VALIDATE_TREE_TABLE(tree, table)

    TDBG_STEP("load_files")
    if (fp64) {
      mat_full_fp64_t* matrix_result = NULL;
      rc = one_off_matrix_T<double,mat_full_fp64_t>(table,tree,unifrac_method,variance_adjust,alpha,bypass_tips,normalize_sample_counts,1,NULL,&matrix_result);
      if (rc==okay) {
         *result = matrix_result->matrix[1];
	 destroy_mat_full_fp64(&matrix_result);
      }
    } else {
      mat_full_fp32_t* matrix_result = NULL;
      rc = one_off_matrix_T<float,mat_full_fp32_t>(table,tree,unifrac_method,variance_adjust,alpha,bypass_tips,normalize_sample_counts,1,NULL,&matrix_result);
      if (rc==okay) {
         *result = matrix_result->matrix[1];
	 destroy_mat_full_fp32(&matrix_result);
      }
    }
    return rc;
}

compute_status one_dense_pair_v3t(unsigned int n_obs, const char ** obs_ids, const double* sample1, const double* sample2,
		                  const opaque_bptree_t* tree_data,
                                  const char* unifrac_method, bool variance_adjust, double alpha,
                                  bool bypass_tips, bool normalize_sample_counts, double* result) {
    if (tree_data==NULL) return tree_missing;
    const su::BPTree &tree = *( (const su::BPTree*) tree_data);
    return one_dense_pair(n_obs,obs_ids,sample1,sample2,
                          tree,
                          unifrac_method,variance_adjust,alpha,bypass_tips,normalize_sample_counts,result);
}

compute_status one_dense_pair_v3(unsigned int n_obs, const char ** obs_ids, const double* sample1, const double* sample2,
		                 const support_bptree_t* tree_data,
                                 const char* unifrac_method, bool variance_adjust, double alpha,
                                 bool bypass_tips, bool normalize_sample_counts, double* result) {
    if (tree_data==NULL) return tree_missing;
    su::BPTree tree(tree_data->structure,
                    tree_data->lengths,
                    tree_data->names,
                    tree_data->n_parens);

    return one_dense_pair(n_obs,obs_ids,sample1,sample2,
                          tree,
                          unifrac_method,variance_adjust,alpha,bypass_tips,normalize_sample_counts,result);
}

// Internal 
inline std::vector<std::string> stringlist_to_vector(const char *stringlist) {
  char *str = strdup(stringlist);
  unsigned int n_els = 1;
  for (const char *p = str; p[0]!=0; ++p) {
    if (p[0]==',') n_els++;
  }

  std::vector<std::string> out(n_els);
  char *prev=str;
  n_els = 0;
  for (char *p = str; p[0]!=0; ++p) {
    if (p[0]==',') {
       p[0] = 0;
       out[n_els] = prev;
       n_els++;
       prev = p+1;
    }
  }
  out[n_els] = prev;

  free(str);
  return out;  
}


#ifndef UNIFRAC_WASM
// Internal: Make sure TReal and real_id match
template<class TReal, class TMat>
inline compute_status compute_permanova_T(const char *grouping_filename, unsigned int n_columns, const char* const* columns,
                                          TMat * result, unsigned int permanova_perms,
                                          TReal *fstats, TReal *pvalues, uint32_t *n_groups) {
     const uint32_t n_samples = result->n_samples;
     uint32_t *grouping = new uint32_t[n_samples];

     indexed_tsv tsv_obj(grouping_filename, n_samples, result->sample_ids);

     // compute each column separately
     for (unsigned int i=0; i<n_columns; i++) {
       try {
         tsv_obj.read_grouping(columns[i], grouping, n_groups[i]);
       } catch(...) {
         delete[] grouping;
         return grouping_missing;
       }

       // seed < 0: the file-based permanova path predates the per-call seed
       // and is documented as concurrency-unsafe in README.md
       su::permanova(result->matrix, n_samples,
                     grouping, permanova_perms,
                     fstats[i], pvalues[i], /*seed*/ -1);
     }
     delete[] grouping;

     return okay;
}

compute_status compute_permanova_fp64(const char *grouping_filename, unsigned int n_columns, const char* *columns,
                                      mat_full_fp64_t * result, unsigned int permanova_perms,
                                      double *fstats, double *pvalues, unsigned int *n_groups) {
  return compute_permanova_T<double,mat_full_fp64_t>(grouping_filename,n_columns,columns,result,permanova_perms,fstats,pvalues,n_groups);
}

compute_status compute_permanova_fp32(const char *grouping_filename, unsigned int n_columns, const char* *columns,
                                      mat_full_fp32_t * result, unsigned int permanova_perms,
                                      float *fstats, float *pvalues, unsigned int *n_groups) {
  return compute_permanova_T<float,mat_full_fp32_t>(grouping_filename,n_columns,columns,result,permanova_perms,fstats,pvalues,n_groups);
}

compute_status unifrac_to_file_v3(const char* biom_filename, const char* tree_filename, const char* out_filename,
                                  const char* unifrac_method, bool variance_adjust, double alpha,
                                  bool bypass_tips, bool normalize_sample_counts, unsigned int n_substeps, const char* format,
                                  unsigned int subsample_depth, bool subsample_with_replacement,
                                  unsigned int pcoa_dims,
                                  unsigned int permanova_perms, const char *grouping_filename, const char *grouping_columns,
                                  const char *mmap_dir)
{
    SETUP_TDBG("unifrac_to_file")

    bool fp64;
    bool save_dist;
    compute_status rc = is_fp64(unifrac_method, format, fp64, save_dist);

    if (rc==okay) {
      if (fp64) {
        mat_full_fp64_t* result = NULL;
        rc = one_off_matrix_v3(biom_filename, tree_filename,
                               unifrac_method, variance_adjust, alpha,
                               bypass_tips, normalize_sample_counts, n_substeps, subsample_depth, subsample_with_replacement, mmap_dir,
                               &result);
        TDBG_STEP("matrix_fp64 computed")

        if (rc==okay) {
          if (permanova_perms>0) {
            typedef const char* Tcstring;

            const auto columns = stringlist_to_vector(grouping_columns);
            const unsigned int n_columns = columns.size();
            Tcstring *columns_c = new Tcstring[n_columns];
            for (unsigned int i=0; i<n_columns; i++)  columns_c[i] = columns[i].c_str();
            double *fstats = new double[n_columns];
            double *pvalues = new double[n_columns];
            uint32_t *n_groups = new uint32_t[n_columns];

            rc = compute_permanova_fp64(grouping_filename,n_columns,columns_c,result,permanova_perms,fstats,pvalues,n_groups);
            TDBG_STEP("permanova_fp64 computed")

            if (rc==okay) {
              std::string stat_method("PERMANOVA");
              Tcstring *stat_methods = new Tcstring[n_columns];
              for (unsigned int i=0; i<n_columns; i++)  stat_methods[i] = stat_method.c_str();
              std::string stat_name("pseudo-F");
              Tcstring *stat_names = new Tcstring[n_columns];
              for (unsigned int i=0; i<n_columns; i++)  stat_names[i] = stat_name.c_str();
              uint32_t *nperm_arr = new uint32_t[n_columns];
              for (unsigned int i=0; i<n_columns; i++)  nperm_arr[i] = permanova_perms;

              IOStatus iostatus = write_mat_from_matrix_hdf5_fp64_v2(out_filename, result, pcoa_dims, save_dist,
                                                                     n_columns, stat_methods, stat_names,
                                                                     fstats, pvalues, nperm_arr,
                                                                     columns_c, n_groups);
              TDBG_STEP("file saved")
              if (iostatus!=write_okay) rc=output_error;
              delete[] nperm_arr;
              delete[] stat_methods;
              delete[] stat_names;
            }
            delete[] n_groups;
            delete[] pvalues;
            delete[] fstats;
            delete[] columns_c;
          } else {
            IOStatus iostatus = write_mat_from_matrix_hdf5_fp64(out_filename, result, pcoa_dims, save_dist);
            TDBG_STEP("file saved")
            if (iostatus!=write_okay) rc=output_error;
          }
          destroy_mat_full_fp64(&result);
        }
      } else {
        mat_full_fp32_t* result = NULL;
        rc = one_off_matrix_fp32_v3(biom_filename, tree_filename,
                                    unifrac_method, variance_adjust, alpha,
                                    bypass_tips, normalize_sample_counts, n_substeps, subsample_depth, subsample_with_replacement, mmap_dir,
                                    &result);
        TDBG_STEP("matrix_fp32 computed")
     
        if (rc==okay) {
          if (permanova_perms>0) {
            typedef const char* Tcstring;

            const auto columns = stringlist_to_vector(grouping_columns);
            const unsigned int n_columns = columns.size();
            Tcstring *columns_c = new Tcstring[n_columns];
            for (unsigned int i=0; i<n_columns; i++)  columns_c[i] = columns[i].c_str();
            float *fstats = new float[n_columns];
            float *pvalues = new float[n_columns];
            uint32_t *n_groups = new uint32_t[n_columns];

            rc = compute_permanova_fp32(grouping_filename,n_columns,columns_c,result,permanova_perms,fstats,pvalues,n_groups);
            TDBG_STEP("permanova_fp32 computed")

            if (rc==okay) {
              std::string stat_method("PERMANOVA");
              Tcstring *stat_methods = new Tcstring[n_columns];
              for (unsigned int i=0; i<n_columns; i++)  stat_methods[i] = stat_method.c_str();
              std::string stat_name("pseudo-F");
              Tcstring *stat_names = new Tcstring[n_columns];
              for (unsigned int i=0; i<n_columns; i++)  stat_names[i] = stat_name.c_str();
              uint32_t *nperm_arr = new uint32_t[n_columns];
              for (unsigned int i=0; i<n_columns; i++)  nperm_arr[i] = permanova_perms;

              IOStatus iostatus = write_mat_from_matrix_hdf5_fp32_v2(out_filename, result, pcoa_dims, save_dist,
                                                                     n_columns, stat_methods, stat_names,
                                                                     fstats, pvalues, nperm_arr,
                                                                     columns_c, n_groups);
              TDBG_STEP("file saved")
              if (iostatus!=write_okay) rc=output_error;

              delete[] nperm_arr;
              delete[] stat_methods;
              delete[] stat_names;
            }
            delete[] n_groups;
            delete[] pvalues;
            delete[] fstats;
            delete[] columns_c;
          } else {
            IOStatus iostatus = write_mat_from_matrix_hdf5_fp32(out_filename, result, pcoa_dims, save_dist);
            TDBG_STEP("file saved")
            if (iostatus!=write_okay) rc=output_error;
          }
          destroy_mat_full_fp32(&result);
        }
      }
    }

    return rc;
}

compute_status unifrac_to_txt_file_v3(const char* biom_filename, const char* tree_filename, const char* out_filename,
                                      const char* unifrac_method, bool variance_adjust, double alpha,
                                      bool bypass_tips, bool normalize_sample_counts, unsigned int n_substeps,
                                      const char *mmap_dir)
{
    SETUP_TDBG("unifrac_to_txt_file")

    bool fp64;
    compute_status rc = is_fp64_method(unifrac_method, fp64);

    if (rc==okay) {
      if (fp64) {
        mat_full_fp64_t* result = NULL;
        rc = one_off_matrix_v3(biom_filename, tree_filename,
                               unifrac_method, variance_adjust, alpha,
                               bypass_tips, normalize_sample_counts, n_substeps,
                               0, true, // not used, just pass disable flags
                               mmap_dir,
                               &result);
        TDBG_STEP("matrix_fp64 computed")

        if (rc==okay) {
          IOStatus iostatus = write_mat_from_matrix(out_filename, result);
          if (iostatus!=write_okay) rc=output_error;
          destroy_mat_full_fp64(&result);
        }
      } else {
        mat_full_fp32_t* result = NULL;
        rc = one_off_matrix_fp32_v3(biom_filename, tree_filename,
                                    unifrac_method, variance_adjust, alpha,
                                    bypass_tips, normalize_sample_counts, n_substeps,
                                    0, true, // not used, just pass disable flags
                                    mmap_dir,
                                    &result);
        TDBG_STEP("matrix_fp32 computed")
     
        if (rc==okay) {
          IOStatus iostatus = write_mat_from_matrix_fp32(out_filename, result);
          if (iostatus!=write_okay) rc=output_error;
          destroy_mat_full_fp32(&result);
        }
      }
    }

    TDBG_STEP("finished")
    return rc;
}

herr_t write_hdf5_string(hid_t output_file_id,const char *dname, const char *str)
{
  // this is the convoluted way to store a string
  // Will use the FORTRAN forma, so we do not depend on null termination
  hid_t filetype_id = H5Tcopy (H5T_FORTRAN_S1);
  H5Tset_size(filetype_id, strlen(str));
  hid_t memtype_id = H5Tcopy (H5T_C_S1);
  H5Tset_size(memtype_id, strlen(str)+1);

  hsize_t  dims[1] = {1};
  hid_t dataspace_id = H5Screate_simple (1, dims, NULL);

  hid_t dataset_id = H5Dcreate(output_file_id,dname, filetype_id, dataspace_id, H5P_DEFAULT, H5P_DEFAULT,
                                H5P_DEFAULT);
  herr_t status = H5Dwrite(dataset_id, memtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, str);

  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);
  H5Tclose(memtype_id);
  H5Tclose(filetype_id);

  return status;
}

// Internal
inline herr_t write_hdf5_stringarray(hid_t output_file_id,
                                     const char *label,
                                     hsize_t n_els, const char * const *els) {
  hsize_t   dims[1];
  dims[0] = n_els;
  hid_t dataspace_id = H5Screate_simple(1, dims, NULL);

  // this is the convoluted way to store an array of strings
  hid_t datatype_id = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype_id,H5T_VARIABLE);

  hid_t dcpl_id = H5Pcreate (H5P_DATASET_CREATE);

  hid_t dataset_id = H5Dcreate1(output_file_id, label , datatype_id, dataspace_id, dcpl_id);

  herr_t status = H5Dwrite(dataset_id, datatype_id, H5S_ALL, H5S_ALL,
                           H5P_DEFAULT, els);

  H5Dclose(dataset_id);
  H5Tclose(datatype_id);
  H5Sclose(dataspace_id);
  H5Pclose(dcpl_id);

  return status;
}

// Internal: Make sure TReal and real_id match
template<class TReal>
inline herr_t write_hdf5_array(hid_t output_file_id, hid_t real_id,
                               const char *label,
                               hsize_t n_els, const TReal *els) {
  hsize_t   dims[1];
  dims[0] = n_els;
  hid_t dataspace_id = H5Screate_simple(1, dims, NULL);

  hid_t dcpl_id = H5Pcreate (H5P_DATASET_CREATE);

  hid_t dataset_id = H5Dcreate2(output_file_id, label, real_id, dataspace_id,
                                H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  herr_t status = H5Dwrite(dataset_id, real_id, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           els);

  H5Pclose(dcpl_id);
  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);

  return status;
}

// Internal: Make sure TReal and real_id match
template<class TReal>
inline herr_t write_hdf5_array2D(hid_t output_file_id, hid_t real_id,
                                 const char *label,
                                 hsize_t dim1, hsize_t dim2, const TReal *els) {
  hsize_t   dims[2];
  dims[0] = dim1;
  dims[1] = dim2;
  hid_t dataspace_id = H5Screate_simple(2, dims, NULL);

  hid_t dcpl_id = H5Pcreate (H5P_DATASET_CREATE);

  hid_t dataset_id = H5Dcreate2(output_file_id, label ,real_id, dataspace_id,
                                H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
  herr_t status = H5Dwrite(dataset_id, real_id, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           els);

  H5Pclose(dcpl_id);
  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);

  return status;
}

template<class TReal>
inline IOStatus append_hdf5_pcoa(hid_t output_file_id, hid_t real_id, unsigned int pcoa_dims, unsigned int n_samples,
                                 const char *eigenvalues_key, const char * samples_key, const char * proportion_explained_key,
                                 const TReal * eigenvalues, const TReal * samples,  const TReal * proportion_explained) {
     // save the eigenvalues
     {
       herr_t status = write_hdf5_array<TReal>(output_file_id,real_id,
                                               eigenvalues_key, pcoa_dims, eigenvalues);
       if (status<0)  return write_error;
     }

     // save the proportion_explained
     {
       herr_t status = write_hdf5_array<TReal>(output_file_id,real_id,
                                               proportion_explained_key, pcoa_dims, proportion_explained);
       if (status<0)  return write_error;
     }

     // save the samples
     {
       herr_t status = write_hdf5_array2D<TReal>(output_file_id,real_id,
                                                 samples_key, n_samples, pcoa_dims, samples);
       if (status<0)  return write_error;
     }

  return write_okay;
}

namespace su {

template<class TReal, class TMat>
class WriteHDF5Multi {
protected:
   hid_t real_id;
   std::string fname;
   unsigned int pcoa_dims;
   hid_t output_file_id;
   unsigned int n_results;
public:
   WriteHDF5Multi(hid_t _real_id, const char* output_filename, unsigned int _pcoa_dims) 
   : real_id(_real_id)
   , fname(output_filename)
   , pcoa_dims(_pcoa_dims)
   , output_file_id(H5Fcreate(output_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT))
   , n_results(0)
   {
      if (output_file_id<0) throw "File open failed";

      // simple header
      if (write_hdf5_string(output_file_id,"format","SUJK")<0) {
         H5Fclose (output_file_id);
         throw "Header write failed";
      }
      if (write_hdf5_string(output_file_id,"version","2023.03")<0) {
         H5Fclose (output_file_id);
         throw "Header write failed";
      }
      if (pcoa_dims>0) {
        if (write_hdf5_string(output_file_id,"pcoa_method","FSVD")<0) {
           H5Fclose (output_file_id);
           throw "Header write failed";
        }
      }
   }

   virtual ~WriteHDF5Multi()
   {
      H5Fclose (output_file_id);
   }

   // Note: May destroy the content of result->matrix
   void add_result(TMat * result, bool save_dist) {
      SETUP_TDBG("WriteHDF5Multi_add_result")

      char fmtstr[64];
      const auto n_samples = result->n_samples;

      if (n_results==0) {
         // save the ids once, they do not change
         {
           herr_t status = write_hdf5_stringarray(output_file_id,
                               "order", n_samples, result->sample_ids);
           if (status<0) throw "Order write failed";
         }
      }

      if (save_dist) {
         // save the matrix
         {
           snprintf(fmtstr,63,"matrix:%d",n_results);
           herr_t status = write_hdf5_array2D<TReal>(output_file_id,real_id,
                            fmtstr, n_samples, n_samples, result->matrix);
           if (status<0) throw "Matrix write failed";
         }
         TDBG_STEP("matrix saved")
      }

      if (pcoa_dims>0) {
         // compute pcoa and save it in the file
         // use inplace variant to keep memory use in check; we don't need matrix anymore
         TReal * eigenvalues;
         TReal * samples;
         TReal * proportion_explained;

         // seed < 0: this path predates the per-call seed and is documented
         // as concurrency-unsafe in README.md
         su::pcoa_inplace(result->matrix, n_samples, pcoa_dims, eigenvalues, samples, proportion_explained, /*seed*/ -1);
         TDBG_STEP("pcoa computed")

         char fmtstr2[64];
         char fmtstr3[64];
         snprintf(fmtstr,63,"pcoa_eigvals:%d",n_results);
         snprintf(fmtstr2,63,"pcoa_samples:%d",n_results);
         snprintf(fmtstr3,63,"pcoa_proportion_explained:%d",n_results);
         IOStatus rc = append_hdf5_pcoa(output_file_id, real_id, pcoa_dims, n_samples,
                             fmtstr, fmtstr2, fmtstr3,
                             eigenvalues, samples,  proportion_explained);
         free(eigenvalues);
         free(proportion_explained);
         free(samples);
         if (rc!=write_okay) throw "PCOA write failed";
         TDBG_STEP("pcoa saved")
      } // if pcoa

      n_results++;
   }

   void write_stats(unsigned int           stat_n_vals,
                    const char* const    * stat_method_arr, const char* const  * stat_name_arr,
                    const TReal          * stat_val_arr,    const TReal        * stat_pval_arr, const uint32_t  * stat_perm_count_arr,
                    const char* const    * stat_group_name_arr, const uint32_t * stat_group_count_arr) {
     SETUP_TDBG("WriteHDF5Multi_write_stats")
     herr_t status = write_hdf5_stringarray(output_file_id,
                         "stat_methods", stat_n_vals, stat_method_arr);
     if (status>=0) {
       status = write_hdf5_stringarray(output_file_id,
                         "stat_test_names", stat_n_vals, stat_name_arr);
     }
     if (status>=0) {
       status = write_hdf5_stringarray(output_file_id,
                         "stat_grouping_names", stat_n_vals, stat_group_name_arr);
     }
     if (status>=0) {
       status = write_hdf5_array<uint32_t>(output_file_id, H5T_STD_U32LE,
                         "stat_n_groups", stat_n_vals, stat_group_count_arr);
     }
     if (status>=0) {
       status = write_hdf5_array<TReal>(output_file_id,real_id,
                         "stat_values", stat_n_vals, stat_val_arr);
     }
     if (status>=0) {
       status = write_hdf5_array<TReal>(output_file_id,real_id,
                         "stat_pvalues", stat_n_vals, stat_pval_arr);
     }
     if (status>=0) {
       status = write_hdf5_array<uint32_t>(output_file_id, H5T_STD_U32LE,
                         "stat_n_permutations", stat_n_vals, stat_perm_count_arr);
     }
     
     // check status after cleanup, for simplicity
     if (status<0) throw "Stats write failed";
     TDBG_STEP("stats saved")
   }

};

class WriteHDF5MultiFP64 : public WriteHDF5Multi<double,mat_full_fp64_t> {
public:
   WriteHDF5MultiFP64(const char* output_filename, unsigned int _pcoa_dims)
   : WriteHDF5Multi(H5T_IEEE_F64LE, output_filename, _pcoa_dims) {}
};

class WriteHDF5MultiFP32 : public WriteHDF5Multi<float,mat_full_fp32_t> {
public:
   WriteHDF5MultiFP32(const char* output_filename, unsigned int _pcoa_dims)
   : WriteHDF5Multi(H5T_IEEE_F32LE, output_filename, _pcoa_dims) {}
};

} // end namespace

template<class TReal, class TMat>
compute_status unifrac_multi_to_file_T(hid_t real_id, const bool save_dist,
                                        const char* biom_filename, const char* tree_filename, const char* out_filename,
                                        const char* unifrac_method, bool variance_adjust, double alpha,
                                        bool bypass_tips, bool normalize_sample_counts, unsigned int nsubsteps, const char* format,
                                        unsigned int n_subsamples, unsigned int subsample_depth, bool subsample_with_replacement,
                                        unsigned int pcoa_dims,
                                        unsigned int permanova_perms, const char *grouping_filename, const char *grouping_columns,
                                        const char *mmap_dir)
{
    compute_status rc = okay;
    SETUP_TDBG("unifrac_multi_to_file")

    if (subsample_depth<1) {
      fprintf(stderr, "ERROR: subsampling depth cannot be 0.\n");
      return table_empty;
    }

    CHECK_FILE(biom_filename, table_missing)
    CHECK_FILE(tree_filename, tree_missing)
    PARSE_TREE_TABLE(tree_filename, biom_filename)
    TDBG_STEP("load_files")

    typedef const char* Tcstring;
    std::vector<std::string> columns;
    Tcstring *columns_c = NULL;

    // permanova values
    TReal *pm_fstats = NULL;
    TReal *pm_pvalues = NULL;
    uint32_t *pm_n_groups = NULL;
    std::vector<std::string> pm_columns;
    Tcstring *pm_columns_c = NULL;

    if (permanova_perms>0) {
         columns = stringlist_to_vector(grouping_columns);
         const unsigned int n_columns = columns.size();
         const unsigned int n_els = n_columns*n_subsamples;

         columns_c = new Tcstring[n_columns];
         for (unsigned int i=0; i<n_columns; i++)  columns_c[i] = columns[i].c_str();

         pm_fstats = new TReal[n_els];
         pm_pvalues = new TReal[n_els];
         pm_n_groups = new uint32_t[n_els];
         pm_columns.resize(n_els);
         pm_columns_c = new Tcstring[n_els];
    }

    try {
      su::WriteHDF5Multi<TReal,TMat> h5obj(real_id, out_filename,pcoa_dims);

      for (unsigned int i=0; i<n_subsamples; i++) {
        su::skbio_biom_subsampled table_subsampled(table, subsample_with_replacement, subsample_depth);
        if ((table_subsampled.n_samples==0) || (table_subsampled.n_obs==0)) {
           rc = table_empty;
           break;
        }
        TDBG_STEP("subsampled")

        TMat* result = NULL;
        rc = one_off_matrix_T<TReal,TMat>(table_subsampled,tree,unifrac_method,variance_adjust,alpha,bypass_tips,normalize_sample_counts,nsubsteps,mmap_dir,&result);
        if (rc!=okay) break;
        TDBG_STEP("matrix computed")

        if (permanova_perms>0) {
            const unsigned int n_columns = columns.size();
            TReal *fstats = pm_fstats+(i*n_columns);
            TReal *pvalues = pm_pvalues+(i*n_columns);
            uint32_t *n_groups = pm_n_groups+(i*n_columns);

            rc = compute_permanova_T<TReal,TMat>(grouping_filename,n_columns,columns_c,result,permanova_perms,fstats,pvalues,n_groups);
            if (rc!=okay) break;
            TDBG_STEP("permanova computed")

            char fmtstr[32];
            snprintf(fmtstr,31,":%d",i);
            for (unsigned int j=0; j<n_columns;j ++) {
              const unsigned int idx = i*n_columns + j;

              pm_columns[idx] = columns[j] + fmtstr;
              pm_columns_c[idx] = pm_columns[idx].c_str();
            }
        }
        h5obj.add_result(result, save_dist);
        destroy_mat_full_T<TMat,TReal>(&result);
      } // for i
      if ((rc==okay)&&(permanova_perms>0)) {
              const unsigned int n_columns = columns.size();
              const unsigned int n_els = n_columns*n_subsamples;
              std::string stat_method("PERMANOVA");
              Tcstring *stat_methods = new Tcstring[n_els];
              for (unsigned int i=0; i<n_els; i++)  stat_methods[i] = stat_method.c_str();
              std::string stat_name("pseudo-F");
              Tcstring *stat_names = new Tcstring[n_els];
              for (unsigned int i=0; i<n_els; i++)  stat_names[i] = stat_name.c_str();
              uint32_t *nperm_arr = new uint32_t[n_els];
              for (unsigned int i=0; i<n_els; i++)  nperm_arr[i] = permanova_perms;

              h5obj.write_stats(n_els, stat_methods, stat_names,
                                pm_fstats, pm_pvalues, nperm_arr,
                                pm_columns_c, pm_n_groups);

              delete[] nperm_arr;
              delete[] stat_names;
              delete[] stat_methods;
      }
    } catch (...) {
       // the only one throwing should be h5obj
       rc = output_error;
    }

    TDBG_STEP("finished")

    if (pm_columns_c!=NULL) delete[] pm_columns_c;
    if (pm_n_groups!=NULL)  delete[] pm_n_groups;
    if (pm_pvalues!=NULL)   delete[] pm_pvalues;
    if (pm_fstats!=NULL)    delete[] pm_fstats;
    if (columns_c!=NULL)    delete[] columns_c;
    return rc;
}

compute_status unifrac_multi_to_file_v3(const char* biom_filename, const char* tree_filename, const char* out_filename,
                                        const char* unifrac_method, bool variance_adjust, double alpha,
                                        bool bypass_tips, bool normalize_sample_counts, unsigned int nsubsteps, const char* format,
                                        unsigned int n_subsamples, unsigned int subsample_depth, bool subsample_with_replacement,
                                        unsigned int pcoa_dims,
                                        unsigned int permanova_perms, const char *grouping_filename, const char *grouping_columns,
                                        const char *mmap_dir)
{
    bool fp64;
    bool save_dist;
    compute_status rc = is_fp64(unifrac_method, format, fp64, save_dist);

    if (rc!=okay) {
      return rc;
    }

    if (fp64) {
      rc = unifrac_multi_to_file_T<double,mat_full_fp64_t>(H5T_IEEE_F64LE, save_dist,
                                     biom_filename, tree_filename, out_filename,
                                     unifrac_method, variance_adjust, alpha,
                                     bypass_tips, normalize_sample_counts, nsubsteps, format,
                                     n_subsamples, subsample_depth, subsample_with_replacement,
                                     pcoa_dims, permanova_perms, grouping_filename, grouping_columns,
                                     mmap_dir);
   } else {
      rc = unifrac_multi_to_file_T<float,mat_full_fp32_t>(H5T_IEEE_F32LE, save_dist,
                                     biom_filename, tree_filename, out_filename,
                                     unifrac_method, variance_adjust, alpha,
                                     bypass_tips, normalize_sample_counts, nsubsteps, format,
                                     n_subsamples, subsample_depth, subsample_with_replacement,
                                     pcoa_dims, permanova_perms, grouping_filename, grouping_columns,
                                     mmap_dir);
   }

   return rc;
}
#endif // UNIFRAC_WASM (compute_permanova / unifrac_to_file / HDF5 helpers / unifrac_multi)

IOStatus write_mat(const char* output_filename, mat_t* result) {
    std::ofstream output;
    SETUP_TDBG("write_mat")

    output.open(output_filename);

    uint64_t comb_N = su::comb_2(result->n_samples);
    uint64_t comb_N_minus = 0;
    double v;

    for(unsigned int i = 0; i < result->n_samples; i++)
        output << "\t" << result->sample_ids[i];
    output << std::endl;
    TDBG_STEP("header saved")

    for(unsigned int i = 0; i < result->n_samples; i++) {
        output << result->sample_ids[i];
        for(unsigned int j = 0; j < result->n_samples; j++) {
            if(i < j) { // upper triangle
                comb_N_minus = su::comb_2(result->n_samples - i);
                v = result->condensed_form[comb_N - comb_N_minus + (j - i - 1)];
            } else if (i > j) { // lower triangle
                comb_N_minus = su::comb_2(result->n_samples - j);
                v = result->condensed_form[comb_N - comb_N_minus + (i - j - 1)];
            } else {
                v = 0.0;
            }
            output << std::setprecision(16) << "\t" << v;
        }
        output << "\n";
    }
    output.close();
    TDBG_STEP("file saved")

    return write_okay;
}

/*
 * Note: std::ofstream is too slow for large files
 *       Use lower-level routines and line caching for speed
 */
template<typename TMat>
inline IOStatus write_mat_from_matrix_txt_T(const char* filename, TMat* result) {
    SETUP_TDBG("write_mat_from_matrix_txt")

    IOStatus outrc = write_okay;

    FILE *output = fopen(filename,"w");
    if (output==NULL) return write_error;

    const unsigned int n_samples = result->n_samples;

    size_t max_id_len = 0;
    size_t head_len = 24*n_samples; // don't expect ids to be more than 20 chars(+tab) in len on average
    char *head_str = (char*) malloc(head_len);
    {
      // create header string in memory first
      char *str_tail_p = head_str;
      for(unsigned int i = 0; i < n_samples; i++) {
        const char* sample_id = result->sample_ids[i];
        size_t id_len =strlen(sample_id);
        {
          const size_t tail_offset = str_tail_p-head_str;
          if ((tail_offset+id_len+4)>head_len) {
            // not enough space, double the buffer
            head_len *=2;
            head_str = (char*) realloc(head_str,head_len);
            str_tail_p = head_str + tail_offset; // everything moves
          }
        }
        str_tail_p[0] = '\t';
        str_tail_p++;
        memcpy(str_tail_p,sample_id,id_len);
        str_tail_p+=id_len;
        if (id_len>max_id_len) max_id_len = id_len;
      }
      str_tail_p[0] = '\n';
      str_tail_p++;
      str_tail_p[0] = 0; // just to be paranoid 
      // now write to file
      if (fwrite(head_str, str_tail_p-head_str,1,output)!=1) return write_error;
    }
    TDBG_STEP("header saved")

    {
#if defined(_OPENMP)
      // allow all threads to create their lines in parallel
      const unsigned int max_threads = omp_get_max_threads();
      const unsigned int use_threads = std::min(8u,max_threads); // no benefit beyond 8 threads, since we have the IO serialization
#else
      const unsigned int max_threads = 1;
#endif
      // max number of chars per element
      constexpr size_t fp_precision = std::is_same<TMat,mat_full_fp32_t>::value ? 12 : 20;

      // pre-allocate a line buffers, so we don't get constant re-allocations
      std::vector<char*>  line_str_vect;  // actual buffers
      std::vector<size_t> line_len_vect;  // size of the buffer
      line_str_vect.resize(max_threads);
      line_len_vect.resize(max_threads);
      line_str_vect[0] = head_str; // reuse the header buffer
      line_len_vect[0] = head_len;
      for(unsigned int t=1; t<max_threads; t++) {
        // all others lazy initialize in loop
        line_str_vect[t] = NULL;
        line_len_vect[t] = 0;
      }

#if defined(_OPENMP)
#pragma omp parallel for ordered schedule(static,1) num_threads(use_threads)
#endif
      for(unsigned int i = 0; i < n_samples; i++) {
#if defined(_OPENMP)
        const unsigned int my_thread = omp_get_thread_num();
#else
        const unsigned int my_thread = 0;
#endif
        auto &line_str = line_str_vect[my_thread];
        auto &line_len = line_len_vect[my_thread];

        if (line_str==NULL) {
          // first time initialization
          line_len = max_id_len+(fp_precision+4)*n_samples;
          line_str = (char*) malloc(line_len);
        }

        // create line string in memory first
        char *str_tail_p = line_str;
        {  // we know we have enough buffer for at least the header
          const char* sample_id = result->sample_ids[i];
          size_t id_len =strlen(sample_id);
          memcpy(str_tail_p,sample_id,id_len);
          str_tail_p+=id_len;
        }
        for(unsigned int j = 0; j < n_samples; j++) {
          {
            const size_t tail_offset = str_tail_p-line_str;
            if ((tail_offset+fp_precision+4)>line_len) {
              // not enough space, double the buffer
              line_len *=2;
              line_str = (char*) realloc(line_str,line_len);
              str_tail_p = line_str + tail_offset; // everything moves
            }
          }
          str_tail_p[0] = '\t';
          str_tail_p++;
          std::to_chars_result tcrc = std::to_chars(str_tail_p, line_str+line_len, result->matrix[i*n_samples+j], std::chars_format::fixed);
          str_tail_p = tcrc.ptr;
        }
        str_tail_p[0] = '\n';
        str_tail_p++;
        str_tail_p[0] = 0; // just to be paranoid 
#if defined(_OPENMP)
#pragma omp ordered
        // compute lines independently, but then write them sequentially
#endif
        {
          if (fwrite(line_str, str_tail_p-line_str,1,output)!=1) outrc = write_error; // cannot exit from OMP, will check outside
        }
      }
      // free all buffers, inlcudeing head_str (==line_str_vect[0])
      for(unsigned int t=0; t<max_threads; t++) {
	if (line_str_vect[t]!=NULL) free(line_str_vect[t]);
      }
    }

    fclose(output);
    if (outrc!=write_okay) return outrc;


    TDBG_STEP("file saved")

    return write_okay;
}

IOStatus write_mat_from_matrix(const char* filename, mat_full_fp64_t* result) {
    return write_mat_from_matrix_txt_T(filename, result);
}

IOStatus write_mat_from_matrix_fp32(const char* filename, mat_full_fp32_t* result) {
    return write_mat_from_matrix_txt_T(filename, result);
}

#ifndef UNIFRAC_WASM
// Internal: Make sure TReal and real_id match
template<class TReal, class TMat>
inline IOStatus write_mat_from_matrix_hdf5_T(const char* output_filename, TMat * result, hid_t real_id,
                                             unsigned int pcoa_dims, bool save_dist,
                                             unsigned int           stat_n_vals,
                                             const char* const    * stat_method_arr, const char* const  * stat_name_arr,
                                             const TReal          * stat_val_arr,    const TReal        * stat_pval_arr, const uint32_t  * stat_perm_count_arr,
                                             const char* const    * stat_group_name_arr, const uint32_t * stat_group_count_arr) {
   SETUP_TDBG("write_mat_from_matrix")
   /* Create a new file using default properties. */
   hid_t output_file_id = H5Fcreate(output_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
   if (output_file_id<0) return write_error;

   const auto n_samples = result->n_samples;

   // simple header
   if (write_hdf5_string(output_file_id,"format","BDSM")<0) {
       H5Fclose (output_file_id);
       return write_error;
   }
   if (write_hdf5_string(output_file_id,"version","2020.12")<0) {
       H5Fclose (output_file_id);
       return write_error;
   }

   // save the ids
   {
     herr_t status = write_hdf5_stringarray(output_file_id,
                         "order", n_samples, result->sample_ids);
     if (status<0) {
       H5Fclose (output_file_id);
       return write_error;
     }
   }
   TDBG_STEP("header saved")

   // save the matrix
   if (save_dist) {
     herr_t status = write_hdf5_array2D<TReal>(output_file_id,real_id,
                         "matrix", n_samples, n_samples, result->matrix);
     if (status<0) {
       H5Fclose (output_file_id);
       return write_error;
     }
     TDBG_STEP("matrix saved")
   }

   if (pcoa_dims>0) {
     // compute pcoa and save it in the file
     // use inplace variant to keep memory use in check; we don't need matrix anymore
     TReal * eigenvalues;
     TReal * samples;
     TReal * proportion_explained;

     // see the note on the other pcoa_inplace call site
     su::pcoa_inplace(result->matrix, n_samples, pcoa_dims, eigenvalues, samples, proportion_explained, /*seed*/ -1);
     TDBG_STEP("pcoa computed")


     IOStatus rc = write_okay;
     if (write_hdf5_string(output_file_id,"pcoa_method","FSVD")<0) {
       rc = write_error;
     } else {
       rc = append_hdf5_pcoa(output_file_id, real_id, pcoa_dims, n_samples,
                             "pcoa_eigvals", "pcoa_samples", "pcoa_proportion_explained",
                             eigenvalues, samples,  proportion_explained);
     }
     free(eigenvalues);
     free(proportion_explained);
     free(samples);
     if (rc!=write_okay) {
       H5Fclose (output_file_id);
       return rc;
     }
     TDBG_STEP("pcoa saved")
   } // if pcoa

   if (stat_n_vals>0) {
     herr_t status = write_hdf5_stringarray(output_file_id,
                         "stat_methods", stat_n_vals, stat_method_arr);
     if (status>=0) {
       status = write_hdf5_stringarray(output_file_id,
                         "stat_test_names", stat_n_vals, stat_name_arr);
     }
     if (status>=0) {
       status = write_hdf5_stringarray(output_file_id,
                         "stat_grouping_names", stat_n_vals, stat_group_name_arr);
     }
     if (status>=0) {
       status = write_hdf5_array<uint32_t>(output_file_id, H5T_STD_U32LE,
                         "stat_n_groups", stat_n_vals, stat_group_count_arr);
     }
     if (status>=0) {
       status = write_hdf5_array<TReal>(output_file_id,real_id,
                         "stat_values", stat_n_vals, stat_val_arr);
     }
     if (status>=0) {
       status = write_hdf5_array<TReal>(output_file_id,real_id,
                         "stat_pvalues", stat_n_vals, stat_pval_arr);
     }
     if (status>=0) {
       status = write_hdf5_array<uint32_t>(output_file_id, H5T_STD_U32LE,
                         "stat_n_permutations", stat_n_vals, stat_perm_count_arr);
     }
     
     // check status after cleanup, for simplicity
     if (status<0) {
       H5Fclose (output_file_id);
       return write_error;
     }
     TDBG_STEP("stats saved")
   } // if stat_n_vals

   H5Fclose (output_file_id);
   return write_okay;
}

// Internal: Make sure TReal and real_id match
// Note: Deprecated, for backwards compatibility only
template<class TReal, class TMat>
inline IOStatus write_mat_hdf5_T(const char* output_filename, mat_t* result,hid_t real_id, unsigned int pcoa_dims, bool save_dist) {
     // compute the matrix
     TMat mat_full;
     mat_full.n_samples = result->n_samples;

     const uint64_t n_samples = result->n_samples;
     mat_full.flags = 0;
     mat_full.matrix = (TReal*) malloc(n_samples*n_samples*sizeof(TReal));
     if (mat_full.matrix==NULL) {
       return open_error; // we don't have a better error code
     }

     mat_full.sample_ids = result->sample_ids; // just link

     condensed_form_to_matrix_T(result->condensed_form, n_samples, mat_full.matrix);
     IOStatus err =  write_mat_from_matrix_hdf5_T<TReal,TMat>(output_filename, &mat_full, real_id, pcoa_dims, save_dist,
                             0,NULL,NULL,NULL,NULL,NULL,NULL,NULL);

     free(mat_full.matrix);
     return err;
}

// Note: Deprecated, for backwards compatibility only
IOStatus write_mat_hdf5_fp64(const char* output_filename, mat_t* result, unsigned int pcoa_dims, int save_dist) {
  return write_mat_hdf5_T<double,mat_full_fp64_t>(output_filename,result,H5T_IEEE_F64LE,pcoa_dims,save_dist);
}

// Note: Deprecated, for backwards compatibility only
IOStatus write_mat_hdf5_fp32(const char* output_filename, mat_t* result, unsigned int pcoa_dims, int save_dist) {
  return write_mat_hdf5_T<float,mat_full_fp32_t>(output_filename,result,H5T_IEEE_F32LE,pcoa_dims,save_dist);
}

IOStatus write_mat_from_matrix_hdf5_fp64_v2(const char* output_filename, mat_full_fp64_t* result,
                                            unsigned int pcoa_dims, int save_dist,
                                            unsigned int stat_n_vals,
                                            const char*  *stat_method_arr,     const char*        *stat_name_arr,
                                            const double *stat_val_arr,        const double       *stat_pval_arr, const unsigned int *stat_perm_count_arr,
                                            const char*  *stat_group_name_arr, const unsigned int *stat_group_count_arr) {
  return write_mat_from_matrix_hdf5_T<double,mat_full_fp64_t>(output_filename,result,H5T_IEEE_F64LE,pcoa_dims,save_dist,
                        stat_n_vals,stat_method_arr,stat_name_arr,stat_val_arr,stat_pval_arr,stat_perm_count_arr,stat_group_name_arr,stat_group_count_arr);
}

IOStatus write_mat_from_matrix_hdf5_fp32_v2(const char* output_filename, mat_full_fp32_t* result,
                                            unsigned int pcoa_dims, int save_dist,
                                            unsigned int stat_n_vals,
                                            const char*  *stat_method_arr,     const char*        *stat_name_arr,
                                            const float  *stat_val_arr,        const float        *stat_pval_arr, const unsigned int *stat_perm_count_arr,
                                            const char*  *stat_group_name_arr, const unsigned int *stat_group_count_arr) {
  return write_mat_from_matrix_hdf5_T<float,mat_full_fp32_t>(output_filename,result,H5T_IEEE_F32LE,pcoa_dims,save_dist,
                        stat_n_vals,stat_method_arr,stat_name_arr,stat_val_arr,stat_pval_arr,stat_perm_count_arr,stat_group_name_arr,stat_group_count_arr);
}
#endif // UNIFRAC_WASM (write_mat_from_matrix_hdf5 family)

IOStatus write_vec(const char* output_filename, r_vec* result) {
    std::ofstream output;
    output.open(output_filename);

    // write sample ids in first column of file and faith's pd in second column
    output << "#SampleID\tfaith_pd" << std::endl;
    for(unsigned int i = 0; i < result->n_samples; i++) {
        output << result->sample_ids[i];
        output << std::setprecision(16) << "\t" << result->values[i];
        output << std::endl;
    }
    output.close();

    return write_okay;
}

#ifndef UNIFRAC_WASM
IOStatus write_partial(const char* output_filename, const partial_mat_t* result) {
    int fd = open(output_filename, O_WRONLY | O_CREAT | O_TRUNC,  S_IRUSR |  S_IWUSR );
    if (fd==-1) return write_error;

    int cnt = -1;

    uint32_t n_stripes = result->stripe_stop - result->stripe_start;

    uint32_t sample_id_length = 0;
    for(unsigned int i = 0; i < result->n_samples; i++) {
        sample_id_length += strlen(result->sample_ids[i])+1;
    }

    {
      char * const samples_buf = (char *)malloc(sample_id_length);
 
      char *samples_ptr = samples_buf;

      /* sample IDs */
      for(unsigned int i = 0; i < result->n_samples; i++) {
          uint32_t length = strlen(result->sample_ids[i])+1;
          memcpy(samples_ptr,result->sample_ids[i],length);
          samples_ptr+= length;
      }

      int max_compressed = LZ4_compressBound(sample_id_length);
      char * const cmp_buf = (char *)malloc(max_compressed);

      int sample_id_length_compressed = LZ4_compress_default(samples_buf,cmp_buf,sample_id_length,max_compressed);
      if (sample_id_length_compressed<1)  {close(fd); return open_error;}

      uint32_t header[8];
      header[0] = PARTIAL_MAGIC_V2;
      header[1] = result->n_samples;
      header[2] = n_stripes;
      header[3] = result->stripe_start;
      header[4] = result->stripe_total;
      header[5] = result->is_upper_triangle;
      header[6] = sample_id_length;
      header[7] = sample_id_length_compressed;

      cnt=write(fd,header, 8 * sizeof(uint32_t));
      if (cnt<1)  {close(fd); return write_error;}

      cnt=write(fd,cmp_buf, sample_id_length_compressed);
      if (cnt<1)  {close(fd); return write_error;}

      free(cmp_buf);
      free(samples_buf);
    }

    {
      int max_compressed = LZ4_compressBound(sizeof(double) * result->n_samples);
      char * const cmp_buf_raw = (char *)malloc(max_compressed+sizeof(uint32_t));
      char * const cmp_buf = cmp_buf_raw + sizeof(uint32_t);

      /* stripe information */
      for(unsigned int i = 0; i < n_stripes; i++) {
        int cmp_size = LZ4_compress_default((const char *) result->stripes[i],cmp_buf,sizeof(double) * result->n_samples,max_compressed);
        if (cmp_size<1)  {close(fd); return open_error;}

        uint32_t *cmp_buf_size_p = (uint32_t *)cmp_buf_raw;
        *cmp_buf_size_p = cmp_size;

        cnt=write(fd, cmp_buf_raw, cmp_size+sizeof(uint32_t));
        if (cnt<1) {return write_error;}
      }

      free(cmp_buf_raw);
    }

    /* footer */
    {
      uint32_t header[1];
      header[0] = PARTIAL_MAGIC_V2;

      cnt=write(fd,header, 1 * sizeof(uint32_t));
      if (cnt<1)  {close(fd); return open_error;}
    }

    close(fd);

    return write_okay;
}

IOStatus _is_partial_file(const char* input_filename) {
    int fd = open(input_filename, O_RDONLY );
    if (fd==-1) return open_error;

    uint32_t header[1];
    int cnt = read(fd,header,sizeof(uint32_t));
    close(fd);

    if (cnt!=sizeof(uint32_t)) return magic_incompatible;
    if ( header[0] != PARTIAL_MAGIC_V2) return magic_incompatible;

    return read_okay;
}

template<class TPMat>
inline IOStatus read_partial_header_fd(int fd, TPMat &result) {
    ssize_t cnt=-1;

    uint32_t header[8];
    cnt = read(fd,header,8*sizeof(uint32_t));
    if (cnt != (8*sizeof(uint32_t))) {return magic_incompatible;}

    if ( header[0] != PARTIAL_MAGIC_V2) {return magic_incompatible;}

    const uint32_t n_samples = header[1];
    const uint32_t n_stripes = header[2];
    const uint32_t stripe_start = header[3];
    const uint32_t stripe_total = header[4];
    const bool is_upper_triangle = header[5];

    /* sanity check header */
    if(n_samples <= 0 || n_stripes <= 0 || stripe_total <= 0 || (!is_upper_triangle))
         {return bad_header;}
    if(stripe_total >= n_samples || n_stripes > stripe_total || stripe_start >= stripe_total || stripe_start + n_stripes > stripe_total)
         {return bad_header;}

    /* initialize the partial result structure */
    result.n_samples = n_samples;
    result.sample_ids = (char**)malloc(sizeof(char*) * n_samples);
    result.stripes = (double**)malloc(sizeof(double*) * (n_stripes));
    result.stripe_start = stripe_start;
    result.stripe_stop = stripe_start + n_stripes;
    result.is_upper_triangle = is_upper_triangle;
    result.stripe_total = stripe_total;

    /* load samples */
    {
      const uint32_t sample_id_length = header[6];
      const uint32_t sample_id_length_compressed = header[7];

      /* sanity check header */
      if (sample_id_length<=0 || sample_id_length_compressed <=0)
         { return bad_header;}

      char * const cmp_buf = (char *)malloc(sample_id_length_compressed);
      if (cmp_buf==NULL) { return bad_header;} // no better error code
      cnt = read(fd,cmp_buf,sample_id_length_compressed);
      if (cnt != ssize_t(sample_id_length_compressed)) {free(cmp_buf); return magic_incompatible;}

      char *samples_buf = (char *)malloc(sample_id_length);
      if (samples_buf==NULL) { free(cmp_buf); return bad_header;} // no better error code

      cnt = LZ4_decompress_safe(cmp_buf,samples_buf,sample_id_length_compressed,sample_id_length);
      if (cnt!=ssize_t(sample_id_length)) {free(samples_buf); free(cmp_buf); return magic_incompatible;}

      const char *samples_ptr = samples_buf;

      for(uint32_t i = 0; i < n_samples; i++) {
        uint32_t sample_length = strlen(samples_ptr);
        if ((samples_ptr+sample_length+1)>(samples_buf+sample_id_length)) {free(samples_buf); free(cmp_buf); return magic_incompatible;}

        result.sample_ids[i] = (char*)malloc(sample_length + 1);
        memcpy(result.sample_ids[i],samples_ptr,sample_length + 1);
        samples_ptr += sample_length + 1;
      }
      free(samples_buf);
      free(cmp_buf);
    }

    return read_okay;
}

template<class TPMat>
inline IOStatus read_partial_data_fd(int fd, TPMat &result) {
    ssize_t cnt=-1;

    const uint32_t n_samples = result.n_samples;
    const uint32_t n_stripes = result.stripe_stop-result.stripe_start;

    /* load stripes */
    {
      int max_compressed = LZ4_compressBound(sizeof(double) * n_samples);
      char * const cmp_buf = (char *)malloc(max_compressed+sizeof(uint32_t));
      if (cmp_buf==NULL) { return bad_header;} // no better error code

      uint32_t *cmp_buf_size_p = (uint32_t *)cmp_buf;

      cnt = read(fd,cmp_buf_size_p , sizeof(uint32_t) );
      if (cnt != ssize_t(sizeof(uint32_t)) ) {free(cmp_buf); return magic_incompatible;}

      for(uint32_t i = 0; i < n_stripes; i++) {
        uint32_t cmp_size = *cmp_buf_size_p;

        uint32_t read_size = cmp_size;
        if ( (i+1)<n_stripes ) read_size += sizeof(uint32_t); // last one does not have the cmp_size

        cnt = read(fd,cmp_buf , read_size );
        if (cnt != ssize_t(read_size)) {free(cmp_buf); return magic_incompatible;}

        result.stripes[i] = (double *) malloc(sizeof(double) * n_samples);
        if(result.stripes[i] == NULL) {
            fprintf(stderr, "failed\n");
            exit(1);
        }
        cnt = LZ4_decompress_safe(cmp_buf, (char *) result.stripes[i],cmp_size,sizeof(double) * n_samples);
        if (cnt != ssize_t( sizeof(double) * n_samples ) ) {free(cmp_buf); return magic_incompatible;}

        cmp_buf_size_p = (uint32_t *)(cmp_buf+cmp_size);
      }

      free(cmp_buf);
    }

    return read_okay;
}

template<class TPMat>
inline IOStatus read_partial_one_stripe_fd(int fd, TPMat &result, uint32_t stripe_idx) {
    ssize_t cnt=-1;

    const uint32_t n_samples = result.n_samples;

    /* load stripes */
    {
      int max_compressed = LZ4_compressBound(sizeof(double) * n_samples);
      char * const cmp_buf = (char *)malloc(max_compressed+sizeof(uint32_t));
      if (cmp_buf==NULL) { return bad_header;} // no better error code

      uint32_t *cmp_buf_size_p = (uint32_t *)cmp_buf;

      uint32_t curr_idx = stripe_idx;
      while (result.offsets[curr_idx]==0) --curr_idx; // must start reading from the first known offset

      for (;curr_idx<stripe_idx; curr_idx++) { // now get all the intermediate indexes
        if (lseek(fd, result.offsets[curr_idx], SEEK_SET)!=off_t(result.offsets[curr_idx])) {
           free(cmp_buf); return bad_header;
        }

        cnt = read(fd,cmp_buf_size_p , sizeof(uint32_t) );
        if (cnt != ssize_t(sizeof(uint32_t)) ) {free(cmp_buf); return magic_incompatible;}

        uint32_t cmp_size = *cmp_buf_size_p;
        uint32_t read_size = cmp_size;

        result.offsets[curr_idx+1] = result.offsets[curr_idx] + sizeof(uint32_t) + read_size;
      }

      // =======================
      // ready to read my stripe

      if (lseek(fd, result.offsets[stripe_idx], SEEK_SET)!=off_t(result.offsets[stripe_idx])) {
         free(cmp_buf); return bad_header;
      }

      cnt = read(fd,cmp_buf_size_p , sizeof(uint32_t) );
      if (cnt != ssize_t(sizeof(uint32_t)) ) {free(cmp_buf); return magic_incompatible;}

      {
        uint32_t cmp_size = *cmp_buf_size_p;

        uint32_t read_size = cmp_size;

        cnt = read(fd,cmp_buf , read_size );
        if (cnt != ssize_t(read_size)) {free(cmp_buf); return magic_incompatible;}

        result.stripes[stripe_idx] = (double *) malloc(sizeof(double) * n_samples);
        if(result.stripes[stripe_idx] == NULL) {
            fprintf(stderr, "failed\n");
            exit(1);
        }
        cnt = LZ4_decompress_safe(cmp_buf, (char *) result.stripes[stripe_idx],cmp_size,sizeof(double) * n_samples);
        if (cnt != ssize_t( sizeof(double) * n_samples ) ) {free(cmp_buf); return magic_incompatible;}
      }

      free(cmp_buf);
    }

    return read_okay;
}

IOStatus read_partial(const char* input_filename, partial_mat_t** result_out) {
    int fd = open(input_filename, O_RDONLY );
    if (fd==-1) return open_error;

    /* initialize the partial result structure */
    partial_mat_t* result = (partial_mat_t*)malloc(sizeof(partial_mat));

    IOStatus sts = magic_incompatible;

    sts = read_partial_header_fd<partial_mat_t>(fd, *result);
    if (sts==read_okay)
       sts = read_partial_data_fd<partial_mat_t>(fd, *result);

    if (sts==read_okay) {
      IOStatus sts = read_okay;
      /* sanity check the footer */
      uint32_t header[1];
      header[0] = 0;
      ssize_t cnt = read(fd,header,sizeof(uint32_t));
      if (cnt != ssize_t(sizeof(uint32_t))) {sts= magic_incompatible;}
    
      if (sts==read_okay) {
        if ( header[0] != PARTIAL_MAGIC_V2) {sts= magic_incompatible;}
      }
    }

    close(fd);

    if (sts==read_okay) {
      (*result_out) = result;
    } else {
      free(result);
      (*result_out) = NULL;
    }
    return sts;
}

IOStatus read_partial_header(const char* input_filename, partial_dyn_mat_t** result_out) {
    int fd = open(input_filename, O_RDONLY );
    if (fd==-1) return open_error;

    /* initialize the partial result structure */
    partial_dyn_mat_t* result = (partial_dyn_mat_t*)malloc(sizeof(partial_dyn_mat));
    {
      IOStatus sts = read_partial_header_fd<partial_dyn_mat_t>(fd, *result);
      if (sts!=read_okay) {free(result); close(fd); return sts;}
    }

    // save the offset of the first stripe
    const uint32_t n_stripes = result->stripe_stop-result->stripe_start;
    result->stripes = (double**) calloc(n_stripes,sizeof(double*));
    result->offsets = (uint64_t*) calloc(n_stripes,sizeof(uint64_t));
    result->offsets[0] = lseek(fd,0,SEEK_CUR);
    
    close(fd);

    result->filename= strdup(input_filename);

    (*result_out) = result;
    return read_okay;
}

IOStatus read_partial_one_stripe(partial_dyn_mat_t* result, uint32_t stripe_idx) {
    if (result->stripes[stripe_idx]!=0) return read_okay; // will not re-read

    int fd = open(result->filename, O_RDONLY );
    if (fd==-1) return open_error;

    IOStatus sts = read_partial_one_stripe_fd<partial_dyn_mat_t>(fd, *result, stripe_idx);

    close(fd);
    return sts;
}


template<class TPMat>
MergeStatus check_partial(const TPMat* const * partial_mats, int n_partials, bool verbose) {
    if(n_partials <= 0) {
        fprintf(stderr, "Zero or less partials.\n");
        exit(EXIT_FAILURE);
    }

    // sanity check
    const uint32_t n_samples = partial_mats[0]->n_samples;
    bool *stripe_map = (bool*)calloc(sizeof(bool), partial_mats[0]->stripe_total);
    uint32_t stripe_count = 0;

    for(int i = 0; i < n_partials; i++) {
        if(partial_mats[i]->n_samples != n_samples) {
            free(stripe_map);
            if (verbose) {
                fprintf(stderr, "Wrong number of samples in %i, %i!=%i\n",
                        i,int(partial_mats[i]->n_samples),int(n_samples));
            }
            return partials_mismatch;
        }

        if(partial_mats[0]->stripe_total != partial_mats[i]->stripe_total) {
            free(stripe_map);
            if (verbose) {
                fprintf(stderr, "Wrong number of stripes in %i, %i!=%i\n",
                        i,int(partial_mats[0]->stripe_total), int(partial_mats[i]->stripe_total));
            }
            return partials_mismatch;
        }
        if(partial_mats[0]->is_upper_triangle != partial_mats[i]->is_upper_triangle) {
            free(stripe_map);
            if (verbose) {
                fprintf(stderr, "Wrong number of is_upper_triangle in %i, %i!=%i\n",
                        i,int(partial_mats[0]->is_upper_triangle),int(partial_mats[i]->is_upper_triangle));
            }
            return square_mismatch;
        }
        for(uint32_t j = 0; j < n_samples; j++) {
            if(strcmp(partial_mats[0]->sample_ids[j], partial_mats[i]->sample_ids[j]) != 0) {
                free(stripe_map);
                if (verbose) {
                    fprintf(stderr, "Wrong number of sample id %i in %i, %s!=%s\n",
                            j,i,partial_mats[0]->sample_ids[j], partial_mats[i]->sample_ids[j]);
                }
                return sample_id_consistency;
            }
        }
        for(uint32_t j = partial_mats[i]->stripe_start; j < partial_mats[i]->stripe_stop; j++) {
            if(stripe_map[j]) {
                if (verbose) {
                    fprintf(stderr, "Overlap in %i vs %i\n",
                            i,j);
                }
                free(stripe_map);
                return stripes_overlap;
            }
            stripe_map[j] = true;
            stripe_count += 1;
        }
    }
    free(stripe_map);

    if(stripe_count != partial_mats[0]->stripe_total) {
        if (verbose) {
            fprintf(stderr, "Insufficient number of stripes found, %i!=%i\n",
                    int(stripe_count), int(partial_mats[0]->stripe_total));
        }
        return incomplete_stripe_set;
    }

    return merge_okay;
}

MergeStatus validate_partial(const partial_dyn_mat_t* const * partial_mats, int n_partials) {
    return check_partial(partial_mats, n_partials, true);
}

// Will keep only the strictly necessary stripes in memory... reading just in time
class PartialStripes : public su::ManagedStripes {
        private:
           const uint32_t n_partials;
           mutable partial_dyn_mat_t* * partial_mats; // link only, not owned

           static bool in_range(const partial_dyn_mat_t &partial_mat, uint32_t stripe) {
             return (stripe>=partial_mat.stripe_start) && (stripe<partial_mat.stripe_stop);
           }

           uint32_t find_partial_idx(uint32_t stripe) const {
              for (uint32_t i=0; i<n_partials; i++) {
                if (in_range(*(partial_mats[i]),stripe)) return i;
              }
              return 0; // should never get here
           }
        public:
           PartialStripes(uint32_t _n_partials, partial_dyn_mat_t* * _partial_mats)
           : n_partials(_n_partials)
           , partial_mats(_partial_mats)
           {}

           virtual const double *get_stripe(uint32_t stripe) const {
              uint32_t pidx = find_partial_idx(stripe);
              partial_dyn_mat_t * const partial_mat = partial_mats[pidx];
              uint32_t sidx = stripe-partial_mat->stripe_start;

              if (partial_mat->stripes[sidx]==NULL) {
                  read_partial_one_stripe(partial_mat,sidx);
                  // ignore any errors, not clear what to do
                  // will just return NULL
              }

              return partial_mat->stripes[sidx];
           }
           virtual void release_stripe(uint32_t stripe) const {
              uint32_t pidx = find_partial_idx(stripe);
              partial_dyn_mat_t * const partial_mat = partial_mats[pidx];
              uint32_t sidx = stripe-partial_mat->stripe_start;

              if (partial_mat->stripes[sidx]!=NULL) {
                 free(partial_mat->stripes[sidx]);
                 partial_mat->stripes[sidx]=NULL;
              }
           }
};

template<class TReal, class TMat>
MergeStatus merge_partial_to_matrix_T(partial_dyn_mat_t* * partial_mats, int n_partials, 
                                      const char *mmap_dir, /* if NULL or "", use malloc */
                                      TMat** result /* out */ ) {
    if (mmap_dir!=NULL) {
     if (mmap_dir[0]==0) mmap_dir = NULL; // easier to have a simple test going on
    }

    MergeStatus err = check_partial(partial_mats, n_partials, false);
    if (err!=merge_okay) return err;

    initialize_mat_full_no_biom_T<TReal,TMat>(*result, partial_mats[0]->sample_ids, partial_mats[0]->n_samples,mmap_dir);

    if ((*result)==NULL) return incomplete_stripe_set;
    if ((*result)->matrix==NULL) return incomplete_stripe_set;
    if ((*result)->sample_ids==NULL) return incomplete_stripe_set;

    PartialStripes ps(n_partials,partial_mats);
    const uint32_t tile_size = (mmap_dir==NULL) ? \
                                  (128/sizeof(TReal)) : /* keep it small for memory access, to fit in chip cache */ \
                                  (4096/sizeof(TReal)); /* make it larger for mmap, as the limiting factor is swapping */
    su::stripes_to_matrix_T<TReal>(ps, partial_mats[0]->n_samples, partial_mats[0]->stripe_total, (*result)->matrix, tile_size);

    return merge_okay;
}

MergeStatus merge_partial_to_matrix(partial_dyn_mat_t* * partial_mats, int n_partials, mat_full_fp64_t** result) {
  return merge_partial_to_matrix_T<double,mat_full_fp64_t>(partial_mats, n_partials, NULL, result);
}

MergeStatus merge_partial_to_matrix_fp32(partial_dyn_mat_t* * partial_mats, int n_partials, mat_full_fp32_t** result) {
  return merge_partial_to_matrix_T<float,mat_full_fp32_t>(partial_mats, n_partials, NULL, result);
}

MergeStatus merge_partial_to_mmap_matrix(partial_dyn_mat_t* * partial_mats, int n_partials, const char *mmap_dir, mat_full_fp64_t** result) {
  return merge_partial_to_matrix_T<double,mat_full_fp64_t>(partial_mats, n_partials, mmap_dir, result);
}

MergeStatus merge_partial_to_mmap_matrix_fp32(partial_dyn_mat_t* * partial_mats, int n_partials, const char *mmap_dir, mat_full_fp32_t** result) {
  return merge_partial_to_matrix_T<float,mat_full_fp32_t>(partial_mats, n_partials, mmap_dir, result);
}

#endif // UNIFRAC_WASM (partial I/O + merge)

// compat versions

#include "api_compat.hpp"

// skbio_alt pass-thoughs


// Find eigen values and vectors
// Based on N. Halko, P.G. Martinsson, Y. Shkolnisky, and M. Tygert.
//     Original Paper: https://arxiv.org/abs/1007.5510
// centered == n x n, must be symmetric, Note: will be used in-place as temp buffer

void find_eigens_fast(const uint32_t n_samples, const uint32_t n_dims, double * centered, double **eigenvalues, double **eigenvectors) {
  su::find_eigens_fast(n_samples, n_dims, centered, *eigenvalues, *eigenvectors);
}

void find_eigens_fast_fp32(const uint32_t n_samples, const uint32_t n_dims, float * centered, float **eigenvalues, float **eigenvectors) {
  su::find_eigens_fast(n_samples, n_dims, centered, *eigenvalues, *eigenvectors);
}

/*
    Perform Principal Coordinate Analysis.

    Principal Coordinate Analysis (PCoA) is a method similar
    to Principal Components Analysis (PCA) with the difference that PCoA
    operates on distance matrices, typically with non-euclidian and thus
    ecologically meaningful distances like UniFrac in microbiome research.

    In ecology, the euclidean distance preserved by Principal
    Component Analysis (PCA) is often not a good choice because it
    deals poorly with double zeros (Species have unimodal
    distributions along environmental gradients, so if a species is
    absent from two sites at the same site, it can't be known if an
    environmental variable is too high in one of them and too low in
    the other, or too low in both, etc. On the other hand, if an
    species is present in two sites, that means that the sites are
    similar.).

    Note that the returned eigenvectors are not normalized to unit length.
*/

// mat       - in, result of unifrac compute
// n_samples - in, size of the matrix (n x n)
// n_dims    - in, Dimensions to reduce the distance matrix to. This number determines how many eigenvectors and eigenvalues will be returned.
// eigenvalues - out, alocated buffer of size n_dims
// samples     - out, alocated buffer of size n_dims x n_samples
// proportion_explained - out, allocated buffer of size n_dims

void pcoa_seeded(const double * mat, const uint32_t n_samples, const uint32_t n_dims, int seed, double * *eigenvalues, double * *samples, double * *proportion_explained) {
  su::pcoa(mat, n_samples, n_dims, *eigenvalues, *samples, *proportion_explained, seed);
}

void pcoa_fp32_seeded(const float * mat, const uint32_t n_samples, const uint32_t n_dims, int seed, float * *eigenvalues, float * *samples, float * *proportion_explained) {
  su::pcoa(mat, n_samples, n_dims, *eigenvalues, *samples, *proportion_explained, seed);
}

void pcoa_mixed_seeded(const double * mat, const uint32_t n_samples, const uint32_t n_dims, int seed, float * *eigenvalues, float * *samples, float * *proportion_explained) {
  su::pcoa(mat, n_samples, n_dims, *eigenvalues, *samples, *proportion_explained, seed);
}

void pcoa(const double * mat, const uint32_t n_samples, const uint32_t n_dims, double * *eigenvalues, double * *samples, double * *proportion_explained) {
  pcoa_seeded(mat, n_samples, n_dims, -1, eigenvalues, samples, proportion_explained);
}

void pcoa_fp32(const float * mat, const uint32_t n_samples, const uint32_t n_dims, float * *eigenvalues, float * *samples, float * *proportion_explained) {
  pcoa_fp32_seeded(mat, n_samples, n_dims, -1, eigenvalues, samples, proportion_explained);
}

void pcoa_mixed(const double * mat, const uint32_t n_samples, const uint32_t n_dims, float * *eigenvalues, float * *samples, float * *proportion_explained) {
  pcoa_mixed_seeded(mat, n_samples, n_dims, -1, eigenvalues, samples, proportion_explained);
}

