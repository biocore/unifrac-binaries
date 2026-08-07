#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <iomanip>
#include <glob.h>
#include "api.hpp"
#include "cmd.hpp"
// Using inlined-header-only funtions
#include "biom.hpp"

enum Format {format_invalid,format_ascii, format_hdf5_fp32, format_hdf5_fp64, format_hdf5_nodist};

void usage() {
    std::cout << "usage: ssu -i <biom> -o <out.dm> -m [METHOD] -t <newick> [-a alpha] [-f]  [--vaw]" << std::endl;
    std::cout << "    [--mode MODE] [--start starting-stripe] [--stop stopping-stripe] [--partial-pattern <glob>]" << std::endl;
    std::cout << "    [--n-partials number_of_partitions] [--report-bare] [--format|-r out-mode]" << std::endl;
    std::cout << "    [--normalize-sample-counts true|false] [--n-substeps n] [--pcoa dims] [--diskbuf path]" << std::endl;
    std::cout << std::endl;
    std::cout << "    -i\t\tThe input BIOM table." << std::endl;
    std::cout << "    -t\t\tThe input phylogeny in newick." << std::endl;
    std::cout << "    -m\t\tThe method, [unweighted | weighted_normalized | weighted_unnormalized | unweighted_unnormalized | generalized |" << std::endl;
    std::cout << "                       unweighted_fp64 | weighted_normalized_fp64 | weighted_unnormalized_fp64 |" << std::endl;
    std::cout << "                       unweighted_unnormalized_fp64 | generalized_fp64 |" << std::endl;
    std::cout << "                       unweighted_fp32 | weighted_normalized_fp32 | weighted_unnormalized_fp32 |" << std::endl;
    std::cout << "                       unweighted_unnormalized_fp32 | generalized_fp32]." << std::endl;
    std::cout << "    -o\t\tThe output distance matrix." << std::endl;
    std::cout << "    -g\t\t[OPTIONAL] The input grouping in TSV." << std::endl;
    std::cout << "    -c\t\t[OPTIONAL] The columns(s) to use for grouping, multiple values comma separated." << std::endl;
    std::cout << "    -a\t\t[OPTIONAL] Generalized UniFrac alpha, default is 1." << std::endl;
    std::cout << "    -f\t\t[OPTIONAL] Bypass tips, reduces compute by about 50%." << std::endl;
    std::cout << "    --vaw\t[OPTIONAL] Variance adjusted, default is to not adjust for variance." << std::endl;
    std::cout << "    --mode\t[OPTIONAL] Mode of operation:" << std::endl;
    std::cout << "    \t\t    one-off : [DEFAULT] compute UniFrac." << std::endl;
    std::cout << "    \t\t    partial : Compute UniFrac over a subset of stripes." << std::endl;
    std::cout << "    \t\t    partial-report : Start and stop suggestions for partial compute." << std::endl;
    std::cout << "    \t\t    merge-partial : Merge partial UniFrac results." << std::endl;
    std::cout << "    \t\t    check-partial : Check partial UniFrac results." << std::endl;
    std::cout << "    \t\t    multi : compute UniFrac multiple times." << std::endl;
    std::cout << "    --start\t[OPTIONAL] If mode==partial, the starting stripe." << std::endl;
    std::cout << "    --stop\t[OPTIONAL] If mode==partial, the stopping stripe." << std::endl;
    std::cout << "    --partial-pattern\t[OPTIONAL] If mode==merge-partial or check-partial, a glob pattern for partial outputs to merge." << std::endl;
    std::cout << "    --n-partials\t[OPTIONAL] If mode==partial-report, the number of partitions to compute." << std::endl;
    std::cout << "    --report-bare\t[OPTIONAL] If mode==partial-report, produce barebones output." << std::endl;
    std::cout << "    --n-substeps\t[OPTIONAL] Internally split the problem in n substeps for reduced memory footprint, default is 1." << std::endl;
    std::cout << "    --normalize-sample-counts\t[OPTIONAL] Should it normalize sample counts?:" << std::endl;
    std::cout << "    \t\t    true  : [DEFAULT] Do normalize, i.e. standard unifrac." << std::endl;
    std::cout << "    \t\t    false : Do not normalize, i.e. absolute quant mode." << std::endl;
    std::cout << "    --format|-r\t[OPTIONAL]  Output format:" << std::endl;
    std::cout << "    \t\t    ascii : Original ASCII format. (default if mode==one-off)" << std::endl;
    std::cout << "    \t\t    hdf5 : HFD5 format.  May be fp32 or fp64, depending on method." << std::endl;
    std::cout << "    \t\t    hdf5_fp32 : HFD5 format, using fp32 precision." << std::endl;
    std::cout << "    \t\t    hdf5_fp64 : HFD5 format, using fp64 precision." << std::endl;
    std::cout << "    \t\t    hdf5_nodist : HFD5 format, no distance matrix. (default if mode==multi)" << std::endl;
    std::cout << "    --subsample-depth\tDepth of subsampling of the input BIOM before computing unifrac (required for mode==multi, optional for one-off)" << std::endl;
    std::cout << "    --subsample-replacement\t[OPTIONAL] Subsample with or without replacement (default is with)" << std::endl;
    std::cout << "    --n-subsamples\t[OPTIONAL] if mode==multi, number of subsampled UniFracs to compute (default: 100)" << std::endl;
    std::cout << "    --permanova\t[OPTIONAL] Number of PERMANOVA permutations to compute (default: 999 with -g, do not compute if 0)" << std::endl;
    std::cout << "    --pcoa\t[OPTIONAL] Number of PCoA dimensions to compute (default: 10, do not compute if 0)" << std::endl;
    std::cout << "    --seed\t[OPTIONAL] Seed to use for initializing the random gnerator" << std::endl;
    std::cout << "    --diskbuf\t[OPTIONAL] Use a disk buffer to reduce memory footprint. Provide path to a fast partition (ideally NVMe)." << std::endl;
    std::cout << "    -n\t\t[OPTIONAL] DEPRECATED, no-op." << std::endl;
    std::cout << std::endl;
    std::cout << "Environment variables: " << std::endl;
    std::cout << "    CPU parallelism is controlled by OMP_NUM_THREADS. If not defined, all detected core will be used." << std::endl;
    std::cout << "    GPU offload can be disabled with UNIFRAC_USE_GPU=N. By default, if a NVIDIA GPU is detected, it will be used." << std::endl;
    std::cout << "    A specific GPU can be selected with ACC_DEVICE_NUM. If not defined, the first GPU will be used." << std::endl;
    std::cout << std::endl;
    std::cout << "Citations: " << std::endl;
    std::cout << "    For UniFrac, please see:" << std::endl;
    std::cout << "        Sfiligoi et al. mSystems 2022; DOI: 10.1128/msystems.00028-22" << std::endl;
    std::cout << "        McDonald et al. Nature Methods 2018; DOI: 10.1038/s41592-018-0187-8" << std::endl;
    std::cout << "        Lozupone and Knight Appl Environ Microbiol 2005; DOI: 10.1128/AEM.71.12.8228-8235.2005" << std::endl;
    std::cout << "        Lozupone et al. Appl Environ Microbiol 2007; DOI: 10.1128/AEM.01996-06" << std::endl;
    std::cout << "        Hamady et al. ISME 2010; DOI: 10.1038/ismej.2009.97" << std::endl;
    std::cout << "        Lozupone et al. ISME 2011; DOI: 10.1038/ismej.2010.133" << std::endl;
    std::cout << "    For Generalized UniFrac, please see: " << std::endl;
    std::cout << "        Chen et al. Bioinformatics 2012; DOI: 10.1093/bioinformatics/bts342" << std::endl;
    std::cout << "    For Variance Adjusted UniFrac, please see: " << std::endl;
    std::cout << "        Chang et al. BMC Bioinformatics 2011; DOI: 10.1186/1471-2105-12-118" << std::endl;
    std::cout << std::endl;
}

/* Indexed by ComputeStatus, so it needs one entry per enumerator. The bound is
 * left off and asserted instead: status_enum.hpp is append-only for ABI, so a
 * status appended without a matching entry here would leave this table short by
 * one and turn any error report for it into an out-of-bounds read.
 */
const char* compute_status_messages[] = {"No error.",
                                          "The tree file cannot be found.", 
                                          "The table file cannot be found.",
                                          "The table file contains an empty table.",
                                          "An unknown method was requested.", 
                                          "Table observation IDs are not a subset of the tree tips. This error can also be triggered if a node name contains a single quote (this is unlikely).",
                                          "Error creating the output.",
                                          "The requested method is not supported.",
                                          "The grouping file cannot be found or does not have the necessary data.",
                                          "Device-resident input and output are not supported."};
static_assert(sizeof(compute_status_messages) / sizeof(compute_status_messages[0])
                  == unsupported_device + 1,
              "compute_status_messages needs one entry per ComputeStatus; "
              "unsupported_device must remain the last enumerator");


// https://stackoverflow.com/questions/8401777/simple-glob-in-c-on-unix-system
inline std::vector<std::string> glob(const std::string& pat){
    using namespace std;
    glob_t glob_result;
    glob(pat.c_str(),GLOB_TILDE,NULL,&glob_result);
    vector<string> ret;
    for(unsigned int i=0;i<glob_result.gl_pathc;++i){
        ret.push_back(string(glob_result.gl_pathv[i]));
    }
    globfree(&glob_result);
    return ret;
}


void err(std::string msg) {
    std::cerr << "ERROR: " << msg << std::endl << std::endl;
    usage();
}

int mode_partial_report(const std::string table_filename, unsigned int npartials, bool bare) {
    if(table_filename.empty()) {
        err("table filename missing");
        return EXIT_FAILURE;
    }
    
    if(npartials < 1) {
        err("--n-partials cannot be < 1");
        exit(EXIT_FAILURE);
    }

    int n_samples = su::biom::load_n_samples(table_filename.c_str());
    int total_stripes = (n_samples + 1) / 2;

    if(!bare) {
        std::cout << "Total samples: " << n_samples << std::endl;
        std::cout << "Total stripes: " << total_stripes << std::endl;
    }

    unsigned int fullchunk = (total_stripes + npartials - 1) / npartials;  // this computes the ceiling
    unsigned int smallchunk = total_stripes / npartials;
    
    unsigned int n_fullbins = total_stripes % npartials;
    if(n_fullbins == 0)
        n_fullbins = npartials;

    unsigned int start = 0;
    unsigned int stop = 0;
    for(unsigned int p = 0; p < npartials; p++) {
        if(p < n_fullbins) {
            stop = start + fullchunk;
            if(bare) 
                std::cout << start << "\t" << stop << std::endl;
            else
                std::cout << "Partition " << p << ", suggested start and stop: " << start << ", " << stop << std::endl;
            start = start + fullchunk;
        } else {
            stop = start + smallchunk;  // stripe end 
            if(bare) 
                std::cout << start << "\t" << stop << std::endl;
            else
                std::cout << "Partition " << p << ", suggested start and stop: " << start << ", " << stop << std::endl;
            start = start + smallchunk;
        }
    }

    return EXIT_SUCCESS;
} 

int mode_merge_partial_fp32(const char * output_filename, Format format_val, unsigned int pcoa_dims,
                            unsigned int permanova_perms, const char *grouping_filename, const char *grouping_columns,
                            size_t partials_size, partial_dyn_mat_t* * partial_mats,
                            const char * mmap_dir) {
    mat_full_fp32_t *result = NULL;

    MergeStatus status = merge_partial_to_mmap_matrix_fp32(partial_mats, partials_size, mmap_dir, &result);

    if(status != merge_okay) {
        std::ostringstream msg;
        msg << "Unable to complete merge; err " << status;
        err(msg.str());
        return EXIT_FAILURE;
    }

    // TODO: Add support for PERMANOVA

    IOStatus iostatus;
    iostatus = write_mat_from_matrix_hdf5_fp32(output_filename, result, pcoa_dims, format_val!=format_hdf5_nodist);
    destroy_mat_full_fp32(&result);
    
    if(iostatus != write_okay) {
        std::ostringstream msg; 
        msg << "Unable to write; err " << iostatus;
        err(msg.str()); 
        return EXIT_FAILURE;
    }   
        
    return EXIT_SUCCESS;
}

int mode_merge_partial_fp64(const char * output_filename, Format format_val, unsigned int pcoa_dims,
                            unsigned int permanova_perms, const char *grouping_filename, const char *grouping_columns,
                            size_t partials_size, partial_dyn_mat_t* * partial_mats,
                            const char * mmap_dir) {
    mat_full_fp64_t *result = NULL;

    MergeStatus status = merge_partial_to_mmap_matrix(partial_mats, partials_size, mmap_dir, &result);

    if(status != merge_okay) {
        std::ostringstream msg;
        msg << "Unable to complete merge; err " << status;
        err(msg.str());
        return EXIT_FAILURE;
    }

    // TODO: Add support for PERMANOVA

    IOStatus iostatus;
    if (format_val!=format_ascii) {
     iostatus = write_mat_from_matrix_hdf5_fp64(output_filename, result, pcoa_dims, format_val!=format_hdf5_nodist);
    } else {
     iostatus = write_mat_from_matrix(output_filename, result);
    }
    destroy_mat_full_fp64(&result);

    if(iostatus != write_okay) {
        std::ostringstream msg;
        msg << "Unable to write; err " << iostatus;
        err(msg.str());
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}


int mode_merge_partial(const std::string &output_filename, Format format_val, unsigned int pcoa_dims,
                       unsigned int permanova_perms, const std::string &grouping_filename, const std::string &grouping_columns,
                       const std::string &partial_pattern,
                       const std::string &mmap_dir) {
    if(output_filename.empty()) {
        err("output filename missing");
        return EXIT_FAILURE;
    }

    if(partial_pattern.empty()) {
        std::string msg("Partial file pattern missing. For instance, if your partial results\n" \
                        "are named 'ssu.unweighted.start0.partial', 'ssu.unweighted.start10.partial', \n" \
                        "etc, then a pattern of 'ssu.unweighted.start*.partial' would make sense");
        err(msg);
        return EXIT_FAILURE;
    }
    
    if((permanova_perms>0) && grouping_filename.empty()) {
        err("grouping filename missing");
        return EXIT_FAILURE;
    }
    
    if((permanova_perms>0) && grouping_columns.empty()) {
        err("grouping columns missing");
        return EXIT_FAILURE;
    }
    
    std::vector<std::string> partials = glob(partial_pattern);
    partial_dyn_mat_t** partial_mats = (partial_dyn_mat_t**)malloc(sizeof(partial_dyn_mat_t*) * partials.size());
    for(size_t i = 0; i < partials.size(); i++) {
        IOStatus io_err = read_partial_header(partials[i].c_str(), &partial_mats[i]);
        if(io_err != read_okay) {
            std::ostringstream msg;
            msg << "Unable to parse file (" << partials[i] << "); err " << io_err;
            err(msg.str());
            return EXIT_FAILURE;
        }
    }

    const char * mmap_dir_c = mmap_dir.empty() ? NULL : mmap_dir.c_str();
    const char * grouping_c = (permanova_perms>0) ? grouping_filename.c_str() : NULL ;
    const char * columns_c = (permanova_perms>0) ? grouping_columns.c_str() : NULL ;

    int status;
    if (format_val==format_hdf5_fp64) {
     status = mode_merge_partial_fp64(output_filename.c_str(), format_val,
                                      pcoa_dims, permanova_perms, grouping_c, columns_c,
                                      partials.size(), partial_mats, mmap_dir_c);
    } else if (format_val==format_hdf5_fp32) {
     status = mode_merge_partial_fp32(output_filename.c_str(), format_val,
                                      pcoa_dims, permanova_perms, grouping_c, columns_c,
                                      partials.size(), partial_mats, mmap_dir_c);
    } else {
     status = mode_merge_partial_fp64(output_filename.c_str(), format_val,
                                      pcoa_dims, permanova_perms, grouping_c, columns_c,
                                      partials.size(), partial_mats, mmap_dir_c);
    }

    for(size_t i = 0; i < partials.size(); i++) {
      destroy_partial_dyn_mat(&partial_mats[i]);
    }

    return status;
}

int mode_check_partial(const std::string &partial_pattern) {
    if(partial_pattern.empty()) {
        std::string msg("Partial file pattern missing. For instance, if your partial results\n" \
                        "are named 'ssu.unweighted.start0.partial', 'ssu.unweighted.start10.partial', \n" \
                        "etc, then a pattern of 'ssu.unweighted.start*.partial' would make sense");
        err(msg);
        return EXIT_FAILURE;
    }

    std::vector<std::string> partials = glob(partial_pattern);
    std::cout << "Processing " << partials.size() << " partial files." << std::endl;


    partial_dyn_mat_t** partial_mats = (partial_dyn_mat_t**)malloc(sizeof(partial_dyn_mat_t*) * partials.size());
    for(size_t i = 0; i < partials.size(); i++) {
        IOStatus io_err = read_partial_header(partials[i].c_str(), &partial_mats[i]);
        if(io_err != read_okay) {
            std::ostringstream msg;
            msg << "Unable to parse file (" << partials[i] << "); err " << io_err;
            err(msg.str());
            return EXIT_FAILURE;
        }
    }

    int status = validate_partial(partial_mats,partials.size());

    if (status==merge_okay) {
        std::cout << "Partials are ready to be merged" << std::endl;
        uint64_t n_samples_64 = partial_mats[0]->n_samples;
        uint64_t msize32 = sizeof(float) * n_samples_64 * n_samples_64;
        uint64_t msize64 = sizeof(double) * n_samples_64 * n_samples_64;

        std::cout << "n_samples = " << n_samples_64 << std::endl;
        std::cout << "Matrix RAM needs: fp32 = " << msize32/(1024.0*1024*1024) << "GB, fp64 = " << msize64/(1024.0*1024*1024) << "GB"<< std::endl;
    }

    for(size_t i = 0; i < partials.size(); i++) {
      destroy_partial_dyn_mat(&partial_mats[i]);
    }

    return status;
}

int mode_partial(std::string table_filename, std::string tree_filename, 
                 std::string output_filename, std::string method_string,
                 bool vaw, double g_unifrac_alpha, bool bypass_tips, bool normalize_sample_counts,
                 unsigned int nsubsteps, int start_stripe, int stop_stripe) {
    if(output_filename.empty()) {
        err("output filename missing");
        return EXIT_FAILURE;
    }

    if(table_filename.empty()) {
        err("table filename missing");
        return EXIT_FAILURE;
    }

    if(tree_filename.empty()) {
        err("tree filename missing");
        return EXIT_FAILURE;
    }
    
    if(method_string.empty()) {
        err("method missing");
        return EXIT_FAILURE;
    }

    if(start_stripe < 0) {
        err("Starting stripe must be >= 0");
        return EXIT_FAILURE;
    }
    if(stop_stripe <= start_stripe) {
        err("In '--mode partial', the stop and start stripes must be specified, and the stop stripe must be > start stripe");
        return EXIT_FAILURE;
    }

    partial_mat_t *result = NULL;
    compute_status status;
    status = partial_v3(table_filename.c_str(), tree_filename.c_str(), method_string.c_str(), 
                        vaw, g_unifrac_alpha, bypass_tips, normalize_sample_counts,
			nsubsteps, start_stripe, stop_stripe, &result);
    if(status != okay || result == NULL) {
        fprintf(stderr, "Compute failed in partial: %s\n", compute_status_messages[status]);
        exit(EXIT_FAILURE);
    }
   
    io_status err = write_partial(output_filename.c_str(), result);
    destroy_partial_mat(&result);

    if(err != write_okay){
        fprintf(stderr, "Write failed: %s\n", err == open_error ? "could not open output" : "unknown error");
        return EXIT_FAILURE;
    } 

    return EXIT_SUCCESS;
}

int mode_one_off(const std::string &table_filename, const std::string &tree_filename, 
                 const std::string &output_filename, const std::string &format_str, Format format_val, 
                 const std::string &method_string, unsigned int subsample_depth, bool subsample_with_replacement, unsigned int pcoa_dims,
                 unsigned int permanova_perms, const std::string &grouping_filename, const std::string &grouping_columns,
                 bool vaw, double g_unifrac_alpha, bool bypass_tips, bool normalize_sample_counts,
                 unsigned int nsubsteps, const std::string &mmap_dir) {
    if(output_filename.empty()) {
        err("output filename missing");
        return EXIT_FAILURE;
    }

    if(table_filename.empty()) {
        err("table filename missing");
        return EXIT_FAILURE;
    }

    if(tree_filename.empty()) {
        err("tree filename missing");
        return EXIT_FAILURE;
    }
    
    if((permanova_perms>0) && grouping_filename.empty()) {
        err("grouping filename missing");
        return EXIT_FAILURE;
    }
    
    if((permanova_perms>0) && grouping_columns.empty()) {
        err("grouping columns missing");
        return EXIT_FAILURE;
    }
    
    if(method_string.empty()) {
        err("method missing");
        return EXIT_FAILURE;
    }

    const char * mmap_dir_c = mmap_dir.empty() ? NULL : mmap_dir.c_str();
    compute_status status = okay;
    if (format_val==format_ascii) {
      mat_t *result = NULL;

      if(subsample_depth>0) {
        err("Subsampling not supported with ASCII output.");
        return EXIT_FAILURE;
      }
      status = unifrac_to_txt_file_v3(table_filename.c_str(), tree_filename.c_str(), output_filename.c_str(),
                                      method_string.c_str(), vaw, g_unifrac_alpha, bypass_tips, normalize_sample_counts, 
				      nsubsteps,
                                      mmap_dir_c);

      if (status != okay) {
        fprintf(stderr, "Compute failed in one_off: %s\n", compute_status_messages[status]);
      }
    } else {
      const char * grouping_c = (permanova_perms>0) ? grouping_filename.c_str() : NULL ;
      const char * columns_c = (permanova_perms>0) ? grouping_columns.c_str() : NULL ;

      status = unifrac_to_file_v3(table_filename.c_str(), tree_filename.c_str(), output_filename.c_str(),
                                  method_string.c_str(), vaw, g_unifrac_alpha, bypass_tips, normalize_sample_counts, 
				  nsubsteps, format_str.c_str(),
                                  subsample_depth, subsample_with_replacement,
                                  pcoa_dims, permanova_perms, grouping_c, columns_c, mmap_dir_c);

      if (status != okay) {
        fprintf(stderr, "Compute failed in one_off: %s\n", compute_status_messages[status]);
      }
    }

    return (status==okay) ? EXIT_SUCCESS : EXIT_FAILURE;
}

int mode_multi(const std::string &table_filename, const std::string &tree_filename, 
               const std::string &output_filename, const std::string &format_str, Format format_val, 
               const std::string &method_string,
               unsigned int n_subsamples, unsigned int subsample_depth, bool subsample_with_replacement,
               unsigned int pcoa_dims,
               unsigned int permanova_perms, const std::string &grouping_filename, const std::string &grouping_columns,
               bool vaw, double g_unifrac_alpha, bool bypass_tips, bool normalize_sample_counts,
               unsigned int nsubsteps, const std::string &mmap_dir) {
    if(output_filename.empty()) {
        err("output filename missing");
        return EXIT_FAILURE;
    }

    if(table_filename.empty()) {
        err("table filename missing");
        return EXIT_FAILURE;
    }

    if(tree_filename.empty()) {
        err("tree filename missing");
        return EXIT_FAILURE;
    }
    
    if((permanova_perms>0) && grouping_filename.empty()) {
        err("grouping filename missing");
        return EXIT_FAILURE;
    }
    
    if((permanova_perms>0) && grouping_columns.empty()) {
        err("grouping columns missing");
        return EXIT_FAILURE;
    }
    
    if(method_string.empty()) {
        err("method missing");
        return EXIT_FAILURE;
    }

    if(subsample_depth<1) {
      err("subsample_depth cannot be 0.");
      return EXIT_FAILURE;
    }

    if(n_subsamples<1) {
      err("n_subsamples cannot be 0.");
      return EXIT_FAILURE;
    }

    if (format_val==format_ascii) {
      err("ASCII format not supported in multi mode");
      return EXIT_FAILURE;
    }

    compute_status status = okay;
    {
      const char * mmap_dir_c = mmap_dir.empty() ? NULL : mmap_dir.c_str();
      const char * grouping_c = (permanova_perms>0) ? grouping_filename.c_str() : NULL ;
      const char * columns_c = (permanova_perms>0) ? grouping_columns.c_str() : NULL ;

      status = unifrac_multi_to_file_v3(table_filename.c_str(), tree_filename.c_str(), output_filename.c_str(),
                                        method_string.c_str(), vaw, g_unifrac_alpha, bypass_tips, normalize_sample_counts,
					nsubsteps, format_str.c_str(),
                                        n_subsamples, subsample_depth, subsample_with_replacement,
                                        pcoa_dims, permanova_perms, grouping_c, columns_c, mmap_dir_c);

      if (status != okay) {
        fprintf(stderr, "Compute failed in multi: %s\n", compute_status_messages[status]);
      }
    }

    return (status==okay) ? EXIT_SUCCESS : EXIT_FAILURE;
}

Format get_format(const std::string &format_string, const std::string &method_string, const std::string &mode_string) {
    Format format_val = format_invalid;
    if (format_string.empty()) {
        if (mode_string!="multi") {
          format_val = format_ascii;
        } else {
          format_val = format_hdf5_nodist;
        }
    } else if (format_string == "ascii") {
        format_val = format_ascii;
    } else if (format_string == "hdf5_fp32") {
        format_val = format_hdf5_fp32;
    } else if (format_string == "hdf5_fp64") {
        format_val = format_hdf5_fp64;
    } else if (format_string == "hdf5_nodist") {
        format_val = format_hdf5_nodist;
    } else if (format_string == "hdf5") {
        if ((method_string=="unweighted_fp64") || (method_string=="weighted_normalized_fp64") || (method_string=="weighted_unnormalized_fp64") || (method_string=="generalized_fp64") || (method_string=="unweighted_unnormalized_fp64"))
           format_val = format_hdf5_fp64;
        else
           format_val = format_hdf5_fp32;
    }

    return format_val;
}

std::string format2str(Format format_val) {
  if (format_val==format_hdf5_nodist) {
    return "hdf5_nodist";
  } else if (format_val==format_hdf5_fp32) {
    return "hdf5_fp32";
  } else if (format_val==format_hdf5_fp64) {
    return "hdf5_fp64";
  } else if (format_val==format_ascii) {
    return "ascii";
  } 
  return "invalid";
}

int main(int argc, char **argv){
    InputParser input(argc, argv);
    if(input.cmdOptionExists("-h") || input.cmdOptionExists("--help") || argc == 1) {
        usage();
        return EXIT_SUCCESS;
    }

    unsigned int nsubsteps;
    std::string table_filename = input.getCmdOption("-i");
    std::string tree_filename = input.getCmdOption("-t");
    std::string output_filename = input.getCmdOption("-o");
    std::string method_string = input.getCmdOption("-m");
    // deprecated, but we still want to support it, even as a no-op
    std::string nold_arg = input.getCmdOption("-n");
    std::string grouping_filename = input.getCmdOption("-g");
    std::string grouping_columns = input.getCmdOption("-c");
    std::string nsubsteps_arg = input.getCmdOption("--n-substeps");
    std::string gunifrac_arg = input.getCmdOption("-a");
    std::string mode_arg = input.getCmdOption("--mode");
    std::string start_arg = input.getCmdOption("--start");
    std::string stop_arg = input.getCmdOption("--stop");
    std::string partial_pattern = input.getCmdOption("--partial-pattern");
    std::string npartials = input.getCmdOption("--n-partials");
    std::string normsc_arg = input.getCmdOption("--normalize-sample-counts");
    std::string report_bare = input.getCmdOption("--report-bare");
    std::string format_arg = input.getCmdOption("--format");
    std::string sformat_arg = input.getCmdOption("-r");
    std::string pcoa_arg = input.getCmdOption("--pcoa");
    std::string permanova_arg = input.getCmdOption("--permanova");
    std::string seed_arg = input.getCmdOption("--seed");
    std::string subsample_depth_arg = input.getCmdOption("--subsample-depth");
    std::string subsample_replacement_arg = input.getCmdOption("--subsample-replacement");
    std::string n_subsamples_arg = input.getCmdOption("--n-subsamples");
    std::string diskbuf_arg = input.getCmdOption("--diskbuf");

    if(nsubsteps_arg.empty()) {
        nsubsteps = 1;
    } else {
        nsubsteps = atoi(nsubsteps_arg.c_str());
    }
    
    bool vaw = input.cmdOptionExists("--vaw"); 
    bool bare = input.cmdOptionExists("--report-bare"); 
    bool bypass_tips = input.cmdOptionExists("-f");
    double g_unifrac_alpha;

    if(gunifrac_arg.empty()) {
        g_unifrac_alpha = 1.0;
    } else {
        g_unifrac_alpha = atof(gunifrac_arg.c_str());
    }

    int start_stripe;
    if(start_arg.empty()) 
        start_stripe = 0;
    else
        start_stripe = atoi(start_arg.c_str());

    int stop_stripe;
    if(stop_arg.empty()) 
        stop_stripe = 0;
    else
        stop_stripe = atoi(stop_arg.c_str());

    int n_partials;
    if(npartials.empty()) 
        n_partials = 1;
    else
        n_partials = atoi(npartials.c_str());

    if(n_partials<1) {
        err("--n-partials cannot be < 1");
        return EXIT_FAILURE;
    }
     if(n_partials>1000000000) {
        err("--n-partials cannot be > 1G");
        return EXIT_FAILURE;
    }

    bool normalize_sample_counts;
    if(normsc_arg.empty()) {
        normalize_sample_counts = true;
    } else if (normsc_arg=="true") {
        normalize_sample_counts = true;
    } else if (normsc_arg=="false") {
        normalize_sample_counts = false;
    } else {
        err("Invalid normalize-sample-counts, must be true or false");
        return EXIT_FAILURE;
    }

    Format format_val = format_invalid;
    if(!format_arg.empty()) {
      format_val = get_format(format_arg,method_string,mode_arg);
    } else {
      format_val = get_format(sformat_arg,method_string,mode_arg);
      format_arg=sformat_arg; // easier to use a single variable
    }
    if(format_val==format_invalid) {
        err("Invalid format, must be one of ascii|hdf5|hdf5_fp32|hdf5_fp64|hdf5_nodist");
        return EXIT_FAILURE;
    }

    unsigned int pcoa_dims;
    if(pcoa_arg.empty())
        pcoa_dims = 10;
    else
        pcoa_dims = atoi(pcoa_arg.c_str());

    unsigned int permanova_perms;
    if(permanova_arg.empty()) {
        if (grouping_filename.empty() || grouping_columns.empty()) {
            // cannot compute permanova without the grouping file and columns
            permanova_perms = 0;
        } else {
            permanova_perms = 999;
        }
    } else {
        permanova_perms = atoi(permanova_arg.c_str());
    }

    if(!seed_arg.empty()) {
         ssu_set_random_seed(atoi(seed_arg.c_str()));
    }

    unsigned int subsample_depth = 0;
    if(subsample_depth_arg.empty()) {
        if(mode_arg == "multi" || mode_arg == "multiple") {
           err("--subsample-depth required in multi mode.");
           return EXIT_FAILURE;
        }
    } else {
        subsample_depth = atoi(subsample_depth_arg.c_str());
    }

    bool subsample_without_replacement = false;
    if(!subsample_replacement_arg.empty()) {
        if (subsample_replacement_arg == "with") {
           subsample_without_replacement = false;
        } else if (subsample_replacement_arg == "without") {
           subsample_without_replacement = true;
        } else {
           err("Invalid --subsample-replacement argument, must be 'with' or 'without'.");
           return EXIT_FAILURE;
        }
    }

    unsigned int n_subsamples = 0;
    if(n_subsamples_arg.empty()) {
        n_subsamples = 100;
    } else {
        if(mode_arg == "multi" || mode_arg == "multiple") {
           n_subsamples = atoi(n_subsamples_arg.c_str());
        } else {
           err("--n-subsamples only allowed in multi mode.");
           return EXIT_FAILURE;
        }
    }

    if(mode_arg.empty() || mode_arg == "one-off")
        return mode_one_off(table_filename, tree_filename, output_filename,  format2str(format_val), format_val, method_string,
                            subsample_depth, !subsample_without_replacement,
                            pcoa_dims, permanova_perms, grouping_filename, grouping_columns,
                            vaw, g_unifrac_alpha, bypass_tips, normalize_sample_counts, nsubsteps, diskbuf_arg);
    else if(mode_arg == "partial") {
        if (subsample_depth>0) {
          err("Cannot subsample in partial mode.");
          return EXIT_FAILURE;
        }
        return mode_partial(table_filename, tree_filename, output_filename, method_string, vaw, g_unifrac_alpha,
			    bypass_tips, normalize_sample_counts, nsubsteps, start_stripe, stop_stripe);
    } else if(mode_arg == "merge-partial")
        return mode_merge_partial(output_filename, format_val,
                                  pcoa_dims, permanova_perms, grouping_filename, grouping_columns,
                                  partial_pattern, diskbuf_arg);
    else if(mode_arg == "check-partial")
        return mode_check_partial(partial_pattern);
    else if(mode_arg == "partial-report")
        return mode_partial_report(table_filename, uint32_t(n_partials), bare);
    else if(mode_arg == "multi" || mode_arg == "multiple")
        return mode_multi(table_filename, tree_filename, output_filename, format2str(format_val), format_val, method_string,
                            n_subsamples,subsample_depth, !subsample_without_replacement,
                            pcoa_dims, permanova_perms, grouping_filename, grouping_columns,
                            vaw, g_unifrac_alpha, bypass_tips, normalize_sample_counts, nsubsteps, diskbuf_arg);
    else 
        err("Unknown mode. Valid options are: one-off, partial, merge-partial, check-partial, partial-report, multi");

    return EXIT_SUCCESS;
}

