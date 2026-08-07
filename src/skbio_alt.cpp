/*
 * Classes, methods and unction that provide skbio-like unctionality
 */

#include "skbio_alt.hpp"
#include <stdlib.h> 

#include <random>
#include <algorithm>

#include <scikit-bio-binaries/util.h>
#include <scikit-bio-binaries/ordination.h>
#include <scikit-bio-binaries/distance.h>

static std::mt19937 myRandomGenerator;

static constexpr int ACC_CPU=0;
static constexpr int ACC_NV=1;
static constexpr int ACC_AMD=2;

#ifdef UNIFRAC_ENABLE_ACC_NV
#define SUCMP_NM  su_acc_nv
#include "unifrac_accapi.hpp"
#undef SUCMP_NM
#endif

#ifdef UNIFRAC_ENABLE_ACC_AMD
#define SUCMP_NM  su_acc_amd
#include "unifrac_accapi.hpp"
#undef SUCMP_NM
#endif


void su::set_random_seed(uint32_t new_seed) {
  myRandomGenerator.seed(new_seed);
  // propagate to the dependency, too
  auto new_seed_skbb = myRandomGenerator();
  skbb_set_random_seed(new_seed_skbb);
}

// test only once, then use persistent value
static int skbio_use_acc = -1;

inline void skbio_check_acc() {
 if (skbio_use_acc!=-1) return; // keep the cached version

 bool print_info = false;

 if (const char* env_p = std::getenv("UNIFRAC_GPU_INFO")) {
   print_info = true;
   std::string env_s(env_p);
   if ((env_s=="NO") || (env_s=="N") || (env_s=="no") || (env_s=="n") ||
       (env_s=="NEVER") || (env_s=="never")) {
     print_info = false;
   }
   // if our library wants reporting, force it at dependency level, too
   if(print_info) skbb_set_acc_reporting_flag(print_info);
 }

 int detected_acc = ACC_CPU;
#if defined(UNIFRAC_ENABLE_ACC_NV)
 bool detected_nv_acc = su_acc_nv::acc_found_gpu();
 if (print_info) {
   if (detected_nv_acc) {
     printf("INFO (skbio_alt): NVIDIA GPU detected\n");
   } else {
     printf("INFO (skbio_alt): NVIDIA GPU not detected\n");
   }
 }
 if ((detected_acc==ACC_CPU) && detected_nv_acc) {
   detected_acc = ACC_NV;
   if (const char* env_p = std::getenv("UNIFRAC_USE_NVIDIA_GPU")) {
     std::string env_s(env_p);
     if ((env_s=="NO") || (env_s=="N") || (env_s=="no") || (env_s=="n") ||
         (env_s=="NEVER") || (env_s=="never")) {
       if (print_info) printf("INFO (skbio_alt): NVIDIA GPU was detected but use explicitly disabled\n");
       detected_acc = ACC_CPU;
     }
   }
 }
#endif
#if defined(UNIFRAC_ENABLE_ACC_AMD)
 bool detected_amd_acc = su_acc_amd::acc_found_gpu();
 if (print_info) {
   if (detected_amd_acc) {
     printf("INFO (skbio_alt): AMD GPU detected\n");
   } else {
     printf("INFO (skbio_alt): AMD GPU not detected\n");
   }
 }
 if ((detected_acc==ACC_CPU) && detected_amd_acc) {
   detected_acc = ACC_AMD;
   if (const char* env_p = std::getenv("UNIFRAC_USE_AMD_GPU")) {
     std::string env_s(env_p);
     if ((env_s=="NO") || (env_s=="N") || (env_s=="no") || (env_s=="n") ||
         (env_s=="NEVER") || (env_s=="never")) {
       if (print_info) printf("INFO (skbio_alt): AMD GPU was detected but use explicitly disabled\n");
       detected_acc = ACC_CPU;
     }
   }
 }
#endif

 if (const char* env_p = std::getenv("UNIFRAC_USE_GPU")) {
   std::string env_s(env_p);
   if ((env_s=="NO") || (env_s=="N") || (env_s=="no") || (env_s=="n") ||
       (env_s=="NEVER") || (env_s=="never")) {
     if (detected_acc!=ACC_CPU) {
        if (print_info) printf("INFO (skbio_alt): GPU was detected but use explicitly disabled\n");
         detected_acc = ACC_CPU;
     }
   }
 }

 if (print_info) {
   if (detected_acc == ACC_CPU) {
     printf("INFO (skbio_alt): Using CPU only (not GPU)\n");
#if defined(UNIFRAC_ENABLE_ACC_NV)
   } else if (detected_acc == ACC_NV) {
     printf("INFO (skbio_alt): Using NVIDIA GPU if feasible\n");
#endif
#if defined(UNIFRAC_ENABLE_ACC_AMD)
   } else if (detected_acc == ACC_AMD) {
     printf("INFO (skbio_alt): Using AMD GPU if feasible \n");
#endif
   }
 }

 if (detected_acc == ACC_CPU) {
   // if our library should not use a GPU, force it at dependency level, too
   if(skbb_set_acc_mode(SKBB_ACC_CPU)){} // ignore eny errors, it is justy advisory
 }

 // we can assume int is atomic
 skbio_use_acc = detected_acc;
}

//
// ======================= PCoA ========================
//

// === Helper functions

void su::mat_to_centered(const double * mat, const uint32_t n_samples, double * centered) {
  skbio_check_acc();
  skbb_center_distance_matrix_fp64(n_samples, mat, centered);
}

void su::mat_to_centered(const float  * mat, const uint32_t n_samples, float  * centered) {
  skbio_check_acc();
  skbb_center_distance_matrix_fp32(n_samples, mat, centered);
}

void su::mat_to_centered(const double * mat, const uint32_t n_samples, float  * centered) {
  skbio_check_acc();
  skbb_center_distance_matrix_fp64_to_fp32(n_samples, mat, centered);
}

void su::find_eigens_fast(const uint32_t n_samples, const uint32_t n_dims, double * centered, double * &eigenvalues, double * &eigenvectors) {
  eigenvalues = (double *) malloc(sizeof(double)*n_dims);
  eigenvectors = (double *) malloc((sizeof(double)*n_dims)*n_samples);
  skbio_check_acc();
  skbb_fsvd_inplace_fp64(n_samples, centered, n_dims, -1, eigenvalues, eigenvectors);
}

void su::find_eigens_fast(const uint32_t n_samples, const uint32_t n_dims, float * centered, float * &eigenvalues, float * &eigenvectors) {

  eigenvalues = (float *) malloc(sizeof(float)*n_dims);
  eigenvectors = (float *) malloc((sizeof(float)*n_dims)*n_samples);
  skbio_check_acc();
  skbb_fsvd_inplace_fp32(n_samples, centered, n_dims, -1, eigenvalues, eigenvectors);
}

// ======================= PCoA proper ========================

void su::pcoa(const double * mat, const uint32_t n_samples, const uint32_t n_dims, double * &eigenvalues, double * &samples, double * &proportion_explained, int seed) {
  eigenvalues = (double *) malloc(sizeof(double)*n_dims);
  samples = (double *) malloc((sizeof(double)*n_dims)*n_samples);
  proportion_explained = (double *) malloc(sizeof(double)*n_dims);
  skbio_check_acc();
  skbb_pcoa_fsvd_fp64(n_samples, mat, n_dims, seed, eigenvalues, samples, proportion_explained);
}

void su::pcoa(const float  * mat, const uint32_t n_samples, const uint32_t n_dims, float  * &eigenvalues, float  * &samples, float  * &proportion_explained, int seed) {
  eigenvalues = (float *) malloc(sizeof(float)*n_dims);
  samples = (float *) malloc((sizeof(float)*n_dims)*n_samples);
  proportion_explained = (float *) malloc(sizeof(float)*n_dims);
  skbio_check_acc();
  skbb_pcoa_fsvd_fp32(n_samples, mat, n_dims, seed, eigenvalues, samples, proportion_explained);
}

void su::pcoa(const double * mat, const uint32_t n_samples, const uint32_t n_dims, float  * &eigenvalues, float  * &samples, float  * &proportion_explained, int seed) {
  eigenvalues = (float *) malloc(sizeof(float)*n_dims);
  samples = (float *) malloc((sizeof(float)*n_dims)*n_samples);
  proportion_explained = (float *) malloc(sizeof(float)*n_dims);
  skbio_check_acc();
  skbb_pcoa_fsvd_fp64_to_fp32(n_samples, mat, n_dims, seed, eigenvalues, samples, proportion_explained);
}

void su::pcoa_inplace(double * mat, const uint32_t n_samples, const uint32_t n_dims, double * &eigenvalues, double * &samples, double * &proportion_explained, int seed) {
  eigenvalues = (double *) malloc(sizeof(double)*n_dims);
  samples = (double *) malloc((sizeof(double)*n_dims)*n_samples);
  proportion_explained = (double *) malloc(sizeof(double)*n_dims);
  skbio_check_acc();
  skbb_pcoa_fsvd_inplace_fp64(n_samples, mat, n_dims, seed, eigenvalues, samples, proportion_explained);
}

void su::pcoa_inplace(float  * mat, const uint32_t n_samples, const uint32_t n_dims, float  * &eigenvalues, float  * &samples, float  * &proportion_explained, int seed) {
  eigenvalues = (float *) malloc(sizeof(float)*n_dims);
  samples = (float *) malloc((sizeof(float)*n_dims)*n_samples);
  proportion_explained = (float *) malloc(sizeof(float)*n_dims);
  skbio_check_acc();
  skbb_pcoa_fsvd_inplace_fp32(n_samples, mat, n_dims, seed, eigenvalues, samples, proportion_explained);
}

//
// ======================= permanova ========================
//

void su::permanova(const double * mat, unsigned int n_dims,
                   const uint32_t *grouping,
                   unsigned int n_perm,
                   double &fstat_out, double &pvalue_out,
                   int seed) {
  skbio_check_acc();
  skbb_permanova_fp64(n_dims, mat, grouping, n_perm, seed, &fstat_out, &pvalue_out);
}

void su::permanova(const float * mat, unsigned int n_dims,
                   const uint32_t *grouping,
                   unsigned int n_perm,
                   float &fstat_out, float &pvalue_out,
                   int seed) {
  skbio_check_acc();
  skbb_permanova_fp32(n_dims, mat, grouping, n_perm, seed, &fstat_out, &pvalue_out);
}

// ======================= skbio_biom_subsampled  ================================


su::skbio_biom_subsampled::skbio_biom_subsampled(const biom_inmem &parent, const bool w_replacement, const uint32_t n)
 : su::biom_subsampled(parent, w_replacement, n, uint32_t(myRandomGenerator()))
{}

