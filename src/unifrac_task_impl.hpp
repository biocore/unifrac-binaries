/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2025, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

#include <math.h> 
#include <algorithm> 
#include "unifrac_task_noclass.hpp"
#include <cstdlib>

#if defined(OMPGPU)

#include <omp.h>

// Biggest vector size supported by the GPU
#define SUCMP_ACC_MAXVEC 2048

#elif defined(_OPENACC)

#include <openacc.h>

// Biggest vector size supported by the GPU
#define SUCMP_ACC_MAXVEC 2048

#endif

// Use a moderate sized step, a few cache lines
#define STEP_SIZE(TFloat) (64*4/sizeof(TFloat))

static inline bool acc_found_gpu_T() {
#if defined(OMPGPU)
  return omp_get_num_devices() > 0;
#elif defined(_OPENACC)
  return acc_get_device_type() != acc_device_host;
#else
  return false;
#endif
}

// is the implementation async, and need the alt structures?
static inline bool acc_need_alt_T() {
#if defined(_OPENACC) || defined(OMPGPU)
   return true;
#else
   return false;
#endif
}

static inline void acc_wait_T() {
#if defined(OMPGPU)
    // TODO: Change if we ever implement async in OMPGPU
#elif defined(_OPENACC)
#pragma acc wait
#endif
}

template<class TNum>
static inline void acc_create_buf_T(
		TNum *buf,
		uint64_t start, uint64_t end) {
#if defined(OMPGPU)
#pragma omp target enter data map(alloc:buf[start:end])
#elif defined(_OPENACC)
#pragma acc enter data create(buf[start:end])
#endif
}

template<class TNum>
static inline void acc_copyin_buf_T(
		TNum *buf,
		uint64_t start, uint64_t end) {
#if defined(OMPGPU)
#pragma omp target enter data map(to:buf[start:end])
#elif defined(_OPENACC)
#pragma acc enter data copyin(buf[start:end])
#endif    
}

template<class TNum>
static inline void acc_update_device_T(
		TNum *buf,
		uint64_t start, uint64_t end) {
#if defined(OMPGPU)
#pragma omp target update to(buf[start:end])
#elif defined(_OPENACC)
#pragma acc update device(buf[start:end])
#endif
}

template<class TNum>
static inline void acc_copyout_buf_T(
		TNum *buf,
		uint64_t start, uint64_t end) {
#if defined(OMPGPU)
#pragma omp target exit data map(from:buf[start:end])
#elif defined(_OPENACC)
#pragma acc exit data copyout(buf[start:end])
#endif
}

template<class TNum>
static inline void acc_destroy_buf_T(
		TNum *buf,
		uint64_t start, uint64_t end) {
#if defined(OMPGPU)
#pragma omp target exit data map(delete:buf[start:end])
#elif defined(_OPENACC)
#pragma acc exit data delete(buf[start:end])
#endif
}

template<class TFloat>
static inline void compute_stripes_totals_T(
		TFloat * const __restrict__ dm_stripes_buf,
		const TFloat * const __restrict__ dm_stripes_total_buf,
		const uint64_t bufels) { 
#if defined(OMPGPU)
#pragma omp target teams distribute parallel for simd default(shared)
#elif defined(_OPENACC)
#pragma acc parallel loop gang vector present(dm_stripes_buf,dm_stripes_total_buf)
#endif
   for(uint64_t idx=0; idx< bufels; idx++)
       dm_stripes_buf[idx]=dm_stripes_buf[idx]/dm_stripes_total_buf[idx];
}


// check for zero values and pre-compute single column sums
template<class TFloat>
static inline void WeightedZerosAndSums(
                      bool   * const __restrict__ zcheck,
                      TFloat * const __restrict__ sums,
                      const TFloat * const __restrict__ embedded_proportions,
                      const TFloat * __restrict__ lengths,
                      const uint64_t embs_stripe,
                      const unsigned int filled_embs,
                      const uint64_t n_samples,
                      const uint64_t n_samples_r) {
#if !(defined(_OPENACC) || defined(OMPGPU))
// CPU version uses transposed embedded_proportions
#pragma omp parallel for default(shared)
    for(uint64_t k=0; k<n_samples; k++) {
            bool all_zeros=true;
            const uint64_t offset = embs_stripe*k;
            TFloat my_sum = 0.0;

            for (uint64_t emb=0; emb<filled_embs; emb++) {
                TFloat u1 = embedded_proportions[offset + emb];
                my_sum += u1*lengths[emb];
                all_zeros = all_zeros && (u1==0.0);
            }

            sums[k]     = my_sum;
            zcheck[k] = all_zeros;
    }
#else
// GPU version uses straight embedded_proportions
#if defined(OMPGPU)
#pragma omp target teams distribute parallel for simd default(shared)
#else
#pragma acc parallel loop gang vector present(embedded_proportions,lengths,zcheck,sums)
#endif
    for(uint64_t k=0; k<n_samples; k++) {
            bool all_zeros=true;
            TFloat my_sum = 0.0;

#if !defined(OMPGPU) && defined(_OPENACC)
#pragma acc loop seq
#endif
            for (uint64_t emb=0; emb<filled_embs; emb++) {
                const uint64_t offset = n_samples_r * emb;

                TFloat u1 = embedded_proportions[offset + k];
                my_sum += u1*lengths[emb];
                all_zeros = all_zeros && (u1==0.0);
            }

            sums[k]     = my_sum;
            zcheck[k] = all_zeros;
    }
#endif
}

// Single step in computing Weighted part of Unifrac
template<class TFloat>
static inline TFloat WeightedVal1(
                      const TFloat * const __restrict__ embedded_proportions,
                      const TFloat * __restrict__ lengths,
                      const uint64_t embs_stripe,
                      const unsigned int filled_embs,
                      const uint64_t n_samples_r,
                      const uint64_t k,
                      const uint64_t l1) {
            TFloat my_stripe = 0.0;
#if !(defined(_OPENACC) || defined(OMPGPU))
// CPU version uses transposed embedded_proportions
            const uint64_t offset_k = embs_stripe*k;
            const uint64_t offset_l = embs_stripe*l1;

#pragma omp simd reduction(+:my_stripe)
            for (uint64_t emb=0; emb<filled_embs; emb++) {
                TFloat u1 = embedded_proportions[offset_k + emb];
                TFloat v1 = embedded_proportions[offset_l + emb];
                TFloat diff1 = u1 - v1;
                TFloat length = lengths[emb];

                my_stripe     += fabs(diff1) * length;
            } // for emb

#else
// GPU sticks with the straight layout
            for (uint64_t emb=0; emb<filled_embs; emb++) {
                const uint64_t offset = n_samples_r * emb;

                TFloat u1 = embedded_proportions[offset + k];
                TFloat v1 = embedded_proportions[offset + l1];
                TFloat diff1 = u1 - v1;
                TFloat length = lengths[emb];

                my_stripe     += fabs(diff1) * length;
            } // for emb
#endif

            return my_stripe;
}


// Single step in computing NormalizedWeighted Unifrac
template<class TFloat>
static inline void UnnormalizedWeighted1(
                      TFloat * const __restrict__ dm_stripes_buf,
                      const bool   * const __restrict__ zcheck,
                      const TFloat * const __restrict__ sums,
                      const TFloat * const __restrict__ embedded_proportions,
                      const TFloat * __restrict__ lengths,
                      const uint64_t embs_stripe,
                      const unsigned int filled_embs,
                      const uint64_t idx,
                      const uint64_t n_samples_r,
                      const uint64_t k,
                      const uint64_t l1) {
       const bool allzero_k = zcheck[k];
       const bool allzero_l1 = zcheck[l1];

       if (allzero_k && allzero_l1) {
         // nothing to do, would have to add 0
       } else {
          TFloat my_stripe;

          if (allzero_k || allzero_l1) {
            // one side has all zeros
            // we can use the distributed property, and use the pre-computed values

            const uint64_t ridx = (allzero_k) ? l1 : // if (nonzero_l1) ridx=l1 // fabs(k-l1), with k==0
                                                k;   // if (nonzero_k)  ridx=k  // fabs(k-l1), with l1==0

            // keep reads in the same place to maximize GPU warp performance
            my_stripe = sums[ridx];
          } else {
            // both sides non zero, use the explicit but slow approach
            my_stripe = WeightedVal1(embedded_proportions, lengths,
                                     embs_stripe, filled_embs, n_samples_r,
                                     k, l1);
          }

          TFloat * const __restrict__ dm_stripe = dm_stripes_buf+idx;
          //TFloat *dm_stripe = dm_stripes[stripe];
          
          // keep all writes in a single place, to maximize GPU warp performance
          dm_stripe[k]       += my_stripe;
       }
}

#if !(defined(_OPENACC) || defined(OMPGPU))
// Vectorized step in computing UnnormalizedWeighted Unifrac
template<class TFloat>
static inline void UnnormalizedWeighted4(
                      TFloat * const __restrict__ dm_stripes_buf,
                      const bool   * const __restrict__ zcheck,
                      const TFloat * const __restrict__ sums,
                      const TFloat * const __restrict__ embedded_proportions,
                      const TFloat * __restrict__ lengths,
                      const uint64_t embs_stripe,
                      const unsigned int filled_embs,
                      const uint64_t idx,
                      const uint64_t n_samples_r,
                      const uint64_t ks,
                      const uint64_t ls) {
       const uint32_t z_k = ((const uint32_t *)(zcheck+ks))[0];
       const uint32_t z_l = ((const uint32_t *)(zcheck+ls))[0];
       const bool allzero_k = z_k==0x01010101;
       const bool allzero_l = z_l==0x01010101;

       if (allzero_k && allzero_l) {
         // nothing to do, would have to add 0
       } else if (allzero_k || allzero_l) {
          TFloat * const __restrict__ dm_stripe = dm_stripes_buf+idx;
          //TFloat *dm_stripe = dm_stripes[stripe];

          if (allzero_k) {
             // one side has all zeros, use distributed property
             // if (nonzero_l1) ridx=l1 // fabs(k-l1), with k==0
             const TFloat sum_l0 = sums[ls];
             const TFloat sum_l1 = sums[ls+1];
             const TFloat sum_l2 = sums[ls+2];
             const TFloat sum_l3 = sums[ls+3];
             dm_stripe[ks]   += sum_l0;
             dm_stripe[ks+1] += sum_l1;
             dm_stripe[ks+2] += sum_l2;
             dm_stripe[ks+3] += sum_l3;
          } else {
             // one side has all zeros, use distributed property
             // if (nonzero_k)  ridx=k  // fabs(k-l1), with l1==0
             const TFloat sum_k0 = sums[ks];
             const TFloat sum_k1 = sums[ks+1];
             const TFloat sum_k2 = sums[ks+2];
             const TFloat sum_k3 = sums[ks+3];
             dm_stripe[ks]   += sum_k0;
             dm_stripe[ks+1] += sum_k1;
             dm_stripe[ks+2] += sum_k2;
             dm_stripe[ks+3] += sum_k3;
          }
       } else {
         // only a few have both sides partially non zero, use the fine grained compute
         for (uint64_t i=0; i<4; i++) {
            UnnormalizedWeighted1<TFloat>(
                                   dm_stripes_buf,
                                   zcheck, sums, embedded_proportions, lengths,
                                   embs_stripe,filled_embs,idx, n_samples_r,
                                   ks+i, ls+i);
         }
       } // (allzero_k && allzero_l)
}

template<class TFloat>
static inline void UnnormalizedWeighted8(
                      TFloat * const __restrict__ dm_stripes_buf,
                      const bool   * const __restrict__ zcheck,
                      const TFloat * const __restrict__ sums,
                      const TFloat * const __restrict__ embedded_proportions,
                      const TFloat * __restrict__ lengths,
                      const uint64_t embs_stripe,
                      const unsigned int filled_embs,
                      const uint64_t idx,
                      const uint64_t n_samples_r,
                      const uint64_t ks,
                      const uint64_t ls) {
       const uint64_t z_k = ((const uint64_t *)(zcheck+ks))[0];
       const uint64_t z_l = ((const uint64_t *)(zcheck+ls))[0];
       const bool allzero_k = z_k==0x0101010101010101;
       const bool allzero_l = z_l==0x0101010101010101;

       if (allzero_k && allzero_l) {
         // nothing to do, would have to add 0
       } else if (allzero_k || allzero_l) {
          TFloat * const __restrict__ dm_stripe = dm_stripes_buf+idx;
          //TFloat *dm_stripe = dm_stripes[stripe];

          if (allzero_k) {
             // one side has all zeros, use distributed property
             // if (nonzero_l1) ridx=l1 // fabs(k-l1), with k==0
             const TFloat sum_l0 = sums[ls];
             const TFloat sum_l1 = sums[ls+1];
             const TFloat sum_l2 = sums[ls+2];
             const TFloat sum_l3 = sums[ls+3];
             const TFloat sum_l4 = sums[ls+4];
             const TFloat sum_l5 = sums[ls+5];
             const TFloat sum_l6 = sums[ls+6];
             const TFloat sum_l7 = sums[ls+7];
             dm_stripe[ks]   += sum_l0;
             dm_stripe[ks+1] += sum_l1;
             dm_stripe[ks+2] += sum_l2;
             dm_stripe[ks+3] += sum_l3;
             dm_stripe[ks+4] += sum_l4;
             dm_stripe[ks+5] += sum_l5;
             dm_stripe[ks+6] += sum_l6;
             dm_stripe[ks+7] += sum_l7;
          } else {
             // one side has all zeros, use distributed property
             // if (nonzero_k)  ridx=k  // fabs(k-l1), with l1==0
             const TFloat sum_k0 = sums[ks];
             const TFloat sum_k1 = sums[ks+1];
             const TFloat sum_k2 = sums[ks+2];
             const TFloat sum_k3 = sums[ks+3];
             const TFloat sum_k4 = sums[ks+4];
             const TFloat sum_k5 = sums[ks+5];
             const TFloat sum_k6 = sums[ks+6];
             const TFloat sum_k7 = sums[ks+7];
             dm_stripe[ks]   += sum_k0;
             dm_stripe[ks+1] += sum_k1;
             dm_stripe[ks+2] += sum_k2;
             dm_stripe[ks+3] += sum_k3;
             dm_stripe[ks+4] += sum_k4;
             dm_stripe[ks+5] += sum_k5;
             dm_stripe[ks+6] += sum_k6;
             dm_stripe[ks+7] += sum_k7;
          } 
       } else {
         // only a few have both sides partially non zero, use the finer grained compute
         for (uint64_t i=0; i<8; i+=4) {
            UnnormalizedWeighted4<TFloat>(
                                   dm_stripes_buf,
                                   zcheck, sums, embedded_proportions, lengths,
                                   embs_stripe, filled_embs,idx, n_samples_r,
                                   ks+i, ls+i);
         }
       } // (allzero_k && allzero_l)
}
#endif

template<class TFloat>
static inline void run_UnnormalizedWeightedTask_T(
		const uint64_t embs_stripe, const unsigned int filled_embs,
		const uint64_t start_idx, const uint64_t stop_idx,
		const uint64_t n_samples, const uint64_t n_samples_r,
		const TFloat * const __restrict__ lengths,
		const TFloat * const __restrict__ embedded_proportions,
		TFloat * const __restrict__ dm_stripes_buf,
		bool * const __restrict__ zcheck,
		TFloat * const __restrict__ sums) {

    constexpr uint64_t step_size = STEP_SIZE(TFloat);
    const uint64_t sample_steps = (n_samples+(step_size-1))/step_size; // round up

    // check for zero values and pre-compute single column sums
    WeightedZerosAndSums(zcheck, sums,
                         embedded_proportions, lengths,
                         embs_stripe, filled_embs, n_samples, n_samples_r);

    // now do the real compute
#if defined(_OPENACC) || defined(OMPGPU)
#if defined(OMPGPU)
    // TODO: Explore async omp target
#pragma omp target teams distribute parallel for simd collapse(3) default(shared)
#else
#pragma acc parallel loop gang vector collapse(3) present(embedded_proportions,dm_stripes_buf,lengths,zcheck,sums) async
#endif
    for(uint64_t sk = 0; sk < sample_steps ; sk++) {
     for(uint64_t stripe = start_idx; stripe < stop_idx; stripe++) {
      // SIMT-based GPU work great one at a time (HW will deal with parallelism)
      for(uint64_t ik = 0; ik < step_size ; ik++) {
       const uint64_t k = sk*step_size + ik;

       if (k>=n_samples) continue; // past the limit

       const uint64_t l1 = (k + stripe + 1)%n_samples; // wraparound
       const uint64_t idx = (stripe-start_idx) * n_samples_r;

       UnnormalizedWeighted1<TFloat>(
                                   dm_stripes_buf,
                                   zcheck, sums, embedded_proportions, lengths,
                                   embs_stripe, filled_embs,idx, n_samples_r,
                                   k, l1);
      } // for ik
     } // for stripe
    } // for sk
#else
    // tilling helps with better cache reuse without the need of multiple cores
    const uint64_t stripe_steps = ((stop_idx-start_idx)+(step_size-1))/step_size; // round up

    // use dynamic scheduling due to non-homogeneity in the loop
    // Use a moderate block to prevent trashing but still have some cache reuse
#pragma omp parallel for collapse(2) schedule(dynamic,step_size) default(shared)
    for(uint64_t ss = 0; ss < stripe_steps ; ss++) {
     for(uint64_t sk = 0; sk < sample_steps ; sk++) {
      // tile to maximize cache reuse
      for(uint64_t is = 0; is < step_size ; is++) {
       const uint64_t stripe = start_idx+ss*step_size + is;
       if (stripe<stop_idx) { // else past limit

      // SIMD-based CPUs need help with vectorization
      const uint64_t idx = (stripe-start_idx) * n_samples_r;
      uint64_t ks = sk*step_size;
      const uint64_t kmax = std::min(ks+step_size,n_samples);
      uint64_t ls = (ks + stripe + 1)%n_samples; // wraparound

      while( ((ks+8) <= kmax) && ((n_samples-ls)>=8) ) {
       UnnormalizedWeighted8<TFloat>(
                                   dm_stripes_buf,
                                   zcheck, sums, embedded_proportions, lengths,
                                   embs_stripe, filled_embs,idx, n_samples_r,
                                   ks, ls);
       ks+=8;
       ls = (ls + 8)%n_samples; // wraparound
      } // for ks+=8

      while( ((ks+4) <= kmax) && ((n_samples-ls)>=4) ) {
       UnnormalizedWeighted4<TFloat>(
                                   dm_stripes_buf,
                                   zcheck, sums, embedded_proportions, lengths,
                                   embs_stripe, filled_embs,idx, n_samples_r,
                                   ks, ls);
       ks+=4;
       ls = (ls + 4)%n_samples; // wraparound
      } // for ks+=4

      // deal with any leftovers in serial way
      for( ; ks < kmax; ks++ ) {
       UnnormalizedWeighted1<TFloat>(
                                   dm_stripes_buf,
                                   zcheck, sums, embedded_proportions, lengths,
                                   embs_stripe, filled_embs,idx, n_samples_r,
                                   ks, ls);
       ls = (ls + 1)%n_samples; // wraparound
      } // for ks

       } // if stripe
      } // for is
     } // for sk
    } // for ss
#endif
}

template<class TFloat>
static inline void run_VawUnnormalizedWeightedTask_T(
		const unsigned int filled_embs,
		const uint64_t start_idx, const uint64_t stop_idx,
		const uint64_t n_samples, const uint64_t n_samples_r,
		const TFloat * const __restrict__ lengths,
		const TFloat * const __restrict__ embedded_proportions,
		const TFloat * const __restrict__ embedded_counts,
		const TFloat * const __restrict__ sample_total_counts,
		TFloat * const __restrict__ dm_stripes_buf) {

    constexpr uint64_t step_size = STEP_SIZE(TFloat);
    const uint64_t sample_steps = (n_samples+(step_size-1))/step_size; // round up

    // point of thread
#if defined(_OPENACC) || defined(OMPGPU)
    // Use as big vector size as we can, to maximize cache line reuse
#if defined(OMPGPU)
    // TODO: Explore async omp target
#pragma omp target teams distribute parallel for simd collapse(3) simdlen( SUCMP_ACC_MAXVEC ) default(shared)
#else
#pragma acc parallel loop collapse(3) vector_length( SUCMP_ACC_MAXVEC ) present(embedded_proportions,embedded_counts,sample_total_counts,dm_stripes_buf,lengths) async
#endif

#else
#pragma omp parallel for default(shared) schedule(dynamic,1)
#endif
    for(uint64_t sk = 0; sk < sample_steps ; sk++) {
      for(uint64_t stripe = start_idx; stripe < stop_idx; stripe++) {
        for(uint64_t ik = 0; ik < step_size ; ik++) {
            const uint64_t k = sk*step_size + ik;
            const uint64_t idx = (stripe-start_idx) * n_samples_r;
            TFloat * const __restrict__ dm_stripe = dm_stripes_buf+idx;
            //TFloat *dm_stripe = dm_stripes[stripe];

            if (k>=n_samples) continue; // past the limit

            const uint64_t l1 = (k + stripe + 1)%n_samples; // wraparound

            const TFloat m = sample_total_counts[k] + sample_total_counts[l1];

            TFloat my_stripe = dm_stripe[k];

#if !defined(OMPGPU) && defined(_OPENACC)
#pragma acc loop seq
#endif
            for (uint64_t emb=0; emb<filled_embs; emb++) {
                const uint64_t offset = n_samples_r*emb;

                TFloat mi = embedded_counts[offset + k] + embedded_counts[offset + l1];
                TFloat vaw = sqrt(mi * (m - mi));

                if(vaw > 0) {
                  TFloat u1 = embedded_proportions[offset + k];
                  TFloat v1 = embedded_proportions[offset + l1];
                  TFloat diff1 = fabs(u1 - v1);
                  TFloat length = lengths[emb];

                   my_stripe += (diff1 * length) / vaw;
                }
            }

            dm_stripe[k]     = my_stripe;
        }

      }
    }
}

// Single step in computing NormalizedWeighted Unifrac
template<class TFloat>
static inline void NormalizedWeighted1(
                      TFloat * const __restrict__ dm_stripes_buf,
                      TFloat * const __restrict__ dm_stripes_total_buf,
                      const bool   * const __restrict__ zcheck,
                      const TFloat * const __restrict__ sums,
                      const TFloat * const __restrict__ embedded_proportions,
                      const TFloat * __restrict__ lengths,
                      const uint64_t embs_stripe,
                      const unsigned int filled_embs,
                      const uint64_t idx,
                      const uint64_t n_samples_r,
                      const uint64_t k,
                      const uint64_t l1) {
       const bool allzero_k = zcheck[k];
       const bool allzero_l1 = zcheck[l1];

       if (allzero_k && allzero_l1) {
         // nothing to do, would have to add 0
       } else {
          TFloat * const __restrict__ dm_stripe = dm_stripes_buf+idx;
          TFloat * const __restrict__ dm_stripe_total = dm_stripes_total_buf+idx;
          //TFloat *dm_stripe = dm_stripes[stripe];
          //TFloat *dm_stripe_total = dm_stripes_total[stripe];

          const TFloat sum_k = sums[k];
          const TFloat sum_l = sums[l1];

          // the totals can always use the distributed property
          dm_stripe_total[k] += sum_k + sum_l;
   
          TFloat my_stripe;

          if (allzero_k || allzero_l1) {
            // one side has all zeros
            // we can use the distributed property, and use the pre-computed values

            my_stripe = (allzero_k) ? sum_l : // if (nonzero_l1) ridx=l1 // fabs(k-l1), with k==0
                                      sum_k;   // if (nonzero_k)  ridx=k  // fabs(k-l1), with l1==0
          } else {
            // both sides non zero, use the explicit but slow approach
            my_stripe = WeightedVal1(embedded_proportions, lengths,
                                     embs_stripe, filled_embs, n_samples_r,
                                     k, l1);
          }

          // keep all writes in a single place, to maximize GPU warp performance
          dm_stripe[k]       += my_stripe;
       }
}

#if !(defined(_OPENACC) || defined(OMPGPU))
// Vectorized step in computing NormalizedWeighted Unifrac
template<class TFloat>
static inline void NormalizedWeighted4(
                      TFloat * const __restrict__ dm_stripes_buf,
                      TFloat * const __restrict__ dm_stripes_total_buf,
                      const bool   * const __restrict__ zcheck,
                      const TFloat * const __restrict__ sums,
                      const TFloat * const __restrict__ embedded_proportions,
                      const TFloat * __restrict__ lengths,
                      const uint64_t embs_stripe,
                      const unsigned int filled_embs,
                      const uint64_t idx,
                      const uint64_t n_samples_r,
                      const uint64_t ks,
                      const uint64_t ls) {
       const uint32_t z_k = ((const uint32_t *)(zcheck+ks))[0];
       const uint32_t z_l = ((const uint32_t *)(zcheck+ls))[0];
       const bool allzero_k = z_k==0x01010101;
       const bool allzero_l = z_l==0x01010101;

       if (allzero_k && allzero_l) {
         // nothing to do, would have to add 0
       } else {
          const TFloat sum_k0 = sums[ks];
          const TFloat sum_k1 = sums[ks+1];
          const TFloat sum_k2 = sums[ks+2];
          const TFloat sum_k3 = sums[ks+3];
          const TFloat sum_l0 = sums[ls];
          const TFloat sum_l1 = sums[ls+1];
          const TFloat sum_l2 = sums[ls+2];
          const TFloat sum_l3 = sums[ls+3];

          const TFloat sum_kl0 = sum_k0 + sum_l0;
          const TFloat sum_kl1 = sum_k1 + sum_l1;
          const TFloat sum_kl2 = sum_k2 + sum_l2;
          const TFloat sum_kl3 = sum_k3 + sum_l3;

          TFloat * const __restrict__ dm_stripe = dm_stripes_buf+idx;
          //TFloat *dm_stripe = dm_stripes[stripe];

          if (allzero_k) {
             // one side has all zeros, use distributed property
             // if (nonzero_l1) ridx=l1 // fabs(k-l1), with k==0
             dm_stripe[ks]   += sum_l0;
             dm_stripe[ks+1] += sum_l1;
             dm_stripe[ks+2] += sum_l2;
             dm_stripe[ks+3] += sum_l3;
          } else if (allzero_l) {
             // one side has all zeros, use distributed property
             // if (nonzero_k)  ridx=k  // fabs(k-l1), with l1==0
             dm_stripe[ks]   += sum_k0;
             dm_stripe[ks+1] += sum_k1;
             dm_stripe[ks+2] += sum_k2;
             dm_stripe[ks+3] += sum_k3;
          }

          // the totals can always use the distributed property
          {
             TFloat * const __restrict__ dm_stripe_total = dm_stripes_total_buf+idx;
             //TFloat *dm_stripe_total = dm_stripes_total[stripe];
             dm_stripe_total[ks]   += sum_kl0;
             dm_stripe_total[ks+1] += sum_kl1;
             dm_stripe_total[ks+2] += sum_kl2;
             dm_stripe_total[ks+3] += sum_kl3;
          }

          if (allzero_k||allzero_l) {
             // already done
          } else {
             // Use UnnormalizedWeighted since we already computed dm_stripe_total
             for (uint64_t i=0; i<4; i++) {
                UnnormalizedWeighted1<TFloat>(
                                   dm_stripes_buf,
                                   zcheck, sums, embedded_proportions, lengths,
                                   embs_stripe,filled_embs,idx, n_samples_r,
                                   ks+i, ls+i);
             }
          }
       } // (allzero_k && allzero_l)
}

template<class TFloat>
static inline void NormalizedWeighted8(
                      TFloat * const __restrict__ dm_stripes_buf,
                      TFloat * const __restrict__ dm_stripes_total_buf,
                      const bool   * const __restrict__ zcheck,
                      const TFloat * const __restrict__ sums,
                      const TFloat * const __restrict__ embedded_proportions,
                      const TFloat * __restrict__ lengths,
                      const uint64_t embs_stripe,
                      const unsigned int filled_embs,
                      const uint64_t idx,
                      const uint64_t n_samples_r,
                      const uint64_t ks,
                      const uint64_t ls) {
       const uint64_t z_k = ((const uint64_t *)(zcheck+ks))[0];
       const uint64_t z_l = ((const uint64_t *)(zcheck+ls))[0];
       const bool allzero_k = z_k==0x0101010101010101;
       const bool allzero_l = z_l==0x0101010101010101;

       if (allzero_k && allzero_l) {
         // nothing to do, would have to add 0
       } else {
          const TFloat sum_k0 = sums[ks];
          const TFloat sum_k1 = sums[ks+1];
          const TFloat sum_k2 = sums[ks+2];
          const TFloat sum_k3 = sums[ks+3];
          const TFloat sum_k4 = sums[ks+4];
          const TFloat sum_k5 = sums[ks+5];
          const TFloat sum_k6 = sums[ks+6];
          const TFloat sum_k7 = sums[ks+7];
          const TFloat sum_l0 = sums[ls];
          const TFloat sum_l1 = sums[ls+1];
          const TFloat sum_l2 = sums[ls+2];
          const TFloat sum_l3 = sums[ls+3];
          const TFloat sum_l4 = sums[ls+4];
          const TFloat sum_l5 = sums[ls+5];
          const TFloat sum_l6 = sums[ls+6];
          const TFloat sum_l7 = sums[ls+7];

          const TFloat sum_kl0 = sum_k0 + sum_l0;
          const TFloat sum_kl1 = sum_k1 + sum_l1;
          const TFloat sum_kl2 = sum_k2 + sum_l2;
          const TFloat sum_kl3 = sum_k3 + sum_l3;
          const TFloat sum_kl4 = sum_k4 + sum_l4;
          const TFloat sum_kl5 = sum_k5 + sum_l5;
          const TFloat sum_kl6 = sum_k6 + sum_l6;
          const TFloat sum_kl7 = sum_k7 + sum_l7;

          TFloat * const __restrict__ dm_stripe = dm_stripes_buf+idx;
          //TFloat *dm_stripe = dm_stripes[stripe];

          if (allzero_k) {
             // one side has all zeros, use distributed property
             // if (nonzero_l1) ridx=l1 // fabs(k-l1), with k==0
             dm_stripe[ks]   += sum_l0;
             dm_stripe[ks+1] += sum_l1;
             dm_stripe[ks+2] += sum_l2;
             dm_stripe[ks+3] += sum_l3;
             dm_stripe[ks+4] += sum_l4;
             dm_stripe[ks+5] += sum_l5;
             dm_stripe[ks+6] += sum_l6;
             dm_stripe[ks+7] += sum_l7;
          } else if (allzero_l) {
             // one side has all zeros, use distributed property
             // if (nonzero_k)  ridx=k  // fabs(k-l1), with l1==0
             dm_stripe[ks]   += sum_k0;
             dm_stripe[ks+1] += sum_k1;
             dm_stripe[ks+2] += sum_k2;
             dm_stripe[ks+3] += sum_k3;
             dm_stripe[ks+4] += sum_k4;
             dm_stripe[ks+5] += sum_k5;
             dm_stripe[ks+6] += sum_k6;
             dm_stripe[ks+7] += sum_k7;
          }

          // the totals can always use the distributed property
          {
             TFloat * const __restrict__ dm_stripe_total = dm_stripes_total_buf+idx;
             //TFloat *dm_stripe_total = dm_stripes_total[stripe];
             dm_stripe_total[ks]   += sum_kl0;
             dm_stripe_total[ks+1] += sum_kl1;
             dm_stripe_total[ks+2] += sum_kl2;
             dm_stripe_total[ks+3] += sum_kl3;
             dm_stripe_total[ks+4] += sum_kl4;
             dm_stripe_total[ks+5] += sum_kl5;
             dm_stripe_total[ks+6] += sum_kl6;
             dm_stripe_total[ks+7] += sum_kl7;
          }

          if (allzero_k||allzero_l) {
             // already done
          } else {
             // Use UnnormalizedWeighted since we already computed dm_stripe_total
             for (uint64_t i=0; i<8; i+=4) {
                UnnormalizedWeighted4<TFloat>(
                                   dm_stripes_buf,
                                   zcheck, sums, embedded_proportions, lengths,
                                   embs_stripe,filled_embs,idx, n_samples_r,
                                   ks+i, ls+i);
             }
          }
       } // (allzero_k && allzero_l)
}
#endif

template<class TFloat>
static inline void run_NormalizedWeightedTask_T(
                const uint64_t embs_stripe,
		const unsigned int filled_embs,
		const uint64_t start_idx, const uint64_t stop_idx,
		const uint64_t n_samples, const uint64_t n_samples_r,
		const TFloat * const __restrict__ lengths,
		const TFloat * const __restrict__ embedded_proportions,
		TFloat * const __restrict__ dm_stripes_buf,
		TFloat * const __restrict__ dm_stripes_total_buf,
		bool * const __restrict__ zcheck,
		TFloat * const __restrict__ sums) {

    constexpr uint64_t step_size = STEP_SIZE(TFloat);
    const uint64_t sample_steps = (n_samples+(step_size-1))/step_size; // round up

    // check for zero values and pre-compute single column sums
    WeightedZerosAndSums(zcheck, sums,
                         embedded_proportions, lengths,
                         embs_stripe,filled_embs, n_samples, n_samples_r);

    // point of thread
#if defined(_OPENACC) || defined(OMPGPU)
#if defined(OMPGPU)
    // TODO: Explore async omp target
#pragma omp target teams distribute parallel for simd collapse(3) default(shared)
#else
#pragma acc parallel loop gang vector collapse(3) present(embedded_proportions,dm_stripes_buf,dm_stripes_total_buf,lengths,zcheck,sums) async
#endif
    for(uint64_t sk = 0; sk < sample_steps ; sk++) {
     for(uint64_t stripe = start_idx; stripe < stop_idx; stripe++) {
      // SIMT-based GPU work great one at a time (HW will deal with parallelism)
      for(uint64_t ik = 0; ik < step_size ; ik++) {
       const uint64_t k = sk*step_size + ik;

       if (k>=n_samples) continue; // past the limit

       const uint64_t l1 = (k + stripe + 1)%n_samples; // wraparound
       const uint64_t idx = (stripe-start_idx) * n_samples_r;

       NormalizedWeighted1<TFloat>(
                                   dm_stripes_buf,dm_stripes_total_buf,
                                   zcheck, sums, embedded_proportions, lengths,
                                   embs_stripe,filled_embs,idx, n_samples_r,
                                   k, l1);
      } // for ik
     } // for stripe
    } // for sk
#else
    // tilling helps with better cache reuse without the need of multiple cores
    const uint64_t stripe_steps = ((stop_idx-start_idx)+(step_size-1))/step_size; // round up

    // use dynamic scheduling due to non-homogeneity in the loop
    // Use a moderate block to prevent trashing but still have some cache reuse
#pragma omp parallel for collapse(2) schedule(dynamic,step_size) default(shared)
    for(uint64_t ss = 0; ss < stripe_steps ; ss++) {
     for(uint64_t sk = 0; sk < sample_steps ; sk++) {
      // tile to maximize cache reuse
      for(uint64_t is = 0; is < step_size ; is++) {
       const uint64_t stripe = start_idx+ss*step_size + is;
       if (stripe<stop_idx) { // else past limit

      // SIMD-based CPUs need help with vectorization
      const uint64_t idx = (stripe-start_idx) * n_samples_r;
      uint64_t ks = sk*step_size;
      const uint64_t kmax = std::min(ks+step_size,n_samples);
      uint64_t ls = (ks + stripe + 1)%n_samples; // wraparound

      while( ((ks+8) <= kmax) && ((n_samples-ls)>=8) ) {
       NormalizedWeighted8<TFloat>(
                                   dm_stripes_buf,dm_stripes_total_buf,
                                   zcheck, sums, embedded_proportions, lengths,
                                   embs_stripe,filled_embs,idx, n_samples_r,
                                   ks, ls);
       ks+=8;
       ls = (ls + 8)%n_samples; // wraparound
      } // for ks+=8

      while( ((ks+4) <= kmax) && ((n_samples-ls)>=4) ) {
       NormalizedWeighted4<TFloat>(
                                   dm_stripes_buf,dm_stripes_total_buf,
                                   zcheck, sums, embedded_proportions, lengths,
                                   embs_stripe,filled_embs,idx, n_samples_r,
                                   ks, ls);
       ks+=4;
       ls = (ls + 4)%n_samples; // wraparound
      } // for ks+=4

      // deal with any leftovers in serial way
      for( ; ks < kmax; ks++ ) {
       NormalizedWeighted1<TFloat>(
                                   dm_stripes_buf,dm_stripes_total_buf,
                                   zcheck, sums, embedded_proportions, lengths,
                                   embs_stripe,filled_embs,idx, n_samples_r,
                                   ks, ls);
       ls = (ls + 1)%n_samples; // wraparound
      } // for ks

       } // if stripe
      } // for is
     } // for sk
    } // for ss
#endif
}

template<class TFloat>
static inline void run_VawNormalizedWeightedTask_T(
		const unsigned int filled_embs,
		const uint64_t start_idx, const uint64_t stop_idx,
		const uint64_t n_samples, const uint64_t n_samples_r,
		const TFloat * const __restrict__ lengths,
		const TFloat * const __restrict__ embedded_proportions,
		const TFloat * const __restrict__ embedded_counts,
		const TFloat * const __restrict__ sample_total_counts,
		TFloat * const __restrict__ dm_stripes_buf,
		TFloat * const __restrict__ dm_stripes_total_buf) {

    constexpr uint64_t step_size = STEP_SIZE(TFloat);
    const uint64_t sample_steps = (n_samples+(step_size-1))/step_size; // round up

    // point of thread
#if defined(_OPENACC) || defined(OMPGPU)
    // Use as big vector size as we can, to maximize cache line reuse
#if defined(OMPGPU)
    // TODO: Explore async omp target
#pragma omp target teams distribute parallel for simd collapse(3) simdlen( SUCMP_ACC_MAXVEC ) default(shared)
#else
#pragma acc parallel loop collapse(3) vector_length( SUCMP_ACC_MAXVEC ) present(embedded_proportions,embedded_counts,sample_total_counts,dm_stripes_buf,dm_stripes_total_buf,lengths) async
#endif

#else
#pragma omp parallel for schedule(dynamic,1) default(shared)
#endif
    for(uint64_t sk = 0; sk < sample_steps ; sk++) {
      for(uint64_t stripe = start_idx; stripe < stop_idx; stripe++) {
        for(uint64_t ik = 0; ik < step_size ; ik++) {
            const uint64_t k = sk*step_size + ik;
            const uint64_t idx = (stripe-start_idx) * n_samples_r;
            TFloat * const __restrict__ dm_stripe = dm_stripes_buf+idx;
            TFloat * const __restrict__ dm_stripe_total = dm_stripes_total_buf+idx;
            //TFloat *dm_stripe = dm_stripes[stripe];
            //TFloat *dm_stripe_total = dm_stripes_total[stripe];

            if (k>=n_samples) continue; // past the limit

            const uint64_t l1 = (k + stripe + 1)%n_samples; // wraparound

            const TFloat m = sample_total_counts[k] + sample_total_counts[l1];

            TFloat my_stripe = dm_stripe[k];
            TFloat my_stripe_total = dm_stripe_total[k];

#if !defined(OMPGPU) && defined(_OPENACC)
#pragma acc loop seq
#endif
            for (uint64_t emb=0; emb<filled_embs; emb++) {
                const uint64_t offset = n_samples_r * emb;

                TFloat mi = embedded_counts[offset + k] + embedded_counts[offset + l1];
                TFloat vaw = sqrt(mi * (m - mi));

                if(vaw > 0) {
                  TFloat u1 = embedded_proportions[offset + k];
                  TFloat v1 = embedded_proportions[offset + l1];
                  TFloat diff1 = fabs(u1 - v1);
                  TFloat length = lengths[emb];

                  my_stripe += (diff1 * length) / vaw;
                  my_stripe_total += ((u1 + v1) * length) / vaw;
                }
            }

            dm_stripe[k]     = my_stripe;
            dm_stripe_total[k]     = my_stripe_total;

        }

      }
    }
}

template<class TFloat>
static inline void run_GeneralizedTask_T(
                const uint64_t embs_stripe,
		const unsigned int filled_embs,
		const uint64_t start_idx, const uint64_t stop_idx,
		const uint64_t n_samples, const uint64_t n_samples_r,
		const TFloat * const __restrict__ lengths,
		const TFloat * const __restrict__ embedded_proportions,
		TFloat * const __restrict__ dm_stripes_buf,
		TFloat * const __restrict__ dm_stripes_total_buf,
		const TFloat g_unifrac_alpha) {
    constexpr uint64_t step_size = STEP_SIZE(TFloat);
    const uint64_t sample_steps = (n_samples+(step_size-1))/step_size; // round up

#if !(defined(_OPENACC) || defined(OMPGPU))
    // CPU version uses transposed embedded_proportions

    // point of thread
#pragma omp parallel for schedule(dynamic,1) default(shared)
    for(uint64_t sk = 0; sk < sample_steps ; sk++) {
      for(uint64_t stripe = start_idx; stripe < stop_idx; stripe++) {
        for(uint64_t ik = 0; ik < step_size ; ik++) {
            const uint64_t k = sk*step_size + ik;
            const uint64_t idx = (stripe-start_idx) * n_samples_r;
            TFloat * const __restrict__ dm_stripe = dm_stripes_buf+idx;
            TFloat * const __restrict__ dm_stripe_total = dm_stripes_total_buf+idx;
            //TFloat *dm_stripe = dm_stripes[stripe];
            //TFloat *dm_stripe_total = dm_stripes_total[stripe];

            if (k>=n_samples) continue; // past the limit

            const uint64_t l1 = (k + stripe + 1)%n_samples; // wraparound

            TFloat my_stripe = dm_stripe[k];
            TFloat my_stripe_total = dm_stripe_total[k];

            const uint64_t offset_k = embs_stripe*k;
            const uint64_t offset_l = embs_stripe*l1;
            for (uint64_t emb=0; emb<filled_embs; emb++) {
                TFloat u1 = embedded_proportions[offset_k + emb];
                TFloat v1 = embedded_proportions[offset_l + emb];
                TFloat sum1 = u1 + v1;

                if(sum1 != 0.0) { 
                   TFloat length = lengths[emb];
                   TFloat diff1 = fabs(u1 - v1);
                   TFloat sum_pow1 = pow(sum1, g_unifrac_alpha) * length; 

                   my_stripe += sum_pow1 * (diff1 / sum1); 
                   my_stripe_total += sum_pow1; 
                }
            }

            dm_stripe[k]     = my_stripe;
            dm_stripe_total[k]     = my_stripe_total;

        }

      }
    }
#else
    // GPU version uses straight embedded_proportions

    // point of thread
    // Use as big vector size as we can, to maximize cache line reuse
#if defined(OMPGPU)
    // TODO: Explore async omp target
#pragma omp target teams distribute parallel for simd collapse(3) simdlen( SUCMP_ACC_MAXVEC ) default(shared)
#else
#pragma acc parallel loop collapse(3) vector_length( SUCMP_ACC_MAXVEC ) present(embedded_proportions,dm_stripes_buf,dm_stripes_total_buf,lengths) async
#endif
    for(uint64_t sk = 0; sk < sample_steps ; sk++) {
      for(uint64_t stripe = start_idx; stripe < stop_idx; stripe++) {
        for(uint64_t ik = 0; ik < step_size ; ik++) {
            const uint64_t k = sk*step_size + ik;
            const uint64_t idx = (stripe-start_idx) * n_samples_r;
            TFloat * const __restrict__ dm_stripe = dm_stripes_buf+idx;
            TFloat * const __restrict__ dm_stripe_total = dm_stripes_total_buf+idx;
            //TFloat *dm_stripe = dm_stripes[stripe];
            //TFloat *dm_stripe_total = dm_stripes_total[stripe];

            if (k>=n_samples) continue; // past the limit

            const uint64_t l1 = (k + stripe + 1)%n_samples; // wraparound

            TFloat my_stripe = dm_stripe[k];
            TFloat my_stripe_total = dm_stripe_total[k];

#if !defined(OMPGPU) && defined(_OPENACC)
#pragma acc loop seq
#endif
            for (uint64_t emb=0; emb<filled_embs; emb++) {
                const uint64_t offset = n_samples_r * emb;

                TFloat u1 = embedded_proportions[offset + k];
                TFloat v1 = embedded_proportions[offset + l1];
                TFloat sum1 = u1 + v1;

                if(sum1 != 0.0) { 
                   TFloat length = lengths[emb];
                   TFloat diff1 = fabs(u1 - v1);
                   TFloat sum_pow1 = pow(sum1, g_unifrac_alpha) * length; 

                   my_stripe += sum_pow1 * (diff1 / sum1); 
                   my_stripe_total += sum_pow1; 
                }
            }

            dm_stripe[k]     = my_stripe;
            dm_stripe_total[k]     = my_stripe_total;

        }

      }
    }
#endif
}

template<class TFloat>
static inline void run_VawGeneralizedTask_T(
		const unsigned int filled_embs,
		const uint64_t start_idx, const uint64_t stop_idx, const uint64_t n_samples, const uint64_t n_samples_r,
		const TFloat * const __restrict__ lengths,
		const TFloat * const __restrict__ embedded_proportions,
		const TFloat * const __restrict__ embedded_counts,
		const TFloat * const __restrict__ sample_total_counts ,
		TFloat * const __restrict__ dm_stripes_buf,
		TFloat * const __restrict__ dm_stripes_total_buf,
		const TFloat g_unifrac_alpha) {

    constexpr uint64_t step_size = STEP_SIZE(TFloat);
    const uint64_t sample_steps = (n_samples+(step_size-1))/step_size; // round up
    // quick hack, to be finished

    // point of thread
#if defined(_OPENACC) || defined(OMPGPU)
    // Use as big vector size as we can, to maximize cache line reuse
#if defined(OMPGPU)
    // TODO: Explore async omp target
#pragma omp target teams distribute parallel for simd collapse(3) simdlen( SUCMP_ACC_MAXVEC ) default(shared)
#else
#pragma acc parallel loop collapse(3) vector_length( SUCMP_ACC_MAXVEC ) present(embedded_proportions,embedded_counts,sample_total_counts,dm_stripes_buf,dm_stripes_total_buf,lengths) async
#endif

#else
#pragma omp parallel for schedule(dynamic,1) default(shared)
#endif
    for(uint64_t sk = 0; sk < sample_steps ; sk++) {
      for(uint64_t stripe = start_idx; stripe < stop_idx; stripe++) {
        for(uint64_t ik = 0; ik < step_size ; ik++) {
            const uint64_t k = sk*step_size + ik;
            const uint64_t idx = (stripe-start_idx) * n_samples_r;
            TFloat * const __restrict__ dm_stripe = dm_stripes_buf+idx;
            TFloat * const __restrict__ dm_stripe_total = dm_stripes_total_buf+idx;
            //TFloat *dm_stripe = dm_stripes[stripe];
            //TFloat *dm_stripe_total = dm_stripes_total[stripe];

            if (k>=n_samples) continue; // past the limit

            const uint64_t l1 = (k + stripe + 1)%n_samples; // wraparound

            const TFloat m = sample_total_counts[k] + sample_total_counts[l1];

            TFloat my_stripe = dm_stripe[k];
            TFloat my_stripe_total = dm_stripe_total[k];

#if !defined(OMPGPU) && defined(_OPENACC)
#pragma acc loop seq
#endif
            for (uint64_t emb=0; emb<filled_embs; emb++) {
                const uint64_t offset = n_samples_r * emb;

                TFloat mi = embedded_counts[offset + k] + embedded_counts[offset + l1];
                TFloat vaw = sqrt(mi * (m - mi));

                if(vaw > 0) {
                  TFloat u1 = embedded_proportions[offset + k];
                  TFloat v1 = embedded_proportions[offset + l1];
                  TFloat length = lengths[emb];

                  TFloat sum1 = (u1 + v1) / vaw;
                  TFloat sub1 = fabs(u1 - v1) / vaw;
                  TFloat sum_pow1 = pow(sum1, g_unifrac_alpha) * length;

                  my_stripe += sum_pow1 * (sub1 / sum1);
                  my_stripe_total += sum_pow1;
                }
            }

            dm_stripe[k]     = my_stripe;
            dm_stripe_total[k]     = my_stripe_total;

        }
      }
    }
}

// Single step in computing Unweighted Unifrac

template<class TFloat>
static inline bool UnweightedOneSide(
                      bool   * const __restrict__ zcheck,
                      TFloat * const  __restrict__ stripe_sums,
                      const TFloat * const   __restrict__ sums,
                      const uint64_t * const __restrict__ embedded_proportions,
                      const uint64_t embs_stripe,
                      const unsigned int filled_embs_els_round,
                      const uint64_t n_samples_r,
                      const uint64_t kl) {
            bool all_zeros=true;
            TFloat my_stripe = 0.0;

#if !(defined(_OPENACC) || defined(OMPGPU))
            // CPUs use transposed layout
            const uint64_t offset = embs_stripe*kl;
#endif
            for (uint64_t emb_el=0; emb_el<filled_embs_els_round; emb_el++) {
                const TFloat * __restrict__ psum = &(sums[emb_el*0x800]);

#if !(defined(_OPENACC) || defined(OMPGPU))
		// CPUs use transposed layout
                uint64_t o1 = embedded_proportions[offset + emb_el];
#else
                const uint64_t offset = n_samples_r * emb_el;
                uint64_t o1 = embedded_proportions[offset + kl];
#endif

                if (o1==0) {  // zeros are prevalent
                    // nothing to do
                } else {
                    all_zeros = false;
                    // Use the pre-computed sums
                    // Each range of 8 lengths has already been pre-computed and stored in psum
                    // Since embedded_proportions packed format is in 64-bit format for performance reasons
                    //    we need to add the 8 sums using the four 8-bits for addressing inside psum

#if !(defined(_OPENACC) || defined(OMPGPU))
                    // CPU/SIMD faster if we check for partial compute
		    TFloat esum = 0.0;
		    for (int i=0; i<8; i++) {
		      uint8_t o1_8 = (uint8_t)(o1);
                      if (o1_8!=0) {
                        my_stripe  += psum[(uint8_t)(o1_8)];
		      }
		      o1 = o1 >> 8;
		      psum += 0x100;
		    }
#else
                    // GPU/SIMT faster if we just go ahead will full at all times
                    TFloat esum      = psum[       (uint8_t)(o1)       ] + 
                                       psum[0x100+((uint8_t)(o1 >>  8))] +
                                       psum[0x200+((uint8_t)(o1 >> 16))] +
                                       psum[0x300+((uint8_t)(o1 >> 24))] +
                                       psum[0x400+((uint8_t)(o1 >> 32))] +
                                       psum[0x500+((uint8_t)(o1 >> 40))] +
                                       psum[0x600+((uint8_t)(o1 >> 48))] +
                                       psum[0x700+((uint8_t)(o1 >> 56))];
#endif
                    my_stripe       += esum;
                }
            }

            stripe_sums[kl] = my_stripe;
	    return all_zeros;

}

#if defined(_OPENACC) || defined(OMPGPU)

template<class TFloat, bool compute_total>
static inline void Unweighted1(
                      TFloat * const __restrict__ dm_stripes_buf,
                      TFloat * const __restrict__ dm_stripes_total_buf,
                      const bool   * const __restrict__ zcheck,
                      const TFloat * const   __restrict__ stripe_sums,
                      const TFloat * const   __restrict__ sums,
                      const uint64_t * const __restrict__ embedded_proportions,
                      const unsigned int filled_embs_els_round,
                      const uint64_t idx,
                      const uint64_t n_samples_r,
                      const uint64_t k,
                      const uint64_t l1) {
       const bool allzero_k = zcheck[k];
       const bool allzero_l1 = zcheck[l1];

       if (allzero_k && allzero_l1) {
          // nothing to do, would have to add 0
       } else {
          TFloat * const __restrict__ dm_stripe = dm_stripes_buf+idx;
          TFloat * const __restrict__ dm_stripe_total = compute_total ? dm_stripes_total_buf+idx : NULL;
          //TFloat *dm_stripe = dm_stripes[stripe];
          //TFloat *dm_stripe_total = dm_stripes_total[stripe];

          bool did_update = false;
          TFloat my_stripe = 0.0;
          TFloat my_stripe_total = 0.0;

          if (allzero_k || allzero_l1) {
            // with one side zero, | and ^ are no-ops
            const uint64_t kl = (allzero_k) ? l1 : k; // only use the non-sero onea
            my_stripe = stripe_sums[kl];
            did_update = (my_stripe!=0.0);
            if constexpr (compute_total) {
            my_stripe_total = my_stripe;
	    }
          } else {
            // we need both sides
            for (uint64_t emb_el=0; emb_el<filled_embs_els_round; emb_el++) {
                const uint64_t offset = n_samples_r * emb_el;
                const TFloat * __restrict__ psum = &(sums[emb_el*0x800]);

                uint64_t u1 = embedded_proportions[offset + k];
                uint64_t v1 = embedded_proportions[offset + l1];
                uint64_t o1 = u1 | v1;

                if (o1!=0) {  // zeros are prevalent
                    did_update=true;
                    uint64_t x1 = u1 ^ v1;

                    // Use the pre-computed sums
                    // Each range of 8 lengths has already been pre-computed and stored in psum
                    // Since embedded_proportions packed format is in 64-bit format for performance reasons
                    //    we need to add the 8 sums using the four 8-bits for addressing inside psum

                    if constexpr (compute_total) {
                    my_stripe_total += psum[       (uint8_t)(o1)       ] + 
                                       psum[0x100+((uint8_t)(o1 >>  8))] +
                                       psum[0x200+((uint8_t)(o1 >> 16))] +
                                       psum[0x300+((uint8_t)(o1 >> 24))] +
                                       psum[0x400+((uint8_t)(o1 >> 32))] +
                                       psum[0x500+((uint8_t)(o1 >> 40))] +
                                       psum[0x600+((uint8_t)(o1 >> 48))] +
                                       psum[0x700+((uint8_t)(o1 >> 56))];
                    }
                    my_stripe       += psum[       (uint8_t)(x1)       ] + 
                                       psum[0x100+((uint8_t)(x1 >>  8))] +
                                       psum[0x200+((uint8_t)(x1 >> 16))] +
                                       psum[0x300+((uint8_t)(x1 >> 24))] +
                                       psum[0x400+((uint8_t)(x1 >> 32))] +
                                       psum[0x500+((uint8_t)(x1 >> 40))] +
                                       psum[0x600+((uint8_t)(x1 >> 48))] +
                                       psum[0x700+((uint8_t)(x1 >> 56))];
                }
            }
          }

          if (did_update) {
            dm_stripe[k]       += my_stripe;
            if constexpr (compute_total) {
            dm_stripe_total[k] += my_stripe_total;
	    }
          }
       } // (allzero_k && allzero_l1)

}
#else

template<class TFloat, bool compute_total>
static inline void Unweighted1(
                      TFloat * const __restrict__ dm_stripes_buf,
                      TFloat * const __restrict__ dm_stripes_total_buf,
                      const bool   * const __restrict__ zcheck,
                      const TFloat * const   __restrict__ stripe_sums,
                      const TFloat * const   __restrict__ sums,
                      const uint64_t * const __restrict__ embedded_proportions,
                      const uint64_t embs_stripe,
                      const unsigned int filled_embs_els_round,
                      const uint64_t idx,
                      const uint64_t n_samples_r,
                      const uint64_t k,
                      const uint64_t l1) {
       const bool allzero_k = zcheck[k];
       const bool allzero_l1 = zcheck[l1];

       if (allzero_k && allzero_l1) {
          // nothing to do, would have to add 0
       } else {
          TFloat * const __restrict__ dm_stripe = dm_stripes_buf+idx;
          TFloat * const __restrict__ dm_stripe_total = compute_total ? dm_stripes_total_buf+idx : NULL;
          //TFloat *dm_stripe = dm_stripes[stripe];
          //TFloat *dm_stripe_total = dm_stripes_total[stripe];

          bool did_update = false;
          TFloat my_stripe = 0.0;
          TFloat my_stripe_total = 0.0;

          if (allzero_k || allzero_l1) {
            // with one side zero, | and ^ are no-ops
            const uint64_t kl = (allzero_k) ? l1 : k; // only use the non-sero onea
            my_stripe = stripe_sums[kl];
            if constexpr (compute_total) {
            my_stripe_total = my_stripe;
            }
            did_update = (my_stripe!=0.0);
          } else {
            // we need both sides
            const uint64_t offset_k = embs_stripe*k;
            const uint64_t offset_l = embs_stripe*l1;
            for (uint64_t emb_el=0; emb_el<filled_embs_els_round; emb_el++) {
                const TFloat * __restrict__ psum = &(sums[emb_el*0x800]);

                uint64_t u1 = embedded_proportions[offset_k + emb_el];
                uint64_t v1 = embedded_proportions[offset_l + emb_el];
                uint64_t o1 = u1 | v1;
                uint64_t x1 = u1 ^ v1;

                if (o1==0) {  // zeros are prevalent
                    // nothing to do
                } else {
                    did_update=true;

                    // Use the pre-computed sums
                    // Each range of 8 lengths has already been pre-computed and stored in psum
                    // Since embedded_proportions packed format is in 64-bit format for performance reasons
                    //    we need to add the 8 sums using the four 8-bits for addressing inside psum
		    for (int i=0; i<8; i++) {
		      uint8_t o1_8 = (uint8_t)(o1);
                      if (o1_8!=0) {
                        if constexpr (compute_total) my_stripe_total += psum[o1_8];
                        my_stripe       += psum[(uint8_t)(x1)];
		      }
		      o1 = o1 >> 8;
		      x1 = x1 >> 8;
		      psum += 0x100;
		    }
                }
            }
          } // (allzero_k || allzero_l1)

          if (did_update) {
              dm_stripe[k]       += my_stripe;
              if constexpr (compute_total) {
              dm_stripe_total[k] += my_stripe_total;
              }
          }
       } // (allzero_k && allzero_l1)

}
#endif

// check for zero values
template<class TFloat, class T>
static inline uint32_t UnweightedZerosAndSums(
                      bool   * const __restrict__ zcheck,
		      uint32_t* const __restrict__ idxs,
                      TFloat * const  __restrict__ stripe_sums,
                      const TFloat * const   __restrict__ el_sums,
                      const T * const __restrict__ embedded_proportions,
                      const uint64_t embs_stripe,
                      const unsigned int filled_embs_els_round,
                      const uint32_t n_samples,
                      const uint64_t n_samples_r) {
    uint32_t n_true_idxs = 0;
#if defined(OMPGPU)
#pragma omp target teams distribute parallel for simd reduction(+:n_true_idxs) default(shared)
#elif defined(_OPENACC)
#pragma acc parallel loop gang vector reduction(+:n_true_idxs) present(embedded_proportions,zcheck,el_sums,stripe_sums)
#else
#pragma omp parallel for reduction(+:n_true_idxs) default(shared)
#endif
    for(uint32_t k=0; k<n_samples; k++) {
            bool all_zeros = UnweightedOneSide(
                      zcheck, stripe_sums,
                      el_sums, embedded_proportions,
                      embs_stripe, filled_embs_els_round, n_samples_r,
                      k);
            zcheck[k] = all_zeros;
	    if (all_zeros) n_true_idxs++;
    }

#if defined(OMPGPU) || defined(_OPENACC)
    // create index of k, first all of those with zcheck true, then all false
    // equivalent to stable_sort, but knowing in advance n_true_idxs
#if defined(OMPGPU)
#pragma omp target teams distribute parallel for simd default(shared)
#elif defined(_OPENACC)
#pragma acc parallel loop gang vector present(idxs,zcheck)
#endif
    for (int b=0; b<2; b++) {
      const bool mytest = (b==0);
      uint32_t icurr =  mytest ? 0 : n_true_idxs; 
      for(uint32_t k=0; k<uint32_t(n_samples); k++) {
        if (zcheck[k]==mytest) {
          idxs[icurr] = k;
          icurr++;;
        }
      }
    }
#else
    // CPU version does not use idxs
#endif

    return n_true_idxs;
}

template<class TFloat>
static inline void run_UnweightedTask_T(
		const uint64_t embs_stripe,
		const unsigned int filled_embs,
		const uint64_t start_idx, const uint64_t stop_idx,
		const uint64_t n_samples, const uint64_t n_samples_r,
		const TFloat * const __restrict__ lengths,
		const uint64_t * const __restrict__ embedded_proportions,
		TFloat * const __restrict__ dm_stripes_buf,
		TFloat * const __restrict__ dm_stripes_total_buf,
		TFloat * const __restrict__ sums,
		bool   * const __restrict__ zcheck,
		uint32_t* const __restrict__ idxs,
		TFloat * const __restrict__ stripe_sums) {
    static constexpr bool compute_total = true;

    constexpr uint64_t step_size = STEP_SIZE(TFloat);
    const uint64_t sample_steps = (n_samples+(step_size-1))/step_size; // round up

    const uint64_t filled_embs_els = filled_embs/64;
    const uint64_t filled_embs_rem = filled_embs%64; 

    const uint64_t filled_embs_els_round = (filled_embs+63)/64;

    // pre-compute sums of length elements, since they are likely to be accessed many times
    // We will use a 8-bit map, to keep it small enough to keep in L1 cache
#if defined(OMPGPU)
    // TODO: Explore async omp target
#pragma omp target teams distribute parallel for collapse(2) default(shared)
#elif defined(_OPENACC)
#pragma acc parallel loop collapse(2) gang present(lengths,sums) async
#else 
#pragma omp parallel for default(shared)
#endif
    for (uint64_t emb_el=0; emb_el<filled_embs_els; emb_el++) {
       for (uint64_t sub8=0; sub8<8; sub8++) {
          const uint64_t emb8 = emb_el*8+sub8;
          TFloat * __restrict__ psum = &(sums[emb8<<8]);
          const TFloat * __restrict__ pl   = &(lengths[emb8*8]);

#if defined(OMPGPU)
#pragma omp simd
#elif defined(_OPENACC)
#pragma acc loop vector
#endif
          // compute all the combinations for this block (8-bits total)
          // psum[0] = 0.0   // +0*pl[0]+0*pl[1]+0*pl[2]+...
          // psum[1] = pl[0] // +0*pl[1]+0*pl[2]+...
          // psum[2] = pl[1] // +0*pl[0]+0*pl[2]+
          // psum[2] = pl[0] + pl[1]
          // ...
          // psum[255] = pl[1] +.. + pl[7] // + 0*pl[0]
          // psum[255] = pl[0] +pl[1] +.. + pl[7]
          for (uint64_t b8_i=0; b8_i<0x100; b8_i++) {
             psum[b8_i] = (((b8_i >> 0) & 1) * pl[0]) + (((b8_i >> 1) & 1) * pl[1]) + 
                          (((b8_i >> 2) & 1) * pl[2]) + (((b8_i >> 3) & 1) * pl[3]) +
                          (((b8_i >> 4) & 1) * pl[4]) + (((b8_i >> 5) & 1) * pl[5]) +
                          (((b8_i >> 6) & 1) * pl[6]) + (((b8_i >> 7) & 1) * pl[7]);
          }
       }
    }
    if (filled_embs_rem>0) { // add also the overflow elements
       const uint64_t emb_el=filled_embs_els;
#if defined(OMPGPU)
#pragma omp target teams distribute parallel for default(shared)
#elif defined(_OPENACC)
#pragma acc parallel loop gang present(lengths,sums) async
#else
       // no advantage of OMP, too small
#endif
       for (uint64_t sub8=0; sub8<8; sub8++) {
          // we are summing we have enough buffer in sums
          const uint64_t emb8 = emb_el*8+sub8;
          TFloat * __restrict__ psum = &(sums[emb8<<8]);

#if defined(OMPGPU)
#pragma omp simd
#elif defined(_OPENACC)
#pragma acc loop vector
#endif
          // compute all the combinations for this block, set to 0 any past the limit
          // as above
          for (uint64_t b8_i=0; b8_i<0x100; b8_i++) {
             TFloat val= 0;
             for (uint64_t li=(emb8*8); li<filled_embs; li++) {
               val += ((b8_i >>  (li-(emb8*8))) & 1) * lengths[li];
             }
             psum[b8_i] = val;
          }
        }
    }

#if !defined(OMPGPU) && defined(_OPENACC)
#pragma acc wait
#endif
    // check for zero values and compute stripe sums
    uint32_t n_true_idxs = 
	    UnweightedZerosAndSums(zcheck, idxs, stripe_sums,
                           sums, embedded_proportions,
                           embs_stripe, filled_embs_els_round, n_samples, n_samples_r);


    // point of thread
#if defined(_OPENACC) || defined(OMPGPU)

    // With idxs ordered, keep elements with same zcheck_k together
    // This minimizes warp divergence, which is important given drastically different cost of some paths
    // It does potentially disrupt memory coalescence, but the major memory cost is access to su, which is pseudo-random anyway
    // (Note that sorting is not affective for the Weighted methods, since they have a much more regular memory access pattern)
#if defined(OMPGPU)
    // TODO: Explore async omp target
#pragma omp target teams distribute parallel for simd collapse(3) default(shared)
#else
#pragma acc parallel loop collapse(3) gang vector present(embedded_proportions,dm_stripes_buf,dm_stripes_total_buf,sums,zcheck,idxs,stripe_sums) async
#endif
    for(uint64_t sk = 0; sk < sample_steps ; sk++) {
      for(uint64_t stripe = start_idx; stripe < stop_idx; stripe++) {
        for(uint64_t ik = 0; ik < step_size ; ik++) {
          const uint64_t k_idx = sk*step_size + ik;
          if (k_idx<n_samples) { // else past the limit
	    const uint64_t k = idxs[k_idx];
            const uint64_t idx = (stripe-start_idx) * n_samples_r;
            const uint64_t l1 = (k + stripe + 1)%n_samples; // wraparound

            Unweighted1<TFloat,compute_total>(
                                dm_stripes_buf,dm_stripes_total_buf,
                                zcheck, stripe_sums,
                                sums, embedded_proportions,
                                filled_embs_els_round,idx, n_samples_r,
                                k, l1);
	  } // if k_idx
        } // for ik
      } // for stripe
    } // for sk
#else
    // tilling helps with better cache reuse without the need of multiple cores
    const uint64_t stripe_steps = ((stop_idx-start_idx)+(step_size-1))/step_size; // round up

    // use dynamic scheduling due to non-homogeneity in the loop
    // Use a moderate block to prevent trashing but still have some cache reuse
#pragma omp parallel for collapse(2) schedule(dynamic,step_size) default(shared)
    for(uint64_t ss = 0; ss < stripe_steps ; ss++) {
     for(uint64_t sk = 0; sk < sample_steps ; sk++) {
       // tile to maximize cache reuse
       for(uint64_t is = 0; is < step_size ; is++) {
        const uint64_t stripe = start_idx+ss*step_size + is;
        if (stripe<stop_idx) { // esle past limit}
         for(uint64_t ik = 0; ik < step_size ; ik++) {
           const uint64_t k = sk*step_size + ik;
           if (k<n_samples) { // elsepast the limit
            const uint64_t idx = (stripe-start_idx) * n_samples_r;
            const uint64_t l1 = (k + stripe + 1)%n_samples; // wraparound

            Unweighted1<TFloat,compute_total>(
                                dm_stripes_buf,dm_stripes_total_buf,
                                zcheck, stripe_sums,
                                sums, embedded_proportions,
                                embs_stripe,filled_embs_els_round,idx, n_samples_r,
                                k, l1);
           } // if k
         } // for ik
        } // if stripe
       } // for is
      } // for sk
    } // for ss
#endif
}

template<class TFloat>
static inline void run_UnnormalizedUnweightedTask_T(
		const uint64_t embs_stripe,
		const unsigned int filled_embs,
		const uint64_t start_idx, const uint64_t stop_idx,
		const uint64_t n_samples, const uint64_t n_samples_r,
		const TFloat * const __restrict__ lengths,
		const uint64_t * const __restrict__ embedded_proportions,
		TFloat * const __restrict__ dm_stripes_buf,
		TFloat * const __restrict__ sums,
		bool   * const __restrict__ zcheck,
		uint32_t* const __restrict__ idxs,
		TFloat * const __restrict__ stripe_sums) {
    static constexpr bool compute_total = false;

    constexpr uint64_t step_size = STEP_SIZE(TFloat);
    const uint64_t sample_steps = (n_samples+(step_size-1))/step_size; // round up

    const uint64_t filled_embs_els = filled_embs/64;
    const uint64_t filled_embs_rem = filled_embs%64; 

    const uint64_t filled_embs_els_round = (filled_embs+63)/64;

    // pre-compute sums of length elements, since they are likely to be accessed many times
    // We will use a 8-bit map, to keep it small enough to keep in L1 cache
#if defined(OMPGPU)
    // TODO: Explore async omp target
#pragma omp target teams distribute parallel for collapse(2) default(shared)
#elif defined(_OPENACC)
#pragma acc parallel loop collapse(2) gang present(lengths,sums) async
#else 
#pragma omp parallel for default(shared)
#endif
    for (uint64_t emb_el=0; emb_el<filled_embs_els; emb_el++) {
       for (uint64_t sub8=0; sub8<8; sub8++) {
          const uint64_t emb8 = emb_el*8+sub8;
          TFloat * __restrict__ psum = &(sums[emb8<<8]);
          const TFloat * __restrict__ pl   = &(lengths[emb8*8]);

#if defined(OMPGPU)
#pragma omp simd
#elif defined(_OPENACC)
#pragma acc loop vector
#endif
          // compute all the combinations for this block (8-bits total)
          // psum[0] = 0.0   // +0*pl[0]+0*pl[1]+0*pl[2]+...
          // psum[1] = pl[0] // +0*pl[1]+0*pl[2]+...
          // psum[2] = pl[1] // +0*pl[0]+0*pl[2]+
          // psum[2] = pl[0] + pl[1]
          // ...
          // psum[255] = pl[1] +.. + pl[7] // + 0*pl[0]
          // psum[255] = pl[0] +pl[1] +.. + pl[7]
          for (uint64_t b8_i=0; b8_i<0x100; b8_i++) {
             psum[b8_i] = (((b8_i >> 0) & 1) * pl[0]) + (((b8_i >> 1) & 1) * pl[1]) + 
                          (((b8_i >> 2) & 1) * pl[2]) + (((b8_i >> 3) & 1) * pl[3]) +
                          (((b8_i >> 4) & 1) * pl[4]) + (((b8_i >> 5) & 1) * pl[5]) +
                          (((b8_i >> 6) & 1) * pl[6]) + (((b8_i >> 7) & 1) * pl[7]);
          }
       }
    }
    if (filled_embs_rem>0) { // add also the overflow elements
       const uint64_t emb_el=filled_embs_els;
#if defined(OMPGPU)
#pragma omp target teams distribute parallel for default(shared)
#elif defined(_OPENACC)
#pragma acc parallel loop gang present(lengths,sums) async
#else
       // no advantage of OMP, too small
#endif
       for (uint64_t sub8=0; sub8<8; sub8++) {
          // we are summing we have enough buffer in sums
          const uint64_t emb8 = emb_el*8+sub8;
          TFloat * __restrict__ psum = &(sums[emb8<<8]);

#if defined(OMPGPU)
#pragma omp simd
#elif defined(_OPENACC)
#pragma acc loop vector
#endif
          // compute all the combinations for this block, set to 0 any past the limit
          // as above
          for (uint64_t b8_i=0; b8_i<0x100; b8_i++) {
             TFloat val= 0;
             for (uint64_t li=(emb8*8); li<filled_embs; li++) {
               val += ((b8_i >>  (li-(emb8*8))) & 1) * lengths[li];
             }
             psum[b8_i] = val;
          }
        }
    }

#if !defined(OMPGPU) && defined(_OPENACC)
#pragma acc wait
#endif
    // check for zero values and compute stripe sums
    uint32_t n_true_idxs = 
            UnweightedZerosAndSums(zcheck, idxs, stripe_sums,
                           sums, embedded_proportions,
                           embs_stripe, filled_embs_els_round, n_samples, n_samples_r);

    // point of thread
#if defined(_OPENACC) || defined(OMPGPU)
    // With idxs ordered, keep elements with same zcheck_k together
    // This minimizes warp divergence, which is important given drastically different cost of some paths
    // It does potentially disrupt memory coalescence, but the major memory cost is access to su, which is pseudo-random anyway
    // (Note that sorting is not affective for the Weighted methods, since they have a much more regular memory access pattern)
#if defined(OMPGPU)
    // TODO: Explore async omp target
#pragma omp target teams distribute parallel for simd collapse(3) default(shared)
#else
#pragma acc parallel loop collapse(3) gang vector present(embedded_proportions,dm_stripes_buf,sums,zcheck,idxs,stripe_sums) async
#endif
    for(uint64_t sk = 0; sk < sample_steps ; sk++) {
      for(uint64_t stripe = start_idx; stripe < stop_idx; stripe++) {
        for(uint64_t ik = 0; ik < step_size ; ik++) {
          const uint64_t k_idx = sk*step_size + ik;
          if (k_idx<n_samples) { // else past the limit
	    const uint64_t k = idxs[k_idx];
            const uint64_t idx = (stripe-start_idx) * n_samples_r;
            const uint64_t l1 = (k + stripe + 1)%n_samples; // wraparound

            Unweighted1<TFloat,compute_total>(
                                dm_stripes_buf,NULL,
                                zcheck, stripe_sums,
                                sums, embedded_proportions,
                                filled_embs_els_round,idx, n_samples_r,
                                k, l1);
	  } // if k_idx
        } // for ik
      } // for stripe
    } // for sk
#else
    // tilling helps with better cache reuse without the need of multiple cores
    const uint64_t stripe_steps = ((stop_idx-start_idx)+(step_size-1))/step_size; // round up

    // use dynamic scheduling due to non-homogeneity in the loop
    // Use a moderate block to prevent trashing but still have some cache reuse
#pragma omp parallel for collapse(2) schedule(dynamic,step_size) default(shared)
    for(uint64_t ss = 0; ss < stripe_steps ; ss++) {
     for(uint64_t sk = 0; sk < sample_steps ; sk++) {
       // tile to maximize cache reuse
       for(uint64_t is = 0; is < step_size ; is++) {
        const uint64_t stripe = start_idx+ss*step_size + is;
        if (stripe<stop_idx) { // esle past limit}
         for(uint64_t ik = 0; ik < step_size ; ik++) {
           const uint64_t k = sk*step_size + ik;
           if (k<n_samples) { // elsepast the limit
            const uint64_t idx = (stripe-start_idx) * n_samples_r;
            const uint64_t l1 = (k + stripe + 1)%n_samples; // wraparound

            Unweighted1<TFloat,compute_total>(
                                dm_stripes_buf,NULL,
                                zcheck, stripe_sums,
                                sums, embedded_proportions,
                                embs_stripe,filled_embs_els_round,idx, n_samples_r,
                                k, l1);
           } // if k
         } // for ik
        } // if stripe
       } // for is
      } // for sk
    } // for ss
#endif
}

template<class TFloat>
static inline void run_VawUnweightedTask_T(
		const unsigned int filled_embs,
		const uint64_t start_idx, const uint64_t stop_idx,
		const uint64_t n_samples, const uint64_t n_samples_r,
		const TFloat * const __restrict__ lengths,
		const uint32_t * const __restrict__ embedded_proportions,
		const TFloat  * const __restrict__ embedded_counts,
		const TFloat  * const __restrict__ sample_total_counts,
		TFloat * const __restrict__ dm_stripes_buf,
		TFloat * const __restrict__ dm_stripes_total_buf) {

    constexpr uint64_t step_size = STEP_SIZE(TFloat);
    const uint64_t sample_steps = (n_samples+(step_size-1))/step_size; // round up

    const uint64_t filled_embs_els = (filled_embs+31)/32; // round up

    // point of thread
#if defined(_OPENACC) || defined(OMPGPU)
    // Use as big vector size as we can, to maximize cache line reuse
#if defined(OMPGPU)
    // TODO: Explore async omp target
#pragma omp target teams distribute parallel for simd collapse(3) simdlen( SUCMP_ACC_MAXVEC ) default(shared)
#else
#pragma acc parallel loop collapse(3) vector_length( SUCMP_ACC_MAXVEC ) present(embedded_proportions,embedded_counts,sample_total_counts,dm_stripes_buf,dm_stripes_total_buf,lengths) async
#endif

#else
#pragma omp parallel for schedule(dynamic,1) default(shared)
#endif
    for(uint64_t sk = 0; sk < sample_steps ; sk++) {
      for(uint64_t stripe = start_idx; stripe < stop_idx; stripe++) {
        for(uint64_t ik = 0; ik < step_size ; ik++) {
            const uint64_t k = sk*step_size + ik;
            const uint64_t idx = (stripe-start_idx) * n_samples_r;
            TFloat * const __restrict__ dm_stripe = dm_stripes_buf+idx;
            TFloat * const __restrict__ dm_stripe_total = dm_stripes_total_buf+idx;
            //TFloat *dm_stripe = dm_stripes[stripe];
            //TFloat *dm_stripe_total = dm_stripes_total[stripe];

            if (k>=n_samples) continue; // past the limit

            const uint64_t l1 = (k + stripe + 1)%n_samples; // wraparound

            TFloat my_stripe = dm_stripe[k];
            TFloat my_stripe_total = dm_stripe_total[k];

            const TFloat m = sample_total_counts[k] + sample_total_counts[l1];

#if !defined(OMPGPU) && defined(_OPENACC)
#pragma acc loop seq
#endif
            for (uint64_t emb_el=0; emb_el<filled_embs_els; emb_el++) {
                const uint64_t offset_p = n_samples_r * emb_el;
                uint32_t u1 = embedded_proportions[offset_p + k];
                uint32_t v1 = embedded_proportions[offset_p + l1];
                uint32_t x1 = u1 ^ v1;
                uint32_t o1 = u1 | v1;

                // embedded_proporions is packed

#if !defined(OMPGPU) && defined(_OPENACC)
#pragma acc loop seq
#endif
                for (uint64_t ei=0; ei<32; ei++) {
                   uint64_t emb=emb_el*32+ei;
                   if (emb<filled_embs) {
                     const uint64_t offset_c = n_samples_r * emb;
                     TFloat mi = embedded_counts[offset_c + k] + embedded_counts[offset_c + l1];
                     TFloat vaw = sqrt(mi * (m - mi));

                     // embedded_counts is not packed

                     if(vaw > 0) {
                       TFloat length = lengths[emb];
                       TFloat lv1 = length / vaw;

                       my_stripe +=       ((x1 >> ei) & 1)*lv1;
                       my_stripe_total += ((o1 >> ei) & 1)*lv1;
                     }
                   }
                }
            }

            dm_stripe[k]     = my_stripe;
            dm_stripe_total[k]     = my_stripe_total;

        }

      }
    }
}

template<class TFloat>
static inline void run_VawUnnormalizedUnweightedTask_T(
		const unsigned int filled_embs,
		const uint64_t start_idx, const uint64_t stop_idx,
		const uint64_t n_samples, const uint64_t n_samples_r,
		const TFloat * const __restrict__ lengths,
		const uint32_t * const __restrict__ embedded_proportions,
		const TFloat  * const __restrict__ embedded_counts,
		const TFloat  * const __restrict__ sample_total_counts,
		TFloat * const __restrict__ dm_stripes_buf) {

    constexpr uint64_t step_size = STEP_SIZE(TFloat);
    const uint64_t sample_steps = (n_samples+(step_size-1))/step_size; // round up

    const uint64_t filled_embs_els = (filled_embs+31)/32; // round up

    // point of thread
#if defined(_OPENACC) || defined(OMPGPU)
    // Use as big vector size as we can, to maximize cache line reuse
#if defined(OMPGPU)
    // TODO: Explore async omp target
#pragma omp target teams distribute parallel for simd collapse(3) simdlen( SUCMP_ACC_MAXVEC ) default(shared)
#else
#pragma acc parallel loop collapse(3) vector_length( SUCMP_ACC_MAXVEC ) present(embedded_proportions,embedded_counts,sample_total_counts,dm_stripes_buf,lengths) async
#endif

#else
#pragma omp parallel for schedule(dynamic,1) default(shared)
#endif
    for(uint64_t sk = 0; sk < sample_steps ; sk++) {
      for(uint64_t stripe = start_idx; stripe < stop_idx; stripe++) {
        for(uint64_t ik = 0; ik < step_size ; ik++) {
            const uint64_t k = sk*step_size + ik;
            const uint64_t idx = (stripe-start_idx) * n_samples_r;
            TFloat * const __restrict__ dm_stripe = dm_stripes_buf+idx;
            //TFloat *dm_stripe = dm_stripes[stripe];

            if (k>=n_samples) continue; // past the limit

            const uint64_t l1 = (k + stripe + 1)%n_samples; // wraparound

            TFloat my_stripe = dm_stripe[k];

            const TFloat m = sample_total_counts[k] + sample_total_counts[l1];

#if !defined(OMPGPU) && defined(_OPENACC)
#pragma acc loop seq
#endif
            for (uint64_t emb_el=0; emb_el<filled_embs_els; emb_el++) {
                const uint64_t offset_p = n_samples_r * emb_el;
                uint32_t u1 = embedded_proportions[offset_p + k];
                uint32_t v1 = embedded_proportions[offset_p + l1];
                uint32_t x1 = u1 ^ v1;

                // embedded_proportions is packed

#if !defined(OMPGPU) && defined(_OPENACC)
#pragma acc loop seq
#endif
                for (uint64_t ei=0; ei<32; ei++) {
                   uint64_t emb=emb_el*32+ei;
                   if (emb<filled_embs) {
                     const uint64_t offset_c = n_samples_r * emb;
                     TFloat mi = embedded_counts[offset_c + k] + embedded_counts[offset_c + l1];
                     TFloat vaw = sqrt(mi * (m - mi));

                     // embedded_counts is not packed

                     if(vaw > 0) {
                       TFloat length = lengths[emb];
                       TFloat lv1 = length / vaw;

                       my_stripe +=       ((x1 >> ei) & 1)*lv1;
                     }
                   }
                }
            }

            dm_stripe[k]     = my_stripe;

        }

      }
    }
}

