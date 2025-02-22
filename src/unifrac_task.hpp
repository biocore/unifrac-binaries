/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2025, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

#ifndef __UNIFRAC_TASKS
#define __UNIFRAC_TASKS 1

#include "unifrac_task_noclass.hpp"
#include "task_parameters.hpp"
#include <math.h>
#include <vector>
#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef SUCMP_NM
/* create a default */
#define SUCMP_NM su_cpu
#endif

// Helper defs to select on SUCMP_NM in preprocessor
#define su_cpu_SUCMP_ID 11
#define su_acc_SUCMP_ID 12
#define su_ompgpu_SUCMP_ID 13

// from https://stackoverflow.com/questions/2335888/how-to-compare-strings-in-c-conditional-preprocessor-directives
#define SUCMP_ID(U) SUCMP_ID_(U)
#define SUCMP_ID_(U) U##_SUCMP_ID

namespace SUCMP_NM {

    // Note: This adds a copy, which is suboptimal
    //       But was the easiest way to get a contiguous buffer
    //       And it does allow for fp32 compute, when desired
    template<class TFloat>
    class UnifracTaskVector {
    private:
      std::vector<double*> &dm_stripes;
      const su::task_parameters* const task_p;

    public:
      const unsigned int start_idx;
      const unsigned int stop_idx;
      const unsigned int n_samples;
      const uint64_t  n_samples_r;
      const uint64_t  bufels;
      TFloat* const buf;

      UnifracTaskVector(std::vector<double*> &_dm_stripes, const su::task_parameters* _task_p)
      : dm_stripes(_dm_stripes), task_p(_task_p)
      , start_idx(task_p->start), stop_idx(task_p->stop), n_samples(task_p->n_samples)
      , n_samples_r(((n_samples + 64-1)/64)*64) // round up to 64 elements (2kbit/4kbits)
      , bufels(n_samples_r * (stop_idx-start_idx))
      , buf((dm_stripes[start_idx]==NULL) ? NULL : (TFloat*) malloc(sizeof(TFloat) * bufels)) // dm_stripes could be null, in which case keep it null
      {
        // keep local copies to avoid the need for *this in the GPU
        const uint64_t  ibufels = bufels;
        TFloat* const ibuf = buf;
        if (ibuf != NULL) {
          for(uint64_t stripe=start_idx; stripe < stop_idx; stripe++) {
             double * dm_stripe = dm_stripes[stripe];
             TFloat * buf_stripe = ibuf+buf_idx(stripe);
             for(uint64_t j=0; j<n_samples; j++) {
                // Note: We could probably just initialize to zero
                buf_stripe[j] = dm_stripe[j];
             }
             for(uint64_t j=n_samples; j<n_samples_r; j++) {
                // Avoid NaNs
                buf_stripe[j] = 0.0;
             }
           }
           acc_copyin_buf(ibuf,0,ibufels);
        }
      }

      UnifracTaskVector(const UnifracTaskVector<TFloat>& ) = delete;
      UnifracTaskVector<TFloat>& operator= (const UnifracTaskVector<TFloat>&) = delete;

      TFloat * operator[](unsigned int idx) { return buf+buf_idx(idx);}
      const TFloat * operator[](unsigned int idx) const { return buf+buf_idx(idx);}


      ~UnifracTaskVector()
      {
        // keep local copies to avoid the need for *this in the GPU
        const uint64_t  ibufels = bufels;
        TFloat* const ibuf = buf;
        if (ibuf != NULL) {
           acc_copyout_buf(ibuf,0,ibufels);
          for(uint64_t stripe=start_idx; stripe < stop_idx; stripe++) {
             double * dm_stripe = dm_stripes[stripe];
             TFloat * buf_stripe = ibuf+buf_idx(stripe);
             for(uint64_t j=0; j<n_samples; j++) {
              dm_stripe[j] = buf_stripe[j];
             }
          }
          free(buf);
        }
      }

    private:
      UnifracTaskVector() = delete;
      UnifracTaskVector operator=(const UnifracTaskVector&other) const = delete;

      uint64_t buf_idx(uint64_t idx) const { return ((idx-start_idx)*n_samples_r);}
    };

    // Base task class to be shared by all tasks
    template<class TFloat, class TEmb>
    class UnifracTaskBase {
      public:
        UnifracTaskVector<TFloat> dm_stripes;
        UnifracTaskVector<TFloat> dm_stripes_total;

        const su::task_parameters* task_p;

        const unsigned int max_embs;
        const uint64_t embsize;
        TFloat * lengths;
       private:
        TEmb * my_embedded_proportions;
        // alternate buffer only needed in async environments, like openacc
        TEmb * my_embedded_proportions_alt; // used as temp
        bool use_alt_emb;

       public:

        UnifracTaskBase(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total, unsigned int _max_embs, const su::task_parameters* _task_p)
        : dm_stripes(_dm_stripes,_task_p), dm_stripes_total(_dm_stripes_total,_task_p), task_p(_task_p)
        , max_embs(_max_embs)
        , embsize(get_embedded_bsize(dm_stripes.n_samples_r,_max_embs))
        , lengths( (TFloat *) malloc(sizeof(TFloat) * _max_embs))
        , my_embedded_proportions((TEmb *) malloc(sizeof(TEmb)*embsize))
        , my_embedded_proportions_alt((TEmb *) NULL)
        , use_alt_emb(false)
        {
            acc_create_buf(lengths,0,max_embs);
            acc_create_buf(my_embedded_proportions,0,embsize);
            if (acc_need_alt()) {
               my_embedded_proportions_alt = (TEmb *) malloc(sizeof(TEmb)*embsize);
               acc_create_buf(my_embedded_proportions_alt,0,embsize);
            }
        }

        UnifracTaskBase(const UnifracTaskBase<TFloat,TEmb>& ) = delete;
        UnifracTaskBase<TFloat,TEmb>& operator= (const UnifracTaskBase<TFloat,TEmb>&) = delete;

        virtual ~UnifracTaskBase()
        {
           if (my_embedded_proportions_alt!=NULL) {
              acc_destroy_buf(my_embedded_proportions_alt,0,embsize);
              free(my_embedded_proportions_alt);
           }

           acc_destroy_buf(my_embedded_proportions,0,embsize);
           acc_destroy_buf(lengths,0,max_embs);

           free(my_embedded_proportions);
           free(lengths);
        }

        TEmb * get_embedded_proportions() {return use_alt_emb ? my_embedded_proportions_alt : my_embedded_proportions;}
        void  set_alt_embedded_proportions() {if (my_embedded_proportions_alt!=NULL) use_alt_emb = !use_alt_emb; /*else , noop */}

        void sync_embedded_proportions(unsigned int filled_embs)
        {
          acc_update_device(this->get_embedded_proportions(),
                            0,dm_stripes.n_samples_r * get_emb_els(filled_embs));
        }

        void sync_lengths(unsigned int filled_embs)
        {
           acc_wait();
           acc_update_device(this->lengths, 0, filled_embs);
        }

        static unsigned int get_emb_els(unsigned int max_embs);

        static uint64_t get_embedded_bsize(const uint64_t  n_samples_r, unsigned int max_embs) {
          const uint64_t bsize = n_samples_r * get_emb_els(max_embs);
          return bsize;
        }

        void embed_proportions_range(const TFloat* __restrict__ in, unsigned int start, unsigned int end, unsigned int emb);
        void embed_proportions(const TFloat* __restrict__ in, unsigned int emb) {embed_proportions_range(in,0,dm_stripes.n_samples,emb);}

        void wait_completion() {
          acc_wait();

        }

        void compute_totals() {
          compute_stripes_totals(this->dm_stripes.buf, this->dm_stripes_total.buf, this->dm_stripes.bufels);
	}

        //
        // ===== Internal, do not use directly =======
        //

        // Just copy from one buffer to another
        // May convert between fp formats in the process (if TOut!=double)
        template<class TOut> void embed_proportions_range_straight(TOut* __restrict__ out, const TFloat* __restrict__ in, unsigned int start, unsigned int end, unsigned int emb) const
        {
          const unsigned int n_samples  = dm_stripes.n_samples;
          const uint64_t n_samples_r  = dm_stripes.n_samples_r;
          const uint64_t offset = emb * n_samples_r;

          for(unsigned int i = start; i < end; i++) {
            out[offset + i] = in[i-start];
          }

          if (end==n_samples) {
            // avoid NaNs
            for(unsigned int i = n_samples; i < n_samples_r; i++) {
              out[offset + i] = 0.0;
            }
          }
        }

        // Just copy from one buffer to another, transpose for faster compute later
        // May convert between fp formats in the process (if TOut!=double)
        template<class TOut> void embed_proportions_range_transpose(TOut* __restrict__ out, const TFloat* __restrict__ in, unsigned int start, unsigned int end, unsigned int emb) const
        {
          const unsigned int n_samples  = dm_stripes.n_samples;
          const uint64_t n_samples_r  = dm_stripes.n_samples_r;
          const uint64_t istripe = max_embs;

          for(unsigned int i = start; i < end; i++) {
            out[i*istripe + emb] = in[i-start];
          }

          if (end==n_samples) {
            // avoid NaNs
            for(unsigned int i = n_samples; i < n_samples_r; i++) {
              out[i*istripe + emb] = 0.0;
            }
          }
        }

        // packed bool
        // Compute (in[:]>0) on each element, and store only the boolean bit.
        // The output values are stored in a multi-byte format, one bit per emb index,
        //    so it will likely take multiple passes to store all the values
        //
        // Note: assumes we are processing emb in increasing order, starting from 0
        template<class TOut> void embed_proportions_range_bool(TOut* __restrict__ out, const TFloat* __restrict__ in, unsigned int start, unsigned int end, unsigned int emb) const
        {
          const unsigned int n_packed = sizeof(TOut)*8;// e.g. 32 for unit32_t
          const unsigned int n_samples  = dm_stripes.n_samples;
          const uint64_t n_samples_r  = dm_stripes.n_samples_r;
          // The output values are stored in a multi-byte format, one bit per emb index
          // Compute the element to store the bit into, as well as whichbit in that element 
          unsigned int emb_block = emb/n_packed; // beginning of the element  block
          unsigned int emb_bit = emb%n_packed;   // bit inside the elements
          const uint64_t offset = emb_block * n_samples_r;

          if  (emb_bit==0) {
            // assign for emb_bit==0, so it clears the other bits
            // assumes we processing emb in increasing order, starting from 0
            for(unsigned int i = start; i < end; i++) {            
              out[offset + i] = (in[i-start] > 0);
            }

            if (end==n_samples) {
              // avoid NaNs
              for(unsigned int i = n_samples; i < n_samples_r; i++) {
                out[offset + i] = 0;
              }
            }
          } else {
            // just update my bit
            for(unsigned int i = start; i < end; i++) {
              out[offset + i] |= (TOut(in[i-start] > 0) << emb_bit);
            }

            // the rest of the els are already OK
          }
        }
        
	// packed bool, transposed
	template<class TOut> void embed_proportions_range_transp_bool(TOut* __restrict__ out, const TFloat* __restrict__ in, unsigned int start, unsigned int end, unsigned int emb) const
        {
          const unsigned int n_packed = sizeof(TOut)*8;// e.g. 32 for unit32_t
          const unsigned int n_samples  = dm_stripes.n_samples;
          const uint64_t n_samples_r  = dm_stripes.n_samples_r;
          // The output values are stored in a multi-byte format, one bit per emb index
          // Compute the element to store the bit into, as well as whichbit in that element 
          const uint64_t istripe = get_emb_els(max_embs);
          unsigned int emb_block = emb/n_packed; // beginning of the element  block
          unsigned int emb_bit = emb%n_packed;   // bit inside the elements

          if  (emb_bit==0) {
            // assign for emb_bit==0, so it clears the other bits
            // assumes we processing emb in increasing order, starting from 0
            for(unsigned int i = start; i < end; i++) {            
              out[i*istripe + emb_block] = (in[i-start] > 0);
            }

            if (end==n_samples) {
              // avoid NaNs
              for(unsigned int i = n_samples; i < n_samples_r; i++) {
                out[i*istripe + emb_block] = 0;
              }
            }
          } else {
            // just update my bit
            for(unsigned int i = start; i < end; i++) {
              out[i*istripe + emb_block] |= (TOut(in[i-start] > 0) << emb_bit);
            }

            // the rest of the els are already OK
          }
        }
    };

#if SUCMP_ID(SUCMP_NM)==su_cpu_SUCMP_ID
    // CPUs better at fine-grained logic, so keeping emb contiguous speeds things up
    // transpose embeded_proportions
    template<> inline void UnifracTaskBase<double,double>::embed_proportions_range(const double* __restrict__ in, unsigned int start, unsigned int end, unsigned int emb) {embed_proportions_range_transpose(get_embedded_proportions(),in,start,end,emb);}
    template<> inline void UnifracTaskBase<double,float>::embed_proportions_range(const double* __restrict__ in, unsigned int start, unsigned int end, unsigned int emb) {embed_proportions_range_transpose(get_embedded_proportions(),in,start,end,emb);}
    template<> inline void UnifracTaskBase<float,float>::embed_proportions_range(const float* __restrict__ in, unsigned int start, unsigned int end, unsigned int emb) {embed_proportions_range_transpose(get_embedded_proportions(),in,start,end,emb);}
#else
    // GPUs need parallelism (while relying on HW masking), so better to keep samples together
    // straight embeded_proportions
    template<> inline void UnifracTaskBase<double,double>::embed_proportions_range(const double* __restrict__ in, unsigned int start, unsigned int end, unsigned int emb) {embed_proportions_range_straight(get_embedded_proportions(),in,start,end,emb);}
    template<> inline void UnifracTaskBase<double,float>::embed_proportions_range(const double* __restrict__ in, unsigned int start, unsigned int end, unsigned int emb) {embed_proportions_range_straight(get_embedded_proportions(),in,start,end,emb);}
    template<> inline void UnifracTaskBase<float,float>::embed_proportions_range(const float* __restrict__ in, unsigned int start, unsigned int end, unsigned int emb) {embed_proportions_range_straight(get_embedded_proportions(),in,start,end,emb);}
#endif
    template<> inline  unsigned int UnifracTaskBase<double,double>::get_emb_els(unsigned int max_embs) {return max_embs;}
    template<> inline  unsigned int UnifracTaskBase<double,float>::get_emb_els(unsigned int max_embs) {return max_embs;}
    template<> inline  unsigned int UnifracTaskBase<float,float>::get_emb_els(unsigned int max_embs) {return max_embs;}

#if SUCMP_ID(SUCMP_NM)==su_cpu_SUCMP_ID
    // CPUs better at fine-grained logic, so keeping emb contiguous speeds things up
    // transposed packed bool embeded_proportions
    template<> inline void UnifracTaskBase<double,uint32_t>::embed_proportions_range(const double* __restrict__ in, unsigned int start, unsigned int end, unsigned int emb) {embed_proportions_range_transp_bool(get_embedded_proportions(),in,start,end,emb);}
    template<> inline void UnifracTaskBase<float,uint32_t>::embed_proportions_range(const float* __restrict__ in, unsigned int start, unsigned int end, unsigned int emb) {embed_proportions_range_transp_bool(get_embedded_proportions(),in,start,end,emb);}

    template<> inline void UnifracTaskBase<double,uint64_t>::embed_proportions_range(const double* __restrict__ in, unsigned int start, unsigned int end, unsigned int emb) {embed_proportions_range_transp_bool(get_embedded_proportions(),in,start,end,emb);}
    template<> inline void UnifracTaskBase<float,uint64_t>::embed_proportions_range(const float* __restrict__ in, unsigned int start, unsigned int end, unsigned int emb) {embed_proportions_range_transp_bool(get_embedded_proportions(),in,start,end,emb);}
#else
    // GPUs need parallelism (while relying on HW masking), so better to keep samples together
    //packed bool embeded_proportions
    template<> inline void UnifracTaskBase<double,uint32_t>::embed_proportions_range(const double* __restrict__ in, unsigned int start, unsigned int end, unsigned int emb) {embed_proportions_range_bool(get_embedded_proportions(),in,start,end,emb);}
    template<> inline void UnifracTaskBase<float,uint32_t>::embed_proportions_range(const float* __restrict__ in, unsigned int start, unsigned int end, unsigned int emb) {embed_proportions_range_bool(get_embedded_proportions(),in,start,end,emb);}

    template<> inline void UnifracTaskBase<double,uint64_t>::embed_proportions_range(const double* __restrict__ in, unsigned int start, unsigned int end, unsigned int emb) {embed_proportions_range_bool(get_embedded_proportions(),in,start,end,emb);}
    template<> inline void UnifracTaskBase<float,uint64_t>::embed_proportions_range(const float* __restrict__ in, unsigned int start, unsigned int end, unsigned int emb) {embed_proportions_range_bool(get_embedded_proportions(),in,start,end,emb);}
#endif
    template<> inline  unsigned int UnifracTaskBase<double,uint32_t>::get_emb_els(unsigned int max_embs) {return (max_embs+31)/32;}
    template<> inline  unsigned int UnifracTaskBase<float,uint32_t>::get_emb_els(unsigned int max_embs) {return (max_embs+31)/32;}
    template<> inline  unsigned int UnifracTaskBase<double,uint64_t>::get_emb_els(unsigned int max_embs) {return (max_embs+63)/64;}
    template<> inline  unsigned int UnifracTaskBase<float,uint64_t>::get_emb_els(unsigned int max_embs) {return (max_embs+63)/64;}

    /* void unifrac tasks
     *
     * all methods utilize the same function signature. that signature is as follows:
     *
     * dm_stripes vector<double> the stripes of the distance matrix being accumulated 
     *      into for unique branch length
     * dm_stripes vector<double> the stripes of the distance matrix being accumulated 
     *      into for total branch length (e.g., to normalize unweighted unifrac)
     * embedded_proportions <double*> the proportions vector for a sample, or rather
     *      the counts vector normalized to 1. this vector is embedded as it is 
     *      duplicated: if A, B and C are proportions for features A, B, and C, the
     *      vector will look like [A B C A B C].
     * length <double> the branch length of the current node to its parent.
     * task_p <task_parameters*> task specific parameters.
     */

    template<class TFloat, class TEmb>
    class UnifracTask : public UnifracTaskBase<TFloat,TEmb> {
      public:

        UnifracTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total, unsigned int _max_embs, const su::task_parameters* _task_p)
        : UnifracTaskBase<TFloat,TEmb>(_dm_stripes, _dm_stripes_total, _max_embs, _task_p) {}

       UnifracTask(const UnifracTask<TFloat,TEmb>& ) = delete;
       UnifracTask<TFloat,TEmb>& operator= (const UnifracTask<TFloat,TEmb>&) = delete;

       virtual ~UnifracTask() {}

       virtual void run(unsigned int filled_embs) = 0;

      protected:
#if SUCMP_ID(SUCMP_NM)==su_cpu_SUCMP_ID
       // size for contiguous access
       static constexpr unsigned int RECOMMENDED_MAX_EMBS_STRAIGHT = 512-16; // optimize for 32k CPU L1 cache... -16 to avoid cache line trashing
       // Must be multiple of 64, since it it bit packed
       // The EMB is 64x smaller, and still msut fit in L1 cache
       // But must also fit sums in L2 cache, since access is pseudo-random (and too big for L1)
       static constexpr unsigned int RECOMMENDED_MAX_EMBS_BOOL = 64*32;  // Optimize for 1M L2 cache... (using 1/4 of it)
#else
       // size for vectorized access
       static constexpr unsigned int RECOMMENDED_MAX_EMBS_STRAIGHT = 128-16; // optimize for GPU L1 cache... -16 to avoid cache line trashing
       // Must be multiple of 64, since it it bit packed
       // The EMB is 64x smaller, and still msut fit in L1 cache
       // But must also fit sums in L2 cache, since access is pseudo-random (and too big for L1)
       static constexpr unsigned int RECOMMENDED_MAX_EMBS_BOOL = 64*32;
#endif

    };


    template<class TFloat>
    class UnifracUnnormalizedWeightedTask : public UnifracTask<TFloat,TFloat> {
      public:
        static constexpr unsigned int RECOMMENDED_MAX_EMBS = UnifracTask<TFloat,TFloat>::RECOMMENDED_MAX_EMBS_STRAIGHT;

        UnifracUnnormalizedWeightedTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total, unsigned int _max_embs, const su::task_parameters* _task_p)
        : UnifracTask<TFloat,TFloat>(_dm_stripes,_dm_stripes_total,_max_embs,_task_p)
        {
          const unsigned int n_samples = this->task_p->n_samples;

          zcheck = (bool*) malloc(sizeof(bool) * n_samples);
          sums = (TFloat*) malloc(sizeof(TFloat) * n_samples);
          acc_create_buf(zcheck, 0, n_samples);
          acc_create_buf(sums, 0 , n_samples);
        }

        UnifracUnnormalizedWeightedTask(const UnifracUnnormalizedWeightedTask<TFloat>& ) = delete;
        UnifracUnnormalizedWeightedTask<TFloat>& operator= (const UnifracUnnormalizedWeightedTask<TFloat>&) = delete;

        virtual ~UnifracUnnormalizedWeightedTask()
        {
          const unsigned int n_samples = this->task_p->n_samples;

          acc_destroy_buf(sums, 0 , n_samples);
          acc_destroy_buf(zcheck, 0, n_samples);

          free(sums);
          free(zcheck);
        }

        virtual void run(unsigned int filled_embs) {_run(filled_embs);}

        void _run(unsigned int filled_embs) {
           run_UnnormalizedWeightedTask(
			  this->max_embs, filled_embs,
			  this->task_p->start, this->task_p->stop, this->task_p->n_samples, this->dm_stripes.n_samples_r,
			  this->lengths,  this->get_embedded_proportions(), this->dm_stripes.buf,
			  this->zcheck, this->sums);

           // next iteration will use the alternative space
           this->set_alt_embedded_proportions();
	}
      protected:
        // temp buffers
        bool     *zcheck;
        TFloat   *sums;
    };
    template<class TFloat>
    class UnifracNormalizedWeightedTask : public UnifracTask<TFloat,TFloat> {
      public:
        static constexpr unsigned int RECOMMENDED_MAX_EMBS = UnifracTask<TFloat,TFloat>::RECOMMENDED_MAX_EMBS_STRAIGHT;

        UnifracNormalizedWeightedTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total, unsigned int _max_embs, const su::task_parameters* _task_p)
        : UnifracTask<TFloat,TFloat>(_dm_stripes,_dm_stripes_total,_max_embs,_task_p)
        {
          const unsigned int n_samples = this->task_p->n_samples;

          zcheck = (bool*) malloc(sizeof(bool) * n_samples);
          sums = (TFloat*) malloc(sizeof(TFloat) * n_samples);
          acc_create_buf(zcheck, 0, n_samples);
          acc_create_buf(sums, 0 , n_samples);
        }

        UnifracNormalizedWeightedTask(const UnifracNormalizedWeightedTask<TFloat>& ) = delete;
        UnifracNormalizedWeightedTask<TFloat>& operator= (const UnifracNormalizedWeightedTask<TFloat>&) = delete;

        virtual ~UnifracNormalizedWeightedTask()
        {
          const unsigned int n_samples = this->task_p->n_samples;

          acc_destroy_buf(sums, 0 , n_samples);
          acc_destroy_buf(zcheck, 0, n_samples);

          free(sums);
          free(zcheck);
        }

        virtual void run(unsigned int filled_embs) {_run(filled_embs);}

        void _run(unsigned int filled_embs) {
          run_NormalizedWeightedTask(
			    this->max_embs, filled_embs,
			    this->task_p->start, this->task_p->stop, this->task_p->n_samples, this->dm_stripes.n_samples_r,
			    this->lengths, this->get_embedded_proportions(),
			    this->dm_stripes.buf, this->dm_stripes_total.buf,
			    this->zcheck, this->sums);

          // next iteration will use the alternative space
          this->set_alt_embedded_proportions();
	}
      protected:
        // temp buffers
        bool     *zcheck;
        TFloat   *sums;
    };
    template<class TFloat>
    class UnifracCommonUnweightedTask : public UnifracTask<TFloat,uint64_t> {
      public:
        static constexpr unsigned int RECOMMENDED_MAX_EMBS = UnifracTask<TFloat,uint64_t>::RECOMMENDED_MAX_EMBS_BOOL;

        // Note: _max_emb MUST be multiple of 64
        UnifracCommonUnweightedTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total, unsigned int _max_embs, const su::task_parameters* _task_p)
        : UnifracTask<TFloat, uint64_t>(_dm_stripes,_dm_stripes_total,_max_embs,_task_p) 
        {
          const unsigned int n_samples = this->task_p->n_samples;
          const unsigned int bsize = _max_embs*(0x400/32);
          zcheck = (bool*) malloc(sizeof(bool) * n_samples);
          stripe_sums = (TFloat*) malloc(sizeof(TFloat) *  n_samples);
          sums = (TFloat*) malloc(sizeof(TFloat) * bsize);

          acc_create_buf(zcheck, 0, n_samples);
          acc_create_buf(stripe_sums, 0 , n_samples);
          acc_create_buf(sums, 0 , bsize);
#if SUCMP_ID(SUCMP_NM)==su_cpu_SUCMP_ID
	  idxs = NULL; // not used in the CPU code
#else
          idxs = (uint32_t*) malloc(sizeof(uint32_t) * n_samples);
          acc_create_buf(idxs, 0, n_samples);
#endif
        }

        virtual ~UnifracCommonUnweightedTask()
        {
          const unsigned int n_samples = this->task_p->n_samples;
          const unsigned int bsize = this->max_embs*(0x400/32);

#if SUCMP_ID(SUCMP_NM)==su_cpu_SUCMP_ID
	  // not allocated in the CPU code
#else
          acc_destroy_buf(idxs, 0, n_samples);
          free(idxs);
#endif
          acc_destroy_buf(sums, 0 , bsize);
          acc_destroy_buf(stripe_sums, 0 , n_samples);
          acc_destroy_buf(zcheck, 0, n_samples);

          free(sums);
          free(stripe_sums);
          free(zcheck);
        }

      protected:
        // temp buffers
        TFloat *sums;
        bool     *zcheck;
        uint32_t *idxs;  // assuming n_samples si really a uint32_t number
        TFloat   *stripe_sums;
    };
    template<class TFloat>
    class UnifracUnweightedTask : public UnifracCommonUnweightedTask<TFloat> {
      public:

        // Note: _max_emb MUST be multiple of 64
        UnifracUnweightedTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total, unsigned int _max_embs, const su::task_parameters* _task_p)
        : UnifracCommonUnweightedTask<TFloat>(_dm_stripes,_dm_stripes_total,_max_embs,_task_p) {}

        UnifracUnweightedTask(const UnifracUnweightedTask<TFloat>& ) = delete;
        UnifracUnweightedTask<TFloat>& operator= (const UnifracUnweightedTask<TFloat>&) = delete;

        virtual void run(unsigned int filled_embs) {_run(filled_embs);}

        void _run(unsigned int filled_embs) {
          run_UnweightedTask(
			 this->get_emb_els(this->max_embs), filled_embs,
			 this->task_p->start, this->task_p->stop, this->task_p->n_samples, this->dm_stripes.n_samples_r,
			 this->lengths, this->get_embedded_proportions(), this->dm_stripes.buf, this->dm_stripes_total.buf,
			 this->sums, this->zcheck, this->idxs, this->stripe_sums);

          // next iteration will use the alternative space
          this->set_alt_embedded_proportions();
	}
    };
    template<class TFloat>
    class UnifracUnnormalizedUnweightedTask : public UnifracCommonUnweightedTask<TFloat> {
      public:

        // Note: _max_emb MUST be multiple of 64
        UnifracUnnormalizedUnweightedTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total, unsigned int _max_embs, const su::task_parameters* _task_p)
        : UnifracCommonUnweightedTask<TFloat>(_dm_stripes,_dm_stripes_total,_max_embs,_task_p) {}

        UnifracUnnormalizedUnweightedTask(const UnifracUnnormalizedUnweightedTask<TFloat>& ) = delete;
        UnifracUnnormalizedUnweightedTask<TFloat>& operator= (const UnifracUnnormalizedUnweightedTask<TFloat>&) = delete;

        virtual void run(unsigned int filled_embs) {_run(filled_embs);}

        void _run(unsigned int filled_embs) {
          run_UnnormalizedUnweightedTask(
			  this->get_emb_els(this->max_embs), filled_embs,
			  this->task_p->start, this->task_p->stop, this->task_p->n_samples, this->dm_stripes.n_samples_r,
			  this->lengths, this->get_embedded_proportions(),
			  this->dm_stripes.buf,
			  this->sums, this->zcheck, this->idxs, this->stripe_sums);

          // next iteration will use the alternative space
          this->set_alt_embedded_proportions();
	}
    };

    template<class TFloat>
    class UnifracGeneralizedTask : public UnifracTask<TFloat,TFloat> {
      public:
        static constexpr unsigned int RECOMMENDED_MAX_EMBS = UnifracTask<TFloat,TFloat>::RECOMMENDED_MAX_EMBS_STRAIGHT;

        UnifracGeneralizedTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total, unsigned int _max_embs, const su::task_parameters* _task_p)
        : UnifracTask<TFloat,TFloat>(_dm_stripes,_dm_stripes_total,_max_embs,_task_p) {}

        UnifracGeneralizedTask(const UnifracGeneralizedTask<TFloat>& ) = delete;
        UnifracGeneralizedTask<TFloat>& operator= (const UnifracGeneralizedTask<TFloat>&) = delete;

        virtual void run(unsigned int filled_embs) {_run(filled_embs);}

        void _run(unsigned int filled_embs) {
          run_GeneralizedTask(
			  this->max_embs, filled_embs,
			  this->task_p->start, this->task_p->stop, this->task_p->n_samples, this->dm_stripes.n_samples_r,
			  this->lengths, this->get_embedded_proportions(),
			  this->dm_stripes.buf, this->dm_stripes_total.buf,
			  (TFloat) this->task_p->g_unifrac_alpha);

          // next iteration will use the alternative space
          this->set_alt_embedded_proportions();
	}
    };

    /* void unifrac_vaw tasks
     *
     * all methods utilize the same function signature. that signature is as follows:
     *
     * dm_stripes vector<double> the stripes of the distance matrix being accumulated 
     *      into for unique branch length
     * dm_stripes vector<double> the stripes of the distance matrix being accumulated 
     *      into for total branch length (e.g., to normalize unweighted unifrac)
     * embedded_proportions <double*> the proportions vector for a sample, or rather
     *      the counts vector normalized to 1. this vector is embedded as it is 
     *      duplicated: if A, B and C are proportions for features A, B, and C, the
     *      vector will look like [A B C A B C].
     * embedded_counts <double*> the counts vector embedded in the same way and order as
     *      embedded_proportions. the values of this array are unnormalized feature 
     *      counts for the subtree.
     * sample_total_counts <double*> the total unnormalized feature counts for all samples
     *      embedded in the same way and order as embedded_proportions.
     * length <double> the branch length of the current node to its parent.
     * task_p <task_parameters*> task specific parameters.
     */
    template<class TFloat, class TEmb>
    class UnifracVawTask : public UnifracTaskBase<TFloat,TEmb> {
      public:
        TFloat * const embedded_counts;
        const TFloat * const sample_total_counts;

        static constexpr unsigned int RECOMMENDED_MAX_EMBS = 128;

        UnifracVawTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total, 
                    const double * _sample_counts,
                    unsigned int _max_embs, const su::task_parameters* _task_p)
        : UnifracTaskBase<TFloat,TEmb>(_dm_stripes, _dm_stripes_total, _max_embs, _task_p)
        , embedded_counts((TFloat *) malloc(sizeof(TFloat)*this->embsize))
        , sample_total_counts(initialize_sample_counts(this->dm_stripes.n_samples, this->dm_stripes.n_samples_r, _task_p, _sample_counts))
        {
          acc_create_buf(embedded_counts, 0, this->embsize);
          acc_copyin_buf(const_cast<TFloat *>(sample_total_counts), 0 , this->dm_stripes.n_samples_r); // const after the contructor
        }

       UnifracVawTask(const UnifracVawTask<TFloat,TEmb>& ) = delete;
       UnifracVawTask<TFloat,TEmb>& operator= (const UnifracVawTask<TFloat,TEmb>&) = delete;

       virtual ~UnifracVawTask() 
       {
          acc_destroy_buf(const_cast<TFloat *>(sample_total_counts), 0 , this->dm_stripes.n_samples_r);
          acc_destroy_buf(embedded_counts, 0, this->embsize);

          free(const_cast<TFloat *>(sample_total_counts)); // while const for the life of this, not const past its lifetime
          free(embedded_counts);
       }

       static TFloat* initialize_sample_counts(uint64_t n_samples, uint64_t n_samples_r, const su::task_parameters* task_p, const double sample_counts[]) {
          TFloat * counts = (TFloat *) malloc(sizeof(TFloat) * n_samples_r);
          for(unsigned int i = 0; i < n_samples; i++) {
            counts[i] = sample_counts[i];
          }
          // avoid NaNs
          for(unsigned int i = n_samples; i < n_samples_r; i++) {
            counts[i] = 0.0;
          }

          return counts;
       }

       void sync_embedded_counts(unsigned int filled_embs)
       {
          acc_update_device(this->embedded_counts, 0, this->dm_stripes.n_samples_r * filled_embs);
       }

       void sync_embedded(unsigned int filled_embs) { this->sync_embedded_proportions(filled_embs); this->sync_embedded_counts(filled_embs);}

        void embed_range(const TFloat* __restrict__ in_proportions, const TFloat* __restrict__ in_counts, unsigned int start, unsigned int end, unsigned int emb);
        void embed(const TFloat* __restrict__ in_proportions, const double* __restrict__ in_counts, unsigned int emb) { embed_range(in_proportions,in_counts,0,this->dm_stripes.n_samples,emb);}

       virtual void run(unsigned int filled_embs) = 0;
    };

    // straight embeded_proportions
    template<> inline void UnifracVawTask<double,double>::embed_range(const double* __restrict__ in_proportions, const double* __restrict__ in_counts, unsigned int start, unsigned int end, unsigned int emb) {
          this->embed_proportions_range_straight(this->get_embedded_proportions(),in_proportions,start,end,emb);
          this->embed_proportions_range_straight(this->embedded_counts,in_counts,start,end,emb);
    }
    template<> inline void UnifracVawTask<double,float>::embed_range(const double* __restrict__ in_proportions, const double* __restrict__ in_counts, unsigned int start, unsigned int end, unsigned int emb) {
          this->embed_proportions_range_straight(this->get_embedded_proportions(),in_proportions,start,end,emb);
          this->embed_proportions_range_straight(this->embedded_counts,in_counts,start,end,emb);
    }
    template<> inline void UnifracVawTask<float,float>::embed_range(const float* __restrict__ in_proportions, const float* __restrict__ in_counts, unsigned int start, unsigned int end, unsigned int emb) {
          this->embed_proportions_range_straight(this->get_embedded_proportions(),in_proportions,start,end,emb);
          this->embed_proportions_range_straight(this->embedded_counts,in_counts,start,end,emb);
    }

    //packed bool embeded_proportions
    template<> inline void UnifracVawTask<double,uint32_t>::embed_range(const double* __restrict__ in_proportions, const double* __restrict__ in_counts, unsigned int start, unsigned int end, unsigned int emb) {
          this->embed_proportions_range_bool(this->get_embedded_proportions(),in_proportions,start,end,emb);
          this->embed_proportions_range_straight(this->embedded_counts,in_counts,start,end,emb);
    }
    template<> inline void UnifracVawTask<double,uint64_t>::embed_range(const double* __restrict__ in_proportions, const double* __restrict__ in_counts, unsigned int start, unsigned int end, unsigned int emb) {
          this->embed_proportions_range_bool(this->get_embedded_proportions(),in_proportions,start,end,emb);
          this->embed_proportions_range_straight(this->embedded_counts,in_counts,start,end,emb);
    }
    template<> inline void UnifracVawTask<float,uint32_t>::embed_range(const float* __restrict__ in_proportions, const float* __restrict__ in_counts, unsigned int start, unsigned int end, unsigned int emb) {
          this->embed_proportions_range_bool(this->get_embedded_proportions(),in_proportions,start,end,emb);
          this->embed_proportions_range_straight(this->embedded_counts,in_counts,start,end,emb);
    }
    template<> inline void UnifracVawTask<float,uint64_t>::embed_range(const float* __restrict__ in_proportions, const float* __restrict__ in_counts, unsigned int start, unsigned int end, unsigned int emb) {
          this->embed_proportions_range_bool(this->get_embedded_proportions(),in_proportions,start,end,emb);
          this->embed_proportions_range_straight(this->embedded_counts,in_counts,start,end,emb);
    }

    template<class TFloat>
    class UnifracVawUnnormalizedWeightedTask : public UnifracVawTask<TFloat,TFloat> {
      public:
        UnifracVawUnnormalizedWeightedTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total, 
                    const double * _sample_counts,
                    unsigned int _max_embs, const su::task_parameters* _task_p)
        : UnifracVawTask<TFloat,TFloat>(_dm_stripes,_dm_stripes_total,_sample_counts,_max_embs,_task_p) {}

        UnifracVawUnnormalizedWeightedTask(const UnifracVawUnnormalizedWeightedTask<TFloat>& ) = delete;
        UnifracVawUnnormalizedWeightedTask<TFloat>& operator= (const UnifracVawUnnormalizedWeightedTask<TFloat>&) = delete;

        virtual void run(unsigned int filled_embs) {_run(filled_embs);}

        void _run(unsigned int filled_embs) {
           run_VawUnnormalizedWeightedTask(
			   filled_embs,
			   this->task_p->start, this->task_p->stop, this->task_p->n_samples, this->dm_stripes.n_samples_r,
			   this->lengths, this->get_embedded_proportions(), this->embedded_counts, this->sample_total_counts,
			   this->dm_stripes.buf);

           // next iteration will use the alternative space
           this->set_alt_embedded_proportions();
	}
    };
    template<class TFloat>
    class UnifracVawNormalizedWeightedTask : public UnifracVawTask<TFloat,TFloat> {
      public:
        UnifracVawNormalizedWeightedTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total, 
                    const double * _sample_counts,
                    unsigned int _max_embs, const su::task_parameters* _task_p)
        : UnifracVawTask<TFloat,TFloat>(_dm_stripes,_dm_stripes_total,_sample_counts,_max_embs,_task_p) {}

        UnifracVawNormalizedWeightedTask(const UnifracVawNormalizedWeightedTask<TFloat>& ) = delete;
        UnifracVawNormalizedWeightedTask<TFloat>& operator= (const UnifracVawNormalizedWeightedTask<TFloat>&) = delete;

        virtual void run(unsigned int filled_embs) {_run(filled_embs);}

        void _run(unsigned int filled_embs) {
          run_VawNormalizedWeightedTask(
			  filled_embs,
			  this->task_p->start, this->task_p->stop, this->task_p->n_samples, this->dm_stripes.n_samples_r,
			  this->lengths, this->get_embedded_proportions(), this->embedded_counts, this->sample_total_counts,
			  this->dm_stripes.buf, this->dm_stripes_total.buf);

          // next iteration will use the alternative space
          this->set_alt_embedded_proportions();
	}
    };
    template<class TFloat>
    class UnifracVawUnweightedTask : public UnifracVawTask<TFloat,uint32_t> {
      public:
        UnifracVawUnweightedTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total, 
                    const double * _sample_counts,
                    unsigned int _max_embs, const su::task_parameters* _task_p)
        : UnifracVawTask<TFloat,uint32_t>(_dm_stripes,_dm_stripes_total,_sample_counts,_max_embs,_task_p) {}

        UnifracVawUnweightedTask(const UnifracVawUnweightedTask<TFloat>& ) = delete;
        UnifracVawUnweightedTask<TFloat>& operator= (const UnifracVawUnweightedTask<TFloat>&) = delete;

        virtual void run(unsigned int filled_embs) {_run(filled_embs);}

        void _run(unsigned int filled_embs) {
          run_VawUnweightedTask(
			  filled_embs,
			  this->task_p->start, this->task_p->stop, this->task_p->n_samples, this->dm_stripes.n_samples_r,
			  this->lengths, this->get_embedded_proportions(), this->embedded_counts, this->sample_total_counts,
			  this->dm_stripes.buf, this->dm_stripes_total.buf);

          // next iteration will use the alternative space
          this->set_alt_embedded_proportions();
	}
    };
    template<class TFloat>
    class UnifracVawUnnormalizedUnweightedTask : public UnifracVawTask<TFloat,uint32_t> {
      public:
        UnifracVawUnnormalizedUnweightedTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total, 
                    const double * _sample_counts,
                    unsigned int _max_embs, const su::task_parameters* _task_p)
        : UnifracVawTask<TFloat,uint32_t>(_dm_stripes,_dm_stripes_total,_sample_counts,_max_embs,_task_p) {}

        UnifracVawUnnormalizedUnweightedTask(const UnifracVawUnweightedTask<TFloat>& ) = delete;
        UnifracVawUnnormalizedUnweightedTask<TFloat>& operator= (const UnifracVawUnweightedTask<TFloat>&) = delete;

        virtual void run(unsigned int filled_embs) {_run(filled_embs);}

        void _run(unsigned int filled_embs) {
          run_VawUnnormalizedUnweightedTask(
			  filled_embs,
			  this->task_p->start, this->task_p->stop, this->task_p->n_samples, this->dm_stripes.n_samples_r,
			  this->lengths, this->get_embedded_proportions(), this->embedded_counts, this->sample_total_counts,
			  this->dm_stripes.buf);

          // next iteration will use the alternative space
          this->set_alt_embedded_proportions();
	}
    };
    template<class TFloat>
    class UnifracVawGeneralizedTask : public UnifracVawTask<TFloat,TFloat> {
      public:
        UnifracVawGeneralizedTask(std::vector<double*> &_dm_stripes, std::vector<double*> &_dm_stripes_total,
                    const double * _sample_counts,
                    unsigned int _max_embs, const su::task_parameters* _task_p)
        : UnifracVawTask<TFloat,TFloat>(_dm_stripes,_dm_stripes_total,_sample_counts,_max_embs,_task_p) {}

        UnifracVawGeneralizedTask(const UnifracVawGeneralizedTask<TFloat>& ) = delete;
        UnifracVawGeneralizedTask<TFloat>& operator= (const UnifracVawGeneralizedTask<TFloat>&) = delete;

        virtual void run(unsigned int filled_embs) {_run(filled_embs);}

        void _run(unsigned int filled_embs) {
          run_VawGeneralizedTask(
			  filled_embs,
			  this->task_p->start, this->task_p->stop, this->task_p->n_samples, this->dm_stripes.n_samples_r,
			  this->lengths, this->get_embedded_proportions(), this->embedded_counts,this->sample_total_counts,
			  this->dm_stripes.buf, this->dm_stripes_total.buf,
			  (TFloat) this->task_p->g_unifrac_alpha);

          // next iteration will use the alternative space
          this->set_alt_embedded_proportions();
	}
    };

}

#endif
