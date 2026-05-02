/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2023-2023, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <random>
#include <vector>
#include <algorithm>
#include "omp_stub.h"

#include "biom_subsampled.hpp"

using namespace su;

linked_sparse_transposed::linked_sparse_transposed(sparse_data &other) 
  : n_obs(other.n_samples)
  , n_samples(other.n_obs) {
   obs_counts_resident = (unsigned int*)calloc(sizeof(unsigned int), n_obs);
   if(obs_counts_resident == NULL) {
        fprintf(stderr, "Failed to allocate %zd bytes; [%s]:%d\n", 
                sizeof(unsigned int) * n_obs, __FILE__, __LINE__);
        exit(EXIT_FAILURE);
   }

   for (uint32_t i=0; i<other.n_obs; i++) {
     const uint32_t cnt = other.obs_counts_resident[i];
     const uint32_t * idxs = other.obs_indices_resident[i];
     for (uint32_t j=0; j<cnt; j++) {
       obs_counts_resident[idxs[j]]++;
     }
   }

   obs_data_resident = (double***)malloc(sizeof(double**) * n_obs);
   if(obs_data_resident == NULL) {
        fprintf(stderr, "Failed to allocate %zd bytes; [%s]:%d\n", 
                sizeof(double**) * n_obs, __FILE__, __LINE__);
        exit(EXIT_FAILURE);
   }
   uint32_t _max_count = 0;
   for (uint32_t i=0; i<n_obs; i++) {
     const uint32_t cnt = obs_counts_resident[i];
     _max_count = std::max(_max_count,cnt);
     obs_data_resident[i] = (double**)malloc(sizeof(double*) * cnt);
     if(obs_data_resident[i] == NULL) {
        fprintf(stderr, "Failed to allocate %zd bytes; [%s]:%d\n", 
                sizeof(double*) * cnt, __FILE__, __LINE__);
        exit(EXIT_FAILURE);
     }
   }
   max_count = _max_count;

   for (uint32_t i=0; i<other.n_obs; i++) {
     const uint32_t cnt = other.obs_counts_resident[i];
     const uint32_t *idxs = other.obs_indices_resident[i];
     double *data = other.obs_data_resident[i];
     for (uint32_t j=0; j<cnt; j++) {
       uint32_t local_i = idxs[j];
       obs_data_resident[local_i][0] = &(data[j]);
       obs_data_resident[local_i]++;
     }
   }

   // we moved the pointers in place, rewind
   for (uint32_t i=0; i<n_obs; i++) {
      const uint32_t cnt = obs_counts_resident[i];
      obs_data_resident[i] -=  cnt;
   }
}

linked_sparse_transposed::~linked_sparse_transposed() {
    if(obs_data_resident != NULL) {
        for (uint32_t i=0; i<n_obs; i++) free(obs_data_resident[i]);
        free(obs_data_resident);
        obs_data_resident = NULL;
    }
    if(obs_counts_resident != NULL) {
       free(obs_counts_resident);
        obs_counts_resident = NULL;
    }
}

namespace su {
  // Equivalent to iterator over np.repeat
  // https://github.com/biocore/biom-format/blob/b0e71a00ecb349a6f5f1ca64a23d71f380ddc19c/biom/_subsample.pyx#LL64C24-L64C55
  class WeightedSampleIterator
  {
  public:
    // While we do not implememnt the whole random_access_iterator interface
    // we want the implementations to use operator- and that requires random
    using iterator_category = std::random_access_iterator_tag;
    using difference_type   = int64_t;
    using value_type        = uint32_t;
    using pointer           = const uint32_t*;
    using reference         = const uint32_t&;

    WeightedSampleIterator(uint64_t *_data_in, uint32_t _idx, uint64_t _cnt)
    : data_in(_data_in)
    , idx(_idx)
    , cnt(_cnt)
    {}

    reference operator*() const { return idx; }
    pointer operator->() const { return &idx; }

    WeightedSampleIterator& operator++()
    {  
       cnt++;
       if (cnt>=data_in[idx]) {
         cnt = 0;
         idx++;
       }
       return *this;
    }

    WeightedSampleIterator operator++(int) { WeightedSampleIterator tmp = *this; ++(*this); return tmp; }

    friend bool operator== (const WeightedSampleIterator& a, const WeightedSampleIterator& b)
    {
       return (a.data_in == b.data_in) && (a.idx == b.idx) && (a.cnt==b.cnt);
    };

    friend bool operator!= (const WeightedSampleIterator& a, const WeightedSampleIterator& b)
    {
       return !((a.data_in == b.data_in) && (a.idx == b.idx) && (a.cnt==b.cnt));
    };

    friend int64_t operator-(const WeightedSampleIterator& b, const WeightedSampleIterator& a)
    {
       int64_t diff = 0;
       //assert(a.data_in == b.data_in);
       //assert(a.idx <= b.idx);
       //assert((a.idx > b.idx) || (a.cnt<=b.cnt));

       for (uint32_t i = a.idx; i<b.idx; i++) {
          diff += a.data_in[i];
       }

       return diff + b.cnt - a.cnt;
    };

  private:

    uint64_t *data_in;
    uint32_t idx; // index of data_in
    uint64_t cnt; // how deep in data_in[idx] are we (must be < data_in[idx])
  };

  class WeightedSample
  {
  public:
    WeightedSample(uint32_t _max_count, uint32_t _n, uint32_t random_seed)
    : max_count(_max_count)
    , n(_n)
    , generator(random_seed)
    , current_count(0)
    , data(max_count)
    , sample_out(n)
    , data_out(max_count)
    {}

   void do_sample(unsigned int length, double* *data_arr) {
        for (unsigned int j=0; j<length; j++) data_out[j] = 0;

        // note: We are assuming length>=n
        //      Enforced by the caller (via filtering)
        assign(length,data_arr);
        std::sample(begin(), end(),
                    sample_out.begin(), n,
                    generator);

        for (uint32_t j=0; j<n; j++) data_out[sample_out[j]]++;

        for (unsigned int j=0; j<length; j++) *(data_arr[j]) = data_out[j];
    }

private:
    void assign(uint32_t length, double * *data_arr) {
       current_count = length;
       for (uint32_t j=0; j<length; j++) data[j] = *(data_arr[j]);
    }

    WeightedSampleIterator begin() { return WeightedSampleIterator(data.data(),0,0); }
    WeightedSampleIterator end()   { return WeightedSampleIterator(data.data(),current_count,0); } // current_count is out of bounds

    uint32_t max_count;
    uint32_t n;
    std::mt19937 generator;

    uint32_t current_count;

    // use persistent buffer to minimize allocation costs
    std::vector<uint64_t> data;  // original values
    std::vector<uint32_t> sample_out;     // random output buffer
    std::vector<uint32_t> data_out; // computed values
  };

  class WeightedSampleWithReplacement
  {
  public:
    WeightedSampleWithReplacement(uint32_t _max_count, uint32_t _n, uint32_t random_seed)
    : max_count(_max_count)
    , n(_n)
    , generator(random_seed)
    , current_count(0)
    , data(max_count)
    , data_out(max_count)
    {}

   void do_sample(unsigned int length, double* *data_arr) {
        assign(length,data_arr);

        // note: We are assuming length>=n
        //      Enforced by the caller (via filtering)
        double *data_in = data.data();
        std::discrete_distribution<uint32_t> multinomial(data_in, data_in+length);
        for (unsigned int j=0; j<length; j++) data_out[j] = 0;
        for (uint32_t j=0; j<n; j++) data_out[multinomial(generator)]++;

        for (unsigned int j=0; j<length; j++) *(data_arr[j]) = data_out[j];
    }

private:
    void assign(uint32_t length, double * *data_arr) {
       current_count = length;
       for (uint32_t j=0; j<length; j++) data[j] = *(data_arr[j]);
    }

    uint32_t max_count;
    uint32_t n;
    std::mt19937 generator;

    uint32_t current_count;

    // use persistent buffer to minimize allocation costs
    std::vector<double> data;  // original values
    std::vector<uint32_t> data_out; // computed values
  };
} // end namespace su


template<class TWork>
inline void linked_sparse_transposed::transposed_subsample(const uint32_t n, const uint32_t random_seed) {
    const uint32_t max_threads = omp_get_max_threads();
    std::mt19937 master_generator(random_seed);

    // use common buffer to minimize allocation cost, but need one per thread
    std::vector<TWork> sample_data_arr;
    for (uint32_t i=0; i<max_threads; i++) sample_data_arr.emplace(sample_data_arr.end(), max_count, n, master_generator());
    
    #pragma omp parallel for
    for (uint32_t i=0; i<n_obs; i++) {
        int my_thread_num = omp_get_thread_num();
        sample_data_arr[my_thread_num].do_sample(obs_counts_resident[i], obs_data_resident[i]);
    }
}

void linked_sparse_transposed::transposed_subsample_without_replacement(const uint32_t n, const uint32_t random_seed) {
    transposed_subsample<WeightedSample>(n, random_seed);
}

void linked_sparse_transposed::transposed_subsample_with_replacement(const uint32_t n, const uint32_t random_seed) {
    transposed_subsample<WeightedSampleWithReplacement>(n, random_seed);
}


// =====================  sparse_data_subsampled  ==========================

void sparse_data_subsampled::subsample_with_replacement(const uint32_t n, const uint32_t random_seed) {
    linked_sparse_transposed transposed(*this);
    transposed.transposed_subsample_with_replacement(n,random_seed);
}

void sparse_data_subsampled::subsample_without_replacement(const uint32_t n, const uint32_t random_seed) {
    linked_sparse_transposed transposed(*this);
    transposed.transposed_subsample_without_replacement(n,random_seed);
}

// =====================  biom_subsampled  ==========================

biom_subsampled::biom_subsampled(const biom_inmem &parent, const bool w_replacement, const uint32_t n, const uint32_t random_seed)
  : biom_inmem(true)
{
   sparse_data_subsampled tmp_obj(parent.get_resident_obj(), parent.get_sample_counts(), n);
   if ((tmp_obj.n_obs==0) || (tmp_obj.n_samples==0)) return; //already everything filtered out

   if (w_replacement) {
     tmp_obj.subsample_with_replacement(n,random_seed);
   } else {
     tmp_obj.subsample_without_replacement(n,random_seed);
   } 
   // Note: We could filter out the zero rows
   // But that's just an optimization and will not be worth it most of the time
   steal_nonzero(parent,tmp_obj);

   /* define a mapping between an ID and its corresponding offset */
   #pragma omp parallel for schedule(static)
   for(int i = 0; i < 3; i++) {
      if(i == 0)
         create_id_index(obs_ids, obs_id_index);
      else if(i == 1) {
         sample_ids.reserve(n_samples);
         const double *parent_sample_counts = parent.get_sample_counts();
         const std::vector<std::string> &parent_sample_ids = parent.get_sample_ids();

         for (uint32_t i=0; i<parent.n_samples; i++) {
            if (parent_sample_counts[i]>=n) {
               sample_ids.push_back(parent_sample_ids[i]);
            }
         }
         create_id_index(sample_ids, sample_id_index);
      } else if(i == 2)
         compute_sample_counts();
   }
}

void biom_subsampled::steal_nonzero(const biom_inmem &parent, sparse_data& subsampled_obj) {
   // initialize data structures
   n_samples = subsampled_obj.n_samples;
   resident_obj.n_samples = subsampled_obj.n_samples;
   resident_obj.n_obs = subsampled_obj.n_obs;
   resident_obj.malloc_resident();
   obs_ids.reserve(parent.n_obs);

   const std::vector<std::string> &parent_obs_ids = parent.get_obs_ids();

   // now do the copy
   n_obs = 0;
   for (uint32_t i=0; i<parent.n_obs; i++) {
     const uint32_t cnt = subsampled_obj.obs_counts_resident[i];
     double *data = subsampled_obj.obs_data_resident[i];
     uint32_t nz = 0; 
     for (uint32_t j=0; j<cnt; j++) if (data[j]>0.0) nz++;

     if (nz>0) {
        // steal non-zero data
        resident_obj.obs_indices_resident[n_obs] = subsampled_obj.steal_indices(i);
        resident_obj.obs_data_resident[n_obs] = subsampled_obj.steal_data(i); 
        resident_obj.obs_counts_resident[n_obs] = cnt;
        obs_ids.push_back(parent_obs_ids[i]);
        n_obs++;
     }
     // else just ignore
   }
   resident_obj.n_obs = n_obs;
   // Note: We could resize the buffersm but it is not worth it
}


