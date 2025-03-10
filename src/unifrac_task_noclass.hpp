/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2025, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

/*
 *
 * This file is used to create the necessary interfaces between
 *   unifrac_task.hpp and unifrac_task_impl.hpp
 * by means of
 *   generate_unifrac_task_noclass.py
 *
 */

#ifndef __UNIFRAC_TASK_NOCLASS
#define __UNIFRAC_TASK_NOCLASS 1

#include <stdint.h>

#ifndef SUCMP_NM
/* create a default */
#define SUCMP_NM su_cpu
#endif

namespace SUCMP_NM {

    // do we have access to a GPU?
    bool acc_found_gpu();

    // is the implementation async, and need the alt structures?
    bool acc_need_alt();

    // wait for the async compute to finish
    void acc_wait();

    // create the equivalent buffer in the device memory space, if partitioned
    // the content in undefined
    template<class T>
    void acc_create_buf(T *buf, uint64_t start, uint64_t end);

    // create the equivalent buffer in the device memory space, if partitioned
    // also copy the buffer over
    template<class T>
    void acc_copyin_buf(T *buf, uint64_t start, uint64_t end);

    // make a copy from host to device buffer, if partitioned
    template<class T>
    void acc_update_device(T *buf, uint64_t start, uint64_t end);

    // make a copy from device to host buffer, if partitioned
    // destroy the equivalent buffer in the device memory space, if partitioned
    template<class T>
    void acc_copyout_buf(T *buf, uint64_t start, uint64_t end);

    // destroy the equivalent buffer in the device memory space, if partitioned
    template<class T>
    void acc_destroy_buf(T *buf, uint64_t start, uint64_t end);

    // compute totals
    template<class TFloat>
    void compute_stripes_totals(
		TFloat * const __restrict__ dm_stripes_buf,
		const TFloat * const __restrict__ dm_stripes_total_buf,
		const uint64_t bufels);

    /* Unifrac tasks
     *
     * All functions utilize the same basic data structures:
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
     */

    // Compute UnnormalizedWeighted step
    template<class TFloat>
    void run_UnnormalizedWeightedTask(
                const uint64_t embs_stripe,
		const unsigned int filled_embs,
		const uint64_t start_idx, const uint64_t stop_idx,
		const uint64_t n_samples, const uint64_t n_samples_r,
		const TFloat * const __restrict__ lengths,
		const TFloat * const __restrict__ embedded_proportions,
		TFloat * const __restrict__ dm_stripes_buf,
		bool * const __restrict__ zcheck,
		TFloat * const __restrict__ sums);

    // Compute NormalizedWeighted step
    template<class TFloat>
    void run_NormalizedWeightedTask(
		const uint64_t embs_stripe,
		const unsigned int filled_embs,
		const uint64_t start_idx, const uint64_t stop_idx,
		const uint64_t n_samples, const uint64_t n_samples_r,
		const TFloat * const __restrict__ lengths,
		const TFloat * const __restrict__ embedded_proportions,
		TFloat * const __restrict__ dm_stripes_buf,
		TFloat * const __restrict__ dm_stripes_total_buf,
		bool * const __restrict__ zcheck,
		TFloat * const __restrict__ sums);

    // Compute Unweighted step
    template<class TFloat>
    void run_UnweightedTask(
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
		TFloat * const __restrict__ stripe_sums);

    // Compute UnnormalizedUnweighted step
    template<class TFloat>
    void run_UnnormalizedUnweightedTask(
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
		TFloat * const __restrict__ stripe_sums);

    // Compute Generalized step
    template<class TFloat>
    void run_GeneralizedTask(
                const uint64_t embs_stripe,
		const unsigned int filled_embs,
		const uint64_t start_idx, const uint64_t stop_idx,
		const uint64_t n_samples, const uint64_t n_samples_r,
		const TFloat * const __restrict__ lengths,
		const TFloat * const __restrict__ embedded_proportions,
		TFloat * const __restrict__ dm_stripes_buf,
		TFloat * const __restrict__ dm_stripes_total_buf,
		const TFloat g_unifrac_alpha);

    /* Unifrac vaw tasks
     *
     * All functions utilize the same basic data structures:
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
     */

    // Compute VawUnnormalizedWeighted step
    template<class TFloat>
    void run_VawUnnormalizedWeightedTask(
		const unsigned int filled_embs,
		const uint64_t start_idx, const uint64_t stop_idx,
		const uint64_t n_samples, const uint64_t n_samples_r,
		const TFloat * const __restrict__ lengths,
		const TFloat * const __restrict__ embedded_proportions,
		const TFloat * const __restrict__ embedded_counts,
		const TFloat * const __restrict__ sample_total_counts,
		TFloat * const __restrict__ dm_stripes_buf);

    // Compute VawNormalizedWeighted step
    template<class TFloat>
    void run_VawNormalizedWeightedTask(
		const unsigned int filled_embs,
		const uint64_t start_idx, const uint64_t stop_idx,
		const uint64_t n_samples, const uint64_t n_samples_r,
		const TFloat * const __restrict__ lengths,
		const TFloat * const __restrict__ embedded_proportions,
		const TFloat * const __restrict__ embedded_counts,
		const TFloat * const __restrict__ sample_total_counts,
		TFloat * const __restrict__ dm_stripes_buf,
		TFloat * const __restrict__ dm_stripes_total_buf);

    // Compute VawUnweighted step
    template<class TFloat>
    void run_VawUnweightedTask(
		const unsigned int filled_embs,
		const uint64_t start_idx, const uint64_t stop_idx,
		const uint64_t n_samples, const uint64_t n_samples_r,
		const TFloat * const __restrict__ lengths,
		const uint32_t * const __restrict__ embedded_proportions,
		const TFloat  * const __restrict__ embedded_counts,
		const TFloat  * const __restrict__ sample_total_counts,
		TFloat * const __restrict__ dm_stripes_buf,
		TFloat * const __restrict__ dm_stripes_total_buf);

    // Compute VawUnnormalizedUnweighted step
    template<class TFloat>
    void run_VawUnnormalizedUnweightedTask(
		const unsigned int filled_embs,
		const uint64_t start_idx, const uint64_t stop_idx,
		const uint64_t n_samples, const uint64_t n_samples_r,
		const TFloat * const __restrict__ lengths,
		const uint32_t * const __restrict__ embedded_proportions,
		const TFloat  * const __restrict__ embedded_counts,
		const TFloat  * const __restrict__ sample_total_counts,
		TFloat * const __restrict__ dm_stripes_buf);

    // Compute Generalized step
    template<class TFloat>
    void run_VawGeneralizedTask(
		const unsigned int filled_embs,
		const uint64_t start_idx, const uint64_t stop_idx, const uint64_t n_samples, const uint64_t n_samples_r,
		const TFloat * const __restrict__ lengths,
		const TFloat * const __restrict__ embedded_proportions,
		const TFloat * const __restrict__ embedded_counts,
		const TFloat * const __restrict__ sample_total_counts ,
		TFloat * const __restrict__ dm_stripes_buf,
		TFloat * const __restrict__ dm_stripes_total_buf,
		const TFloat g_unifrac_alpha);
}

#endif
