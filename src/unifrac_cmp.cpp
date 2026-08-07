/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2025, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

#include "tree.hpp"
#include "biom_interface.hpp"
#include <unordered_map>
#include <cstdlib>
#include <thread>
#include <algorithm>

#include "unifrac_cmp.hpp"
#include "unifrac_internal.hpp"
#include "unifrac_task.hpp"

using namespace SUCMP_NM;

// do we have access to a GPU?
// Just a wrapper, to expose it to the rest of the code
bool SUCMP_NM::found_gpu() {
  return acc_found_gpu();
}

template<class TaskT, class TFloat>
inline void unifracTT(const su::biom_interface &table,
                      const su::BPTree &tree,
                      const bool want_total,
                      std::vector<double*> &dm_stripes,
                      std::vector<double*> &dm_stripes_total,
                      const su::task_parameters* task_p) {
    // no processor affinity whenusing openacc or openmp

    if(table.n_samples != task_p->n_samples) {
        fprintf(stderr, "Task and table n_samples not equal\n");
        exit(EXIT_FAILURE);
    }

    const bool normalize_sample_counts = task_p->normalize_sample_counts;
    su::PropStackMulti<TFloat> propstack_multi(table.n_samples);

    const unsigned int max_emb =  TaskT::RECOMMENDED_MAX_EMBS;

    su::initialize_stripes(std::ref(dm_stripes), std::ref(dm_stripes_total), want_total, task_p);

    TaskT taskObj(std::ref(dm_stripes), std::ref(dm_stripes_total),max_emb,task_p);

    TFloat * const lengths = taskObj.lengths;

        /*
         * The values in the example vectors correspond to index positions of an
         * element in the resulting distance matrix. So, in the example below,
         * the following can be interpreted:
         *
         * [0 1 2]
         * [1 2 3]
         *
         * As comparing the sample for row 0 against the sample for col 1, the
         * sample for row 1 against the sample for col 2, the sample for row 2
         * against the sample for col 3.
         *
         * In other words, we're computing stripes of a distance matrix. In the
         * following example, we're computing over 6 samples requiring 3
         * stripes.
         *
         * A; stripe == 0
         * [0 1 2 3 4 5]
         * [1 2 3 4 5 0]
         *
         * B; stripe == 1
         * [0 1 2 3 4 5]
         * [2 3 4 5 0 1]
         *
         * C; stripe == 2
         * [0 1 2 3 4 5]
         * [3 4 5 0 1 2]
         *
         * The stripes end up computing the following positions in the distance
         * matrix.
         *
         * x A B C x x
         * x x A B C x
         * x x x A B C
         * C x x x A B
         * B C x x x A
         * A B C x x x
         *
         * However, we store those stripes as vectors, ie
         * [ A A A A A A ]
         *
         * We end up performing N / 2 redundant calculations on the last stripe
         * (see C) but that is small over large N.
         */

    unsigned int k = 0; // index in tree
    const unsigned int max_k = (tree.nparens>1) ? ((tree.nparens / 2) - 1) : 0;

    const unsigned int num_prop_chunks = propstack_multi.get_num_stacks();
    while (k<max_k) {
          const unsigned int k_start = k;
          unsigned int filled_emb = 0;

          // chunk the progress to maximize cache reuse
#pragma omp parallel for 
          for (unsigned int ck=0; ck<num_prop_chunks; ck++) {
            su::PropStack<TFloat> &propstack = propstack_multi.get_prop_stack(ck);
            const unsigned int tstart = propstack_multi.get_start(ck);
            const unsigned int tend = propstack_multi.get_end(ck);
            unsigned int my_filled_emb = 0;
            unsigned int my_k=k_start;

            while ((my_filled_emb<max_emb) && (my_k<max_k)) {
              const uint32_t node = tree.postorderselect(my_k);
              my_k++;

              TFloat *node_proportions = propstack.pop(node);
              su::set_proportions_range(node_proportions, tree, node, table, tstart, tend, propstack, normalize_sample_counts);

              if(task_p->bypass_tips && tree.isleaf(node))
                  continue;

              if (ck==0) { // they all do the same thing, so enough for the first to update the global state
                lengths[filled_emb] = tree.lengths[node];
                filled_emb++;
              }
              taskObj.embed_proportions_range(node_proportions, tstart, tend, my_filled_emb);
              my_filled_emb++;
            }
            if (ck==0) { // they all do the same thing, so enough for the first to update the global state
              k=my_k;
            }
          }

          taskObj.sync_embedded_proportions(filled_emb);
          taskObj.sync_lengths(filled_emb);
          taskObj._run(filled_emb);
          filled_emb=0;
    }

    taskObj.wait_completion();

    if(want_total) {
        taskObj.compute_totals();
    }

}

void SUCMP_NM::unifrac(const su::biom_interface &table,
                       const su::BPTree &tree,
                       su::Method unifrac_method,
                       std::vector<double*> &dm_stripes,
                       std::vector<double*> &dm_stripes_total,
                        const su::task_parameters* task_p) {
    switch(unifrac_method) {
        case su::unweighted:
            unifracTT<SUCMP_NM::UnifracUnweightedTask<double>,double>(           table, tree, true,  dm_stripes,dm_stripes_total,task_p);
            break;
        case su::unweighted_unnormalized:
            unifracTT<SUCMP_NM::UnifracUnnormalizedUnweightedTask<double>,double>(table,tree, false, dm_stripes,dm_stripes_total,task_p);
            break;
        case su::weighted_normalized:
            unifracTT<SUCMP_NM::UnifracNormalizedWeightedTask<double>,double>(   table, tree, true,  dm_stripes,dm_stripes_total,task_p);
            break;
        case su::weighted_unnormalized:
            unifracTT<SUCMP_NM::UnifracUnnormalizedWeightedTask<double>,double>( table, tree, false, dm_stripes,dm_stripes_total,task_p);
            break;
        case su::generalized:
            unifracTT<SUCMP_NM::UnifracGeneralizedTask<double>,double>(          table, tree, true,  dm_stripes,dm_stripes_total,task_p);
            break;
        case su::unweighted_fp32:
            unifracTT<SUCMP_NM::UnifracUnweightedTask<float >,float>(            table, tree, true,  dm_stripes,dm_stripes_total,task_p);
            break;
        case su::unweighted_unnormalized_fp32:
            unifracTT<SUCMP_NM::UnifracUnnormalizedUnweightedTask<float >,float>(table, tree, false, dm_stripes,dm_stripes_total,task_p);
            break;
        case su::weighted_normalized_fp32:
            unifracTT<SUCMP_NM::UnifracNormalizedWeightedTask<float >,float>(    table, tree, true,  dm_stripes,dm_stripes_total,task_p);
            break;
        case su::weighted_unnormalized_fp32:
            unifracTT<SUCMP_NM::UnifracUnnormalizedWeightedTask<float >,float>(  table, tree, false, dm_stripes,dm_stripes_total,task_p);
            break;
        case su::generalized_fp32:
            unifracTT<SUCMP_NM::UnifracGeneralizedTask<float >,float>(           table, tree, true,  dm_stripes,dm_stripes_total,task_p);
            break;
        default:
            fprintf(stderr, "Unknown unifrac task\n");
            exit(1);
            break;
    }
}


template<class TaskT, class TFloat>
inline void unifrac_vawTT(const su::biom_interface &table,
                          const su::BPTree &tree,
                          const bool want_total,
                          std::vector<double*> &dm_stripes,
                          std::vector<double*> &dm_stripes_total,
                          const su::task_parameters* task_p) {
    // no processor affinity whenusing openacc or openmp

    if(table.n_samples != task_p->n_samples) {
        fprintf(stderr, "Task and table n_samples not equal\n");
        exit(EXIT_FAILURE);
    }

    const bool normalize_sample_counts = task_p->normalize_sample_counts;
    su::PropStackMulti<TFloat> propstack_multi(table.n_samples);
    su::PropStackMulti<TFloat> countstack_multi(table.n_samples);

    const unsigned int max_emb = TaskT::RECOMMENDED_MAX_EMBS;

    su::initialize_stripes(std::ref(dm_stripes), std::ref(dm_stripes_total), want_total, task_p);

    TaskT taskObj(std::ref(dm_stripes), std::ref(dm_stripes_total), table.get_sample_counts(), max_emb, task_p);

    TFloat * const lengths = taskObj.lengths;

    unsigned int k = 0; // index in tree
    const unsigned int max_k = (tree.nparens>1) ? ((tree.nparens / 2) - 1) : 0;

    const unsigned int num_prop_chunks = propstack_multi.get_num_stacks();
    while (k<max_k) {
          const unsigned int k_start = k;
          unsigned int filled_emb = 0;

          // chunk the progress to maximize cache reuse
#pragma omp parallel for 
          for (unsigned int ck=0; ck<num_prop_chunks; ck++) {
            su::PropStack<TFloat> &propstack = propstack_multi.get_prop_stack(ck);
            su::PropStack<TFloat> &countstack = countstack_multi.get_prop_stack(ck);
            const unsigned int tstart = propstack_multi.get_start(ck);
            const unsigned int tend = propstack_multi.get_end(ck);
            unsigned int my_filled_emb = 0;
            unsigned int my_k=k_start;

            while ((my_filled_emb<max_emb) && (my_k<max_k)) {
              const uint32_t node = tree.postorderselect(my_k);
              my_k++;

              TFloat *node_proportions = propstack.pop(node);
              TFloat *node_counts = countstack.pop(node);

              su::set_proportions_range(node_proportions, tree, node, table, tstart, tend, propstack, normalize_sample_counts);
              su::set_proportions_range(node_counts, tree, node, table, tstart, tend, countstack, false);

              if(task_p->bypass_tips && tree.isleaf(node))
                  continue;

              if (ck==0) { // they all do the same thing, so enough for the first to update the global state
                lengths[filled_emb] = tree.lengths[node];
                filled_emb++;
              }
              taskObj.embed_range(node_proportions, node_counts, tstart, tend, my_filled_emb);
              my_filled_emb++;
            }
            if (ck==0) { // they all do the same thing, so enough for the first to update the global state
              k=my_k;
            }
          }

	  taskObj.sync_lengths(filled_emb);
          taskObj.sync_embedded(filled_emb);
          taskObj._run(filled_emb);
          filled_emb = 0;
    }

    taskObj.wait_completion();

    if(want_total) {
        taskObj.compute_totals();
    }


}

void SUCMP_NM::unifrac_vaw(const su::biom_interface &table,
                           const su::BPTree &tree,
                           su::Method unifrac_method,
                           std::vector<double*> &dm_stripes,
                           std::vector<double*> &dm_stripes_total,
                           const su::task_parameters* task_p) {
    switch(unifrac_method) {
        case su::unweighted:
            unifrac_vawTT<SUCMP_NM::UnifracVawUnweightedTask<double>,double>(           table, tree, true,  dm_stripes,dm_stripes_total,task_p);
            break;
        case su::unweighted_unnormalized:
            unifrac_vawTT<SUCMP_NM::UnifracVawUnnormalizedUnweightedTask<double>,double>(table,tree, false, dm_stripes,dm_stripes_total,task_p);
            break;
        case su::weighted_normalized:
            unifrac_vawTT<SUCMP_NM::UnifracVawNormalizedWeightedTask<double>,double>(   table, tree, true,  dm_stripes,dm_stripes_total,task_p);
            break;
        case su::weighted_unnormalized:
            unifrac_vawTT<SUCMP_NM::UnifracVawUnnormalizedWeightedTask<double>,double>( table, tree, false, dm_stripes,dm_stripes_total,task_p);
            break;
        case su::generalized:
            unifrac_vawTT<SUCMP_NM::UnifracVawGeneralizedTask<double>,double>(          table, tree, true,  dm_stripes,dm_stripes_total,task_p);
            break;
        case su::unweighted_fp32:
            unifrac_vawTT<SUCMP_NM::UnifracVawUnweightedTask<float >,float >(           table, tree, true,  dm_stripes,dm_stripes_total,task_p);
            break;
        case su::unweighted_unnormalized_fp32:
            unifrac_vawTT<SUCMP_NM::UnifracVawUnnormalizedUnweightedTask<float >,float >(table,tree, false, dm_stripes,dm_stripes_total,task_p);
            break;
        case su::weighted_normalized_fp32:
            unifrac_vawTT<SUCMP_NM::UnifracVawNormalizedWeightedTask<float >,float >(   table, tree, true,  dm_stripes,dm_stripes_total,task_p);
            break;
        case su::weighted_unnormalized_fp32:
            unifrac_vawTT<SUCMP_NM::UnifracVawUnnormalizedWeightedTask<float >,float >( table, tree, false, dm_stripes,dm_stripes_total,task_p);
            break;
        case su::generalized_fp32:
            unifrac_vawTT<SUCMP_NM::UnifracVawGeneralizedTask<float >,float >(          table, tree, true,  dm_stripes,dm_stripes_total,task_p);
            break;
        default:
            fprintf(stderr, "Unknown unifrac task\n");
            exit(1);
            break;
    }
}

