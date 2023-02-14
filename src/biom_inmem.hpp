/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2023, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */


#ifndef _UNIFRAC_BIOM_INMEM_H
#define _UNIFRAC_BIOM_INMEM_H

#include <vector>
#include <unordered_map>

#include "biom_interface.hpp"

namespace su {
    class sparse_data {
        public:
            /* default constructor */
            sparse_data(bool _clean_on_destruction);

            /* constructor from compress sparse data
             *
             * @param n_obs number of observations
             * @param index vector of index positions
             * @param indptr vector of indptr positions
             * @param data vector of observation counts
             */
            sparse_data(const uint32_t n_obs,
                        const uint32_t n_samples,
                        uint32_t* index,
                        uint32_t* indptr,
                        double* data);

            /* default destructor */
            virtual ~sparse_data();

            /* modified copy constructor */
            sparse_data(const sparse_data& other, bool _clean_on_destruction);

            /* prevent default copy constructor and operator from being generated */
            sparse_data(const sparse_data& other) = delete;
            sparse_data& operator= (const sparse_data&) = delete;

            /* Helper functions */
            void malloc_resident();
            void free_resident();
            template<class TData> TData *copy_resident_el(unsigned int cnt, const TData *other) const;

        public:  // keep it open for ease of access
            uint32_t n_obs;     // row dimension
            uint32_t n_samples; // column dimension
            bool clean_on_destruction;
            
            uint32_t **obs_indices_resident;
            double **obs_data_resident;
            unsigned int *obs_counts_resident;

    };

    class biom_inmem : public biom_interface {
        public:
            /* default constructor */
            biom_inmem(bool _clean_on_destruction);

            /* constructor from compress sparse data
             *
             * @param obs_ids vector of observation identifiers
             * @param samp_ids vector of sample identifiers
             * @param index vector of index positions
             * @param indptr vector of indptr positions
             * @param data vector of observation counts
             * @param n_obs number of observations
             * @param n_samples number of samples
             * @param nnz number of data points
             */
            biom_inmem(const char* const * obs_ids,
                       const char* const * samp_ids,
                       uint32_t* index,
                       uint32_t* indptr,
                       double* data,
                       const int n_obs,
                       const int n_samples,
                       const int nnz);

            /* default destructor */
            virtual ~biom_inmem();

            /* get a dense vector of observation data
             *
             * @param id The observation ID to fetch
             * @param out An allocated array of at least size n_samples. 
             *      Values of an index position [0, n_samples) which do not
             *      have data will be zero'd.
             */
            void get_obs_data(const std::string &id, double* out) const; 
            void get_obs_data(const std::string &id, float* out) const;

            /* get a dense vector of a range of observation data
             *
             * @param id The observation ID to fetch
             * @param start Initial index
             * @param end   First index past the end
             * @param normalize If set, divide by sample_counts
             * @param out An allocated array of at least size (end-start). First element will corrrectpoint to index start. 
             *      Values of an index position [0, (end-start)) which do not
             *      have data will be zero'd.
             */
            void get_obs_data_range(const std::string &id, unsigned int start, unsigned int end, bool normalize, double* out) const;
            void get_obs_data_range(const std::string &id, unsigned int start, unsigned int end, bool normalize, float* out) const;

            /* getters to local variables */
            virtual const std::vector<std::string> &get_sample_ids() const;
            virtual const std::vector<std::string> &get_obs_ids() const;
            virtual const double *get_sample_counts() const;
        protected:
            sparse_data resident_obj;

            // distilled from the above resident values
            double *sample_counts;

            /* At construction, lookups mapping IDs -> index position within an
             * axis are defined
             */
            std::unordered_map<std::string, uint32_t> obs_id_index;
            std::unordered_map<std::string, uint32_t> sample_id_index;
 
            // cache the IDs contained within the table
            std::vector<std::string> sample_ids;
            std::vector<std::string> obs_ids;

        protected:
            /* Only allow copy constructr by children classes */
            biom_inmem(const biom_inmem& other, bool _clean_on_destruction);

            void compute_sample_counts();


            /* create an index mapping an ID to its corresponding index 
             * position.
             *
             * @param ids A vector of IDs to index
             * @param map A hash table to populate
             */
            void create_id_index(const std::vector<std::string> &ids, 
                                 std::unordered_map<std::string, uint32_t> &map);


            // templatized version
            template<class TFloat> void get_obs_data_TT(const uint32_t idx, TFloat* out) const;
            template<class TFloat> void get_obs_data_TT(const std::string &id, TFloat* out) const;
            template<class TFloat> void get_obs_data_range_TT(const uint32_t idx, unsigned int start, unsigned int end, bool normalize, TFloat* out) const;
            template<class TFloat> void get_obs_data_range_TT(const std::string &id, unsigned int start, unsigned int end, bool normalize, TFloat* out) const;
        public:
            /* prevent default copy constructor and operator from being generated */
            biom_inmem(const biom_inmem& other) = delete;
            biom_inmem& operator= (const biom_inmem&) = delete;
    };
}

#endif /* _UNIFRAC_BIOM_H */
