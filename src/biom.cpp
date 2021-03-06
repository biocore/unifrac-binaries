/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2016-2021, UniFrac development team.
 * All rights reserved.
 *
 * See LICENSE file for more details
 */

#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include "biom.hpp"

using namespace H5;
using namespace su;

/* datasets defined by the BIOM 2.x spec */ 
const std::string OBS_INDPTR = std::string("/observation/matrix/indptr");
const std::string OBS_INDICES = std::string("/observation/matrix/indices");
const std::string OBS_DATA = std::string("/observation/matrix/data");
const std::string OBS_IDS = std::string("/observation/ids");

const std::string SAMPLE_INDPTR = std::string("/sample/matrix/indptr");
const std::string SAMPLE_INDICES = std::string("/sample/matrix/indices");
const std::string SAMPLE_DATA = std::string("/sample/matrix/data");
const std::string SAMPLE_IDS = std::string("/sample/ids");

biom::biom(std::string filename) : has_hdf5_backing(true) {
    file = H5File(filename.c_str(), H5F_ACC_RDONLY);

    /* establish the datasets */
    obs_indices = file.openDataSet(OBS_INDICES.c_str());
    obs_data = file.openDataSet(OBS_DATA.c_str());
    sample_indices = file.openDataSet(SAMPLE_INDICES.c_str());
    sample_data = file.openDataSet(SAMPLE_DATA.c_str());
    
    /* cache IDs and indptr */
    sample_ids = std::vector<std::string>();
    obs_ids = std::vector<std::string>();
    sample_indptr = std::vector<uint32_t>();
    obs_indptr = std::vector<uint32_t>();

    load_ids(OBS_IDS.c_str(), obs_ids);
    load_ids(SAMPLE_IDS.c_str(), sample_ids);
    load_indptr(OBS_INDPTR.c_str(), obs_indptr);    
    load_indptr(SAMPLE_INDPTR.c_str(), sample_indptr);    

    /* cache shape and nnz info */
    n_samples = sample_ids.size();
    n_obs = obs_ids.size();
    set_nnz();

    /* define a mapping between an ID and its corresponding offset */
    obs_id_index = std::unordered_map<std::string, uint32_t>();
    sample_id_index = std::unordered_map<std::string, uint32_t>();

    #pragma omp parallel for schedule(static)
    for(int i = 0; i < 3; i++) {
        if(i == 0)
            create_id_index(obs_ids, obs_id_index);
        else if(i == 1)
            create_id_index(sample_ids, sample_id_index);
        else if(i == 2)
            malloc_resident(n_obs);
    }

    uint32_t *current_indices = NULL;
    double *current_data = NULL;
    for(unsigned int i = 0; i < obs_ids.size(); i++)  {
        std::string id_ = obs_ids[i];
        unsigned int n = get_obs_data_direct(id_, current_indices, current_data);
        obs_counts_resident[i] = n;
        obs_indices_resident[i] = current_indices;
        obs_data_resident[i] = current_data;
    }
    sample_counts = get_sample_counts();
}

biom::~biom() {
    if(has_hdf5_backing) {
        if(obs_indices_resident != NULL && obs_data_resident != NULL) {
            for(unsigned int i = 0; i < n_obs; i++) {
                if(obs_indices_resident[i] != NULL)
                    free(obs_indices_resident[i]);
                if(obs_data_resident[i] != NULL)
                    free(obs_data_resident[i]);
            }
        }
        
        if(obs_indices_resident != NULL)
            free(obs_indices_resident);
        if(obs_data_resident != NULL)
            free(obs_data_resident);
        if(obs_counts_resident != NULL)
            free(obs_counts_resident);
    } 
    // else, it is the responsibility of the entity constructing this object
    // to clean itself up
}

void biom::malloc_resident(uint32_t n_obs) { 
    /* load obs sparse data */
    obs_indices_resident = (uint32_t**)malloc(sizeof(uint32_t**) * n_obs);
    if(obs_indices_resident == NULL) {
        fprintf(stderr, "Failed to allocate %zd bytes; [%s]:%d\n", 
                sizeof(uint32_t**) * n_obs, __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    obs_data_resident = (double**)malloc(sizeof(double**) * n_obs);
    if(obs_data_resident == NULL) {
        fprintf(stderr, "Failed to allocate %zd bytes; [%s]:%d\n", 
                sizeof(double**) * n_obs, __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    obs_counts_resident = (unsigned int*)malloc(sizeof(unsigned int) * n_obs);
    if(obs_counts_resident == NULL) {
        fprintf(stderr, "Failed to allocate %zd bytes; [%s]:%d\n", 
                sizeof(unsigned int) * n_obs, __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
}

biom::biom() : has_hdf5_backing(false) { 
    n_obs = 0;
    malloc_resident(0);
}

// not using const on indices/indptr/data as the pointers are being borrowed
biom::biom(char** obs_ids_in,
           char** samp_ids_in,
           uint32_t* indices,
           uint32_t* indptr,
           double* data,
           const int n_obs,
           const int n_samples,
           const int nnz) : has_hdf5_backing(false) {

    this->nnz = nnz;
    this->n_samples = n_samples;
    this->n_obs = n_obs;

    sample_ids = std::vector<std::string>();
    sample_ids.resize(n_samples);
    obs_ids = std::vector<std::string>();
    obs_ids.resize(n_obs);
   
    #pragma omp parallel for schedule(static)
    for(int x = 0; x < 2; x++) {
        if(x == 0) {
            for(int i = 0; i < n_obs; i++) {
                obs_ids[i] = std::string(obs_ids_in[i]);
            }
        } else {
            for(int i = 0; i < n_samples; i++) {
                sample_ids[i] = std::string(samp_ids_in[i]);
            }
        }
    }

    /* define a mapping between an ID and its corresponding offset */
    obs_id_index = std::unordered_map<std::string, uint32_t>();
    sample_id_index = std::unordered_map<std::string, uint32_t>();

    #pragma omp parallel for schedule(static)
    for(int i = 0; i < 3; i++) {
        if(i == 0)
            create_id_index(obs_ids, obs_id_index);
        else if(i == 1)
            create_id_index(sample_ids, sample_id_index);
        else if(i == 2)
            malloc_resident(n_obs);
    }

    #pragma omp parallel for schedule(static)
    for(unsigned int i = 0; i < n_obs; i++)  {
        int32_t start = indptr[i];
        int32_t end = indptr[i + 1];
        unsigned int count = end - start;

        uint32_t* index_ptr = (indices + start);
        double* data_ptr = (data + start);
        
        obs_indices_resident[i] = index_ptr;
        obs_data_resident[i] = data_ptr;
        obs_counts_resident[i] = count;
    }
    sample_counts = get_sample_counts();
}

void biom::set_nnz() {
    if(!has_hdf5_backing) {
        fprintf(stderr, "Lacks HDF5 backing; [%s]:%d\n", 
                __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }

    // should these be cached?
    DataType dtype = obs_data.getDataType();
    DataSpace dataspace = obs_data.getSpace();

    hsize_t dims[1];
    dataspace.getSimpleExtentDims(dims, NULL);
    nnz = dims[0];
}

void biom::load_ids(const char *path, std::vector<std::string> &ids) {
    if(!has_hdf5_backing) {
        fprintf(stderr, "Lacks HDF5 backing; [%s]:%d\n", 
                __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }

    DataSet ds_ids = file.openDataSet(path);
    DataType dtype = ds_ids.getDataType();
    DataSpace dataspace = ds_ids.getSpace();

    hsize_t dims[1];
    dataspace.getSimpleExtentDims(dims, NULL);

    /* the IDs are a dataset of variable length strings */
    char **dataout = (char**)malloc(sizeof(char*) * dims[0]);
    if(dataout == NULL) {
        fprintf(stderr, "Failed to allocate %zd bytes; [%s]:%d\n", 
                sizeof(char*) * dims[0], __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    ds_ids.read((void*)dataout, dtype);

    ids.reserve(dims[0]);
    for(unsigned int i = 0; i < dims[0]; i++) {
        ids.push_back(dataout[i]);
    }
    
    for(unsigned int i = 0; i < dims[0]; i++) 
        free(dataout[i]);
    free(dataout);
}

void biom::load_indptr(const char *path, std::vector<uint32_t> &indptr) {
    if(!has_hdf5_backing) {
        fprintf(stderr, "Lacks HDF5 backing; [%s]:%d\n", 
                __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }

    DataSet ds = file.openDataSet(path);
    DataType dtype = ds.getDataType();
    DataSpace dataspace = ds.getSpace();

    hsize_t dims[1];
    dataspace.getSimpleExtentDims(dims, NULL);
    
    uint32_t *dataout = (uint32_t*)malloc(sizeof(uint32_t) * dims[0]);
    if(dataout == NULL) {
        fprintf(stderr, "Failed to allocate %zd bytes; [%s]:%d\n", 
                sizeof(uint32_t) * dims[0], __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    ds.read((void*)dataout, dtype);

    indptr.reserve(dims[0]);
    for(unsigned int i = 0; i < dims[0]; i++)
        indptr.push_back(dataout[i]);
    free(dataout);
}

void biom::create_id_index(const std::vector<std::string> &ids, 
                           std::unordered_map<std::string, uint32_t> &map) {
    uint32_t count = 0;
    map.reserve(ids.size());
    for(auto i = ids.begin(); i != ids.end(); i++, count++) {
        map[*i] = count;
    }
}

unsigned int biom::get_obs_data_direct(const std::string &id, uint32_t *& current_indices_out, double *& current_data_out) {
    if(!has_hdf5_backing) {
        fprintf(stderr, "Lacks HDF5 backing; [%s]:%d\n", 
                __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }

    uint32_t idx = obs_id_index.at(id);
    uint32_t start = obs_indptr[idx];
    uint32_t end = obs_indptr[idx + 1];
    
    hsize_t count[1] = {end - start};
    hsize_t offset[1] = {start};

    DataType indices_dtype = obs_indices.getDataType();
    DataType data_dtype = obs_data.getDataType();

    DataSpace indices_dataspace = obs_indices.getSpace();
    DataSpace data_dataspace = obs_data.getSpace();
    
    DataSpace indices_memspace(1, count, NULL);
    DataSpace data_memspace(1, count, NULL);

    indices_dataspace.selectHyperslab(H5S_SELECT_SET, count, offset); 
    data_dataspace.selectHyperslab(H5S_SELECT_SET, count, offset); 

    current_indices_out = (uint32_t*)malloc(sizeof(uint32_t) * count[0]);
    if(current_indices_out == NULL) {
        fprintf(stderr, "Failed to allocate %zd bytes; [%s]:%d\n", 
                sizeof(uint32_t) * count[0], __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    current_data_out = (double*)malloc(sizeof(double) * count[0]);
    if(current_data_out == NULL) {
        fprintf(stderr, "Failed to allocate %zd bytes; [%s]:%d\n", 
                sizeof(double) * count[0], __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }

    obs_indices.read((void*)current_indices_out, indices_dtype, indices_memspace, indices_dataspace);
    obs_data.read((void*)current_data_out, data_dtype, data_memspace, data_dataspace);
    
    return count[0];
}

template<class TFloat>
void biom::get_obs_data_TT(const std::string &id, TFloat* out) const {
    uint32_t idx = obs_id_index.at(id);
    unsigned int count = obs_counts_resident[idx];
    const uint32_t * const indices = obs_indices_resident[idx];
    const double * const data = obs_data_resident[idx];

    // reset our output buffer
    for(unsigned int i = 0; i < n_samples; i++)
        out[i] = 0.0;
    
    for(unsigned int i = 0; i < count; i++) {
        out[indices[i]] = data[i];
    }
}

void biom::get_obs_data(const std::string &id, double* out) const {
  biom::get_obs_data_TT(id,out);
}

void biom::get_obs_data(const std::string &id, float* out) const {
  biom::get_obs_data_TT(id,out);
}


// note: out is supposed to be fully filled, i.e. out[start:end]
template<class TFloat>
void biom::get_obs_data_range_TT(const std::string &id, unsigned int start, unsigned int end, bool normalize, TFloat* out) const {
    uint32_t idx = obs_id_index.at(id);
    unsigned int count = obs_counts_resident[idx];
    const uint32_t * const indices = obs_indices_resident[idx];
    const double * const data = obs_data_resident[idx];

    // reset our output buffer
    for(unsigned int i = start; i < end; i++)
        out[i-start] = 0.0;

    if (normalize) {
      for(unsigned int i = 0; i < count; i++) {
        const int32_t j = indices[i];
        if ((j>=start)&&(j<end)) { 
          out[j-start] = data[i]/sample_counts[j];
        }
      }
    } else {
      for(unsigned int i = 0; i < count; i++) {
        const uint32_t j = indices[i];
        if ((j>=start)&&(j<end)) {
          out[j-start] = data[i];
        }
      }
    }
}

void biom::get_obs_data_range(const std::string &id, unsigned int start, unsigned int end, bool normalize, double* out) const {
  biom::get_obs_data_range_TT(id,start,end,normalize,out);
}

void biom::get_obs_data_range(const std::string &id, unsigned int start, unsigned int end, bool normalize, float* out) const {
  biom::get_obs_data_range_TT(id,start,end,normalize,out);
}

unsigned int biom::get_sample_data_direct(const std::string &id, uint32_t *& current_indices_out, double *& current_data_out) {
    if(!has_hdf5_backing) {
        fprintf(stderr, "Lacks HDF5 backing; [%s]:%d\n", 
                __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }

    uint32_t idx = sample_id_index.at(id);
    uint32_t start = sample_indptr[idx];
    uint32_t end = sample_indptr[idx + 1];

    hsize_t count[1] = {end - start};
    hsize_t offset[1] = {start};

    DataType indices_dtype = sample_indices.getDataType();
    DataType data_dtype = sample_data.getDataType();

    DataSpace indices_dataspace = sample_indices.getSpace();
    DataSpace data_dataspace = sample_data.getSpace();
    
    DataSpace indices_memspace(1, count, NULL);
    DataSpace data_memspace(1, count, NULL);

    indices_dataspace.selectHyperslab(H5S_SELECT_SET, count, offset); 
    data_dataspace.selectHyperslab(H5S_SELECT_SET, count, offset); 

    current_indices_out = (uint32_t*)malloc(sizeof(uint32_t) * count[0]);
    if(current_indices_out == NULL) {
        fprintf(stderr, "Failed to allocate %zd bytes; [%s]:%d\n", 
                sizeof(uint32_t) * count[0], __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    current_data_out = (double*)malloc(sizeof(double) * count[0]);
    if(current_data_out == NULL) {
        fprintf(stderr, "Failed to allocate %zd bytes; [%s]:%d\n", 
                sizeof(double) * count[0], __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }

    sample_indices.read((void*)current_indices_out, indices_dtype, indices_memspace, indices_dataspace);
    sample_data.read((void*)current_data_out, data_dtype, data_memspace, data_dataspace);

    return count[0];
}

double* biom::get_sample_counts() {
    double *sample_counts = (double*)calloc(sizeof(double), n_samples);

    for(unsigned int i = 0; i < n_obs; i++) {
        unsigned int count = obs_counts_resident[i];
        uint32_t *indices = obs_indices_resident[i];
        double *data = obs_data_resident[i];
        for(unsigned int j = 0; j < count; j++) {
            uint32_t index = indices[j];
            double datum = data[j];
            sample_counts[index] += datum;
        }
    }
    return sample_counts;
}
