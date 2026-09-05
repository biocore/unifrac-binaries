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

biom::biom(std::string filename) 
  : biom_inmem(true)
  , has_hdf5_backing(true) {
    file = H5File(filename.c_str(), H5F_ACC_RDONLY);

    /* establish the datasets */
    obs_indices = file.openDataSet(OBS_INDICES);
    obs_data = file.openDataSet(OBS_DATA);
    sample_indices = file.openDataSet(SAMPLE_INDICES);
    sample_data = file.openDataSet(SAMPLE_DATA);
    
    /* cache IDs and indptr */
    load_ids(OBS_IDS, obs_ids);
    load_ids(SAMPLE_IDS, sample_ids);
    load_indptr(OBS_INDPTR, obs_indptr);    
    load_indptr(SAMPLE_INDPTR, sample_indptr);    

    /* cache shape and nnz info */
    n_samples = sample_ids.size();
    n_obs = obs_ids.size();
    set_nnz();

    /* define a mapping between an ID and its corresponding offset */
    #pragma omp parallel for schedule(static)
    for(int i = 0; i < 3; i++) {
        if(i == 0)
            create_id_index(obs_ids, obs_id_index);
        else if(i == 1)
            create_id_index(sample_ids, sample_id_index);
        else if(i == 2) {
            resident_obj.n_obs = n_obs;
            resident_obj.n_samples = n_samples;
            resident_obj.malloc_resident();
        }
    }

    uint32_t *current_indices = NULL;
    double *current_data = NULL;
    for(unsigned int i = 0; i < obs_ids.size(); i++)  {
        std::string id_ = obs_ids[i];
        unsigned int n = get_obs_data_direct(id_, current_indices, current_data);
        resident_obj.obs_counts_resident[i] = n;
        resident_obj.obs_indices_resident[i] = current_indices;
        resident_obj.obs_data_resident[i] = current_data;
    }
    compute_sample_counts();
}

biom::~biom() {}

biom::biom() 
  : biom_inmem(false)
  , has_hdf5_backing(false)
  , nnz(0) { 
    resident_obj.malloc_resident();
}

// not using const on indices/indptr/data as the pointers are being borrowed
// Note: Pass through, for historical reasons
biom::biom(const char* const * obs_ids_in,
          const  char* const * samp_ids_in,
           uint32_t* indices,
           uint32_t* indptr,
           double* data,
           const int n_obs,
           const int n_samples,
           const int _nnz) 
  : biom_inmem(obs_ids_in,samp_ids_in,indices,indptr,data,n_obs,n_samples)
  , has_hdf5_backing(false)
  , nnz(_nnz)
{}

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

    DataType fdtype = ds_ids.getDataType();
    hid_t dtype_id = H5Tget_native_type(fdtype.getId(), H5T_DIR_ASCEND);
    DataType dtype(dtype_id);

    DataSpace dataspace = ds_ids.getSpace();

    hsize_t dims[1];
    dataspace.getSimpleExtentDims(dims, NULL);

    /* the IDs are a dataset of variable length strings */
    char **dataout = (char**)malloc(sizeof(char*) * dims[0]);
    if(dataout == NULL) {
        fprintf(stderr, "Failed to allocate %zd bytes; [%s]:%d\n", 
                long(sizeof(char*) * dims[0]), __FILE__, __LINE__);
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

    DataType fdtype = ds.getDataType();
    hid_t dtype_id = H5Tget_native_type(fdtype.getId(), H5T_DIR_ASCEND);
    DataType dtype(dtype_id);

    DataSpace dataspace = ds.getSpace();

    hsize_t dims[1];
    dataspace.getSimpleExtentDims(dims, NULL);
    
    uint32_t *dataout = (uint32_t*)malloc(sizeof(uint32_t) * dims[0]);
    if(dataout == NULL) {
        fprintf(stderr, "Failed to allocate %zd bytes; [%s]:%d\n", 
                long(sizeof(uint32_t) * dims[0]), __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    ds.read((void*)dataout, dtype);

    indptr.reserve(dims[0]);
    for(unsigned int i = 0; i < dims[0]; i++)
        indptr.push_back(dataout[i]);
    free(dataout);
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

    DataType indices_fdtype = obs_indices.getDataType();
    hid_t indices_dtype_id = H5Tget_native_type(indices_fdtype.getId(), H5T_DIR_ASCEND);
    DataType indices_dtype(indices_dtype_id);

    DataType data_fdtype = obs_data.getDataType();
    hid_t data_dtype_id = H5Tget_native_type(data_fdtype.getId(), H5T_DIR_ASCEND);
    DataType data_dtype(data_dtype_id);

    DataSpace indices_dataspace = obs_indices.getSpace();
    DataSpace data_dataspace = obs_data.getSpace();
    
    DataSpace indices_memspace(1, count, NULL);
    DataSpace data_memspace(1, count, NULL);

    indices_dataspace.selectHyperslab(H5S_SELECT_SET, count, offset); 
    data_dataspace.selectHyperslab(H5S_SELECT_SET, count, offset); 

    current_indices_out = (uint32_t*)malloc(sizeof(uint32_t) * count[0]);
    if(current_indices_out == NULL) {
        fprintf(stderr, "Failed to allocate %zd bytes; [%s]:%d\n", 
                long(sizeof(uint32_t) * count[0]), __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    current_data_out = (double*)malloc(sizeof(double) * count[0]);
    if(current_data_out == NULL) {
        fprintf(stderr, "Failed to allocate %zd bytes; [%s]:%d\n", 
                long(sizeof(double) * count[0]), __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }

    obs_indices.read((void*)current_indices_out, indices_dtype, indices_memspace, indices_dataspace);
    obs_data.read((void*)current_data_out, data_dtype, data_memspace, data_dataspace);
    
    return count[0];
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

    DataType indices_fdtype = sample_indices.getDataType();
    hid_t indices_dtype_id = H5Tget_native_type(indices_fdtype.getId(), H5T_DIR_ASCEND);
    DataType indices_dtype(indices_dtype_id);

    DataType data_fdtype = sample_data.getDataType();
    hid_t data_dtype_id = H5Tget_native_type(data_fdtype.getId(), H5T_DIR_ASCEND);
    DataType data_dtype(data_dtype_id);

    DataSpace indices_dataspace = sample_indices.getSpace();
    DataSpace data_dataspace = sample_data.getSpace();
    
    DataSpace indices_memspace(1, count, NULL);
    DataSpace data_memspace(1, count, NULL);

    indices_dataspace.selectHyperslab(H5S_SELECT_SET, count, offset); 
    data_dataspace.selectHyperslab(H5S_SELECT_SET, count, offset); 

    current_indices_out = (uint32_t*)malloc(sizeof(uint32_t) * count[0]);
    if(current_indices_out == NULL) {
        fprintf(stderr, "Failed to allocate %zd bytes; [%s]:%d\n", 
                long(sizeof(uint32_t) * count[0]), __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }
    current_data_out = (double*)malloc(sizeof(double) * count[0]);
    if(current_data_out == NULL) {
        fprintf(stderr, "Failed to allocate %zd bytes; [%s]:%d\n", 
                long(sizeof(double) * count[0]), __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }

    sample_indices.read((void*)current_indices_out, indices_dtype, indices_memspace, indices_dataspace);
    sample_data.read((void*)current_data_out, data_dtype, data_memspace, data_dataspace);

    return count[0];
}

