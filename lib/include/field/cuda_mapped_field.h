// ======================================================================================
// Name                :    GeoClasses : Generic Geophysical Flow Modelling Framework
// Description         :    This code pack provides a generic framework for developing 
//                          Geophysical CFD software.
// ======================================================================================
// Version             :    0.1 
// Author              :    Xilin Xia (PhD candidate in Newcastle University)
// Create Time         :    2014/10/04
// Update Time         :    2015/10/09
// ======================================================================================
// Copyright @ Xilin Xia 2015 . All rights reserved.
// ======================================================================================

/*!
 \file cuda_mapped_field.h
 \brief Header file for cuda mapped field class 

 \version 0.1
 \author xilin xia
*/ 


#ifndef CUDA_MAPPED_FIELD_H
#define CUDA_MAPPED_FIELD_H

#include "mapped_field.h"
#include <memory>
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include "cuda_arrays.h"
#include "cuda_mesh_fv.h"
#include "Flag.h"
//#include "cuda_boundary.h"
#include "cuda_kernel_launch_parameters.h"
#include <thrust/device_vector.h>

namespace GC{

  template <typename T>
  __global__ void update_boundary_values_kernel(ShortDualHandle* boundary2opposite_handles, unsigned int boundary2opposite_handles_length, ShortDualHandle* cell_neighbours, Flag* cell_halffacets, unsigned int cell_halffacets_length, Vector* halffacet_normal, T* data, ShortTripleFlag* boundary_type, T* boundary_value, T* boundary_source, unsigned int boundary_source_length, unsigned int* cursors, Scalar* time_series, Scalar current_t){

    unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
    while (index < boundary2opposite_handles_length){
      auto opposite_cell = boundary2opposite_handles[index];
      Flag cell_id = opposite_cell.get_global_id();
      Flag halffacet_cnt = opposite_cell.get_local_id();
      Flag halffacet_id = cell_halffacets[halffacet_cnt*cell_halffacets_length + index];
      auto in_value = data[cell_id];
      auto normal = halffacet_normal[halffacet_id];
      auto _boundary_type = boundary_type[index];
      Flag primary_type = _boundary_type.getx();
      Flag secondary_type = _boundary_type.gety();
      Flag source_id = _boundary_type.getz();
      if (0 == primary_type){ //cannot evolve
//        std::cout << "Type is 0, can not be updated!" << std::endl;
      }
      else if (1 == primary_type){ //fixed gradient 

      }
      else if (2 == primary_type){ //zero gradient or wall
        switch (secondary_type){
        case 0:
          boundary_value[index] = in_value;
          break;
        case 1:
          boundary_value[index] = in_value;
          break;
        case 2:
          boundary_value[index] = -1.0*in_value;
          break;
        case 3:
          Vector3 _in_value = in_value;
          Vector3 _normal = uni(normal);
          Vector3 shearx = uni(perpend(_normal));
          Vector3 sheary = uni(cross(_normal, shearx));
          Scalar normal_value = -1.0*dot(_in_value, _normal);
          Scalar shearx_value = dot(_in_value, shearx);
          Scalar sheary_value = dot(_in_value, sheary);
          boundary_value[index] = normal_value*_normal + shearx_value*shearx + sheary_value*sheary;
          break;
        }
      }
      else if (3 == primary_type){ //fixed value
        switch (secondary_type){
        case 0:
          unsigned int id = cursors[source_id];
          auto source_value_first = boundary_source[(id - 1)*boundary_source_length + source_id];
          auto source_value_second = boundary_source[id*boundary_source_length + source_id];
          auto time_first = time_series[(id - 1)*boundary_source_length + source_id];
          auto time_second = time_series[id*boundary_source_length + source_id];
          if (current_t <= time_second){
            boundary_value[index] = source_value_first + (source_value_second - source_value_first)*(current_t - time_first) / (time_second - time_first);
          }
          else{
            boundary_value[index] = source_value_second;
          }
          break;
        }

      }
      else if (4 == primary_type){//exchanged from neighbouring domains
        boundary_value[index] = in_value;
      }
      index += blockDim.x * gridDim.x;
    }

  }

  template <typename T>
  __global__ void update_data_values_kernel(Flag* region_mask, T* data, unsigned int data_length, T* data_source_current){
    unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
    while (index < data_length){
      auto source_id = region_mask[index];
      data[index] = data_source_current[source_id];
      index += blockDim.x * gridDim.x;
    }
  }

  ///This class implements field data mapped on finite volume mesh by CUDA
  template <typename T, MAPPING_MODES C>
  class cuFvMappedField{
    public:
      ///constructor initialize by field on host
      cuFvMappedField(fvMappedField<T,C>& field_host, std::shared_ptr<cuUnstructuredFvMesh> _mesh) : mesh(_mesh){
        data.initialize_from_host(field_host.data_begin(), field_host.data_end());
        boundary.initialize_from_host(field_host.boundary_type_begin(), field_host.boundary_type_end());
        boundary_value.initialize_from_host(field_host.boundary_value_begin(), field_host.boundary_value_end());
        time_series.initialize_from_host(field_host.time_series_begin(), field_host.time_series_end());
        boundary_source.initialize_from_host(field_host.boundary_source_begin(), field_host.boundary_source_end());
        region_mask.initialize_from_host(field_host.region_mask_begin(), field_host.region_mask_end());
        data_time_series.initialize_from_host(field_host.data_time_series_begin(), field_host.data_time_series_end());
        data_source.initialize_from_host(field_host.data_source_begin(), field_host.data_source_end());
        auto boundary_source_size = boundary_source.length();
        cursors.resize(boundary_source_size); ///there are maximum 100 inflow points
        for (unsigned int i = 0; i < boundary_source_size; i++){
          (cursors.host_ptr())[i] = 1;
        }
        cursors.sync_by_host();
        auto data_source_size = data_source.length();
        data_source_current.resize(data_source_size);
        data_source_cursors.resize(data_source_size);
        for (unsigned int i = 0; i < data_source_size; i++){
          data_source_cursors[i] = 1;
        }
      }
      ///constructor initialize by mesh only, data will be zero
      cuFvMappedField(std::shared_ptr<cuUnstructuredFvMesh> _mesh) : mesh(_mesh){
        switch (C){
        case on_cell: //on cell
          data.resize(mesh->cell_neighbours.length());
          boundary.resize(mesh->boundary2opposite_handles.size()); ///to do: boundary type information on either cell or half facet may not be needed
          boundary_value.resize(mesh->boundary2opposite_handles.size());
          break;
        case on_halffacet: //on halffacet
          data.resize(mesh->halffacet_normal_directions.size());
          boundary.resize(mesh->boundary2opposite_handles.size());
          break;
        case on_vertex: //on vertex
          break;
        }
      }
      ///copy constructor, it has two modes: full and partial, default value is full
      cuFvMappedField(const cuFvMappedField& other, COPY_MODES c_m = full) : mesh(other.mesh){
        if(c_m == full){
          data = other.data;
          boundary = other.boundary;
          boundary_value = other.boundary_value;
          time_series = other.time_series;
          boundary_source = other.boundary_source;
          data_time_series = other.data_time_series;
          data_source = other.data_source;
          region_mask = other.region_mask;
          data_source_cursors = other.data_source_cursors;
          data_source_current = other.data_source_current;
          cursors = other.cursors;
          current_t = other.current_t;
        }else{
          data.resize(other.data.size());
          boundary = other.boundary;
          boundary_value = other.boundary_value;
          time_series = other.time_series;
          boundary_source = other.boundary_source;
          data_time_series = other.data_time_series;
          data_source = other.data_source;
          region_mask = other.region_mask;
          data_source_cursors = other.data_source_cursors;
          data_source_current = other.data_source_current;
          cursors = other.cursors;
          current_t = other.current_t;
        }
      }

    ///update the boundary values
      void update_boundary_values(){
        //Firstly update the cursors
        for (unsigned int i = 0; i < time_series.length(); i++){
          unsigned int max_time_points = (time_series.dims_host_ptr())[i] - 1;
          unsigned int old_cursor = (cursors.host_ptr())[i];
          Scalar old_cursor_time = (time_series.host_ptr())[old_cursor*time_series.length() + i];
          while (old_cursor_time < current_t && old_cursor < max_time_points){
            (cursors.host_ptr())[i] += 1;
            old_cursor = (cursors.host_ptr())[i];
            old_cursor_time = (time_series.host_ptr())[old_cursor*time_series.length() + i];
          }
          cursors.sync_by_host();
        }

        //secondly update the boundary values
        update_boundary_values_kernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(
          mesh->boundary2opposite_handles.dev_ptr(), 
          mesh->boundary2opposite_handles.size(), 
          mesh->cell_neighbours.dev_ptr(), 
          mesh->cell_halffacets.dev_ptr(), 
          mesh->cell_halffacets.length(),         
          mesh->halffacet_normal_directions.dev_ptr(), 
          data.dev_ptr(), 
          boundary.dev_ptr(), 
          boundary_value.dev_ptr(), 
          boundary_source.dev_ptr(), 
          boundary_source.length(),
          cursors.dev_ptr(), 
          time_series.dev_ptr(), 
          current_t);

      }

      void update_data_values(){
        //Firstly update the cursors
        for (unsigned int i = 0; i < data_time_series.length(); i++){
          unsigned int max_time_points = (data_time_series.dims_host_ptr())[i] - 1;
          unsigned int old_cursor = data_source_cursors[i];
          Scalar old_cursor_time = (data_time_series.host_ptr())[old_cursor*data_time_series.length() + i];
          while (old_cursor_time < current_t && old_cursor < max_time_points){
            data_source_cursors[i] += 1;
            old_cursor = data_source_cursors[i];
            old_cursor_time = (data_time_series.host_ptr())[old_cursor*data_time_series.length() + i];
          }
          unsigned int id = data_source_cursors[i];
          unsigned int length = data_time_series.length();
          auto source_value_first = (data_source.host_ptr())[(id - 1)*length + i];
          auto source_value_second = (data_source.host_ptr())[id*length + i];
          auto time_first = (data_time_series.host_ptr())[(id - 1)*length + i];
          auto time_second = (data_time_series.host_ptr())[id*length + i];
          if (current_t <= time_second){
            (data_source_current.host_ptr())[i] = source_value_first + (source_value_second - source_value_first)*(current_t - time_first) / (time_second - time_first);
          }
          else{
            (data_source_current.host_ptr())[i] = source_value_second;
          }
        }
        data_source_current.sync_by_host();
        update_data_values_kernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(region_mask.dev_ptr(), data.dev_ptr(), data.size(), data_source_current.dev_ptr());
      }

      void update_time(const Scalar& t, const Scalar& dt){
        current_t = t + dt;
      }

      void update_data_source(const char* path, const char* field_name){
        std::vector< std::vector<Scalar> > data_time_series_input;
        std::vector< std::vector<T> > data_source_input;
        for (unsigned int i = 0;; ++i){
          std::string filename;
          std::ostringstream file_id;
          file_id << i;
          filename = std::string(path) + std::string(field_name) + "_source_" + file_id.str() + ".dat";
          std::ifstream input;
          input.open(filename);
          if (!input){
            break;
          }
          else{
            data_time_series_input.resize(i + 1);
            data_source_input.resize(i + 1);
            std::istringstream stream;
            std::string line;
            Scalar t;
            Vector3 value;
            while (!input.eof()){
              getline(input, line);
              stream.str(line);
              stream >> t >> value;
              data_time_series_input[i].push_back(t);
              data_source_input[i].push_back(value);
              stream.str(std::string());
              stream.clear();		//clear the istringstream
            }
          }
        }
        data_time_series.initialize_from_host(data_time_series_input.begin(), data_time_series_input.end());
        data_source.initialize_from_host(data_source_input.begin(), data_source_input.end());
        update_data_values();
      }

      void update_boundary_source(const char* path, const char* field_name){
        std::vector< std::vector<Scalar> > time_series_input;
        std::vector< std::vector<T> > boundary_source_input;
        for (unsigned int i = 0;; ++i){
          std::string filename;
          std::ostringstream file_id;
          file_id << i;
          filename = std::string(path) + std::string(field_name) + "_BC_" + file_id.str() + ".dat";
          std::ifstream input;
          input.open(filename);
          if (!input){
            break;
          }
          else{
            time_series_input.resize(i + 1);
            boundary_source_input.resize(i + 1);
            std::istringstream stream;
            std::string line;
            Scalar t;
            Vector3 value;
            while (!input.eof()){
              getline(input, line);
              stream.str(line);
              stream >> t >> value;
              time_series_input[i].push_back(t);
              boundary_source_input[i].push_back(value);
              stream.str(std::string());
              stream.clear();		//clear the istringstream
            }
          }
        }
        time_series.initialize_from_host(time_series_input.begin(), time_series_input.end());
        boundary_source.initialize_from_host(boundary_source_input.begin(), boundary_source_input.end());
        update_boundary_values();
      }

    public:
      std::shared_ptr<cuUnstructuredFvMesh> mesh;
    public:
      cuArray<T> data;
      cuArray<ShortTripleFlag> boundary;
      cuArray<T> boundary_value;
      cu2dArray<Scalar> time_series;
      cuArray<Flag> region_mask;
      cu2dArray<T> boundary_source;
      cu2dArray<Scalar> data_time_series;
      cu2dArray<T> data_source;
      thrust::device_vector<T> received_upper; //resize when called by syncronizer
      thrust::device_vector<T> received_lower;
    protected:
      Scalar current_t;
      cuArray<unsigned int> cursors;
      cuArray<T> data_source_current;
      std::vector<unsigned int> data_source_cursors;
  };

}//--end of namespace GC----------

#endif

