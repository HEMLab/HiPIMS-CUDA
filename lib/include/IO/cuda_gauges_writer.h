// ======================================================================================
// Name                :    High-Performance Integrated Modelling System
// Description         :    This code pack provides a generic framework for developing 
//                          Geophysical CFD software. Legacy name: GeoClasses
// ======================================================================================
// Version             :    1.0.1 
// Author              :    Xilin Xia
// Create Time         :    2014/10/04
// Update Time         :    2020/04/26
// ======================================================================================
// LICENCE: GPLv3 
// ======================================================================================

/*!
\file cuda_gauges_writer.h
\brief Header file for cuda simple field writer class

*/

#ifndef CUDA_GAUGES_WRITER_H
#define CUDA_GAUGES_WRITER_H

#include <vector>
#include <fstream>
#include <sstream>
#include <memory>
#include <cuda_runtime_api.h>
#include <cuda.h>
#include "Flag.h"
#include "Scalar.h"
#include "Vector.h"
#include "mesh_interface.h"
#include "mapped_field.h"
#include "cuda_mapped_field.h"
#include "cuda_kernel_launch_parameters.h"

namespace GC{
  
  template <typename T, MAPPING_MODES C>
  class cuGaugesWriter{
    
  public:
    
    cuGaugesWriter() = default;
    cuGaugesWriter(const fvMeshQueries& _mesh, const cuFvMappedField<T, C>& phi, const char* gauge_pos_filename, const char* output_filename) :mesh(_mesh){
      std::ifstream input;
      input.open(gauge_pos_filename);
      if (!input){
        return;
      }
      else{
        std::istringstream stream;
        std::string line;
        Vector3 value;
        while (!input.eof()){
          getline(input, line);
          stream.str(line);
          stream >> value;
          gauge_postions.push_back(value);
          gauge_cell_ids.push_back(0);
          gauge_distance2cell.push_back(1e35);
          stream.str(std::string());
          stream.clear();		//clear the istringstream
        }
        for (auto cell_postions_iter = mesh.Cell.Centroid.begin(); cell_postions_iter < mesh.Cell.Centroid.end(); ++cell_postions_iter){
          Flag index = cell_postions_iter - mesh.Cell.Centroid.begin();
          for (Flag gauge_id = 0; gauge_id < gauge_postions.size(); ++gauge_id){
            Vector distance_vec = gauge_postions[gauge_id] - *cell_postions_iter;
            Scalar distance = norm(distance_vec);
            if (distance < gauge_distance2cell[gauge_id]){
              gauge_distance2cell[gauge_id] = distance;
              gauge_cell_ids[gauge_id] = index;
            }
          }
        }
        output_file.open(output_filename);
        phi_ptr = phi.data.dev_ptr();
      }
    }

    void write(Scalar t){
      if (gauge_cell_ids.size() > 0){
        output_file << t;
        for (Flag index = 0; index < gauge_cell_ids.size(); ++index){
          T phi_tmp;
          Flag cell_id = gauge_cell_ids[index];
          checkCuda(cudaMemcpy(&phi_tmp, phi_ptr + cell_id, sizeof(T), cudaMemcpyDeviceToHost));
          output_file << " " << phi_tmp;
        }
        output_file << std::endl;
      }
      
    }

  private:
    
    T* phi_ptr;
    fvMeshQueries mesh;
    std::vector<Vector> gauge_postions;
    std::vector<Scalar> gauge_distance2cell;
    std::vector<Flag> gauge_cell_ids;
    std::ofstream output_file;

  };

}



#endif