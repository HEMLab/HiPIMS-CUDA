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
\file cuda_backup_writer.h
\brief Header file for cuda simple field writer class

*/

#ifndef CUDA_BACKUP_WRITER_H
#define CUDA_BACKUP_WRITER_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdio>
#include "Flag.h"
#include "Vector.h"
#include "Scalar.h"
#include "cuda_mapped_field.h"

namespace GC{

  template <class T>
  void cuBackupWriter(cuFvMappedField<T, on_cell>& phi, const char* filename, Scalar t, const char* directory = "output/"){
    std::ofstream fout;
    std::ostringstream file_id, out_time;
    out_time << t;
    file_id.fill('0');
    file_id.width(3);
    std::string name = std::string(directory) + std::string(filename) + "_" + out_time.str() + ".dat";
    fout.open(name.c_str());
    if (!fout){
      std::cout << "Unable to create output file." << std::endl;
    }
    phi.data.sync();
    phi.boundary.sync();
    phi.mesh->boundary2opposite_handles.sync();
    fout << "$Element Number" << std::endl;
    fout << phi.data.size() << std::endl;
    fout << "$Element_id  Value" << std::endl;
    fout.flags(std::ios::scientific);
    auto size = phi.data.size();
    auto data_array = phi.data.host_ptr();
    auto boundary_size = phi.boundary.size();
    auto boundary_array = phi.boundary.host_ptr();
    auto boundary_cell_handles = phi.mesh->boundary2opposite_handles.host_ptr();
    for (unsigned int i = 0; i < size; ++i){
      fout << i << "  " << std::setprecision(20) << data_array[i] << std::endl;
    }
    std::map<Flag, ShortTripleFlag> boundary_type;
    for (unsigned int i = 0; i < boundary_size; ++i){
      unsigned int cell_id = boundary_cell_handles[i].get_global_id();
      boundary_type.insert({ cell_id, boundary_array[i] });
    }
    fout << "$Boundary Numbers" << std::endl;
    fout << boundary_type.size() << std::endl;
    fout << "$Element_id Boundary_type" << std::endl;
    for (auto cellid_boundary : boundary_type){
      fout << cellid_boundary.first << "  " << cellid_boundary.second << std::endl;
    }
    fout.close();
  }

}

#endif
