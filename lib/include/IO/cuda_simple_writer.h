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
  \file cuda_simple_writer.h
  \brief Header file for cuda simple field writer class

*/

#ifndef CUDA_SIMPLE_WRITER_H
#define CUDA_SIMPLE_WRITER_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include "Flag.h"
#include "Vector.h"
#include "Scalar.h"
#include "cuda_mapped_field.h"

namespace GC{

  template <class T>
  void cuSimpleWriter(cuFvMappedField<T,on_cell>& phi, const char* filename, int cnt, Scalar t){
    std::ofstream fout;
    std::ostringstream file_id, out_time;
    out_time << t;
    file_id.fill('0');
    file_id.width(3);
    file_id << cnt;
    std::string name = "output/" + std::string(filename) + "_" + file_id.str() + "_" + out_time.str() + ".dat";
    fout.open(name.c_str());
    if (!fout){
      std::cout << "Unable to create output file." << std::endl;
    }
    fout.flags(std::ios::scientific);
    phi.data.sync();
    auto size = phi.data.size();
    auto positions_array = phi.mesh->cell_centre_positions.host_ptr();
    auto data_array = phi.data.host_ptr();
    for (unsigned int i = 0; i < size; ++i){
		fout << positions_array[i] << "  " << data_array[i] << std::endl;
    }
    fout.close();
  }

  template <class T>
  void cuSimpleWriterLowPrecision(cuFvMappedField<T, on_cell>& phi, const char* filename, Scalar t, const char* directory = "output/"){
	  std::ofstream fout;
	  std::ostringstream file_id, out_time,out_string;
	  out_time << t;
	  std::string name = std::string(directory) + std::string(filename) + "_" + out_time.str() + ".dat";
	  fout.open(name.c_str());
	  if (!fout){
		  std::cout << "Unable to create output file." << std::endl;
	  }
    out_string.flags(std::ios::scientific);
	  phi.data.sync();
	  auto size = phi.data.size();
	  auto positions_array = phi.mesh->cell_centre_positions.host_ptr();
	  auto data_array = phi.data.host_ptr();
	  for (unsigned int i = 0; i < size; ++i){
		  //fout.precision(3);
      out_string << std::fixed << std::setprecision(4) << positions_array[i] << "  " << std::scientific << std::setprecision(6) << data_array[i] << std::endl;
	  }
    std::string s = out_string.str();
    fout.write(s.c_str(), strlen(s.c_str()));
	  fout.close();
  }

  template <class T>
  void cuSimpleWriterLowPrecisionNoCoordinates(cuFvMappedField<T, on_cell>& phi, const char* filename, Scalar t){
    std::ofstream fout;
    std::ostringstream file_id, out_time;
    out_time << t;
    std::string name = "output/" + std::string(filename) + "_" + out_time.str() + ".dat";
    fout.open(name.c_str());
    if (!fout){
      std::cout << "Unable to create output file." << std::endl;
    }
    fout.flags(std::ios::scientific);
    phi.data.sync();
    auto size = phi.data.size();
    auto positions_array = phi.mesh->cell_centre_positions.host_ptr();
    auto data_array = phi.data.host_ptr();
    for (unsigned int i = 0; i < size; ++i){
      //fout.precision(3);
      fout << std::scientific << std::setprecision(6) << data_array[i] << std::endl;
    }
    fout.close();
  }

  template <class T>
  void cuSimpleWriterHighPrecision(cuFvMappedField<T, on_cell>& phi, const char* filename, Scalar t){
	  std::ofstream fout;
	  std::ostringstream file_id, out_time;
	  out_time << t;
	  std::string name = "output/" + std::string(filename) + "_" + out_time.str() + ".dat";
	  fout.open(name.c_str());
	  if (!fout){
      std::cout << "Unable to create output file." << std::endl;
	  }
	  fout.flags(std::ios::scientific);
	  phi.data.sync();
	  auto size = phi.data.size();
	  auto positions_array = phi.mesh->cell_centre_positions.host_ptr();
	  auto data_array = phi.data.host_ptr();
	  for (unsigned int i = 0; i < size; ++i){
      fout << std::scientific << std::setprecision(20) << positions_array[i] << "  " << std::setprecision(20) << data_array[i] << std::endl;
	  }
	  fout.close();
  }

}//--end namespace GC


#endif
