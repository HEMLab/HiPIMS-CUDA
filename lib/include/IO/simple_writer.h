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
  \file simple_writer.h
  \brief Header file for simple field writer class

*/

#ifndef SIMPLE_WRITER_H
#define SIMPLE_WRITER_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "Flag.h"
#include "Vector.h"
#include "Scalar.h"
#include "mapped_field.h"

namespace GC{

  template <class T>
  void simpleWriter(fvMappedField<T,on_cell>& phi, const char* filename, int cnt, Scalar t){
    std::ofstream fout;
    std::ostringstream file_id, out_time;
    out_time << t;
    file_id.fill('0');
    file_id.width(3);
    file_id << cnt;
    std::string name = std::string(filename) + "_" + file_id.str() + "_" +out_time.str() + ".dat";
    fout.open(name.c_str());
    fout.flags(std::ios::scientific);
    auto cell_iter = phi.mesh.Cell.Centroid.begin();
    for (auto phi_iter = phi.data_begin(); phi_iter < phi.data_end(); ++phi_iter, ++cell_iter){
      fout << *cell_iter << "  " << *phi_iter << std::endl;
    }
    fout.close();
  }

}//--end namespace GC


#endif
