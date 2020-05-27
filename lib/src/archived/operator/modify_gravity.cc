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
 \file modify_gravity.cc
 \Source code file for modify gravity function 

 \version 0.1
 \author xilin xia
*/ 

#include "modify_gravity.h"

namespace GC{

  namespace fv{

    fvScalarFieldOnCell& modifyGravity::operator() (fvScalarFieldOnCell& gravity, fvVectorFieldOnCell& z_gradient){
      buffer.initialize_by_field(gravity);
      auto gravity_begin = gravity.data_begin();
      auto modified_gravity_begin = buffer.data_begin();
      auto z_gradient_begin = z_gradient.data_begin();
      for (auto modified_gravity_iter = buffer.data_begin(); modified_gravity_iter < buffer.data_end(); ++modified_gravity_iter){
        auto id = modified_gravity_iter - modified_gravity_begin;
        auto gradient_ = *(z_gradient_begin + id);
        auto gravity_ = *(gravity_begin + id);
        *modified_gravity_iter = gravity_/(1.0 + dot(gradient_, gradient_));
      }
      return buffer;
    }


  }//end of namespace fv

}//end of namespace GC
