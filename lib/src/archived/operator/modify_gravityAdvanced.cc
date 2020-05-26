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
 \file modify_gravityAdvanced.cc
 \Source code file for advanced modify gravity function 

 \version 0.1
 \author xilin xia
*/ 

#include "modify_gravityAdvanced.h"
#include <algorithm>

namespace GC{

  namespace fv{

    fvScalarFieldOnCell& modifyGravityAdvanced::operator() (fvScalarFieldOnCell& gravity, fvVectorFieldOnCell& z_gradient, fvVectorFieldOnCell& h_gradient){
      buffer.initialize_by_field(gravity);
      auto gravity_begin = gravity.data_begin();
      auto modified_gravity_begin = buffer.data_begin();
      auto z_gradient_begin = z_gradient.data_begin();
      auto h_gradient_begin = h_gradient.data_begin();
      for (auto modified_gravity_iter = buffer.data_begin(); modified_gravity_iter < buffer.data_end(); ++modified_gravity_iter){
        auto id = modified_gravity_iter - modified_gravity_begin;
        auto z_gradient_ = *(z_gradient_begin + id);
        auto h_gradient_ = *(h_gradient_begin + id);
        auto gravity_ = *(gravity_begin + id);
        auto new_gravity = gravity_/(1.0 + dot(z_gradient_, z_gradient_ + h_gradient_));
        new_gravity = std::max(0.0, std::min(gravity_,new_gravity));
        *modified_gravity_iter = new_gravity;
      }
      return buffer;
    }

  }

}
