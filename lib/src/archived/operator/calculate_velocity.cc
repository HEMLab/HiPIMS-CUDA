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
\file calculate_velocity.cc
\brief Source file for calculating velocity from discharge

\version 0.1
\author xilin xia

*/

#include "calculate_velocity.h"

namespace GC{

  namespace fv{

    fvMappedField<Vector, on_cell>& clacVelocity::operator()(fvMappedField<Scalar, on_cell>& h, fvMappedField<Vector, on_cell>& hU){
      vector_buffer.initialize_by_field(hU);
      auto h_begin = h.data_begin();
      auto hU_begin = hU.data_begin();
      auto u_begin = vector_buffer.data_begin();

      for(auto u_iter = u_begin; u_iter < vector_buffer.data_end(); ++u_iter){
        auto id = u_iter - u_begin;
        auto _h = *(h_begin + id);
        auto _hU = *(hU_begin + id);
        if (_h != 0.0){
          *u_iter = _hU/_h;
        }else{
          *u_iter = Vector(0.0);
        }
      }

      return vector_buffer;
    }

  }
}
