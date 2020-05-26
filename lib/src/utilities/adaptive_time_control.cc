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
\file adaptive_time_control.cc
\brief Source file for adaptive time control class

\version 0.1
\author xilin xia

*/

#include "adaptive_time_control.h"
#include <algorithm>

namespace GC{

  void adaptiveTimeControl2D::updateByCFL(fvScalarFieldOnCell gravity, fvScalarFieldOnCell h, fvVectorFieldOnCell hU){
    auto gravity_begin = gravity.data_begin();
    auto h_begin = h.data_begin();
    auto hU_begin = hU.data_begin();
    auto cell_volume_begin = h.mesh.Cell.Volume.begin();

    Scalar dt = 1e10;

    for(auto h_iter = h_begin; h_iter < h.data_end(); ++h_iter){
      auto id = h_iter - h_begin;
      auto gravity_ = *(gravity_begin + id);
      auto h_ = *(h_iter);
      auto hU_ = *(hU_begin + id);
      auto c_sound = sqrt(gravity_*h_);
      if (h_ > 1e-10){
        Vector u_ = hU_/h_;
        Scalar volume_ = *(cell_volume_begin + id);
        Scalar dl = sqrt(volume_);
        Scalar dt1 = dl/(c_sound + fabs(u_.x));
        Scalar dt2 = dl/(c_sound + fabs(u_.y));
        Scalar dt_tmp = Courant_*std::min(dt1,dt2);
        dt = std::min(dt, dt_tmp);
      }
    }

    dt_ = dt;

  }

}
