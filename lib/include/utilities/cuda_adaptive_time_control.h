// ======================================================================================
// Name                :    GeoClasses : Generic Geophysical Flow Modelling Framework
// Description         :    This code pack provides a generic framework for developing 
//                          Geophysical CFD software.
// ======================================================================================
// Version             :    0.1 
// Author              :    Xilin Xia (PhD candidate in Newcastle University)
// Create Time         :    2014/10/04
// Update Time         :    2015/10/25
// ======================================================================================
// Copyright @ Xilin Xia 2015 . All rights reserved.
// ======================================================================================

/*!
\file cuda_adaptive_time_control.h
\brief Header file for adaptive time control class

\version 0.1
\author xilin xia

*/

#ifndef CUDA_ADAPTIVE_TIME_CONTROL_H
#define CUDA_ADAPTIVE_TIME_CONTROL_H

#include "Scalar.h"
#include <memory>
#include "time_control.h"
#include "cuda_mapped_field.h"

namespace GC{

  class cuAdaptiveTimeControl2D{
  public:
    void updateByCFL(cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Vector2, on_cell>& hU);
  public:
    cuAdaptiveTimeControl2D(Scalar _dt, Scalar _end, Scalar _Courant = 0.5, Scalar _current = 0.0) :
      dt_(_dt), end_(_end), Courant_(_Courant), current_(_current){
      is_beginning = true;
    }
    bool is_end();
    Scalar dt();
    Scalar current();
    void forward();
    void set_dt(Scalar dt);
  private:
    Scalar dt_;
    const Scalar end_;
    const Scalar Courant_;
    Scalar current_;
    bool is_beginning;
    std::shared_ptr<cuFvMappedField<Scalar, on_cell> > dt_array_ptr;
  };

}


#endif
