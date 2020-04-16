// ======================================================================================
// Name                :    GeoClasses : Generic Geophysical Flow Modelling Framework
// Description         :    This code pack provides a generic framework for developing 
//                          Geophysical CFD software.
// ======================================================================================
// Version             :    0.1 
// Author              :    Xilin Xia (PhD candidate in Newcastle University)
// Create Time         :    2014/10/04
// Update Time         :    2012/10/29
// ======================================================================================
// Copyright @ Xilin Xia 2014 . All rights reserved.
// ======================================================================================

/*!
\file time_control.h
\brief Header file for time control class

\version 0.1
\author xilin xia

*/

#ifndef TIME_CONTROL_H
#define TIME_CONTROL_H

#include "Scalar.h"

namespace GC{

  class timeControl{
  public:
    timeControl(Scalar _dt, Scalar _end, Scalar _Courant = 0.5):
      dt_(_dt), end_(_end), Courant_(_Courant), current_(0.0){}
    bool is_end();
    Scalar dt();
    Scalar current();
    void forward();
  protected:
    Scalar dt_;
    const Scalar end_;
    const Scalar Courant_;
    Scalar current_;
  };

}


#endif
