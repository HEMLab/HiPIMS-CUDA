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
\file adaptive_time_control.h
\brief Header file for adaptive time control class

\version 0.1
\author xilin xia

*/

#ifndef ADAPTIVE_TIME_CONTROL_H
#define ADAPTIVE_TIME_CONTROL_H

#include "Scalar.h"
#include "time_control.h"
#include "mapped_field.h"

namespace GC{

  class adaptiveTimeControl2D:public timeControl{
  public:
    adaptiveTimeControl2D(Scalar _dt, Scalar _end, Scalar _Courant = 0.5):timeControl(_dt, _end, _Courant){}
    void updateByCFL(fvScalarFieldOnCell gravity, fvScalarFieldOnCell h, fvVectorFieldOnCell hU);
  };

}


#endif
