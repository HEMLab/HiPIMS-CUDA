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
\file friction.h
\brief Header file for Mohr Columb friction formula

\version 0.1
\author xilin xia

*/

#ifndef MOHR_COLUMB_H
#define MOHR_COLUMB_H

#include "Scalar.h"
#include "Vector.h"

namespace GC{

  namespace fv{

    Vector MCFriction(Scalar& g, Scalar& miu, Scalar& h, Vector& hU, Vector& a);

  }//end of namespace fv

}//end of namespace GC

#endif

