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

