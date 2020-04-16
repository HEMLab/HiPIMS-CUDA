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
\flie boundary.h
\brief Header file for boundary functions

\version 0.1
\author xilin xia
*/

#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <vector>
#include "Scalar.h"
#include "Vector.h"
#include "Flag.h"

namespace GC{

  ///This function returns boundary values according to internal values
  Vector3 GetBoundary(const ShortTripleFlag& boundary_type,
                      const Vector3& in_value,
                      const Vector3& normal);


}

#endif