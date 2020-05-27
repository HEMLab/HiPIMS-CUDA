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