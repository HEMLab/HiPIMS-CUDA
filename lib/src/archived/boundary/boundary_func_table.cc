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
\flie boundary_func_table.cc
\brief Source file for boundary functions

\version 0.1
\author xilin xia
*/

#include "boundary_func_table.h"

namespace GC{

  Vector3 ZeroGradientScalar(const Vector3& in_value, const Vector3& _normal);
  Vector3 ZeroGradientVector(const Vector3& in_value, const Vector3& _normal);
  Vector3 WallNonSlip(const Vector3& in_value, const Vector3& _normal);
  Vector3 WallSlip(const Vector3& in_value, const Vector3& _normal);

  
  namespace BoundaryFuncTable{

    const std::map < Flag, ZeroGradientBoundaryFunc > ZeroGradientAndWall =

    std::map < Flag, ZeroGradientBoundaryFunc > {
        {0, ZeroGradientScalar},
        {1, ZeroGradientVector},
        {2, WallNonSlip},
        {3, WallSlip}
    };

  }
}