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
\flie boundary_func_table.h
\brief Header file for boundary functions

\version 0.1
\author xilin xia
*/

#ifndef BOUNDARY_FUNC_TABLE_H
#define BOUNDARY_FUNC_TABLE_H

#include <map>
#include <functional>
#include "Flag.h"
#include "Vector.h"

namespace GC{

  namespace BoundaryFuncTable{

    typedef std::function<Vector3(const Vector3& in_value, const Vector3& normal)> ZeroGradientBoundaryFunc;

    extern const std::map < Flag, ZeroGradientBoundaryFunc > ZeroGradientAndWall;

  }

}


#endif