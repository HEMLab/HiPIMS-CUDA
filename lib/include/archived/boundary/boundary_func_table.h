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