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
\flie boundary_funcs.cc
\brief Source file for boundary functions

\version 0.1
\author xilin xia
*/

#include "Scalar.h"
#include "Vector.h"

namespace GC{


  //The returning values here are all Vector3, they will be automatically transformed
  Vector3 ZeroGradientScalar(const Vector3& in_value, const Vector3& _normal){
    return Vector3(in_value.x);
  }

  Vector3 ZeroGradientVector(const Vector3& in_value, const Vector3& _normal){
    return Vector3(in_value);
  }

  Vector3 WallNonSlip(const Vector3& in_value, const Vector3& _normal){
    return Vector3(-1.0*in_value);
  }

  Vector3 WallSlip(const Vector3& in_value, const Vector3& _normal){
    Vector3 normal = uni(_normal);
    Vector3 shearx = uni(perpend(normal));
    Vector3 sheary = uni(cross(normal, shearx));
    Scalar normal_value = -1.0*dot(in_value, normal);
    Scalar shearx_value = dot(in_value, shearx);
    Scalar sheary_value = dot(in_value, sheary);
    return normal_value*normal + shearx_value*shearx + sheary_value*sheary;
  }


}