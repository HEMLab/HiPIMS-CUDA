// ======================================================================================
// Name                :    GeoClasses : Generic Geophysical Flow Modelling Framework
// Description         :    This code pack provides a generic framework for developing 
//                          Geophysical CFD software.
// ======================================================================================
// Version             :    0.1 
// Author              :    Xilin Xia (PhD candidate in Newcastle University)
// Create Time         :    2014/10/04
// Update Time         :    2015/10/29
// ======================================================================================
// Copyright @ Xilin Xia 2015 . All rights reserved.
// ======================================================================================

/*!
\flie boundary_funcs.h
\brief Header file for boundary functions

\version 0.1
\author xilin xia
*/

#ifndef BOUNDARY_FUNCS_H
#define BOUNDARY_FUNCS_H

namespace GC{

  Vector3 ZeroGradientScalar(const Vector3& in_value, const Vector3& _normal);
  Vector3 ZeroGradientVector(const Vector3& in_value, const Vector3& _normal);
  Vector3 WallNonSlip(const Vector3& in_value, const Vector3& _normal);
  Vector3 WallSlip(const Vector3& in_value, const Vector3& _normal);

}//end of namespace GC

#endif