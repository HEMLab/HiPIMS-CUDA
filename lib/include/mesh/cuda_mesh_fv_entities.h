// ======================================================================================
// Name                :    GeoClasses : Generic Geophysical Flow Modelling Framework
// Description         :    This code pack provides a generic framework for developing 
//                          Geophysical CFD software.
// ======================================================================================
// Version             :    0.2 
// Author              :    Xilin Xia (PhD candidate in Newcastle University)
// Create Time         :    2014/10/04
// Update Time         :    2015/10/08
// ======================================================================================
// Copyright @ Xilin Xia 2015 . All rights reserved.
// ======================================================================================

/*!
 \file cuda_mesh_fv_entities.h
 \brief Header file for basic mesh entities for cuda 

*/

#ifndef CUDA_MESH_FV_ENTITIES_H
#define CUDA_MESH_FV_ENTITIES_H

#include "mesh_fv_entities.h"
#include "cuda_arrays.h"

namespace GC{

  typedef cuArray<ShortDualHandle> cuShortDualHandleList;
  typedef cuArray<Flag> cuHandleList;
  typedef cuArray<Scalar> cuScalarList;
  typedef cuArray<Vector> cuVectorList;
  typedef cu2dArray<ShortDualHandle> cu2dShortDualHandleList;
  typedef cu2dArray<Flag> cu2dHandleList;

} //--end of namespace GC------------


#endif
