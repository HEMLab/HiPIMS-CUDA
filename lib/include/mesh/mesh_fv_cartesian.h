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
\file mesh_fv_cartesian.h
\brief Header file for finite volume mesh class with reduced memory consumption

\version 0.1
\authot xilin xia
*/

#ifndef MESH_FV_CARTESIAN_H  //header file protector
#define MESH_FV_CARTESIAN_H

#include "mesh_fv.h"

namespace GC{

  ///this class implements data structure for the unstructured finite volume mesh with reduced memory consumption
  class CartesianFvMesh :public unstructuredFvMesh{
  public:
    CartesianFvMesh(const char* filename);
  private:
    const char* type;
  };

}


#endif