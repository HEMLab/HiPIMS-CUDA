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