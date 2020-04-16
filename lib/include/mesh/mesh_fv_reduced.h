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
\file mesh_fv_reduced.h
\brief Header file for finite volume mesh class with reduced memory consumption

\version 0.1
\authot xilin xia
*/

#ifndef MESH_FV_REDUCED_H  //header file protector
#define MESH_FV_REDUCED_H

#include "mesh_fv.h"

namespace GC{

  ///this class implements data structure for the unstructured finite volume mesh with reduced memory consumption
  class unstructuredReducedFvMesh:public unstructuredFvMesh{
  public:
    unstructuredReducedFvMesh(const meshReader& mesh_reader);
    unstructuredReducedFvMesh(const meshReader&& mesh_reader);
  protected:
    virtual void BuildEntities() override; /// This function builds edges and faces  
    virtual void BuildTopDown() override; ///This function builds the top-down connectivity, i.e. face2edge, edge2vertex
    virtual void BuildBottomUp() override; ///vice-versa
    virtual void BuildGeometry() override; ///Build Geometry properties
  private:
    const char* type;
  };

}//---end namespace-------------------------------




#endif