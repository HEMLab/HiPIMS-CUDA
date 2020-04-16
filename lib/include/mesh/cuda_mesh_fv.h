// ======================================================================================
// Name                :    GeoClasses : Generic Geophysical Flow Modelling Framework
// Description         :    This code pack provides a generic framework for developing 
//                          Geophysical CFD software.
// ======================================================================================
// Version             :    0.1 
// Author              :    Xilin Xia (PhD candidate in Newcastle University)
// Create Time         :    2014/10/04
// Update Time         :    2015/10/09
// ======================================================================================
// Copyright @ Xilin Xia 2015 . All rights reserved.
// ======================================================================================

/*!
 \file cuda_mesh_fv.h
 \brief Header file for finite volume mesh class in cuda

 \version 0.1
 \authot xilin xia
*/

#ifndef CUDA_MESH_FV_H
#define CUDA_MESH_FV_H

#include "cuda_mesh_fv_entities.h"
#include "mesh_fv_queries.h"

namespace GC{

  ///This class implements mesh data structure for finite volume method with cuda
  class cuUnstructuredFvMesh{
    public:
      cuUnstructuredFvMesh(fvMeshQueries mesh);
    public:
      cu2dShortDualHandleList cell_neighbours;
      cu2dHandleList cell_halffacets;
      cuShortDualHandleList halffacet2cell_handles;
      cuShortDualHandleList boundary2opposite_handles;
      cuVectorList halffacet_centre_positions;
      cuVectorList halffacet_normal_directions;
      cuScalarList halffacet_areas;
      cuVectorList cell_centre_positions;
      cuScalarList cell_volumes;
  };


}//--end of namespace GC

#endif
