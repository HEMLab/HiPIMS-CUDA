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
 \file cuda_mesh_fv.cu
 \brief Source file for cuda mesh fv class

*/

#include "cuda_mesh_fv.h"
#include <algorithm>

namespace GC{

  cuUnstructuredFvMesh::cuUnstructuredFvMesh(fvMeshQueries mesh){

    cell_neighbours.initialize_from_host(mesh.Cell.Neighbours.begin(), mesh.Cell.Neighbours.end());
    cell_halffacets.initialize_from_host(mesh.Cell.HalfFacets.begin(), mesh.Cell.HalfFacets.end());
    halffacet2cell_handles.initialize_from_host(mesh.HalfFacet.Cell.begin(),mesh.HalfFacet.Cell.end());
    boundary2opposite_handles.initialize_from_host(mesh.Boundary.Opposite.begin(),mesh.Boundary.Opposite.end());
    halffacet_centre_positions.initialize_from_host(mesh.HalfFacet.Centroid.begin(), mesh.HalfFacet.Centroid.end()); 
    halffacet_normal_directions.initialize_from_host(mesh.HalfFacet.Normal.begin(), mesh.HalfFacet.Normal.end()); 
    halffacet_areas.initialize_from_host(mesh.HalfFacet.Area.begin(), mesh.HalfFacet.Area.end());
    cell_centre_positions.initialize_from_host(mesh.Cell.Centroid.begin(), mesh.Cell.Centroid.end());
    cell_volumes.initialize_from_host(mesh.Cell.Volume.begin(), mesh.Cell.Volume.end());

  }

}//--end of namespace GC
