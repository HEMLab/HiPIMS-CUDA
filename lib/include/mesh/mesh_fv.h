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
 \file mesh_fv.h
 \brief Header file for finite volume mesh class

 \version 0.1
 \authot xilin xia
*/


#ifndef MESH_FV_H  //header file protector
#define MESH_FV_H

#include <map>
#include <utility>
#include "mesh_basic.h"
#include "mesh_fv_entities.h"

namespace GC{

  ///This class implements mesh data structure for finite volume method
  class unstructuredFvMesh:public basicMesh{
    public:
      ///default constructor - do nothing
      unstructuredFvMesh() = default;
      ///constructor - initialize by reading from file
      unstructuredFvMesh(const meshReader& mesh_reader);
      //--iterators for vertices-------------------------------------------------
      virtual size_t vertex_size();
      //cells incident of
      typedef MultiHandleSet::iterator vertex2cells_iterator;
      typedef MultiHandleSet::const_iterator vertex2cells_const_iterator;
      virtual vertex2cells_iterator vertex2cells_begin();
      virtual vertex2cells_iterator vertex2cells_end();
      //halffacet orgins from
      typedef ShortDualHandleSet::iterator vertex2halffacet_iterator;
      typedef ShortDualHandleSet::const_iterator vertex2halffacet_const_iterator;
      virtual vertex2halffacet_iterator vertex2halffacet_begin();
      virtual vertex2halffacet_iterator vertex2halffacet_end();
      //--end iterators for vertices---------------------------------------------
      
      //--iterators for half facets-------------------------------------
      virtual size_t halffacet_size();
      //centre positions
      typedef vectorProps::iterator halffacet_centre_positions_iterator;
      typedef vectorProps::const_iterator halffacet_centre_positions_const_iterator;
      virtual halffacet_centre_positions_iterator
              halffacet_centre_positions_begin();
      virtual halffacet_centre_positions_iterator
              halffacet_centre_positions_end();  
      //areas
      typedef scalarProps::iterator halffacet_areas_iterator;
      typedef scalarProps::const_iterator halffacet_areas_const_iterator;
      virtual halffacet_areas_iterator halffacet_areas_begin();
      virtual halffacet_areas_iterator halffacet_areas_end();
      //normal directions
      typedef vectorProps::iterator halffacet_normal_directions_iterator;
      typedef vectorProps::const_iterator halffacet_normal_directions_const_iterator;
      virtual halffacet_normal_directions_iterator halffacet_normal_directions_begin();
      virtual halffacet_normal_directions_iterator halffacet_normal_directions_end();
      //face belongs to
      typedef ShortDualHandleSet::iterator halffacet2cell_iterator;
      typedef ShortDualHandleSet::const_iterator halffacet2cell_const_iterator;
      virtual halffacet2cell_iterator halffacet2cell_begin();
      virtual halffacet2cell_iterator halffacet2cell_end();
      //--end iterators for half facets---------------------------------

      //--iterators for cells----------------------------------------
      virtual size_t cell_size();
      typedef vectorProps::iterator cell_centre_positions_iterator;
      typedef vectorProps::const_iterator cell_centre_positions_const_iterator;
      virtual cell_centre_positions_iterator cell_centre_positions_begin();
      virtual cell_centre_positions_iterator cell_centre_positions_end();

      typedef scalarProps::iterator cell_volumes_iterator;
      typedef scalarProps::const_iterator cell_volumes_const_iterator;
      virtual cell_volumes_iterator cell_volumes_begin();
      virtual cell_volumes_iterator cell_volumes_end();
      
      typedef CellNeighbourSet::iterator cell_neighbours_iterator;
      typedef CellNeighbourSet::const_iterator cell_neighbours_const_iterator;
      virtual cell_neighbours_iterator cell_neighbours_begin();
      virtual cell_neighbours_iterator cell_neighbours_end();

      typedef MultiHandleSet::iterator cell_halffacets_iterator;
      typedef MultiHandleSet::const_iterator cell_halffacets_const_iterator;
      virtual cell_halffacets_iterator cell_halffacets_begin();
      virtual cell_halffacets_iterator cell_halffacets_end();
      //--end iterators for faces-------------------------------------

      //--iterators for boundary half facets--------------------------
      virtual size_t boundary_size();
      typedef ShortDualHandleSet::iterator boundary2opposite_iterator;
      typedef ShortDualHandleSet::const_iterator boundary2opposite_const_iterator;
      virtual boundary2opposite_iterator boundary2opposite_begin();
      virtual boundary2opposite_iterator boundary2opposite_end();
      //--end iterators for boundary half facets----------------------

      ///Return the type of this mesh
      virtual const char* type_name(); 
    protected:
      //all private member functions
      virtual void BuildEntities(); /// This function builds edges and faces  
   //   virtual void SortMesh(); /// This function sorts the entire mesh according to space filling curve
      virtual void BuildTopDown(); ///This function builds the top-down connectivity, i.e. face2edge, edge2vertex
      virtual void BuildBottomUp(); ///vice-versa
      virtual void BuildGeometry(); ///Build Geometry properties
      MultiHandle next_order_hffacet(const MultiHandle& hffacet);
      std::map<MultiHandle, ShortDualHandle> halffacet2cell_map;
      std::vector < std::pair<MultiHandle, ShortDualHandle>> halffacetcell_pairs;
    public:
      //all data members
      //--Comprehensive mesh data structure
      ShortDualHandleSet vertex2halffacet_handles; //x - (>0) local index (0)boundary y - face id or boundary edge id
      MultiHandleSet  vertex2cell_handles;
      ShortDualHandleSet halffacet2cell_handles; //x - local index y - face id
      ShortDualHandleSet boundary2opposite_handles;
      CellNeighbourSet cell_neighbours;
      MultiHandleSet cell_halffacets;
      //--Geometric properties
      vectorProps halffacet_centre_positions;
      vectorProps halffacet_normal_directions;
      scalarProps halffacet_areas;
      vectorProps cell_centre_positions;
      scalarProps cell_volumes;      
      //--Sorted indices for synchronization between mesh and attached field data
      HandleSet vertex_indices;
      HandleSet halffacet_indices;
      HandleSet cell_indices;
      //--type
    private:
      const char* type;
  };

}

#endif  //end header file protector
