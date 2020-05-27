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
\file mesh_fv_reduced.cc
\brief Source file for finite volume mesh class

*/

#include "mesh_fv_reduced.h"
#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include "hilbert.h"
#include "element_table.h"


namespace GC{

  unstructuredReducedFvMesh::unstructuredReducedFvMesh(const meshReader& mesh_reader):type("unstructuredReducedFvMesh")
  {
    BuildElementTable(mesh_reader);
    BuildEntities();
  }

  unstructuredReducedFvMesh::unstructuredReducedFvMesh(const meshReader&& mesh_reader) :type("unstructuredReducedFvMesh")
  {
    BuildElementTable(std::move(mesh_reader));
    BuildEntities();
  }

  ///This function builds the comprehensive mesh data structure
  void unstructuredReducedFvMesh::BuildEntities(){
    BuildTopDown();
    BuildBottomUp();
    BuildGeometry();
  }

  void unstructuredReducedFvMesh::BuildTopDown(){
    cell_halffacets.resize(cell_vertices.size());
    Flag id_halffacets = 0;
    Flag id_cells = 0;
    ///This loop builds all the facets and faces
    for (auto& c_it : cell_vertices){
      Flag cell_type = cell_types[id_cells];
      auto cell_property = ElementPropertyTable.find(cell_type)->second; //to do: throw an exception when type value is invalid
      Flag id_local_halffacet = 0;
      for (auto local_halffacet : cell_property.facets){
        MultiHandle halffacet_;
        Flag facet_type = local_halffacet.type_facet;
        for (auto local_vertex_id : local_halffacet.vertices){
          halffacet_.push_back(c_it[local_vertex_id]);
        }
        halffacetcell_pairs.push_back(std::make_pair(halffacet_, ShortDualHandle(id_local_halffacet + 1, id_cells)));
        halffacet2cell_handles.push_back(ShortDualHandle(id_local_halffacet + 1, id_cells));
        id_local_halffacet++;
        cell_halffacets[id_cells].push_back(id_halffacets);
        id_halffacets++;
      }
      id_cells++;
    }

    //sort this vector for binary search
    std::sort(halffacetcell_pairs.begin(), halffacetcell_pairs.end(), [](const std::pair<MultiHandle, ShortDualHandle>& a,
                                                                         const std::pair<MultiHandle, ShortDualHandle>& b)
                                                                                   {return a.first < b.first; });

  }

  ///This function builds other datas
  void unstructuredReducedFvMesh::BuildBottomUp(){
    Flag id_boundary = 0;
    Flag id_cells = 0;
    cell_neighbours.resize(cell_vertices.size());  //request memory for cell_neighbours
    for (auto& it : cell_neighbours){
      it.resize(cell_halffacets[id_cells++].size());
    }

    for (auto& it : halffacetcell_pairs){
      MultiHandle _facet = it.first;
      ShortDualHandle _facet2cell_handle = it.second;
      MultiHandle _opposite_facet = _facet;
      std::vector < std::pair<MultiHandle, ShortDualHandle>>::iterator opposite_it;
      bool opposite_exists = false;
      bool halffacet_found = false;

      //start finding from itself
      do{

        std::pair<MultiHandle, ShortDualHandle> opposit_it_temp(_opposite_facet, 0);
        //check whether it exists
        halffacet_found = std::binary_search(halffacetcell_pairs.begin(), halffacetcell_pairs.end(), opposit_it_temp, [](const std::pair<MultiHandle, ShortDualHandle>& a,
                                                                                                                         const std::pair<MultiHandle, ShortDualHandle>& b)
                                                                                                                                   {return a.first < b.first; });

        //return the iterator
        if (halffacet_found) {
          opposite_it = std::upper_bound(halffacetcell_pairs.begin(), halffacetcell_pairs.end(), opposit_it_temp, [](const std::pair<MultiHandle, ShortDualHandle>& a,
                                                                                                                     const std::pair<MultiHandle, ShortDualHandle>& b)
                                                                                                                               {return a.first < b.first; });
          opposite_it -= 1;

          //opposite exists only when another halffacet that belong to another cell has the same vertices
          opposite_exists = halffacet_found && (it.second.get_global_id() != opposite_it->second.get_global_id());
        }

        if (opposite_exists) break;

        _opposite_facet = next_order_hffacet(_opposite_facet);
      } while (_facet != _opposite_facet);

      if (opposite_exists){ //not a boundary
        if (_facet.front() > _facet.back()){ //otherwise opperate twice
          ShortDualHandle _opposite_facet2cell_handle = opposite_it->second;
          Flag _local_id = _facet2cell_handle.get_local_id();
          Flag _global_id = _facet2cell_handle.get_global_id();
          Flag _opposite_local_id = _opposite_facet2cell_handle.get_local_id();
          Flag _opposite_global_id = _opposite_facet2cell_handle.get_global_id();
          cell_neighbours[_global_id][_local_id] = _opposite_facet2cell_handle;
          cell_neighbours[_opposite_global_id][_opposite_local_id] = _facet2cell_handle;
        }
      }
      else{ // is a boundary
        ShortDualHandle _boundary_handle(0, id_boundary);
        Flag _local_id = _facet2cell_handle.get_local_id();
        Flag _global_id = _facet2cell_handle.get_global_id();
        cell_neighbours[_global_id][_local_id] = _boundary_handle;
        boundary2opposite_handles.push_back(_facet2cell_handle);
        id_boundary++;
      }
    }

    //release the memory of halffacet cell pairs list
    halffacetcell_pairs.clear();
    halffacetcell_pairs.shrink_to_fit();
  }

  void unstructuredReducedFvMesh::BuildGeometry(){
    Flag id_cells = 0;
    ///This loop builds all the facets and faces
    for (auto& c_it : cell_vertices){
      std::vector<Vector3> all_cell_vertices;
      for (auto v_it : c_it){
        all_cell_vertices.push_back(vertex_positions[v_it]);
      }
      Flag cell_type = cell_types[id_cells];
      auto cell_property = ElementPropertyTable.find(cell_type)->second; //to do: throw an exception when type value is invalid
      cell_volumes.push_back(cell_property.volume(all_cell_vertices));
      cell_centre_positions.push_back(cell_property.centroid(all_cell_vertices));
      for (auto local_halffacet : cell_property.facets){
        std::vector<Vector3> all_facet_vertices;
        Flag facet_type = local_halffacet.type_facet;
        auto facet_property = ElementPropertyTable.find(facet_type)->second;
        for (auto local_vertex_id : local_halffacet.vertices){
          all_facet_vertices.push_back(vertex_positions[c_it[local_vertex_id]]);
        }
        halffacet_areas.push_back(facet_property.volume(all_facet_vertices));
        halffacet_normal_directions.push_back(facet_property.normal(all_facet_vertices));
        halffacet_centre_positions.push_back(facet_property.centroid(all_facet_vertices));
      }
      id_cells++;
    }

  }


}
