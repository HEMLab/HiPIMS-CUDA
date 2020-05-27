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
 \file mesh_fv.cc
 \brief Source file for finite volume mesh class

*/

#include "mesh_fv.h"
#include <iostream>
#include <map>
#include <vector>
#include "hilbert.h"
#include "element_table.h"

namespace GC{

  //--initialization----------------------------------------------------------------

  unstructuredFvMesh::unstructuredFvMesh(const meshReader& mesh_reader):type("unstructuredFvMesh") {
    BuildElementTable(mesh_reader);
    BuildEntities();
  }

  ///This function builds the comprehensive mesh data structure
  void unstructuredFvMesh::BuildEntities(){
      BuildTopDown();
      BuildBottomUp();
      BuildGeometry();
  }

  void unstructuredFvMesh::BuildTopDown(){
    vertex2halffacet_handles.resize(vertex_positions.size());
    vertex2cell_handles.resize(vertex_positions.size());
    cell_halffacets.resize(cell_vertices.size());
    Flag id_halffacets = 0;
    Flag id_cells = 0;
    ///This loop builds all the facets and faces
    for(auto& c_it:cell_vertices){
      for (auto v_it : c_it){
        vertex2cell_handles[v_it].push_back(id_cells);
      }
      Flag cell_type = cell_types[id_cells];
      auto cell_property = ElementPropertyTable.find(cell_type)->second; //to do: throw an exception when type value is invalid
      Flag id_local_halffacet = 0;
      for(auto local_halffacet : cell_property.facets){
        MultiHandle halffacet_;
        for(auto local_vertex_id : local_halffacet.vertices){
          halffacet_.push_back(c_it[local_vertex_id]);
        }
        halffacet2cell_map[halffacet_] = ShortDualHandle(id_local_halffacet+1,id_cells); // +1 indicates that it is not a boundary half facet 
        vertex2halffacet_handles[halffacet_.front()] = ShortDualHandle(id_local_halffacet+1,id_cells);
        halffacet2cell_handles.push_back(ShortDualHandle(id_local_halffacet+1,id_cells));
        id_local_halffacet++;
        cell_halffacets[id_cells].push_back(id_halffacets);
        id_halffacets++;
      }
      id_cells++;
    }

  }

  ///This function builds other datas
  void unstructuredFvMesh::BuildBottomUp(){
    Flag id_boundary = 0;
    Flag id_cells = 0;
    cell_neighbours.resize(cell_vertices.size());  //request memory for cell_neighbours
    for(auto& it:cell_neighbours){
      it.resize(cell_halffacets[id_cells++].size());
    }

    for(auto& it : halffacet2cell_map){
      MultiHandle _facet = it.first;
      ShortDualHandle _facet2cell_handle = it.second;
      MultiHandle _opposite_facet = next_order_hffacet(_facet);
      std::map<MultiHandle, ShortDualHandle>::iterator opposite_it;
      bool opposite_exists = false;
      while(_facet != _opposite_facet){
        opposite_exists = (opposite_it = halffacet2cell_map.find(_opposite_facet)) != halffacet2cell_map.end();
        if(opposite_exists) break;
        _opposite_facet = next_order_hffacet(_opposite_facet);
      }
      if(opposite_exists){ //not a boundary
        if(_facet.front() > _facet.back()){ //otherwise opperate twice
          ShortDualHandle _opposite_facet2cell_handle = opposite_it->second;
          Flag _local_id = _facet2cell_handle.get_local_id();
          Flag _global_id = _facet2cell_handle.get_global_id();
          Flag _opposite_local_id = _opposite_facet2cell_handle.get_local_id();
          Flag _opposite_global_id = _opposite_facet2cell_handle.get_global_id();
          cell_neighbours[_global_id][_local_id] = _opposite_facet2cell_handle;
          cell_neighbours[_opposite_global_id][_opposite_local_id] = _facet2cell_handle;
        }
      }else{ // is a boundary
        ShortDualHandle _boundary_handle(0, id_boundary);
        Flag _local_id = _facet2cell_handle.get_local_id();
        Flag _global_id = _facet2cell_handle.get_global_id();
        cell_neighbours[_global_id][_local_id] = _boundary_handle;
        vertex2halffacet_handles[_opposite_facet.front()] = _boundary_handle;
        boundary2opposite_handles.push_back(_facet2cell_handle);
        id_boundary++;
      }
    }
    halffacet2cell_map.clear();
  }

  void unstructuredFvMesh::BuildGeometry(){
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

  MultiHandle unstructuredFvMesh::next_order_hffacet(const MultiHandle& hffacet){
      MultiHandle new_hffacet;
      new_hffacet.reserve(hffacet.size());
      new_hffacet.push_back(hffacet.back());
      for(auto it = hffacet.begin(); it != hffacet.end() - 1; ++it){
        new_hffacet.push_back(*it);
      }
      return new_hffacet;
  }



  //--iterators-------------------------------------------------------------
  //--iterators for vertices-------------------------------------------------
  size_t unstructuredFvMesh::vertex_size(){
    return vertex_positions.size();
  }

  //cells incident of
  unstructuredFvMesh::vertex2cells_iterator
  unstructuredFvMesh::vertex2cells_begin(){
    return vertex2cell_handles.begin();
  }

  unstructuredFvMesh::vertex2cells_iterator
  unstructuredFvMesh::vertex2cells_end(){
    return vertex2cell_handles.end();
  }
  //out going half facets 
  unstructuredFvMesh::vertex2halffacet_iterator
  unstructuredFvMesh::vertex2halffacet_begin(){
    return vertex2halffacet_handles.begin();
  }

  unstructuredFvMesh::vertex2halffacet_iterator
  unstructuredFvMesh::vertex2halffacet_end(){
    return vertex2halffacet_handles.end();
  }

  //--end iterators for vertices---------------------------------------------
      
  //--iterators for half facets-------------------------------------
  size_t unstructuredFvMesh::halffacet_size(){
    return halffacet_centre_positions.size();
  }

  //centre positions
  unstructuredFvMesh::halffacet_centre_positions_iterator
  unstructuredFvMesh::halffacet_centre_positions_begin(){
    return halffacet_centre_positions.begin();
  }

  unstructuredFvMesh::halffacet_centre_positions_iterator
  unstructuredFvMesh::halffacet_centre_positions_end(){
    return halffacet_centre_positions.end();
  }

  //areas
  unstructuredFvMesh::halffacet_areas_iterator 
  unstructuredFvMesh::halffacet_areas_begin(){
    return halffacet_areas.begin();
  }

  unstructuredFvMesh::halffacet_areas_iterator 
  unstructuredFvMesh::halffacet_areas_end(){
    return halffacet_areas.end();
  }

  //normal directions
  unstructuredFvMesh::halffacet_normal_directions_iterator 
  unstructuredFvMesh::halffacet_normal_directions_begin(){
    return halffacet_normal_directions.begin();
  }

  unstructuredFvMesh::halffacet_normal_directions_iterator 
  unstructuredFvMesh::halffacet_normal_directions_end(){
    return halffacet_normal_directions.end();
  }

  //face belongs to
  unstructuredFvMesh::halffacet2cell_iterator 
  unstructuredFvMesh::halffacet2cell_begin(){
    return halffacet2cell_handles.begin();
  }

  unstructuredFvMesh::halffacet2cell_iterator 
  unstructuredFvMesh::halffacet2cell_end(){
    return halffacet2cell_handles.end();
  }

  //--end iterators for half facets---------------------------------

  //--iterators for cells------------------------------------------------
  size_t unstructuredFvMesh::cell_size(){
    return cell_centre_positions.size();
  }

  //positions
  unstructuredFvMesh::cell_centre_positions_iterator
  unstructuredFvMesh::cell_centre_positions_begin(){
    return cell_centre_positions.begin();
  }

  unstructuredFvMesh::cell_centre_positions_iterator
  unstructuredFvMesh::cell_centre_positions_end(){
    return cell_centre_positions.end();
  }

  //volumes
  unstructuredFvMesh::cell_volumes_iterator
  unstructuredFvMesh::cell_volumes_begin(){
    return cell_volumes.begin();
  }

  unstructuredFvMesh::cell_volumes_iterator
  unstructuredFvMesh::cell_volumes_end(){
    return cell_volumes.end();
  }

  //neighbours
  unstructuredFvMesh::cell_neighbours_iterator
  unstructuredFvMesh::cell_neighbours_begin(){
    return cell_neighbours.begin();
  }

  unstructuredFvMesh::cell_neighbours_iterator
  unstructuredFvMesh::cell_neighbours_end(){
    return cell_neighbours.end();
  }

  //half facets
  unstructuredFvMesh::cell_halffacets_iterator
  unstructuredFvMesh::cell_halffacets_begin(){
    return cell_halffacets.begin();
  }

  unstructuredFvMesh::cell_halffacets_iterator
  unstructuredFvMesh::cell_halffacets_end(){
    return cell_halffacets.end();
  }

  const char* unstructuredFvMesh::type_name(){
    return type;
  }

  //--end iterator for cells

  //--iterator for boundary---------------------------------
  size_t unstructuredFvMesh::boundary_size(){
    return boundary2opposite_handles.size();
  }

  unstructuredFvMesh::boundary2opposite_iterator
  unstructuredFvMesh::boundary2opposite_begin(){
    return boundary2opposite_handles.begin();
  }


  unstructuredFvMesh::boundary2opposite_iterator
  unstructuredFvMesh::boundary2opposite_end(){
    return boundary2opposite_handles.end();
  }

}
