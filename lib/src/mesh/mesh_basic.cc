// ======================================================================================
// Name                :    GeoClasses : Generic Geophysical Flow Modelling Framework
// Description         :    This code pack provides a generic framework for developing 
//                          Geophysical CFD software.
// ======================================================================================
// Version             :    0.2 
// Author              :    Xilin Xia (PhD candidate in Newcastle University)
// Create Time         :    2014/10/04
// Update Time         :    2015/10/09
// ======================================================================================
// Copyright @ Xilin Xia 2015 . All rights reserved.
// ======================================================================================

/*!
 \file mesh_basic.cc
 \brief Source file for basic mesh class

*/

#include "mesh_basic.h"
#include <iostream>
#include <utility>

namespace GC{

  //--initialization----------------------------------------------------------------
  basicMesh::basicMesh(const meshReader& mesh_reader) : type("basicMesh"){
    BuildElementTable(mesh_reader);
  }

  basicMesh::basicMesh(const meshReader&& mesh_reader) : type("basicMesh"){
    BuildElementTable(std::move(mesh_reader));
  }

  //--iterators---------------------------------------------------------------------
  basicMesh::vertex_positions_iterator
  basicMesh::vertex_positions_begin(){
    return vertex_positions.begin();
  }

  basicMesh::vertex_positions_iterator
  basicMesh::vertex_positions_end(){
    return vertex_positions.end();
  }

  basicMesh::cell_vertices_iterator
  basicMesh::cell_vertices_begin(){
    return cell_vertices.begin();
  }

  basicMesh::cell_vertices_iterator
  basicMesh::cell_vertices_end(){
    return cell_vertices_end();
  }

  void basicMesh::BuildElementTable(const meshReader& mesh_reader){
    vertex_positions = mesh_reader.vertex_positions;
    cell_vertices = mesh_reader.element_vertices;
    cell_types = mesh_reader.element_types;
  }

  void basicMesh::BuildElementTable(const meshReader&& mesh_reader){
    vertex_positions = std::move(mesh_reader.vertex_positions);
    cell_vertices = std::move(mesh_reader.element_vertices);
    cell_types = std::move(mesh_reader.element_types);
  }

  const char* basicMesh::type_name(){
    return type;
  }

}
