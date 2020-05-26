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
 \file mesh_basic.h
 \brief Header file for basic mesh class

 \version 0.1
 \authot xilin xia
*/


#ifndef MESH_BASIC_H  //header file protector
#define MESH_BASIC_H

#include "mesh_reader.h"
#include "mesh_fv_entities.h"

namespace GC{

  ///This class implements the basic mesh data structure, i.e. the element connectivity table
  class basicMesh{
    public:
      ///default constructor - do nothing
      basicMesh() = default;
      ///constructor - initialize by reading from file
      basicMesh(const meshReader& mesh_reader);
      basicMesh(const meshReader&& mesh_reader);
      //--iterators for vertices-------------------------------------------------
      //positions
      typedef vectorProps::iterator vertex_positions_iterator;
      typedef vectorProps::const_iterator vertex_positions_const_iterator;
      virtual vertex_positions_iterator vertex_positions_begin();
      virtual vertex_positions_iterator vertex_positions_end();        
      //--end iterators for vertices---------------------------------------------
      

      //--iterators for cells----------------------------------------
      typedef CellVerticeSet::iterator cell_vertices_iterator;
      typedef CellVerticeSet::const_iterator cell_vertices_const_iterator;
      virtual cell_vertices_iterator cell_vertices_begin();
      virtual cell_vertices_iterator cell_vertices_end();
      //--end iterators for faces-------------------------------------

      ///Return the type of this mesh
      virtual const char* type_name();
    protected:
      void BuildElementTable(const meshReader& mesh_reader);
      void BuildElementTable(const meshReader&& mesh_reader);
    protected:
      //all data members
      //--Element connectivity table
      vectorProps vertex_positions; 
      CellVerticeSet cell_vertices;
      HandleSet cell_types;
    private:
      //--type
      const char* type;
  };

}

#endif  //end header file protector
