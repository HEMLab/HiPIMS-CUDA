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
  \file mesh_reader.h
  \brief Header file for basic mesh file reader class

*/

#ifndef MESHREADER_H
#define MESHREADER_H

#include "mesh_fv_entities.h"

namespace GC{

  class meshReader{
    public:
      ///constructor
      meshReader() = default;
      //element connectivity table
      vectorProps vertex_positions;
      HandleSet element_types;
      CellVerticeSet element_vertices;
  };

}
#endif
