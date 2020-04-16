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
