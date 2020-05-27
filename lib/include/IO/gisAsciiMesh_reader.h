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
  \file gisAscii_reader.h
  \brief Header file for gis ascii file reader class

*/

#ifndef GIS_ASCII_MESH_READER_H
#define GIS_ASCII_MESH_READER_H

#include "mesh_reader.h"

namespace GC{

  class gisAsciiMeshReader : public meshReader{
    public:
      ///constructor
      gisAsciiMeshReader(const char* filename);
    private:
      ///readin element connectivity table from file
      virtual void readin(const char* filename);
  }; 

}


#endif

