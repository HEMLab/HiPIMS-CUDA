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
  \file gmsh_reader.h
  \brief Header file for gmsh file reader class

*/

#ifndef GMSH_READER_H
#define GMSH_READER_H

#include "mesh_reader.h"

namespace GC{

  class gmshReader : public meshReader{
    public:
      ///constructor
      gmshReader(const char* filename);
    private:  
      ///readin element connectivity table from file
      virtual void readin(const char* filename);
  };

}

#endif
