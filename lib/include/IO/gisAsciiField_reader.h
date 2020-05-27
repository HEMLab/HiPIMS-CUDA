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
  \file field_reader.h
  \brief Header file for basic field file reader class

*/

#include <map>
#include "Flag.h"
#include "Scalar.h"
#include "Vector.h"
#include "field_reader.h"

#ifndef GIS_ASCII_FIELD_READER_H
#define GIS_ASCII_FIELD_READER_H

namespace GC{

  class gisAsciiFieldReader : public fieldReader{
  public:
    gisAsciiFieldReader(const char* filename);
  protected:
    virtual void readin(const char* filename);
  };

}

#endif

