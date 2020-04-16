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

