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
#include <vector>
#include "Flag.h"
#include "Scalar.h"
#include "Vector.h"

#ifndef FIELD_READER_H
#define FIELD_READER_H

namespace GC{

  class fieldReader{
  public:
    fieldReader() = default;
    fieldReader(const char* filename);
    std::map<Flag, Vector3> data;
    std::map<Flag, ShortTripleFlag> boundary_type;
  protected:
    virtual void readin_field(const char* filename);
  };

  class completeFieldReader :public fieldReader{
  public:
    std::vector< std::vector<Scalar> > time_series;
    std::vector< std::vector<Vector3> > boundary_source;
    std::vector< std::vector<Scalar> > data_time_series;
    std::vector< std::vector<Vector3> > data_source;
    std::map<Flag, Flag> region_mask;
    completeFieldReader(const char* path, const char* field_name);
  protected:
    virtual bool readin_boundary_source(const char* filename, const unsigned int cnt);
    virtual bool readin_data_source(const char* filename, const unsigned int cnt);
    virtual bool readin_region_mask(const char* filename);
  };
}

#endif
