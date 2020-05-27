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
  \file element_table.h
  \brief Header file for element property table

  \version 0.1
  \author xilin xia

*/

#ifndef ELEMENT_TABLE_H
#define ELEMENT_TABLE_H


#include <cstdint>
#include <vector>
#include <map>
#include <iostream>
#include "Flag.h"
#include "Vector.h"
#include "Scalar.h"

namespace GC{

  Scalar point_volume(const std::vector<Vector3>& vertices);
  Vector3 point_normal(const std::vector<Vector3>& vertices);
  Vector3 point_centre(const std::vector<Vector3>& vertices);

  Scalar segment_volume(const std::vector<Vector3>& vertices);
  Vector3 segment_normal(const std::vector<Vector3>& vertices);
  Vector3 segment_centre(const std::vector<Vector3>& vertices);

  Scalar triangle_volume(const std::vector<Vector3>& vertices);
  Vector3 triangle_normal(const std::vector<Vector3>& vertices);
  Vector3 triangle_centre(const std::vector<Vector3>& vertices);

  Scalar quadrilateral_volume(const std::vector<Vector3>& vertices);
  Vector3 quadrilateral_normal(const std::vector<Vector3>& vertices);
  Vector3 quadilateral_centre(const std::vector<Vector3>& vertices);

  ///This class defines an element property table



  struct ElementProperty{
    typedef std::vector<short> Vertices;

    struct Facet{
      unsigned short type_facet;
      Vertices vertices;
    };

    unsigned short num_facets;
    std::vector<Facet> facets;
    Scalar(*volume) (const std::vector<Vector3>& vertices);
    Vector3(*normal) (const std::vector<Vector3>& vertices);
    Vector3(*centroid) (const std::vector<Vector3>& vertices);
  };


  extern const std::map<Flag, ElementProperty> ElementPropertyTable;

}

#endif
