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
  \file geometry_func.cc 
  \brief Source file for geometry functions

  \version 0.1
  \author xilin xia

*/

#include <vector>
#include "Scalar.h"
#include "Vector.h"

namespace GC{

  ///This function calculates the volume of a point, it returns 1
  Scalar point_volume(const std::vector<Vector3>& vertices){
    return 1.0;
  }

  ///This function calculates the normal direction of a point, it returns(0.0, 0.0, 0.0)
  Vector3 point_normal(const std::vector<Vector3>& vertices){
    return Vector3(0.0,0.0,0.0);
  }

  ///This function calculates the centre of a point
  Vector3 point_centre(const std::vector<Vector3>& vertices){
    return vertices[0];
  }

  ///This function calculates the volume of a segment, it returns its length
  Scalar segment_volume(const std::vector<Vector3>& vertices){
    return norm(vertices[0] - vertices[1]);
  }

  ///This function calculates the normal direction of a segment
  Vector3 segment_normal(const std::vector<Vector3>& vertices){
    return uni(perpend(vertices[0] - vertices[1]));
  }

  ///This function calculates the centre of a segment
  Vector3 segment_centre(const std::vector<Vector3>& vertices){
    Scalar a = 0.5;
    return a*(vertices[0]+vertices[1]);
  }

  ///This function calculates the volume of a triangle, it returns its area
  Scalar triangle_volume(const std::vector<Vector3>& vertices){
    Vector3 ab = vertices[1] - vertices[0];
    Vector3 bc = vertices[2] - vertices[1];
    return 0.5*norm(cross(ab,bc));
  }

  ///This function calculates the normal direction of a triangle
  Vector3 triangle_normal(const std::vector<Vector3>& vertices){
    Vector3 ab = vertices[1] - vertices[0];
    Vector3 bc = vertices[2] - vertices[1];
    return uni(cross(ab,bc));
  }

  ///This function calculates the centre of a triangle
  Vector3 triangle_centre(const std::vector<Vector3>& vertices){
    Scalar a = 3.0;
    return (vertices[0] + vertices[1] + vertices[2])/a;
  }

  ///This function calculates the volume of a quadrilateral, it return its area
  Scalar quadrilateral_volume(const std::vector<Vector3>& vertices){
    Vector3 ab = vertices[1] - vertices[0];
    Vector3 bc = vertices[2] - vertices[1];
    Vector3 cd = vertices[3] - vertices[2];
    Vector3 da = vertices[0] - vertices[3];
    return 0.5*norm(cross(ab,bc)) + 0.5*norm(cross(cd, da));
  }

  ///This function calculates the normal direction of a quadrilateral
  Vector3 quadrilateral_normal(const std::vector<Vector3>& vertices){
    Vector3 ab = vertices[1] - vertices[0];
    Vector3 bc = vertices[2] - vertices[1];
    return uni(cross(ab,bc));
  }

  ///This function calculates the centre of a quadrilateral 
  Vector3 quadilateral_centre(const std::vector<Vector3>& vertices){
    Vector3 ab = vertices[1] - vertices[0];
    Vector3 bc = vertices[2] - vertices[1];
    Vector3 cd = vertices[3] - vertices[2];
    Vector3 da = vertices[0] - vertices[3];
    Scalar A_abc = 0.5*norm(cross(ab,bc));
    Scalar A_acd = 0.5*norm(cross(cd,da));
    Scalar a = 3.0;
    Vector3 C_abc = (vertices[0] + vertices[1] + vertices[2])/a;
    Vector3 C_acd = (vertices[0] + vertices[2] + vertices[3])/a;
    return (C_abc*A_abc + C_acd*A_acd)/(A_abc + A_acd);
  }

}
