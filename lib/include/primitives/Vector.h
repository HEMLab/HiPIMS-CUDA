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
  \file Vector.h
  \brief Header file for Vector class

  \version 0.1
  \author xilin xia
*/

#ifndef VECTOR_H
#define VECTOR_H

#include "Scalar.h"
#include <cmath>
#include <iostream>
#include "hemi.h"

namespace GC{

  struct Vector1;
  struct Vector2;
  struct Vector3;


  struct Vector3{
    //class for three dimensional vector
      HEMI_DEV_CALLABLE_INLINE_MEMBER
      Vector3(Scalar x = 0, Scalar y = 0, Scalar z = 0)
        :x(x), y(y), z(z){}
      HEMI_DEV_CALLABLE_INLINE_MEMBER
      Vector3(const Vector2& a);
      HEMI_DEV_CALLABLE_INLINE_MEMBER
      Vector3(const Vector1& a);
      HEMI_DEV_CALLABLE_INLINE_MEMBER
      operator Scalar() const{return x;}
      Scalar x;
      Scalar y;
      Scalar z;
  };

  struct Vector2{
    //class for two dimensional vector
      HEMI_DEV_CALLABLE_INLINE_MEMBER
      Vector2(Scalar x = 0, Scalar y = 0):x(x), y(y){};
      HEMI_DEV_CALLABLE_INLINE_MEMBER
      Vector2(const Vector3& a);
      HEMI_DEV_CALLABLE_INLINE_MEMBER
      Vector2(const Vector1& a);
      HEMI_DEV_CALLABLE_INLINE_MEMBER
      operator Scalar() const{return x;}
      Scalar x;
      Scalar y;
  };

  struct Vector1{
    //class for on dimensional vector, i.e. a scalar
      HEMI_DEV_CALLABLE_INLINE_MEMBER
      Vector1(Scalar x = 0):x(x){};
      HEMI_DEV_CALLABLE_INLINE_MEMBER
      Vector1(const Vector3& a);
      HEMI_DEV_CALLABLE_INLINE_MEMBER
      Vector1(const Vector2& a);
      HEMI_DEV_CALLABLE_INLINE_MEMBER
      operator Scalar() const{return x;}
      Scalar x;
  };

  Vector3::Vector3(const Vector2& a):x(a.x), y(a.y), z(0){}
  Vector3::Vector3(const Vector1& a):x(a.x), y(0), z(0){}

  Vector2::Vector2(const Vector3& a):x(a.x), y(a.y){}
  Vector2::Vector2(const Vector1& a):x(a.x), y(0){}

  Vector1::Vector1(const Vector3& a):x(a.x){}
  Vector1::Vector1(const Vector2& a):x(a.x){}

  HEMI_DEV_CALLABLE_INLINE Vector2 operator+ (const Vector2& lhs, const Vector2& rhs){	//vector adding
    return Vector2(lhs.x + rhs.x, lhs.y + rhs.y);
  }

  HEMI_DEV_CALLABLE_INLINE Vector2& operator+= (Vector2& lhs, const Vector2& rhs){
    lhs.x += rhs.x;
    lhs.y += rhs.y;
    return lhs;
  }

  HEMI_DEV_CALLABLE_INLINE Vector2 operator- (const Vector2& lhs, const Vector2& rhs){	//vector reducing
    return Vector2(lhs.x - rhs.x, lhs.y - rhs.y);
  }
  
  HEMI_DEV_CALLABLE_INLINE Vector2 operator- (const Vector2& rhs){	//vector reducing
    return Vector2(- rhs.x, - rhs.y);
  }

  HEMI_DEV_CALLABLE_INLINE Vector2& operator-= (Vector2& lhs, const Vector2& rhs){
    lhs.x -= rhs.x;
    lhs.y -= rhs.y;
    return lhs;
  }

  HEMI_DEV_CALLABLE_INLINE Vector2 operator* (const Scalar& lhs, const Vector2& rhs){	//multiplying
    return Vector2(lhs*rhs.x, lhs*rhs.y);
  }

  HEMI_DEV_CALLABLE_INLINE Vector2 operator* (const Vector2& lhs, const Scalar& rhs){	//multiplying
    return Vector2(lhs.x*rhs, lhs.y*rhs);
  }

  HEMI_DEV_CALLABLE_INLINE Vector2& operator*= (Vector2& lhs, const Scalar& rhs){
    lhs.x *= rhs;
    lhs.y *= rhs;
    return lhs;
  }

  HEMI_DEV_CALLABLE_INLINE Vector2 operator/ (const Vector2& lhs, const Scalar& rhs){	//dividing
    return Vector2(lhs.x/rhs, lhs.y/rhs);
  }

  HEMI_DEV_CALLABLE_INLINE Vector2& operator/= (Vector2& lhs, const Scalar& rhs){
    lhs.x /= rhs;
    lhs.y /= rhs;
    return lhs;
  }

  HEMI_DEV_CALLABLE_INLINE Scalar dot (const Vector2& lhs, const Vector2& rhs){		//dot product
    return lhs.x*rhs.x + lhs.y*rhs.y;
  }

  HEMI_DEV_CALLABLE_INLINE Scalar norm (const Vector2& a){					//modulus
    return sqrt(dot(a,a));
  }

  HEMI_DEV_CALLABLE_INLINE Vector2 uni(const Vector2& a){					//unity vector
    if (norm(a) > 1e-15){
      return Vector2(a.x / norm(a), a.y / norm(a));
    }
    else{
      return Vector2(0.0, 0.0);
    }
  }

  HEMI_DEV_CALLABLE_INLINE Vector2 project (const Vector2& a, const Vector2& b){		//project a vector to another
    return dot(a, uni(b))*uni(b);
  }	

  HEMI_DEV_CALLABLE_INLINE Vector2 perpend (const Vector2& a){				//perpending vector
    return Vector2(-a.y, a.x);
  }

  inline std::istream& operator>> (std::istream& in, Vector2& a){
    in>>a.x>>a.y;
    return in;
  }

  inline std::ostream& operator<< (std::ostream& out, const Vector2& a){
    out<<a.x<<"  "<<a.y;
    return out;
  }


  HEMI_DEV_CALLABLE_INLINE Vector3 operator+ (const Vector3& lhs, const Vector3& rhs){	//vector adding
    return Vector3(lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z);
  }

  HEMI_DEV_CALLABLE_INLINE Vector3& operator+= (Vector3& lhs, const Vector3& rhs){
    lhs.x += rhs.x;
    lhs.y += rhs.y;
    lhs.z += rhs.z;
    return lhs;
  }

  HEMI_DEV_CALLABLE_INLINE Vector3 operator- (const Vector3& lhs, const Vector3& rhs){	//vector reducing
    return Vector3(lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z);
  }

  HEMI_DEV_CALLABLE_INLINE Vector3 operator- (const Vector3& rhs){	//vector reducing
    return Vector3( - rhs.x, - rhs.y, - rhs.z);
  }

  HEMI_DEV_CALLABLE_INLINE Vector3& operator-= (Vector3& lhs, const Vector3& rhs){
    lhs.x -= rhs.x;
    lhs.y -= rhs.y;
    lhs.z -= rhs.z;
    return lhs;
  }

  HEMI_DEV_CALLABLE_INLINE Vector3 operator* (const Scalar& lhs, const Vector3& rhs){	//multiplying
    return Vector3(lhs*rhs.x, lhs*rhs.y, lhs*rhs.z);
  }

  HEMI_DEV_CALLABLE_INLINE Vector3 operator* (const Vector3& lhs, const Scalar& rhs){	//multiplying
    return Vector3(lhs.x*rhs, lhs.y*rhs, lhs.z*rhs);
  }

  HEMI_DEV_CALLABLE_INLINE Vector3& operator*= (Vector3& lhs, const Scalar& rhs){
    lhs.x *= rhs;
    lhs.y *= rhs;
    lhs.z *= rhs;
    return lhs;
  }

  HEMI_DEV_CALLABLE_INLINE Vector3 operator/ (const Vector3& lhs, const Scalar& rhs){	//dividing
    return Vector3(lhs.x/rhs, lhs.y/rhs, lhs.z/rhs);
  }

  HEMI_DEV_CALLABLE_INLINE Vector3& operator/= (Vector3& lhs, const Scalar& rhs){
    lhs.x /= rhs;
    lhs.y /= rhs;
    lhs.z /= rhs;
    return lhs;
  }

  HEMI_DEV_CALLABLE_INLINE Scalar dot (const Vector3& lhs, const Vector3& rhs){		//dot product
    return lhs.x*rhs.x + lhs.y*rhs.y + lhs.z*rhs.z;
  }

  HEMI_DEV_CALLABLE_INLINE Scalar norm (const Vector3& a){					//modulus
    return sqrt(dot(a,a));
  }

  HEMI_DEV_CALLABLE_INLINE Vector3 uni (const Vector3& a){//unity vector
    if (norm(a) > 1e-15){
      return Vector3(a.x / norm(a), a.y / norm(a), a.z / norm(a));
    }
    else{
      return Vector2(0.0, 0.0);
    }
  }

  HEMI_DEV_CALLABLE_INLINE Vector3 project (const Vector3& a, const Vector3& b){		//project a vector to another
    return dot(a, uni(b))*uni(b);
  }	

  HEMI_DEV_CALLABLE_INLINE Vector3 perpend(const Vector3& a){				//perpending vector
    return Vector3(-a.y, a.x);
  }

  HEMI_DEV_CALLABLE_INLINE Vector3 cross (const Vector3& a, const Vector3& b){
   return Vector3(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
  } 

  inline std::istream& operator>> (std::istream& in, Vector3& a){
    in>>a.x>>a.y>>a.z;
    return in;
  }

  inline std::ostream& operator<< (std::ostream& out, const Vector3& a){
    out<<a.x<<"  "<<a.y<<"  "<<a.z;
    return out;
  }

  typedef Vector2 Vector;

}

#endif
