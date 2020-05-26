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
  \file Tensor.h
  \brief Header file for Tensor class

  \version 0.1
  \author xilin xia
*/


#ifndef TENSOR_H
#define TENSOR_H

#include "Scalar.h"
#include "Vector.h"
#include <cmath>
#include <iostream>

namespace GC{

  struct Tensor2;

  struct Tensor2{
    //class for two dimensional second order tensor
    HEMI_DEV_CALLABLE_INLINE_MEMBER
    Tensor2(Scalar _xx = 0, Scalar _xy = 0, Scalar _yx = 0, Scalar _yy = 0):xx(_xx), xy(_xy), yx(_yx), yy(_yy){}
    Scalar xx;
    Scalar xy;
    Scalar yx;
    Scalar yy;

  };


  HEMI_DEV_CALLABLE_INLINE Tensor2 operator+ (const Tensor2& lhs, const Tensor2& rhs){	//vector adding
    return Tensor2(lhs.xx + rhs.xx, lhs.xy + rhs.xy, lhs.yx + rhs.yx, lhs.yy + rhs.yy);
  }

  HEMI_DEV_CALLABLE_INLINE Tensor2& operator+= (Tensor2& lhs, const Tensor2& rhs){
    lhs.xx += rhs.xx;
    lhs.xy += rhs.xy;
    lhs.yx += rhs.yx;
    lhs.yy += rhs.yy;
    return lhs;
  }

  HEMI_DEV_CALLABLE_INLINE Tensor2 operator- (const Tensor2& lhs, const Tensor2& rhs){	//vector reducing
    return Tensor2(lhs.xx - rhs.xx, lhs.xy - rhs.xy, lhs.yx - rhs.yx, lhs.yy - rhs.yy);

  }
  
  HEMI_DEV_CALLABLE_INLINE Tensor2 operator- (const Tensor2& rhs){	//vector reducing
    return Tensor2(- rhs.xx, - rhs.xy, - rhs.yx, - rhs.yy);
  }

  HEMI_DEV_CALLABLE_INLINE Tensor2& operator-= (Tensor2& lhs, const Tensor2& rhs){
    lhs.xx -= rhs.xx;
    lhs.xy -= rhs.xy;
    lhs.yx -= rhs.yx;
    lhs.yy -= rhs.yy;
    return lhs;
  }

  HEMI_DEV_CALLABLE_INLINE Tensor2 operator* (const Scalar& lhs, const Tensor2& rhs){	//multiplying
    return Tensor2(lhs*rhs.xx, lhs*rhs.xy, lhs*rhs.yx, lhs*rhs.yy);
  }

  HEMI_DEV_CALLABLE_INLINE Tensor2 operator* (const Tensor2& lhs, const Scalar& rhs){	//multiplying
    return Tensor2(lhs.xx*rhs, lhs.xy*rhs, lhs.yx*rhs, lhs.yy*rhs);
  }

  HEMI_DEV_CALLABLE_INLINE Tensor2& operator*= (Tensor2& lhs, const Scalar& rhs){
    lhs.xx *= rhs;
    lhs.xy *= rhs;
    lhs.yx *= rhs;
    lhs.yy *= rhs;
    return lhs;
  }

  HEMI_DEV_CALLABLE_INLINE Tensor2 operator/ (const Tensor2& lhs, const Scalar& rhs){	//dividing
    return Tensor2(lhs.xx/rhs, lhs.xy/rhs, lhs.yx/rhs, lhs.yy/rhs);
  }

  HEMI_DEV_CALLABLE_INLINE Tensor2& operator/= (Tensor2& lhs, const Scalar& rhs){
    lhs.xx /= rhs;
    lhs.xy /= rhs;
    lhs.yx /= rhs;
    lhs.yy /= rhs;
    return lhs;
  }

  HEMI_DEV_CALLABLE_INLINE Vector2 dot (const Tensor2& lhs, const Vector2& rhs){		//
    return Vector2(lhs.xx*rhs.x + lhs.xy*rhs.y, lhs.yx*rhs.x + lhs.yy*rhs.y);
  }

  HEMI_DEV_CALLABLE_INLINE Tensor2 product(const Vector2& lhs, const Vector2& rhs){  //Tensor product
    return Tensor2(lhs.x*rhs.x, lhs.x*rhs.y, lhs.y*rhs.x, lhs.y*rhs.y);
  }

  HEMI_DEV_CALLABLE_INLINE Tensor2 inverse(const Tensor2& phi){  //Tensor product
    Scalar det = phi.xx*phi.yy - phi.xy*phi.yx;
    return Tensor2(phi.yy/det, -phi.xy/det,-phi.yx/det, phi.xx/det);
  }

  inline std::istream& operator>> (std::istream& in, Tensor2& a){
    in>>a.xx>>a.xy>>a.yx>>a.yy;
    return in;
  }

  inline std::ostream& operator<< (std::ostream& out, const Tensor2& a){
    out<<a.xx<<"  "<<a.xy<<"  "<<a.yx<<"  "<<a.yy;
    return out;
  }

  typedef Tensor2 Tensor;


}//--end of namespace GC

#endif
