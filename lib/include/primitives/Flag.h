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

// header file for flag class

#ifndef FLAG_H
#define FLAG_H

#include <iostream>
#include "hemi.h"

namespace GC{

  typedef unsigned int Flag;

  class DualFlag{
    public:
      unsigned int x;
      unsigned int y;
  };


  ///This class uses a single unsigned int to store two values: x - short index, y - long index
  class ShortDualFlag{
    protected:
      unsigned int mask;
    public:
      HEMI_DEV_CALLABLE_INLINE_MEMBER 
      ShortDualFlag(unsigned int x = 0, unsigned int y = 0):mask(0){
        setx(x);
        sety(y);
      }

      HEMI_DEV_CALLABLE_INLINE_MEMBER 
      void setx(unsigned int x){
        mask = (mask & (~7)) | x;
      }

      HEMI_DEV_CALLABLE_INLINE_MEMBER 
      unsigned int getx() const{
        return mask & 7;
      }

      HEMI_DEV_CALLABLE_INLINE_MEMBER
      void sety(unsigned int y){
        mask = (mask & 7) | ( y << 3);
      }

      HEMI_DEV_CALLABLE_INLINE_MEMBER 
      unsigned int gety() const{
        return mask >> 3;	
      }
  };

  ///This class uses a single unsigned int to store three values: x - 4bits y - 12bits z - 16bits
  class ShortTripleFlag{
    protected:
      unsigned int mask;
    public:
      HEMI_DEV_CALLABLE_INLINE_MEMBER 
      ShortTripleFlag(unsigned int x = 0, unsigned int y = 0, unsigned int z = 0):mask(0){
        setx(x);
        sety(y);
        setz(z);
      }

      HEMI_DEV_CALLABLE_INLINE_MEMBER 
      void setx(unsigned int x){
        mask = (mask & (~15)) | x;
      }

      HEMI_DEV_CALLABLE_INLINE_MEMBER 
      void sety(unsigned int y){
        mask = (mask & (~(4095<<4))) | (y << 4);
      }

      HEMI_DEV_CALLABLE_INLINE_MEMBER 
      void setz(unsigned int z){
        mask = (mask & 65535) | (z << 16);
      }

      HEMI_DEV_CALLABLE_INLINE_MEMBER 
      unsigned int getx() const{
        return mask & 15;
      }

      HEMI_DEV_CALLABLE_INLINE_MEMBER 
      unsigned int gety() const{
        return (mask >> 4) & 4095;
      }

      HEMI_DEV_CALLABLE_INLINE_MEMBER 
      unsigned int getz() const{
        return mask >> 16;
      }
  };

  inline std::ostream& operator<< (std::ostream& out, const ShortDualFlag& a){
    out<<a.getx()<<"  "<<a.gety();
    return out;
  }

  inline std::ostream& operator<< (std::ostream& out, const ShortTripleFlag& a){
    out<<a.getx()<<"  "<<a.gety()<<"  "<<a.getz();
    return out;
  }

}

#endif
