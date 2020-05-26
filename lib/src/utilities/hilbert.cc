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
 \file hilbert.cc
 \brief Source file for generatiing hilbert key

 \version 0.1
 \authot xilin xia
*/

#include "hilbert.h"

namespace GC{

  ///convert (x,y) to hilbert key
  unsigned int xy2hilbertkey (unsigned int n, unsigned int x, unsigned int y){
    unsigned int rx, ry, s, d=0;
    for (s=n/2; s>0; s/=2) {
      rx = (x & s) > 0;
      ry = (y & s) > 0;
      d += s * s * ((3 * rx) ^ ry);
      hilbert_rot(s, &x, &y, rx, ry);
    }
    return d;
  }
   
  //convert hilbert key to (x,y)
  void hilbertkey2xy(unsigned int n, unsigned int d, unsigned int *x, unsigned int *y){
    unsigned int rx, ry, s, t=d;
    *x = *y = 0;
    for (s=1; s<n; s*=2){
      rx = 1 & (t/2);
      ry = 1 & (t ^ rx);
      hilbert_rot(s, x, y, rx, ry);
      *x += s * rx;
      *y += s * ry;
      t /= 4;
    }
  }
   
  //rotate/flip a quadrant appropriately
  void hilbert_rot(unsigned int n, unsigned int *x, unsigned int *y, unsigned int rx, unsigned int ry){
    if (ry == 0){
      if (rx == 1) {
        *x = n-1 - *x;
        *y = n-1 - *y;
      }

      //Swap x and y
      unsigned int t  = *x;
      *x = *y;
      *y = t;
    }
  }

}
