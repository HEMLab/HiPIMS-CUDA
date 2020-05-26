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
 \file hilbert.h
 \brief Header file for generatiing hilbert key

 \version 0.1
 \authot xilin xia
*/

#ifndef HILBERT_H
#define HILBERT_H

namespace GC{

  void hilbert_rot(unsigned int n, unsigned int *x, unsigned int *y, unsigned int rx, unsigned int ry);
  unsigned int xy2hilbertkey (unsigned int n, unsigned int x, unsigned int y); 
  void hilbertkey2xy(unsigned int n, unsigned int d, unsigned int *x, unsigned int *y);

}

#endif
