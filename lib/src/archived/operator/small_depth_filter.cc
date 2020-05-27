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
\file friction.h
\brief Source file for very small depth filter

\version 0.1
\author xilin xia

*/

#include "small_depth_filter.h"

namespace GC{

  namespace fv{

    void smallDepthFilter::operator()(fvScalarFieldOnCell& h, fvVectorFieldOnCell& hU){
      auto h_begin = h.data_begin();
      auto hU_begin = hU.data_begin();

      for(auto h_iter = h_begin; h_iter < h.data_end(); ++h_iter){
        if(*h_iter < small_value){
          *h_iter = 0.0;
          auto id = h_iter - h_begin;
          *(hU_begin + id) = 0.0;
        }
      }
      
    } 



  }

}
