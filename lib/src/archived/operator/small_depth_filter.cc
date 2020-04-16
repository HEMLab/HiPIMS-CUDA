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
