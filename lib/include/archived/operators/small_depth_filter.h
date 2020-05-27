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
\brief Header file for very small depth filter

\version 0.1
\author xilin xia

*/

#ifndef SMALL_DEPTH_FILTER_H
#define SMALL_DEPTH_FILTER_H

#include "mapped_field.h"

namespace GC{

  namespace fv{

    class smallDepthFilter{
      public:
        smallDepthFilter(Scalar _small_value):small_value(_small_value){};
        void operator() (fvScalarFieldOnCell& h, fvVectorFieldOnCell& hU);

      private:
        const Scalar small_value;
    };

  }


}

#endif

