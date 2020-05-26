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
  \flie gradient.h
  \brief Header file for gradient of an arbitary variable phi

  \version 0.1
  \author xilin xia
*/

#ifndef CALC_VELOCITY_H
#define CALC_VELOCITY_H

#include "mapped_field.h"

namespace GC{

  namespace fv{
    
    ///This is a function class for calculating velocity from discharge
    class clacVelocity{
      public:
        fvMappedField<Vector, on_cell>& operator() (fvMappedField<Scalar, on_cell>& h, fvMappedField<Vector, on_cell>& hU);
      private:
        fvMappedField<Vector, on_cell> vector_buffer;
    };

  }

}


#endif
