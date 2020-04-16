// ======================================================================================
// Name                :    GeoClasses : Generic Geophysical Flow Modelling Framework
// Description         :    This code pack provides a generic framework for developing 
//                          Geophysical CFD software.
// ======================================================================================
// Version             :    0.1 
// Author              :    Xilin Xia (PhD candidate in Newcastle University)
// Create Time         :    2014/10/04
// Update Time         :    2015/10/29
// ======================================================================================
// Copyright @ Xilin Xia 2015 . All rights reserved.
// ======================================================================================

/*!
  \flie gradient.h
  \brief Header file for gradient of an arbitary variable phi

  \version 0.1
  \author xilin xia
*/

#ifndef GRADIENT_H
#define GRADIENT_H

#include "matrix.h"
#include "mapped_field.h"
#include "Scalar.h"
#include "Vector.h"
#include "Tensor.h"
#include "Flag.h"

namespace GC{

  namespace fv{
    
    ///This is a function class for calculating gradient of an arbitary variable phi
    class gradient{
      public:
        fvMappedField<Vector, on_cell>& operator() (fvMappedField<Scalar, on_cell>& phi);
        fvMappedField<Tensor, on_cell>& operator() (fvMappedField<Vector, on_cell>& phi);
      private:
        fvMappedField<Vector, on_cell> vector_buffer;
        fvMappedField<Tensor, on_cell> tensor_buffer;
    };

    void Gradient(fvMappedField<Scalar, on_cell>& phi, fvMappedField<Vector, on_cell>& phi_grad);
    void Gradient(fvMappedField<Vector, on_cell>& phi, fvMappedField<Tensor, on_cell>& phi_grad);

  }

}


#endif
