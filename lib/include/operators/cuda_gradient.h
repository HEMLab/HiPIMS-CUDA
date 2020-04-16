// ======================================================================================
// Name                :    GeoClasses : Generic Geophysical Flow Modelling Framework
// Description         :    This code pack provides a generic framework for developing 
//                          Geophysical CFD software.
// ======================================================================================
// Version             :    0.1 
// Author              :    Xilin Xia (PhD candidate in Newcastle University)
// Create Time         :    2014/10/04
// Update Time         :    2015/10/15
// ======================================================================================
// Copyright @ Xilin Xia 2015 . All rights reserved.
// ======================================================================================

/*!
  \flie cuda_gradient.h
  \brief Header file for gradient of an arbitary variable phi by CUDA

  \version 0.1
  \author xilin xia
*/

#ifndef CUDA_GRADIENT_H
#define CUDA_GRADIENT_H

#include "cuda_mapped_field.h"
#include "Scalar.h"
#include "Vector.h"
#include "Tensor.h"
#include "Flag.h"

namespace GC{

  namespace fv{

    void cuGradient(cuFvMappedField<Scalar, on_cell>& phi, cuFvMappedField<Vector, on_cell>& phi_gradient);
    void cuGradient(cuFvMappedField<Vector, on_cell>& phi, cuFvMappedField<Tensor, on_cell>& phi_gradient);
    void cuLimitedGradientCartesian(cuFvMappedField<Scalar, on_cell>& phi, cuFvMappedField<Vector, on_cell>& phi_gradient);
    void cuLimitedGradientCartesian(cuFvMappedField<Vector, on_cell>& phi, cuFvMappedField<Tensor, on_cell>& phi_gradient);

  }


}//---end of namespace GC--------------

#endif
