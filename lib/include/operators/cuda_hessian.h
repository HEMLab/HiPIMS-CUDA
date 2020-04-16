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
  \flie cuda_hessian.h
  \brief Header file for Hessian of an arbitary variable phi by CUDA

  \version 0.1
  \author xilin xia
*/

#ifndef CUDA_HESSIAN_H
#define CUDA_HESSIAN_H

#include "cuda_mapped_field.h"
#include "Scalar.h"
#include "Vector.h"
#include "Tensor.h"
#include "Flag.h"

namespace GC{

  namespace fv{
    
    ///function for calculating Hessian matrix on a 2D Cartesian grid
    void cuHessianCartesian2D(cuFvMappedField<Scalar, on_cell>& phi, cuFvMappedField<Tensor2, on_cell>& phi_hessian);

  }

}

#endif

