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
  \flie cuda_transport_NSWEs.h
  \brief Header file for advection incorporating transport of non hydrostatic shallow water equations

  \version 0.1
  \author xilin xia
*/

#ifndef CUDA_TRANSPORT_NSWES_H
#define CUDA_TRANSPORT_NSWES_H

#include "cuda_mapped_field.h"
#include "Scalar.h"
#include "Vector.h"
#include "Tensor.h"
#include "Flag.h"

namespace GC{

  namespace fv{

    ///A fast 1st order surface reconstruction method for SWEs on Cartesian grids only 
    void cuTransportNSWEsSRMCartesian(cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Scalar, on_cell>& z, cuFvMappedField<Vector, on_cell>& z_gradient, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Scalar, on_cell> hC, cuFvMappedField<Scalar, on_cell>& h_advection, cuFvMappedField<Vector, on_cell>& hU_advection, cuFvMappedField<Scalar, on_cell>& hC_advection);

  }


}  

#endif
