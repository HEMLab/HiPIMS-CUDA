// ======================================================================================
// Name                :    GeoClasses : Generic Geophysical Flow Modelling Framework
// Description         :    This code pack provides a generic framework for developing 
//                          Geophysical CFD software.
// ======================================================================================
// Version             :    0.1 
// Author              :    Xilin Xia (PhD candidate in Newcastle University)
// Create Time         :    2014/10/04
// Update Time         :    2016/11/14
// ======================================================================================
// Copyright @ Xilin Xia 2016 . All rights reserved.
// ======================================================================================

/*!
\file cuda_erosion_deposition.h
\brief Header file for erosion and deposition rate

\version 0.1
\author xilin xia

*/

#ifndef CUDA_EROSION_DEPOSITION_H
#define CUDA_EROSION_DEPOSITION_H

#include "cuda_mapped_field.h"
#include "Scalar.h"
#include "Vector.h"
#include "Tensor.h"
#include "Flag.h"


namespace GC{

  namespace fv{


    void cuEDMeyerPeterMuller(Scalar dt, cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Scalar, on_cell>& hC, cuFvMappedField<Scalar, on_cell>& manning_coeff, cuFvMappedField<Scalar, on_cell>& ED_rate, Scalar rho_solid, Scalar rho_water, Scalar dim_mid);

    void cuMomentumCorrection(cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Scalar, on_cell>& hC, cuFvMappedField<Vector, on_cell>& C_grad, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Scalar, on_cell>& ED_rate,
      cuFvMappedField<Vector, on_cell>& mom_correction, Scalar rho_solid, Scalar rho_water, Scalar porosity);

    void cuEDTakahashiIversonXia(cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Scalar, on_cell>& hC, cuFvMappedField<Scalar, on_cell>& manning_coeff, cuFvMappedField<Scalar, on_cell>& miu_dynamic, cuFvMappedField<Scalar, on_cell>& miu_static, cuFvMappedField<Scalar, on_cell>& ED_rate, Scalar rho_solid, Scalar rho_water, Scalar porosity, Scalar dim_mid);

    void cuBankCollapse(cuFvMappedField<Scalar, on_cell>& z, cuFvMappedField<Scalar, on_cell>& z_new, Scalar critical_slope);

  }

}

#endif
