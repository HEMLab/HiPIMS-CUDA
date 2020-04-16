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
  \flie cuda_advection_NSWEs.h
  \brief Header file for advection of non hydrostatic shallow water equations

  \version 0.1
  \author xilin xia
*/

#ifndef CUDA_ADVECTION_NSWES_H
#define CUDA_ADVECTION_NSWES_H

#include "cuda_mapped_field.h"
#include "Scalar.h"
#include "Vector.h"
#include "Tensor.h"
#include "Flag.h"


namespace GC{


  namespace fv{

    void cuAdvectionNSWEs2nd(cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Scalar, on_cell>& z, cuFvMappedField<Vector, on_cell>& u, cuFvMappedField<Vector, on_cell>& h_gradient, cuFvMappedField<Vector, on_cell>& eta_gradient, cuFvMappedField<Tensor, on_cell>& u_gradient, cuFvMappedField<Scalar, on_cell>& h_advection, cuFvMappedField<Vector, on_cell>& hU_advection);

    void cuAdvectionNSWEs2ndRobust(cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Scalar, on_cell>& z, cuFvMappedField<Vector, on_cell>& u, cuFvMappedField<Vector, on_cell>& h_gradient, cuFvMappedField<Vector, on_cell>& eta_gradient, cuFvMappedField<Tensor, on_cell>& u_gradient, cuFvMappedField<Scalar, on_cell>& h_advection, cuFvMappedField<Vector, on_cell>& hU_advection);

    void cuAdvectionNSWEs2ndCurv(cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& centrifugal, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Scalar, on_cell>& z, cuFvMappedField<Vector, on_cell>& u, cuFvMappedField<Vector, on_cell>& h_gradient, cuFvMappedField<Vector, on_cell>& eta_gradient, cuFvMappedField<Tensor, on_cell>& u_gradient, cuFvMappedField<Scalar, on_cell>& h_advection, cuFvMappedField<Vector, on_cell>& hU_advection);

    ///This function is an implementation which avoids accumulating velocity
    void cuAdvectionNSWEs2ndRobustCurv(cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& centrifugal, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Scalar, on_cell>& z, cuFvMappedField<Vector, on_cell>& u, cuFvMappedField<Vector, on_cell>& h_gradient, cuFvMappedField<Vector, on_cell>& eta_gradient, cuFvMappedField<Tensor, on_cell>& u_gradient, cuFvMappedField<Scalar, on_cell>& h_advection, cuFvMappedField<Vector, on_cell>& hU_advection);

    ///This function is faster, but it takes two more arguments
    void cuAdvectionNSWEs2ndFast(cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Scalar, on_cell>& z, cuFvMappedField<Vector, on_cell>& u, cuFvMappedField<Vector, on_cell>& h_gradient, cuFvMappedField<Vector, on_cell>& eta_gradient, cuFvMappedField<Tensor, on_cell>& u_gradient, cuFvMappedField<Scalar, on_halffacet>& h_flux, cuFvMappedField<Vector, on_halffacet>& hU_flux, cuFvMappedField<Scalar, on_cell>& h_advection, cuFvMappedField<Vector, on_cell>& hU_advection);

    ///Calculating advaction based on a modified SWEs, e.g. Berthon (2012) Efficient well-balanced hydrostatic upwind schemes for shallow-water equations
    void cuAdvectionMSWEs(cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Scalar, on_cell>& z, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Scalar, on_cell>& h_advection, cuFvMappedField<Vector, on_cell>& hU_advection, cuFvMappedField<Scalar, on_cell>& dt_mass);
  
    ///A fast 1st order hydrostatic reconstruction method for SWEs on Cartesian grids only 
    void cuAdvectionSWEsCartesian(cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Scalar, on_cell>& z, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Scalar, on_cell>& h_advection, cuFvMappedField<Vector, on_cell>& hU_advection);
    
    ///A fast 1st order hydrostatic reconstruction method for SWEs on Cartesian grids only Recording mass flux
    void cuAdvectionSWEsRFCartesian(cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Scalar, on_cell>& z, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Scalar, on_cell>& h_advection, cuFvMappedField<Vector, on_cell>& hU_advection, cuFvMappedField<Vector, on_cell>& hU_EW, cuFvMappedField<Vector, on_cell>& hU_NS);

    ///A fast 1st order hydrostatic reconstruction method for SWEs on Cartesian grids only Asymptotic Preserving    
    void cuAdvectionSWEsAPCartesian(cuFvMappedField<Scalar, on_cell>& manning, cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Scalar, on_cell>& z, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Scalar, on_cell>& h_advection, cuFvMappedField<Vector, on_cell>& hU_advection);

    ///A fast 1st order surface reconstruction method for SWEs on Cartesian grids only 
    void cuAdvectionMSWEsCartesian(cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Scalar, on_cell>& z, cuFvMappedField<Vector, on_cell>& z_gradient, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Scalar, on_cell>& h_advection, cuFvMappedField<Vector, on_cell>& hU_advection);
    
    ///A fast 1st order surface reconstruction method for non-hydrostatic SWEs on Cartesian grids only
    void cuAdvectionNSWEsSRMCartesian(cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& centrifugal, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Scalar, on_cell>& z, cuFvMappedField<Vector, on_cell>& z_gradient, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Scalar, on_cell>& h_advection, cuFvMappedField<Vector, on_cell>& hU_advection);

    ///A fast 2nd order hydrostatic reconstruction method for SWEs on Cartesian grids only 
    void cuAdvectionNSWEs2ndCartesian(cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Scalar, on_cell>& z, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Vector, on_cell>& h_gradient, cuFvMappedField<Vector, on_cell>& eta_gradient, cuFvMappedField<Tensor, on_cell>& hU_gradient, cuFvMappedField<Scalar, on_cell>& h_advection, cuFvMappedField<Vector, on_cell>& hU_advection);

    ///D. L. George's Riemann solver with source terms on Cartesian grids only
    void cuAdvectionSWEsGeorgeCartesian(cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Scalar, on_cell>& z, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Scalar, on_cell>& h_advection, cuFvMappedField<Vector, on_cell>& hU_advection);
    
    ///D. L. George's Riemann solver with source terms on Cartesian grids only Recording mass flux
    void cuAdvectionSWEsRFGeorgeCartesian(cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Scalar, on_cell>& z, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Scalar, on_cell>& h_advection, cuFvMappedField<Vector, on_cell>& hU_advection, cuFvMappedField<Vector, on_cell>& hU_EW, cuFvMappedField<Vector, on_cell>& hU_NS);

  }//end of namespace fv---------------


}//--end of namespace GC---------------


#endif
