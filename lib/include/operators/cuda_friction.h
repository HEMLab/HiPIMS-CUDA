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
\brief Header file for friction operator

\version 0.1
\author xilin xia

*/

#ifndef CUDA_FRICTION_H
#define CUDA_FRICTION_H

#include "cuda_mapped_field.h"
#include "Scalar.h"
#include "Vector.h"
#include "Tensor.h"
#include "Flag.h"

namespace GC{

  namespace fv{

    ///This is a function class for calculation friction force
    void cuFrictionMC(cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& friction_coeff, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Vector, on_cell>& z_grad, cuFvMappedField<Vector, on_cell>& friction_force);

    void cuFrictionMCBalanced(Scalar dt, cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& friction_coeff, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Vector, on_cell>& z_grad, cuFvMappedField<Vector, on_cell>& friction_force);

    void cuFrictionMCWithWall(Scalar dt, cuFvMappedField<Scalar, on_cell>& gravity, Scalar friction_coeff1, Scalar friction_coeff2, Scalar wall_width, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Vector, on_cell>& z_grad, cuFvMappedField<Vector, on_cell>& friction_force);

    void cuFrictionPlastic(Scalar dt, cuFvMappedField<Scalar, on_cell>& retarding_stress, cuFvMappedField<Scalar, on_cell>& rho, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Vector, on_cell>& z_grad, cuFvMappedField<Vector, on_cell>& friction_force);

    void cuFrictionMCPlastic(Scalar dt, cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& miu, cuFvMappedField<Scalar, on_cell>& retarding_stress, cuFvMappedField<Scalar, on_cell>& rho, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Vector, on_cell>& z_grad, cuFvMappedField<Vector, on_cell>& friction_force);

    void cuFrictionMuI(Scalar dt, cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& miu_1, cuFvMappedField<Scalar, on_cell>& miu_2, Scalar beta, Scalar L, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Vector, on_cell>& z_grad, cuFvMappedField<Vector, on_cell>& friction_force);

    void cuFrictionLucas(Scalar dt, cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& miu_1, cuFvMappedField<Scalar, on_cell>& miu_2, Scalar U_w, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Vector, on_cell>& z_grad, cuFvMappedField<Vector, on_cell>& friction_force);

    void cuFrictionManningImplicit(Scalar dt, cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& manning_coeff, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Vector, on_cell>& hU_advection);

    void cuFrictionManningImplicitWithForce(Scalar dt, cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& manning_coeff, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Vector, on_cell>& hU_advection, cuFvMappedField<Vector, on_cell>& friction_force);
    
    void cuFrictionManningExplicit(Scalar dt, cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& manning_coeff, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Vector, on_cell>& hU_advection);

    void cuFrictionManningSplitting(Scalar dt, cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& manning_coeff, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Vector, on_cell>& hU_advection);
  
    void cuFrictionManningSteepImplicit(Scalar dt, cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& manning_coeff, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Vector, on_cell>& hU_advection, cuFvMappedField<Vector, on_cell>& z_grad);
  
    void cuFrictionManningMCImplicit(Scalar dt, Scalar rho_water, Scalar rho_solid, cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& manning_coeff, cuFvMappedField<Scalar, on_cell>& friction_coeff, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Scalar, on_cell>& hC, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Vector, on_cell>& hU_advection);
  }//end of namespace fv

}//end of namespace GC



#endif

