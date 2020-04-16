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
 \file riemann.h
 \brief Header file for riemann solvers for SWEs 

*/

#ifndef RIEMANN_H
#define RIEMANN_H

#include "Scalar.h"
#include "Vector.h"
#include "hemi.h"

namespace GC{

  /// namespace for finite volume method
  namespace fv{

    /// This is a template class for riemann state
    template <typename T> class RiemannState{
      public:
        HEMI_DEV_CALLABLE_INLINE_MEMBER
        RiemannState(const T& _L, const T& _R) : L(_L), R(_R){}
        T L;
        T R;
    };

    typedef RiemannState<Scalar> ScalarRiemannState;
    typedef RiemannState<Vector3> VectorRiemannState;

    /// This is a class for riemann flux 
    class RiemannFluxSWEs{
      public:
        HEMI_DEV_CALLABLE_INLINE_MEMBER
        RiemannFluxSWEs(const Scalar& _h, const Vector2& _q) : h(_h), q(_q){}
        Scalar h;
        Vector2 q;
    };

    ///This function is the hllc riemann solver for shallow water equations
    RiemannFluxSWEs hllcRiemannSolverSWEs(const Scalar& g, const ScalarRiemannState& h,const VectorRiemannState& q);
    ///This function is the modified hllc riemann solver for shallow water equations with non-uniform gravity field
    RiemannFluxSWEs hllcRiemannSolverNonUniGravitySWEs(const ScalarRiemannState& g, const ScalarRiemannState& h,const VectorRiemannState& q);


  } //--end namespace of finite volume method

}

#endif
