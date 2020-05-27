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
  \flie flux_SWEs.h
  \brief Header file for flux in shallow water equations

  \version 0.1
  \author xilin xia
*/

#ifndef FLUX_BALANCED_SWES_H
#define FLUX_BALANCED_SWES_H

#include "mapped_field.h"
#include <functional>
#include "riemann.h"

namespace GC{
  
  ///fvc denotes for finite volume method on Cartesian grids
  namespace fv{
    ///calculation flux in shallow water equations by MUSCL scheme
    class FluxBalancedSWEs_1st{
      public:
        FluxBalancedSWEs_1st(std::function<fv::RiemannFluxSWEs(const Scalar&, 
                                                      const ScalarRiemannState&, 
                                                      const VectorRiemannState&)> _RiemannSolver)
         :RiemannSolver(_RiemannSolver){}
        void operator() (Scalar g, fvScalarFieldOnCell& h, fvVectorFieldOnCell& hU,
                        fvScalarFieldOnCell& z, fvScalarFieldOnCell& h_flux, 
                        fvVectorFieldOnCell& hU_flux);
      private:
        std::function<fv::RiemannFluxSWEs(const Scalar&, 
                                      const ScalarRiemannState&, 
                                      const VectorRiemannState&)> RiemannSolver;
    };


  }

}

#endif

