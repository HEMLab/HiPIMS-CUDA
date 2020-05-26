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
  \flie fluxBalanced_NSWEs.h
  \brief Header file for flux in non hydrostatic shallow water equations

  \version 0.1
  \author xilin xia
*/

#ifndef FLUX_BALANCED_NSWES_H
#define FLUX_BALANCED_NSWES_H

#include "mapped_field.h"
#include <functional>
#include "riemann.h"

namespace GC{
  
  ///fvc denotes for finite volume method on Cartesian grids
  namespace fv{
    ///calculation flux in shallow water equations by MUSCL scheme
    class FluxBalancedNSWEs_1st{
      public:
        FluxBalancedNSWEs_1st(std::function<fv::RiemannFluxSWEs(const Scalar&, 
                                                      const ScalarRiemannState&, 
                                                      const VectorRiemannState&)> _RiemannSolver)
         :RiemannSolver(_RiemannSolver){}
        void operator() (fvScalarFieldOnCell& g, fvScalarFieldOnCell& h, fvVectorFieldOnCell& hU,
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

