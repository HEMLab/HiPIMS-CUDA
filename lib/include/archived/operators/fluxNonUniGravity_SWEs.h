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
  \brief Header file for shallow water equations flux with non uniform gravity field

  \version 0.1
  \author xilin xia
*/


#ifndef FLUX_NONUNI_GRAVITY_SWES_H
#define FLUX_NONUNI_GRAVITY_SWES_H

#include "mapped_field.h"
#include <functional>
#include "riemann.h"

namespace GC{

  namespace fv{

    ///This is a function class for calculating face fluxes of shallow water equations
    class FluxNonUniGravitySWEs_1st{
      public:
        FluxNonUniGravitySWEs_1st(std::function<RiemannFluxSWEs(const ScalarRiemannState&, 
                                                   const ScalarRiemannState&, 
                                                   const VectorRiemannState&)> _RiemannSolver)
         :RiemannSolver(_RiemannSolver){} 
        void operator() (fvScalarFieldOnCell& g,
                         fvScalarFieldOnCell& h, 
                         fvVectorFieldOnCell& hU, 
                         fvScalarFieldOnCell& h_flux,
                         fvVectorFieldOnCell& hU_flux); 

      private:
        std::function<RiemannFluxSWEs(const ScalarRiemannState&, 
                                      const ScalarRiemannState&, 
                                      const VectorRiemannState&)> RiemannSolver;



    };

  }

}

#endif
