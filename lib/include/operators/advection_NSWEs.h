// ======================================================================================
// Name                :    GeoClasses : Generic Geophysical Flow Modelling Framework
// Description         :    This code pack provides a generic framework for developing 
//                          Geophysical CFD software.
// ======================================================================================
// Version             :    0.1 
// Author              :    Xilin Xia (PhD candidate in Newcastle University)
// Create Time         :    2014/10/04
// Update Time         :    2015/10/09
// ======================================================================================
// Copyright @ Xilin Xia 2015 . All rights reserved.
// ======================================================================================

/*!
  \flie advection_NSWEs.h
  \brief Header file for flux in non hydrostatic shallow water equations

  \version 0.1
  \author xilin xia
*/

#ifndef ADVECTION_NSWES_H
#define ADVECTION_NSWES_H

#include "mapped_field.h"
#include <functional>
#include "riemann.h"

namespace GC{
  
  ///fvc denotes for finite volume method on Cartesian grids
  namespace fv{
    ///calculation flux in shallow water equations by MUSCL scheme
    void AdvectionNSWEs_2nd(fvScalarFieldOnCell& g, fvScalarFieldOnCell& h, fvVectorFieldOnCell& hU, fvScalarFieldOnCell& z, fvVectorFieldOnCell& grad_h, fvVectorFieldOnCell& grad_surface, fvTensorFieldOnCell& grad_u, fvScalarFieldOnCell& h_flux,
      fvVectorFieldOnCell& hU_flux);

  }

}

#endif
