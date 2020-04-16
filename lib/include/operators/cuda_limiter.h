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
  \flie cuda_limiter.h
  \brief Header file for gradient limiter of an arbitary variable phi by CUDA

  \version 0.1
  \author xilin xia
*/

#ifndef CUDA_LIMITER_H
#define CUDA_LIMITER_H

#include "cuda_mapped_field.h"


namespace GC{

  namespace fv{
    
    ///This function limits the 2D scalar gradient by solving an LP program (Berger et al 2005)
    void cuGradientLimiter(cuFvMappedField<Scalar, on_cell>& phi, cuFvMappedField<Vector2, on_cell>& phi_gradient);

    ///This function limits the 2D vector gradient by solving an LP program (Berger et al 2005)
    void cuGradientLimiter(cuFvMappedField<Vector2, on_cell>& phi, cuFvMappedField<Tensor2, on_cell>& phi_gradient);

    ///This function limits the 2D scalar gradient on Cartesian grids
    void cuGradientLimiterCartesian(cuFvMappedField<Scalar, on_cell>& phi, cuFvMappedField<Vector2, on_cell>& phi_gradient);

    ///This function limits the 2D vector gradient on Cartesian grids
    void cuGradientLimiterCartesian(cuFvMappedField<Vector2, on_cell>& phi, cuFvMappedField<Tensor2, on_cell>& phi_gradient);

  }//end of namespace fv

}//--end of namespace GC


#endif
