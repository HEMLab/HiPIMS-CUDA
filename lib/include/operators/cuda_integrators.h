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
  \flie cuda_integrators.h
  \brief Header file for integrator template functions
  \version 0.1
  \author xilin xia
*/

#ifndef CUDA_INTEGRATORS_H
#define CUDA_INTEGRATORS_H

#include "mapped_field.h"
#include "cuda_mapped_field.h"
#include "cuda_kernel_launch_parameters.h"
#include "Scalar.h"

namespace GC{


  namespace fv{

    template <typename T>
    __global__ void cuEulerIntegratorKernel(T* phi, T* phi_acc, Scalar delta_t, unsigned int N){
      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      while(index  < N){
        phi[index] = phi[index] + delta_t*phi_acc[index];
        index += blockDim.x * gridDim.x;
      }
    }


    template <typename T, MAPPING_MODES C>
    void cuEulerIntegrator(cuFvMappedField<T, C>& phi, cuFvMappedField<T, C>& phi_acc, Scalar delta_t, Scalar current_t){

      cuEulerIntegratorKernel<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK>>>(phi.data.dev_ptr(), phi_acc.data.dev_ptr(), delta_t, phi.data.size());
      phi.update_time(current_t, delta_t);
      phi.update_boundary_values();

    }

  }//end of namespace fv


}//--end of namespace GC


#endif
