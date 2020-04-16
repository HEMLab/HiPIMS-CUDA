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
  \flie cuda_field_algebra.h
  \brief Header file for field algebra template functions
  \version 0.1
  \author xilin xia
*/


#ifndef CUDA_FIELD_ALGEBRA_H
#define CUDA_FIELD_ALGEBRA_H

#include "mapped_field.h"
#include "cuda_mapped_field.h"
#include "cuda_kernel_launch_parameters.h"


namespace GC{

  namespace fv{

    template <typename T_a, typename T_b, typename F>
    __global__ void cuUnaryKernel(T_a* a, T_b* b, F func, unsigned int N){
      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      while(index  < N){
        b[index] = func(a[index]);
        index += blockDim.x * gridDim.x;
      }
    }
    
    ///Unary field operation by CUDA, output is written to the second argument
    template <typename T_a, typename T_b, MAPPING_MODES C, typename F>
    void cuUnary(cuFvMappedField<T_a, C>& a, cuFvMappedField<T_b, C>& b, F func){
      
      cuUnaryKernel<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK>>>(a.data.dev_ptr(), b.data.dev_ptr(), func, a.data.size());
      cuUnaryKernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(a.boundary_value.dev_ptr(), b.boundary_value.dev_ptr(), func, a.boundary_value.size());


    }

    template <typename T, typename F>
    __global__ void cuUnaryOnKernel(T* a, F func, unsigned int N){
      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      while(index  < N){
        a[index] = func(a[index]);
        index += blockDim.x * gridDim.x;
      }
    }
    
    ///Unary field operation by CUDA, output is written to the first argument, i.e. the original field
    template <typename T, MAPPING_MODES C, typename F>
    void cuUnaryOn(cuFvMappedField<T, C>& a, F func){
      
      cuUnaryOnKernel<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK>>>(a.data.dev_ptr(), func, a.data.size());
      cuUnaryOnKernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(a.boundary_value.dev_ptr(), func, a.boundary_value.size());

    }

    template <typename T_a, typename T_b, typename T_c, typename F>
    __global__ void cuBinaryKernel(T_a* a, T_b* b, T_c* c, F func, unsigned int N){
      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      while(index  < N){
        c[index] = func(a[index], b[index]);
        index += blockDim.x * gridDim.x;
      }
    }

    ///Binary field operation by CUDA, output is written to the third argument
    template <typename T_a, typename T_b, typename T_c, MAPPING_MODES C, typename F>
    void cuBinary(cuFvMappedField<T_a, C>& a, cuFvMappedField<T_b, C>& b, cuFvMappedField<T_c, C>& c, F func){
      
      cuBinaryKernel<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK>>>(a.data.dev_ptr(), b.data.dev_ptr(), c.data.dev_ptr(), func, a.data.size());
      cuBinaryKernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(a.boundary_value.dev_ptr(), b.boundary_value.dev_ptr(), c.boundary_value.dev_ptr(), func, a.boundary_value.size());


    }

    template <typename T_a, typename T_b, typename F>
    __global__ void cuBinaryOnKernel(T_a* a, T_b* b, F func, unsigned int N){
      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      while(index  < N){
        a[index] = func(a[index], b[index]);
        index += blockDim.x * gridDim.x;
      }
    }
    
    ///Binary field operation by CUDA, output is written to the first argument, i.e. the first original field
    template <typename T_a, typename T_b, MAPPING_MODES C, typename F>
    void cuBinaryOn(cuFvMappedField<T_a, C>& a, cuFvMappedField<T_b, C>& b, F func){
      
      cuBinaryOnKernel<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK>>>(a.data.dev_ptr(), b.data.dev_ptr(), func, a.data.size());
      cuBinaryOnKernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(a.boundary_value.dev_ptr(), b.boundary_value.dev_ptr(), func, a.boundary_value.size());


    }


  }//end of namespace fv

}//--end of namespace GC


#endif
