// ======================================================================================
// Name                :    GeoClasses : Generic Geophysical Flow Modelling Framework
// Description         :    This code pack provides a generic framework for developing 
//                          Geophysical CFD software.
// ======================================================================================
// Version             :    0.1 
// Author              :    Xilin Xia (PhD candidate in Newcastle University)
// Create Time         :    2014/10/04
// Update Time         :    2015/10/25
// ======================================================================================
// Copyright @ Xilin Xia 2015 . All rights reserved.
// ======================================================================================

/*!
\file cuda_adaptive_time_control.cu
\brief Source file for adaptive time control class

\version 0.1
\author xilin xia

*/

#include "cuda_adaptive_time_control.h"
#include <thrust/device_ptr.h>
#include <thrust/reduce.h>
#include <algorithm>
#include "cuda_kernel_launch_parameters.h"
#include "hemi_error.h"


namespace GC{

  const int threadsPerBlock = 256;

  __global__ void cuUpdateByCFLReduceKernel1(Scalar* dt_array, unsigned int N){

    __shared__ Scalar cache[threadsPerBlock];
    int index = threadIdx.x + blockIdx.x * blockDim.x; 
    int cacheIndex = threadIdx.x;
    Scalar temp = 3e35;
    while (index < N){
      temp = fmin(temp, dt_array[index]);
      index += blockDim.x * gridDim.x;
    }

    cache[cacheIndex] = temp;

    __syncthreads();

    int i = blockDim.x / 2;
    while (i != 0) {
      if (cacheIndex < i){
        cache[cacheIndex] = fmin(cache[cacheIndex], cache[cacheIndex + i]);
      }
      __syncthreads();
      i /= 2;
    }

    if (cacheIndex == 0){
      dt_array[blockIdx.x] = cache[0];
    }

  }

  __global__ void cuUpdateByCFLReduceKernel2(Scalar* dt_array, unsigned int N){

    __shared__ Scalar cache[threadsPerBlock];
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    int cacheIndex = threadIdx.x;

    if (cacheIndex < N){
      cache[cacheIndex] = dt_array[index];
    }
    else{
      cache[cacheIndex] = 1e35;
    }

    __syncthreads();

    int i = blockDim.x / 2;
    while (i != 0) {
      if (cacheIndex < i){
        cache[cacheIndex] = fmin(cache[cacheIndex], cache[cacheIndex + i]);
      }
      __syncthreads();
      i /= 2;
    }

    if (cacheIndex == 0){
      dt_array[blockIdx.x] = cache[0];
    }

  }


  __global__ void cuUpdateByCFLKernel(Scalar* gravity, Scalar* h, Vector2* hU, Scalar* cell_volumes, Scalar* dt_array, Scalar Courant_, unsigned int N){
    
    unsigned int index = blockDim.x * blockIdx.x + threadIdx.x; 
    while(index < N){
      auto gravity_ = gravity[index];
      auto h_ = h[index];
      auto hU_ = hU[index];
      Vector2 u_ = 0.0;
      if (h_ >= 1e-10){
        u_ = hU_/h_;
        auto c_sound = sqrt(gravity_*h_);
        auto volume_ = cell_volumes[index];
        Scalar dl = sqrt(volume_);
        Scalar dt1 = dl/(c_sound + fabs(u_.x));
        Scalar dt2 = dl/(c_sound + fabs(u_.y));
        dt_array[index] = Courant_*fmin(dt1,dt2);
      }else{
        dt_array[index] = 3e35;
      }
      index += blockDim.x * gridDim.x;
    }

  }


  void cuAdaptiveTimeControl2D::updateByCFL(cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Vector2, on_cell>& hU){

    if(is_beginning){
      dt_array_ptr = std::make_shared<cuFvMappedField<Scalar, on_cell> >(h, partial);
      is_beginning = false;
    }

    cuUpdateByCFLKernel<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK>>>(gravity.data.dev_ptr(), h.data.dev_ptr(), hU.data.dev_ptr(), h.mesh->cell_volumes.dev_ptr(), dt_array_ptr->data.dev_ptr(), Courant_, h.data.size());
    cuUpdateByCFLReduceKernel1 <<<BLOCKS_PER_GRID, threadsPerBlock >>>(dt_array_ptr->data.dev_ptr(),  dt_array_ptr->data.size());
    cuUpdateByCFLReduceKernel2 << <1, threadsPerBlock >> >(dt_array_ptr->data.dev_ptr(), dt_array_ptr->data.size());
//    dt_ = thrust::reduce(thrust::device_ptr <Scalar> (dt_array_ptr->data.dev_ptr()), thrust::device_ptr <Scalar> (dt_array_ptr->data.dev_ptr() + dt_array_ptr->data.size()),(Scalar) 3e35, thrust::minimum<Scalar>());
    checkCuda(cudaMemcpy(&dt_, dt_array_ptr->data.dev_ptr(), sizeof(Scalar), cudaMemcpyDeviceToHost));
    dt_ = std::min(dt_, 60.0); //maximum time step is 1 minite
  }

  bool cuAdaptiveTimeControl2D::is_end(){
    return current_ >= end_;
  }

  Scalar cuAdaptiveTimeControl2D::current(){
    return current_;
  }

  Scalar cuAdaptiveTimeControl2D::dt(){
    return dt_;
  }

  void cuAdaptiveTimeControl2D::forward(){
    dt_ = std::min(dt_, end_ - current_);
    current_ += dt_;
  }

  void cuAdaptiveTimeControl2D::set_dt(Scalar dt){
    dt_ = dt;
  }

}//--end of namespace GC
