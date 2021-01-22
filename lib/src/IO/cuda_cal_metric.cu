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
  \file cuda_cal_metric.cu
  \brief Source file for writing output as arcgis ascii files

*/

#include "cuda_cal_metric.h"

namespace GC{
  __global__ void cuCalDepthDurationKernel(Scalar* h_old, Scalar* h, Scalar* t_hGTx, Scalar x, Scalar dt, unsigned int phi_size){

    unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
    while(index < phi_size){
      Scalar t = t_hGTx[index];
      if(h_old[index] > x &&  h[index] > x){
        t += dt;
      }else{
        t = 0.0;
      }
      t_hGTx[index] = t;
      index += blockDim.x * gridDim.x;
    }

  }

  void cuCalDepthDuration(cuFvMappedField<Scalar, on_cell>& h_old, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Scalar, on_cell>& t_hGTx, Scalar x, Scalar dt){

    cuCalDepthDurationKernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(h_old.data.dev_ptr(),
    h.data.dev_ptr(),
    t_hGTx.data.dev_ptr(),
    x,
    dt,
    h.data.size());

  }

}