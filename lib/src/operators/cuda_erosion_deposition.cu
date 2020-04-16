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
  \flie cuda_erosion_deposition.cu
  \brief Source file for erosion and deposition rate

  \version 0.1
  \author xilin xia
*/

#include "cuda_erosion_deposition.h"
#include "cuda_kernel_launch_parameters.h"


namespace GC{

  namespace fv{

    __global__ void cuEDMeyerPeterMullerKernel(Scalar dt, Scalar* gravity, Scalar* h, Vector* hU, Scalar* hC, Scalar* manning_coeff, Scalar* ED_rate, Scalar rho_solid, Scalar rho_water, Scalar dim_mid, unsigned int size){
      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      while (index < size){
        auto g = gravity[index];
        auto h_ = h[index];
        auto hC_ = hC[index];
        auto hU_ = hU[index];
        auto manning_ = manning_coeff[index];
        Vector2 U = 0.0;
        Scalar C = 0.0;
        Scalar small_value = 1e-10;
        if (h_ > small_value){
          C = hC_ / h_;
          U = hU_ / h_;
        }
        else{
          C = 0.0;
          U = 0.0;
        }
        //Scalar rho_mix = rho_water;
        //Scalar rho_mix = rho_water*(1.0-C) + rho_solid*C;
        Scalar rho_mix = rho_water;
        Scalar tau = 0.0;
        if (h_ > dim_mid){
          tau = rho_mix*g*pow(manning_,2.0)*pow(h_, -1.0/3.0)*dot(U,U);
        }
        else{
          tau = 0.0;
        }
        Scalar tau_cr = (rho_solid - rho_water)*g*dim_mid;
        Scalar shields = tau / tau_cr;
        Scalar qb_sat = 8.0*pow(fmax(shields - 0.047,0.0), 1.5)*sqrt((rho_solid/rho_water-1.0)*g*pow(dim_mid, 3.0));
        Scalar L_sat = 2.67*rho_solid/rho_water*dim_mid*(sqrt(shields/0.047) - 1.0); //This equation comes from the book "Granular Media: Between ..." p 339
        //Scalar visc = 1e-6;
        //Scalar v_settling = sqrt(pow(13.95*visc / dim_mid, 2.0) + 1.09*g*(rho_solid-rho_water)/rho_water*dim_mid) - 13.95*visc / dim_mid;
        //Scalar L_sat = hC_/(1-0.42)*norm(U) / v_settling; //Xia
        //L_sat = fmax(8.0*h_, L_sat);
        L_sat = fmax(5.0*dim_mid, L_sat);
        Scalar ED_rate_ = (qb_sat - norm(hU_)*C) / L_sat;
        Scalar hC_sat = 0.0;
        if (norm(U) > 0.0){
          hC_sat = qb_sat / norm(U);
        }
        else{
          hC_sat = 0.0;
        }
        //hC_sat = fmin(hC_sat, 1.0);
        if (hC_ > hC_sat){
          ED_rate_ = fmax((hC_sat - hC_) / dt, ED_rate_);
        }
        else{
          ED_rate_ = fmin((hC_sat - hC_) / dt, ED_rate_);
        }
        //if (index == 382168){
        //  printf("h %f, hC_sat %f, hC_ %f, ED_rate %f\n", h_, hC_sat, hC_, ED_rate_);
        //}
        ED_rate[index] = ED_rate_;
        index += blockDim.x * gridDim.x;
      }
    }

    void cuEDMeyerPeterMuller(Scalar dt, cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Scalar, on_cell>& hC, cuFvMappedField<Scalar, on_cell>& manning_coeff, cuFvMappedField<Scalar, on_cell>& ED_rate, Scalar rho_solid, Scalar rho_water, Scalar dim_mid){

      auto mesh = h.mesh;

      cuEDMeyerPeterMullerKernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(
        dt,
        gravity.data.dev_ptr(), 
        h.data.dev_ptr(),
        hU.data.dev_ptr(),
        hC.data.dev_ptr(),
        manning_coeff.data.dev_ptr(),
        ED_rate.data.dev_ptr(),
        rho_solid,
        rho_water,
        dim_mid,
        h.data.size()
        );
    }

    __global__ void cuMomentumCorrectionKernel(Scalar* gravity, Scalar* h, Scalar* hC, Vector* C_grad, Vector* hU, Scalar* ED_rate, Vector* mom_correction, Scalar rho_solid, Scalar rho_water, Scalar porosity, unsigned int size){
      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      while (index < size){
        auto g = gravity[index];
        auto h_ = h[index];
        auto hC_ = hC[index];
        auto C_grad_ = C_grad[index];
        auto hU_ = hU[index];
        Scalar small_value = 1e-4;
        Vector2 U = 0.0;
        Scalar C = 0.0;
        if (h_ > small_value){
          C = hC_ / h_;
          U = hU_ / h_;
        }
        else{
          C = 0.0;
          U = 0.0;
        }
        auto ED_rate_ = ED_rate[index];
        Scalar rho_mix = rho_water*(1 - C) + rho_solid*C;
        Scalar rho_0 = rho_water*porosity + rho_solid*(1.0 - porosity);
        Vector2 mom_correction_ = -1.0*(rho_solid - rho_water)*g*pow(h_, 2.0) / rho_mix / 2.0*C_grad_ - (rho_0 - rho_mix) / rho_mix / (1.0 - porosity)*ED_rate_*U;
        mom_correction[index] = mom_correction_;
        index += blockDim.x * gridDim.x;
      }
    }

    void cuMomentumCorrection(cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Scalar, on_cell>& hC, cuFvMappedField<Vector, on_cell>& C_grad, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Scalar, on_cell>& ED_rate, 
      cuFvMappedField<Vector, on_cell>& mom_correction, Scalar rho_solid, Scalar rho_water, Scalar porosity){

      cuMomentumCorrectionKernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(
        gravity.data.dev_ptr(),
        h.data.dev_ptr(), 
        hC.data.dev_ptr(), 
        C_grad.data.dev_ptr(), 
        hU.data.dev_ptr(), 
        ED_rate.data.dev_ptr(), 
        mom_correction.data.dev_ptr(), 
        rho_solid,
        rho_water,
        porosity,
        h.data.size());

    }

    __global__ void cuEDTakahashiIversonXiaKernel(Scalar* gravity, Scalar* h, Vector* hU, Scalar* hC, Scalar* manning_coeff, Scalar* miu_dynamic, Scalar* miu_static, Scalar* ED_rate, Scalar rho_solid, Scalar rho_water, Scalar porosity, Scalar dim_mid, unsigned int size){
      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      while (index < size){
        Scalar g = gravity[index];
        Scalar h_ = h[index];
        Scalar hC_ = hC[index];
        Vector2 hU_ = hU[index];
        Vector2 U_ = 0.0;
        Scalar C_ = 0.0;
        Scalar manning_ = manning_coeff[index];
        Scalar miu_dynamic_ = miu_dynamic[index];
        Scalar miu_static_ = miu_static[index];
        Scalar small_value = 1e-10;
        if (h_ < small_value){
          C_ = 0.0;
          U_ = 0.0;
        }
        else{
          C_ = hC_ / h_;
          U_ = hU_ / h_;
        }
        Scalar rho_mix = rho_water*(1.0 - C_) + rho_solid*C_;
        Scalar tau = 0.0;
        if (h_ > dim_mid){
          tau = rho_mix*g*pow(manning_, 2.0)*pow(h_, -1.0 / 3.0)*dot(U_,U_) + (rho_solid - rho_water)*hC_*g*miu_dynamic_;
        }
        else{
          tau = 0.0;
        }
        Scalar tau_0 = (rho_solid - rho_water)*hC_*g*miu_static_;
        Scalar rho_bed = porosity*rho_water + (1.0 - porosity)*rho_solid;
        Scalar _E_rate = 0.0;
        Scalar _D_rate = 0.0;
        if(tau > tau_0){
          _E_rate = (tau - tau_0)*(1.0 - porosity) / (rho_bed*fmax(norm(U_),0.1));
          _D_rate = 0.0;
        }else{
          _E_rate = 0.0;
          _D_rate = (tau - tau_0)*(1.0 - porosity) / (rho_bed*fmax(norm(U_),0.1));
        }
        Scalar _ED_rate = _E_rate + _D_rate;
        //if(index == 27475){
        //  printf("%f %f %f\n", _ED_rate, hC_, C_);
        //}
        ED_rate[index] = _ED_rate;
        index += blockDim.x * gridDim.x;
      }
    }

    void cuEDTakahashiIversonXia(cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Scalar, on_cell>& hC, cuFvMappedField<Scalar, on_cell>& manning_coeff, cuFvMappedField<Scalar, on_cell>& miu_dynamic, cuFvMappedField<Scalar, on_cell>& miu_static, cuFvMappedField<Scalar, on_cell>& ED_rate, Scalar rho_solid, Scalar rho_water, Scalar porosity, Scalar dim_mid){

      cuEDTakahashiIversonXiaKernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(
        gravity.data.dev_ptr(),
        h.data.dev_ptr(),
        hU.data.dev_ptr(),
        hC.data.dev_ptr(),
        manning_coeff.data.dev_ptr(),
        miu_dynamic.data.dev_ptr(),
        miu_static.data.dev_ptr(),
        ED_rate.data.dev_ptr(),
        rho_solid,
        rho_water,
        porosity,
        dim_mid,
        h.data.size()
        );

    }


    __global__ void cuBankCollapseKernel(Scalar* z, Scalar* z_new, Scalar critical_slope, unsigned int size, ShortDualHandle* cell_neigbours, unsigned int cell_neighbours_length, Scalar* cell_volume){
      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      while (index < size){
        auto z_this = z[index];
        Scalar z_new_ = z_this;
        for (Flag i = 0; i < 4; ++i){
          ShortDualHandle neib = cell_neigbours[i*cell_neighbours_length + index];
          Scalar volume = cell_volume[index];
          Scalar length = sqrt(volume);
          if (!neib.is_boundary()){
            Flag id_neib = neib.get_global_id();
            Scalar z_neib = z[id_neib];
            Scalar real_slope = abs(z_neib - z_this) / length;
            if (real_slope > critical_slope){
              Scalar dz = (real_slope - critical_slope)*length;
              if (z_this > z_neib){
                z_new_ -= 0.125*dz;
              }
              else{
                z_new_ += 0.125*dz;
              }
            }
          }
          else{
            //do nothing for boundary cells
          }
        }
        z_new[index] = z_new_;
        index += blockDim.x * gridDim.x;
      }
    }

    void cuBankCollapse(cuFvMappedField<Scalar, on_cell>& z, cuFvMappedField<Scalar, on_cell>& z_new, Scalar critical_slope){

      auto mesh = z.mesh;

      cuBankCollapseKernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(
        z.data.dev_ptr(),
        z_new.data.dev_ptr(),
        critical_slope,
        z.data.size(),
        mesh->cell_neighbours.dev_ptr(),
        mesh->cell_neighbours.length(),
        mesh->cell_volumes.dev_ptr()
        );

    }

  }//--end of namespace fv-------

}//--end of namespace GC---------
