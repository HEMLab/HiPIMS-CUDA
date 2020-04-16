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
  \flie cuda_transport_NSWEs.cu
  \brief Source file for advection of non hydrostatic shallow water equations

  \version 0.1
  \author xilin xia
*/

#include "cuda_transport_NSWEs.h"
#include "cuda_kernel_launch_parameters.h"
//#include "cuda_boundary.h"
#include "riemann.h"


namespace GC{

  namespace fv{  
  

  __device__ RiemannFluxSWEs cuHLLCRiemannSolverSWEsTransport(const Scalar& g,const ScalarRiemannState& h,const VectorRiemannState& q){
      Scalar h_L = h.L;
      Scalar h_R = h.R;
      Scalar qx_L = q.L.x;
      Scalar qy_L = q.L.y;
      Scalar qx_R = q.R.x;
      Scalar qy_R = q.R.y;
      Scalar u_L, u_R, v_L, v_R;
      Scalar h_small = 1.0e-10;

      if(h_L <= h_small && h_R <= h_small){
        return RiemannFluxSWEs(0.0, Vector2(0.0, 0.0));
      }
      
      if(h_L <= h_small){
        u_L = 0.0;
        v_L = 0.0;
      }else{
        u_L = qx_L/h_L;
        v_L = qy_L/h_L;
      }

      if(h_R <= h_small){
        u_R = 0.0;
        v_R = 0.0;
      }else{
        u_R = qx_R/h_R;
        v_R = qy_R/h_R;
      }
      
      Scalar a_L = sqrt(g*h_L);
      Scalar a_R = sqrt(g*h_R);
      
      Scalar h_star = pow((a_L + a_R)/2.0 + (u_L - u_R)/4.0, 2.0)/g;
      Scalar u_star = (u_L + u_R)/2.0 + a_L - a_R;
      Scalar a_star = sqrt(g*h_star);

      Scalar s_L, s_R;
      if(h_L <= h_small){
        s_L = u_R - 2.0*a_R;
      }else{
        s_L = fmin(u_L - a_L, u_star - a_star);
      }
      
      if(h_R <= h_small){
        s_R = u_L + 2.0*a_L;
      }else{
        s_R = fmax(u_R + a_R, u_star + a_star);
      }

      Scalar s_M = (s_L*h_R*(u_R - s_R) - s_R*h_L*(u_L - s_L))/(h_R*(u_R - s_R) - h_L*(u_L - s_L));
    
      Scalar h_flux_L = qx_L;
      Scalar qx_flux_L = u_L*qx_L + 0.5*g*h_L*h_L;
      Scalar qy_flux_L = u_L*qy_L;

      Scalar h_flux_R = qx_R;
      Scalar qx_flux_R = u_R*qx_R + 0.5*g*h_R*h_R;
      Scalar qy_flux_R = u_R*qy_R;

      Scalar h_flux_M = (s_R*h_flux_L - s_L*h_flux_R + s_L*s_R*(h_R - h_L))/(s_R - s_L);
      Scalar qx_flux_M = (s_R*qx_flux_L - s_L*qx_flux_R + s_L*s_R*(qx_R - qx_L))/(s_R - s_L);

      Scalar h_flux, qx_flux, qy_flux;
      if(0.0 <= s_L){
        h_flux = h_flux_L;
        qx_flux = qx_flux_L;
        qy_flux = qy_flux_L;
      }else if(s_L <= 0.0 && 0.0 <= s_M){
        h_flux = h_flux_M;
        qx_flux = qx_flux_M;
        qy_flux = h_flux_M*v_L; 
      }else if(s_M <= 0.0 && 0.0 <= s_R){
        h_flux = h_flux_M;
        qx_flux = qx_flux_M;
        qy_flux = h_flux_M*v_R; 
      }else{
        h_flux = h_flux_R;
        qx_flux = qx_flux_R;
        qy_flux = qy_flux_R;
      }
      return RiemannFluxSWEs(h_flux, Vector2(qx_flux, qy_flux));
    };


      __global__ void cuTransportNSWEsSRMCartesianKernel(Scalar* gravity, Scalar* h, Scalar* _h_bound, Scalar* z, Scalar* _z_bound, Vector* z_gradient, Vector* hU, Vector* _hU_bound, Scalar* hC, Scalar* _hC_bound, unsigned int phi_size, ShortDualHandle* cell_neigbours, unsigned int cell_neighbours_length, Scalar* cell_volume, Scalar* h_advection, Vector* hU_advection, Scalar* hC_advection){

      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      Scalar h_small = 1e-10;
      Vector2 face_normal[4];
      Vector2 face_shear[4];
      face_normal[0] = Vector2(0.0, -1.0);
      face_normal[1] = Vector2(1.0, 0.0);
      face_normal[2] = Vector2(0.0, 1.0);
      face_normal[3] = Vector2(-1.0, 0.0);
      face_shear[0] = Vector2(1.0, 0.0);
      face_shear[1] = Vector2(0.0, 1.0);
      face_shear[2] = Vector2(-1.0, 0.0);
      face_shear[3] = Vector2(0.0, -1.0);
      while (index < phi_size){
        Scalar g_this = gravity[index];
        Scalar h_this = h[index];
        Scalar z_this = z[index];
        Vector2 _z_gradient_this = z_gradient[index];
        Scalar eta_this = h_this + z_this;
        Vector2 hU_this = hU[index];
        Vector2 u_this = 0.0;
        Scalar hC_this = hC[index];
        Scalar c_this = 0.0;
        if (h_this < h_small){
          u_this = 0.0;
          c_this = 0.0;
        }
        else{
          u_this = hU_this / h_this;
          c_this = hC_this / h_this;
        }
        Scalar volume = cell_volume[index];
        Scalar area = sqrt(volume);
        Scalar _h_advection(0.0);
        Vector2 _hU_advection(0.0, 0.0);
        Scalar _hC_advection(0.0);
        for (Flag i = 0; i < 4; ++i){
          Vector2 normal = face_normal[i];
          Vector2 shear = face_shear[i];
          ShortDualHandle neib = cell_neigbours[i*cell_neighbours_length + index];
          Scalar g_neib = 0.0;
          Scalar h_neib = 0.0;
          Scalar z_neib = 0.0;
          Vector2 hU_neib = 0.0;
          Vector2 _z_gradient_neib = 0.0;
          Scalar hC_neib = 0.0;
          if (!neib.is_boundary()){
            Flag id_neib = neib.get_global_id();
            g_neib = gravity[id_neib];
            h_neib = h[id_neib];
            z_neib = z[id_neib];
            hU_neib = hU[id_neib];
            hC_neib = hC[id_neib];
            _z_gradient_neib = z_gradient[id_neib];
          }
          else{
            Flag id_boundary = neib.get_global_id();
            g_neib = g_this;
            h_neib = _h_bound[id_boundary];
            z_neib = _z_bound[id_boundary];
            //z_neib = z_this; //This is important for accuracy at outlet
            hU_neib = _hU_bound[id_boundary];
            hC_neib = _hC_bound[id_boundary];
          }
          if (h_this < h_small && h_neib < h_small){
            continue;
          }
          Scalar eta_neib = z_neib + h_neib;
          Vector2 u_neib = 0.0;
          Scalar c_neib = 0.0;
          if (h_neib < h_small){
            u_neib = 0.0;
            c_neib = 0.0;
          }
          else{
            u_neib = hU_neib / h_neib;
            c_neib = hC_neib / h_neib;
          }
          Vector2 direction_this = 0.5*normal*area;
          Vector2 direction_neib = -0.5*normal*area;
          Scalar _z_this = z_this + dot(_z_gradient_this, direction_this);
          Scalar _z_neib = z_neib + dot(_z_gradient_neib, direction_neib);
          Scalar z_f = fmax(z_this, z_neib);
          Scalar delta_z = 0.0;
          Scalar dz_clip = _z_neib - _z_this;
          if (neib.is_boundary()){
            dz_clip = 0.0;
          }
          Scalar dz = z_neib - z_this - dz_clip;
          Scalar deta_this = fmax(0.0, fmin(dz, eta_neib - eta_this));
          Scalar deta_neib = fmax(0.0, fmin(-dz, - eta_neib + eta_this));
          Scalar eta_L = eta_this + deta_this;
          Scalar eta_R = eta_neib + deta_neib;
          Scalar h_L = fmax(0.0, eta_L - z_f);
          Scalar h_R = fmax(0.0, eta_R - z_f);
          Vector2 u_L(dot(u_this, normal), dot(u_this, shear));
          Vector2 u_R(dot(u_neib, normal), dot(u_neib, shear));
          Scalar g = 0.5*(g_this + g_neib);
          auto flux = cuHLLCRiemannSolverSWEsTransport(g, ScalarRiemannState(h_L, h_R), VectorRiemannState(h_L*u_L, h_R*u_R));
          Scalar _h_flux = flux.h;
          Scalar _hC_flux = 0.0;
          if (_h_flux >= 0.0){
            _hC_flux = _h_flux*c_this;
          }else{
            _hC_flux = _h_flux*c_neib;
          }
          Vector2 _hU_flux = (flux.q.x*normal + flux.q.y*shear);
          if (h_neib < h_small){
            delta_z = fmax(0.0, z_f - eta_this);
          }
          else{
            delta_z = fmax(0.0, fmin(dz_clip, z_f - eta_this));
          }
          z_f -= delta_z;
          Vector2 _z_flux = 0.5*g*(h_L + h_this)*(z_f - z_this)*normal;
          _h_advection += _h_flux*area / volume;
          _hU_advection += (_hU_flux + _z_flux)*area / volume;
          _hC_advection += _hC_flux*area / volume;
          //if (index == 379938){
            //printf("h_i %f, h_j %f, hU_i %f, hU_j %f, slope source %f, advection %f\n", h_this, h_neib, norm(hU_this), norm(hU_neib), norm(_z_flux), norm(_hU_advection));
          //}
        }
        h_advection[index] = _h_advection;
        hU_advection[index] = _hU_advection;
        hC_advection[index] = _hC_advection;
        //if (index == 379938){
        //  printf("advection %f\n", norm(_hU_advection));
        //}
        __syncthreads();
        index += blockDim.x * gridDim.x;
      }
    }

    void cuTransportNSWEsSRMCartesian(cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Scalar, on_cell>& z, cuFvMappedField<Vector, on_cell>& z_gradient, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Scalar, on_cell> hC, cuFvMappedField<Scalar, on_cell>& h_advection, cuFvMappedField<Vector, on_cell>& hU_advection, cuFvMappedField<Scalar, on_cell>& hC_advection){

      auto mesh = h.mesh;

      cuTransportNSWEsSRMCartesianKernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(gravity.data.dev_ptr(),
        h.data.dev_ptr(),
        h.boundary_value.dev_ptr(),
        z.data.dev_ptr(),
        z.boundary_value.dev_ptr(),
        z_gradient.data.dev_ptr(),
        hU.data.dev_ptr(),
        hU.boundary_value.dev_ptr(),
        hC.data.dev_ptr(),
        hC.boundary_value.dev_ptr(),
        h.data.size(),
        mesh->cell_neighbours.dev_ptr(),
        mesh->cell_neighbours.length(),
        mesh->cell_volumes.dev_ptr(),
        h_advection.data.dev_ptr(),
        hU_advection.data.dev_ptr(),
        hC_advection.data.dev_ptr());
    }  

  }

}
