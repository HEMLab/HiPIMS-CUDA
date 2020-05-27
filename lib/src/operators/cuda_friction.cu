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
\file cuda_friction.cu
\brief Source file for friction operator

\version 0.1
\author xilin xia

*/

#include "cuda_friction.h"
#include "cuda_kernel_launch_parameters.h"

namespace GC{

  namespace fv{

    __device__ Vector MCFriction(Scalar& g, Scalar& miu, Scalar& h, Vector& hU, Vector& z_grad){
      Scalar small_error = 1e-10;

      Vector u, f;
      if (h <= small_error){
        u = 0.0;
      }
      else{
        u = hU / h;
      }

      Scalar phi = sqrt(dot(u, u) + pow(dot(u, z_grad), 2)) / sqrt(1.0 + dot(z_grad, z_grad));
      if (phi <= small_error){
        f = 0.0;
      }
      else{
        f = -1.0*g*h*miu*u / phi;
      }
      return f;
    }

    __global__ void cuFrictionMCKernel(Scalar* grav, Scalar* friction_coeff, Scalar* h, Vector* hU, Vector* z_grad, Vector* friction_force, unsigned int size){

      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      while (index < size){
        auto g = grav[index];
        auto miu = friction_coeff[index];
        auto h_ = h[index];
        auto hU_ = hU[index];
        auto z_grad_ = z_grad[index];
        friction_force[index] = MCFriction(g, miu, h_, hU_, z_grad_);
        index += blockDim.x * gridDim.x;
      }
    }

    void cuFrictionMC(cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& friction_coeff, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Vector, on_cell>& z_grad, cuFvMappedField<Vector, on_cell>& friction_force){

      cuFrictionMCKernel <<<BLOCKS_PER_GRID, THREADS_PER_BLOCK >>>(gravity.data.dev_ptr(), friction_coeff.data.dev_ptr(), h.data.dev_ptr(), hU.data.dev_ptr(), z_grad.data.dev_ptr(), friction_force.data.dev_ptr(), friction_force.data.size());

    }

    __global__ void cuFrictionMCBalancedKernel(Scalar dt, Scalar* grav, Scalar* friction_coeff, Scalar* h, Vector* hU, Vector* z_grad, Vector* friction_force, unsigned int size){

      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      while (index < size){
        auto g = grav[index];
        auto miu = friction_coeff[index];
        auto h_ = h[index];
        auto hU_ = hU[index];
        auto z_grad_ = z_grad[index];
        auto _friction_force = MCFriction(g, miu, h_, hU_, z_grad_);
        if (dot(_friction_force*dt, _friction_force*dt) >= dot(hU_, hU_)){
          _friction_force = -1.0 / dt*hU_;
        }
        friction_force[index] = _friction_force;
        index += blockDim.x * gridDim.x;
     
      }
    }

    void cuFrictionMCBalanced(Scalar dt, cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& friction_coeff, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Vector, on_cell>& z_grad, cuFvMappedField<Vector, on_cell>& friction_force){

      cuFrictionMCBalancedKernel <<<BLOCKS_PER_GRID, THREADS_PER_BLOCK >>>(dt, gravity.data.dev_ptr(), friction_coeff.data.dev_ptr(), h.data.dev_ptr(), hU.data.dev_ptr(), z_grad.data.dev_ptr(), friction_force.data.dev_ptr(), friction_force.data.size());

    }

    __global__ void cuFrictionMCWithWallKernel(Scalar dt, Scalar* grav, Scalar friction_coeff1, Scalar friction_coeff2, Scalar wall_width, Scalar* h, Vector* hU, Vector* z_grad, Vector* friction_force, unsigned int size){

      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      while (index < size){
        auto g = grav[index];
        auto h_ = h[index];
        auto hU_ = hU[index];
        auto z_grad_ = z_grad[index];
        auto miu = friction_coeff1 + friction_coeff2*h_ / wall_width / sqrt(1.0 + dot(z_grad_, z_grad_));
		//if (g >= 8.0*9.81 || g <= 0.2*9.81){
        //  printf("miu %d %f %f %f \n", index, miu, h_, g);
        //}
        auto _friction_force = MCFriction(g, miu, h_, hU_, z_grad_);
        if (dot(_friction_force*dt, _friction_force*dt) >= dot(hU_, hU_)){
          _friction_force = -1.0 / dt*hU_;
        }
        friction_force[index] = _friction_force;
        index += blockDim.x * gridDim.x;

      }
    }

    void cuFrictionMCWithWall(Scalar dt, cuFvMappedField<Scalar, on_cell>& gravity, Scalar friction_coeff1, Scalar friction_coeff2, Scalar wall_width, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Vector, on_cell>& z_grad, cuFvMappedField<Vector, on_cell>& friction_force){

      cuFrictionMCWithWallKernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(dt, gravity.data.dev_ptr(), friction_coeff1, friction_coeff2, wall_width, h.data.dev_ptr(), hU.data.dev_ptr(), z_grad.data.dev_ptr(), friction_force.data.dev_ptr(), friction_force.data.size());

    }

    __device__ Vector PlasticFriction(Scalar& tau, Scalar& rho, Scalar& h, Vector& hU, Vector& z_grad){
      Scalar small_error = 1e-10;

      Vector u, f;
      if (h <= small_error){
        u = 0.0;
      }
      else{
        u = hU / h;
      }

      Scalar vel = sqrt(dot(u, u) + pow(dot(u, z_grad), 2));
      if (vel <= small_error){
        f = 0.0;
      }
      else{
        f = -1.0*tau/rho*u / vel;
      }
      return f;
    }

    __global__ void cuFrictionPlasticKernel(Scalar dt, Scalar* retarding_stress, Scalar* rho, Scalar* h, Vector* hU, Vector* z_grad, Vector* friction_force, unsigned int size){

      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      while (index < size){
        auto tau_ = retarding_stress[index];
        auto rho_ = rho[index];
        auto h_ = h[index];
        auto hU_ = hU[index];
        auto z_grad_ = z_grad[index];
        auto _friction_force = PlasticFriction(tau_, rho_, h_, hU_, z_grad_);
        if (dot(_friction_force*dt, _friction_force*dt) >= dot(hU_, hU_)){
          _friction_force = -1.0 / dt*hU_;
        }
        friction_force[index] = _friction_force;
        index += blockDim.x * gridDim.x;

      }
    }

    void cuFrictionPlastic(Scalar dt, cuFvMappedField<Scalar, on_cell>& retarding_stress, cuFvMappedField<Scalar, on_cell>& rho, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Vector, on_cell>& z_grad, cuFvMappedField<Vector, on_cell>& friction_force){

      cuFrictionPlasticKernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(dt, retarding_stress.data.dev_ptr(), rho.data.dev_ptr(), h.data.dev_ptr(), hU.data.dev_ptr(), z_grad.data.dev_ptr(), friction_force.data.dev_ptr(), friction_force.data.size());

    }


    __device__ Vector MCPlasticFriction(Scalar& g, Scalar& miu, Scalar& tau, Scalar& rho, Scalar& h, Vector& hU, Vector& z_grad){
      Scalar small_error = 1e-10;
      Vector u, f;
      if (h <= small_error){
        u = 0.0;
      }
      else{
        u = hU / h;
      }

      Scalar vel = sqrt(dot(u, u) + pow(dot(u, z_grad), 2));
      Scalar phi = sqrt(1.0 + dot(z_grad, z_grad));
      if (vel <= small_error){
        f = 0.0;
      }
      else{
        f = -1.0*tau*phi / rho*u / vel - 1.0*g*h*miu*u*phi / vel;
      }
      return f;
    }

    __global__ void cuFrictionMCPlasticKernel(Scalar dt, Scalar* gravity, Scalar* miu, Scalar* retarding_stress, Scalar* rho, Scalar* h, Vector* hU, Vector* z_grad, Vector* friction_force, unsigned int size){

      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      while (index < size){
        auto g_ = gravity[index];
        auto miu_ = miu[index];
        auto tau_ = retarding_stress[index];
        auto rho_ = rho[index];
        auto h_ = h[index];
        auto hU_ = hU[index];
        auto z_grad_ = z_grad[index];
        auto _friction_force = MCPlasticFriction(g_, miu_, tau_, rho_, h_, hU_, z_grad_);
        if (dot(_friction_force*dt, _friction_force*dt) >= dot(hU_, hU_)){
          _friction_force = -1.0 / dt*hU_;
        }
        friction_force[index] = _friction_force;
        index += blockDim.x * gridDim.x;

      }
    }


    void cuFrictionMCPlastic(Scalar dt, cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& miu, cuFvMappedField<Scalar, on_cell>& retarding_stress, cuFvMappedField<Scalar, on_cell>& rho, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Vector, on_cell>& z_grad, cuFvMappedField<Vector, on_cell>& friction_force){

      cuFrictionMCPlasticKernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(dt, gravity.data.dev_ptr(), miu.data.dev_ptr(), retarding_stress.data.dev_ptr(), rho.data.dev_ptr(), h.data.dev_ptr(), hU.data.dev_ptr(), z_grad.data.dev_ptr(), friction_force.data.dev_ptr(), friction_force.data.size());

    }

    __device__ Vector MuIFriction(Scalar& g, Scalar& miu_1, Scalar& miu_2, Scalar& beta, Scalar& L, Scalar& h, Vector& hU, Vector& z_grad){
      Scalar small_error = 1e-10;
      Vector u, f;
      if (h <= small_error){
        u = 0.0;
      }
      else{
        u = hU / h;
      }
      Scalar miu = 0.0;
      Scalar vel = sqrt(dot(u, u) + pow(dot(u, z_grad), 2));
      Scalar phi = sqrt(1.0 + dot(z_grad, z_grad));
      if (h <= small_error){
        miu = miu_1;
      }
      else if (vel <= small_error){
        miu = miu_2;
      }
      else{
        Scalar Froude = vel / sqrt(g*h / sqrt(1.0 + dot(z_grad, z_grad)));
        miu = miu_1 + (miu_2 - miu_1) / (beta*h / L / Froude + 1);
      }
      if (vel <= small_error){
        f = 0.0;
      }
      else{
        f = - 1.0*g*h*miu*u*phi / vel;
      }
      return f;
    }

    __global__ void cuFrictionMuIKernel(Scalar dt, Scalar* gravity, Scalar* miu_1, Scalar* miu_2, Scalar beta, Scalar L, Scalar* h, Vector* hU, Vector* z_grad, Vector* friction_force, unsigned int size){

      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      while (index < size){
        auto g_ = gravity[index];
        auto miu_1_ = miu_1[index];
        auto miu_2_ = miu_2[index];
        auto h_ = h[index];
        auto hU_ = hU[index];
        auto z_grad_ = z_grad[index];
        auto _friction_force = MuIFriction(g_, miu_1_, miu_2_, beta, L, h_, hU_, z_grad_);
        if (dot(_friction_force*dt, _friction_force*dt) >= dot(hU_, hU_)){
          _friction_force = -1.0 / dt*hU_;
        }
        friction_force[index] = _friction_force;
        index += blockDim.x * gridDim.x;

      }
    }

    void cuFrictionMuI(Scalar dt, cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& miu_1, cuFvMappedField<Scalar, on_cell>& miu_2, Scalar beta, Scalar L, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Vector, on_cell>& z_grad, cuFvMappedField<Vector, on_cell>& friction_force){

      cuFrictionMuIKernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(dt, gravity.data.dev_ptr(), miu_1.data.dev_ptr(), miu_2.data.dev_ptr(), beta, L, h.data.dev_ptr(), hU.data.dev_ptr(), z_grad.data.dev_ptr(), friction_force.data.dev_ptr(), friction_force.data.size());

    }

    __device__ Vector LucasFriction(Scalar& g, Scalar& miu_1, Scalar& miu_2, Scalar& U_w, Scalar& h, Vector& hU, Vector& z_grad){
      Scalar small_error = 1e-10;
      Vector u, f;
      if (h <= small_error){
        u = 0.0;
      }
      else{
        u = hU / h;
      }

      Scalar vel = sqrt(dot(u, u) + pow(dot(u, z_grad), 2));
      Scalar phi = sqrt(1.0 + dot(z_grad, z_grad));
      Scalar miu = miu_1 + (miu_2 - miu_1) / (1 + vel/U_w);

      if (vel <= small_error){
        f = 0.0;
      }
      else{
        f = -1.0*g*h*miu*u*phi / vel;
      }
      return f;
    }

    __global__ void cuFrictionLucasKernel(Scalar dt, Scalar* gravity, Scalar* miu_1, Scalar* miu_2, Scalar U_w, Scalar* h, Vector* hU, Vector* z_grad, Vector* friction_force, unsigned int size){

      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      while (index < size){
        auto g_ = gravity[index];
        auto miu_1_ = miu_1[index];
        auto miu_2_ = miu_2[index];
        auto h_ = h[index];
        auto hU_ = hU[index];
        auto z_grad_ = z_grad[index];
        auto _friction_force = LucasFriction(g_, miu_1_, miu_2_, U_w, h_, hU_, z_grad_);
        if (dot(_friction_force*dt, _friction_force*dt) >= dot(hU_, hU_)){
          _friction_force = -1.0 / dt*hU_;
        }
        friction_force[index] = _friction_force;
        index += blockDim.x * gridDim.x;

      }
    }

    void cuFrictionLucas(Scalar dt, cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& miu_1, cuFvMappedField<Scalar, on_cell>& miu_2, Scalar U_w, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Vector, on_cell>& z_grad, cuFvMappedField<Vector, on_cell>& friction_force){

        cuFrictionLucasKernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(dt, gravity.data.dev_ptr(), miu_1.data.dev_ptr(), miu_2.data.dev_ptr(), U_w, h.data.dev_ptr(), hU.data.dev_ptr(), z_grad.data.dev_ptr(), friction_force.data.dev_ptr(), friction_force.data.size());

    }

    __device__ Vector2 ManningNewton(Vector2 S_b, Scalar C_f, Scalar dt, Vector2 U){
      Vector2 S_f = -C_f*norm(U)*U;
      if (norm(S_b + S_f) <= 1e-10){ //steady state, return directly
        return U;
      }
      Scalar epsilon = 0.001; //termination criteria
      Vector2 U_k = U;
      Vector2 U_k1 = U;
      unsigned int k = 0;
      while (true){
        Tensor2 inv;
        if (norm(U_k) <= 1e-10){ //Jacobian matrix is 0
          inv = Tensor2(1.0, 0.0, 0.0, 1.0);
        }
        else{
          Scalar J_xx = -C_f*(2.0*U_k.x*U_k.x + U_k.y*U_k.y) / norm(U_k);
          Scalar J_yy = -C_f*(U_k.x*U_k.x + 2.0*U_k.y*U_k.y) / norm(U_k);
          Scalar J_xy = -C_f*(U_k.x*U_k.y) / norm(U_k);
          Tensor2 ijac = Tensor2(1.0 - dt*J_xx, -dt*J_xy, -dt*J_xy, 1.0 - dt*J_yy);
          inv = inverse(ijac);
        }
        Vector2 S = S_b - C_f*norm(U_k)*U_k;
        U_k1 = U_k + dot(inv, dt*S - U_k + U);
        if (norm(U_k1 - U_k) <= epsilon* norm(U_k) || k > 10){
          break;
        }
        U_k = U_k1;
        k++;
      }
      return U_k1;
    }

    __global__ void cuFrictionManningImplicitKernel(Scalar* gravity, Scalar* manning_coeff, Scalar* h, Vector* hU, Vector* hU_advection, Scalar dt, unsigned int size){

      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      while (index < size){
        auto g_ = gravity[index];
        auto manning_ = manning_coeff[index];
        auto h_ = h[index];
        auto hU_ = hU[index];
        auto hU_advection_ = hU_advection[index];
        Scalar h_small = 1e-10;
        if (h_ < h_small){
          hU[index] = hU[index] +dt*hU_advection_;
        }
        else{
          //updating velocity rather than discharge
// Iterative fully implicit scheme
          Scalar C_f = g_*manning_*manning_*pow(h_, -4.0 / 3.0);
          Vector2 U_ = hU_ / h_;
          Vector2 acc_ = hU_advection_ / h_;
          Vector2 U_1 = ManningNewton(acc_, C_f, dt, U_);
          hU[index] = U_1*h_;
// Non-iterative fully-implicit scheme
 /*       Scalar C_f = g_*manning_*manning_*pow(h_, -4.0 / 3.0);
          Vector2 m_ = hU[index] + dt*hU_advection_;
          Vector2 U_ = m_ / h_;
          Scalar vel = norm(U_);
          Scalar alpha = (1.0 - sqrt(1.0 + 4.0*dt*C_f*vel)) / (1e-10 - 2.0*dt*C_f*vel);
          hU[index] = alpha*m_; */
// Point implicit version 1 (Song et al. 2011)
          /* Scalar C_f = g_*manning_*manning_*pow(h_, -4.0 / 3.0);
          Vector2 m_ = hU[index] + dt*hU_advection_;
          Vector2 U_ = hU_ / h_;
          Scalar vel = norm(U_);
          Scalar alpha = 1.0 / (1 + dt*C_f*vel);
          hU[index] = alpha*m_; */
          //Scalar C_f = g_*manning_*manning_*pow(h_, -1.0 / 3.0);
          //Scalar D_x = 0.0, D_y = 0.0;
          //Vector2 U_ = hU_ / h_;
          //if (norm(U_) > 1e-10){
          //  D_x = dt*C_f / (h_)*norm(U_);
          //  D_y = dt*C_f / (h_)*norm(U_);
          //}
          //Vector2 hU_1 = 0.0;
          //hU_1.x = (hU_.x + dt*hU_advection_.x) / (1 + D_x);
          //hU_1.y = (hU_.y + dt*hU_advection_.y) / (1 + D_y);
          //hU[index] = hU_1;
// Point Implicit version 2
/*        Scalar C_f = g_*manning_*manning_*pow(h_, -1.0 / 3.0);
          Scalar D_x = 0.0, D_y = 0.0;
          Vector2 U_ = hU_ / h_;          
          if (norm(U_) > 1e-10){
            D_x = dt*C_f / (h_)*(2.0*U_.x*U_.x + U_.y* U_.y) / norm(U_);
            D_y = dt*C_f / (h_)*(U_.x*U_.x + 2.0*U_.y* U_.y) / norm(U_);
          }
          Vector2 hU_1 = 0.0;
          Scalar S_fx = 0.0;
          Scalar S_fy = 0.0;
          S_fx = -C_f *U_.x*norm(U_);
          S_fy = -C_f *U_.y*norm(U_);
          hU_1.x = (dt*hU_advection_.x + dt*S_fx )/ (1 + D_x);
          hU_1.y = (dt*hU_advection_.y + dt*S_fy )/ (1 + D_y);
          hU[index] = hU_ + hU_1;  */
// Point Implicit Liang
/*          hU_ = hU[index] +dt*hU_advection_;
          Scalar C_f = g_*manning_*manning_*pow(h_, -1.0 / 3.0);
          Scalar D_x = 0.0, D_y = 0.0;
          Vector2 U_ = hU_ / h_;          
          if (norm(U_) > 1e-10){
            D_x = dt*C_f / (h_)*(2.0*U_.x*U_.x + U_.y* U_.y) / norm(U_);
            D_y = dt*C_f / (h_)*(U_.x*U_.x + 2.0*U_.y* U_.y) / norm(U_);
          }
          Scalar S_fx = 0.0;
          Scalar S_fy = 0.0;
          S_fx = -C_f *U_.x*norm(U_);
          S_fy = -C_f *U_.y*norm(U_);
          hU_.x = hU_.x + dt*S_fx/ (1 + D_x);
          hU_.y = hU_.y + dt*S_fy/ (1 + D_y);
          hU[index] = hU_;  */
        }
        index += blockDim.x * gridDim.x;
      }       
    }

    void cuFrictionManningImplicit(Scalar dt, cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& manning_coeff, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Vector, on_cell>& hU_advection){

      cuFrictionManningImplicitKernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(gravity.data.dev_ptr(), manning_coeff.data.dev_ptr(), h.data.dev_ptr(), hU.data.dev_ptr(), hU_advection.data.dev_ptr(), dt, h.data.size());
    
    }

    __global__ void cuFrictionManningImplicitWithForceKernel(Scalar* gravity, Scalar* manning_coeff, Scalar* h, Vector* hU, Vector* hU_advection, Vector* friction_force, Scalar dt, unsigned int size){

      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      while (index < size){
        auto g_ = gravity[index];
        auto manning_ = manning_coeff[index];
        auto h_ = h[index];
        auto hU_ = hU[index];
        auto hU_advection_ = hU_advection[index];
        Scalar h_small = 1e-10;
        if (h_ < h_small){
          hU[index] = hU[index]; // +dt*hU_advection_;
        }
        else{
          //updating velocity rather than discharge
          Scalar C_f = g_*manning_*manning_*pow(h_, -4.0 / 3.0);
          Vector2 U_ = hU_ / h_;
          Vector2 acc_ = hU_advection_ / h_;
          Vector2 U_1 = ManningNewton(acc_, C_f, dt, U_);
          hU[index] = U_1*h_;
        }
        friction_force[index] = (hU[index] - hU_) / dt - hU_advection_;
        //if (index == 379938){
        //  printf("U0x %f U0y %f U1x %f U1y %f frictionx %f frictiony %f\n", hU_.x, hU_.y, hU[index].x, hU[index].y, friction_force[index].x, friction_force[index].y);
        //  printf("advectionx %f advectiony%f\n", hU_advection_.x, hU_advection_.y);
        //}
        index += blockDim.x * gridDim.x;
      }
    }

    void cuFrictionManningImplicitWithForce(Scalar dt, cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& manning_coeff, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Vector, on_cell>& hU_advection, cuFvMappedField<Vector, on_cell>& friction_force){

      cuFrictionManningImplicitWithForceKernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(gravity.data.dev_ptr(), manning_coeff.data.dev_ptr(), h.data.dev_ptr(), hU.data.dev_ptr(), hU_advection.data.dev_ptr(), friction_force.data.dev_ptr(), dt, h.data.size());

    }

    __global__ void cuFrictionManningExplicitKernel(Scalar* gravity, Scalar* manning_coeff, Scalar* h, Vector* hU, Vector* hU_advection, Scalar dt, unsigned int size){

      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      while (index < size){
        auto g_ = gravity[index];
        auto manning_ = manning_coeff[index];
        auto h_ = h[index];
        auto hU_ = hU[index];
        auto hU_advection_ = hU_advection[index];
        Scalar h_small = 1e-10;
        if (h_ < h_small){
          hU[index] = hU[index]; // +dt*hU_advection_;
        }
        else{
          //updating velocity rather than discharge
/*          Scalar C_f = g_*manning_*manning_*pow(h_, -4.0 / 3.0);
          Vector2 U_ = hU_ / h_;
          Vector2 acc_ = hU_advection_ / h_;
          Vector2 U_1 = ManningNewton(acc_, C_f, dt, U_);
          hU[index] = U_1*h_;  */
// Point Implicit
          Scalar C_f = g_*manning_*manning_*pow(h_, -1.0 / 3.0);
          Scalar D_x = 0.0, D_y = 0.0;
          Vector2 U_ = hU_ / h_;          
          if (norm(U_) > 1e-10){
            D_x = dt*C_f / (h_)*(U_.x*U_.x + U_.y* U_.y) / norm(U_);
            D_y = dt*C_f / (h_)*(U_.x*U_.x + U_.y* U_.y) / norm(U_);
          }
          Vector2 hU_1 = 0.0;
          Scalar S_fx = 0.0;
          Scalar S_fy = 0.0;
          S_fx = -C_f *U_.x*norm(U_);
          S_fy = -C_f *U_.y*norm(U_);
          hU_1.x = (hU_.x + dt*hU_advection_.x )/ (1 + D_x);
          hU_1.y = (hU_.y + dt*hU_advection_.y )/ (1 + D_y);
          hU[index] = hU_1;  
// Point Implicit Liang
/*          hU_ = hU[index] +dt*hU_advection_;
          Scalar C_f = g_*manning_*manning_*pow(h_, -1.0 / 3.0);
          Scalar D_x = 0.0, D_y = 0.0;
          Vector2 U_ = hU_ / h_;          
          if (norm(U_) > 1e-10){
            D_x = dt*C_f / (h_)*(2.0*U_.x*U_.x + U_.y* U_.y) / norm(U_);
            D_y = dt*C_f / (h_)*(U_.x*U_.x + 2.0*U_.y* U_.y) / norm(U_);
          }
          Scalar S_fx = 0.0;
          Scalar S_fy = 0.0;
          S_fx = -C_f *U_.x*norm(U_);
          S_fy = -C_f *U_.y*norm(U_);
          hU_.x = hU_.x + dt*S_fx/ (1 + D_x);
          hU_.y = hU_.y + dt*S_fy/ (1 + D_y);
          hU[index] = hU_;  */
        }
        index += blockDim.x * gridDim.x;
      }       
    }

    void cuFrictionManningExplicit(Scalar dt, cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& manning_coeff, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Vector, on_cell>& hU_advection){

      cuFrictionManningExplicitKernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(gravity.data.dev_ptr(), manning_coeff.data.dev_ptr(), h.data.dev_ptr(), hU.data.dev_ptr(), hU_advection.data.dev_ptr(), dt, h.data.size());
    
    }

    __global__ void cuFrictionManningSplittingKernel(Scalar* gravity, Scalar* manning_coeff, Scalar* h, Vector* hU, Vector* hU_advection, Scalar dt, unsigned int size){

      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      while (index < size){
        auto g_ = gravity[index];
        auto manning_ = manning_coeff[index];
        auto h_ = h[index];
        auto hU_ = hU[index];
        auto hU_advection_ = hU_advection[index];
        Scalar h_small = 1e-10;
        if (h_ < h_small){
          hU[index] = hU[index]; // +dt*hU_advection_;
        }
        else{
          //updating velocity rather than discharge
/*          Scalar C_f = g_*manning_*manning_*pow(h_, -4.0 / 3.0);
          Vector2 U_ = hU_ / h_;
          Vector2 acc_ = hU_advection_ / h_;
          Vector2 U_1 = ManningNewton(acc_, C_f, dt, U_);
          hU[index] = U_1*h_;  */
// Point Implicit
/*          Scalar C_f = g_*manning_*manning_*pow(h_, -1.0 / 3.0);
          Scalar D_x = 0.0, D_y = 0.0;
          Vector2 U_ = hU_ / h_;          
          if (norm(U_) > 1e-10){
            D_x = dt*C_f / (h_)*(U_.x*U_.x + U_.y* U_.y) / norm(U_);
            D_y = dt*C_f / (h_)*(U_.x*U_.x + U_.y* U_.y) / norm(U_);
          }
          Vector2 hU_1 = 0.0;
          Scalar S_fx = 0.0;
          Scalar S_fy = 0.0;
          S_fx = -C_f *U_.x*norm(U_);
          S_fy = -C_f *U_.y*norm(U_);
          hU_1.x = (dt*hU_advection_.x + dt*S_fx )/ (1 + D_x);
          hU_1.y = (dt*hU_advection_.y + dt*S_fy )/ (1 + D_y);
          hU[index] = hU_ + hU_1;  */
// Point Implicit Liang
          hU_ = hU[index] +dt*hU_advection_;
          Scalar C_f = g_*manning_*manning_*pow(h_, -1.0 / 3.0);
          Scalar D_x = 0.0, D_y = 0.0;
          Vector2 U_ = hU_ / h_;          
          if (norm(U_) > 1e-10){
            D_x = dt*C_f / (h_)*(2.0*U_.x*U_.x + U_.y* U_.y) / norm(U_);
            D_y = dt*C_f / (h_)*(U_.x*U_.x + 2.0*U_.y* U_.y) / norm(U_);
          }
          Scalar S_fx = 0.0;
          Scalar S_fy = 0.0;
          S_fx = -C_f *U_.x*norm(U_);
          S_fy = -C_f *U_.y*norm(U_);
          hU_.x = hU_.x + dt*S_fx/ (1 + D_x);
          hU_.y = hU_.y + dt*S_fy/ (1 + D_y);
          hU[index] = hU_;  
        }
        index += blockDim.x * gridDim.x;
      }       
    }

    void cuFrictionManningSplitting(Scalar dt, cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& manning_coeff, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Vector, on_cell>& hU_advection){

      cuFrictionManningSplittingKernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(gravity.data.dev_ptr(), manning_coeff.data.dev_ptr(), h.data.dev_ptr(), hU.data.dev_ptr(), hU_advection.data.dev_ptr(), dt, h.data.size());
    
    }

    __device__ Vector2 ManningSteepNewton(Vector2 S_b, Scalar C_f, Scalar dt, Vector2 U, Vector2 z_grad){
      Vector2 S_f = -C_f*norm(U)*U;
      if (norm(S_b + S_f) <= 1e-10){ //steady state, return directly
        return U;
      }
      Scalar epsilon = 0.001; //termination criteria
      Vector2 U_k = U;
      Vector2 U_k1 = U;
      unsigned int k = 0;
      while (true){
        Tensor2 inv;
        Vector3 U3_k(U_k.x, U_k.y, dot(U_k, z_grad));
        if (norm(U3_k) <= 1e-10){ //Jacobian matrix is 0
          inv = Tensor2(1.0, 0.0, 0.0, 1.0);
        }
        else{
          Scalar J_xx = -C_f*(2.0*U_k.x*U_k.x + U_k.y*U_k.y + U3_k.z*U3_k.z + U_k.x*z_grad.x*U3_k.z) / norm(U3_k);
          Scalar J_xy = -C_f*(U_k.x*U_k.y + U_k.x*z_grad.y*U3_k.z) / norm(U3_k);
          Scalar J_yx = -C_f*(U_k.x*U_k.y + U_k.y*z_grad.x*U3_k.z) / norm(U3_k);
          Scalar J_yy = -C_f*(U_k.x*U_k.x + 2.0*U_k.y*U_k.y + U3_k.z*U3_k.z + U_k.y*z_grad.y*U3_k.z) / norm(U3_k);
          Tensor2 ijac = Tensor2(1.0 - dt*J_xx, -dt*J_xy, -dt*J_yx, 1.0 - dt*J_yy);
          inv = inverse(ijac);
        }
        Vector2 S = S_b - C_f*norm(U3_k)*U_k;
        U_k1 = U_k + dot(inv, dt*S - U_k + U);
        if (norm(U_k1 - U_k) <= epsilon* norm(U_k) || k > 10){
          break;
        }
        U_k = U_k1;
        k++;
      }
      return U_k1;
    }

    __global__ void cuFrictionManningSteepImplicitKernel(Scalar* gravity, Scalar* manning_coeff, Scalar* h, Vector* hU, Vector* hU_advection, Vector* z_grad, Scalar dt, unsigned int size){

      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      while (index < size){
        auto g_ = gravity[index];
        auto manning_ = manning_coeff[index];
        auto h_ = h[index];
        auto hU_ = hU[index];
        auto hU_advection_ = hU_advection[index];
        auto z_grad_ = z_grad[index];
        Scalar phi = sqrt(1.0 / (1.0 + dot(z_grad_,z_grad_))); //cosine theta
        Scalar h_small = 1e-10;
        if (h_ < h_small){
          hU[index] = hU[index]; // +dt*hU_advection_;
        }
        else{
          //updating velocity rather than discharge
          Scalar C_f = g_*manning_*manning_*pow(h_, -4.0 / 3.0)*pow(phi, -1.0/3.0);
          Vector2 U_ = hU_ / h_;
          Vector2 acc_ = hU_advection_ / h_;
          Vector2 U_1 = ManningSteepNewton(acc_, C_f, dt, U_, z_grad_);
          hU[index] = U_1*h_;
        }
        index += blockDim.x * gridDim.x;
      }
    }
    
    void cuFrictionManningSteepImplicit(Scalar dt, cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& manning_coeff, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Vector, on_cell>& hU_advection, cuFvMappedField<Vector, on_cell>& z_grad){

      cuFrictionManningSteepImplicitKernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(gravity.data.dev_ptr(), manning_coeff.data.dev_ptr(), h.data.dev_ptr(), hU.data.dev_ptr(), hU_advection.data.dev_ptr(), z_grad.data.dev_ptr(), dt, h.data.size());

    }

    __global__ void cuFrictionManningMCImplicitKernel(Scalar* gravity, Scalar* manning_coeff, Scalar* friction_coeff, Scalar* h, Scalar* hC, Vector* hU, Vector* hU_advection, Scalar rho_water, Scalar rho_solid, Scalar dt, unsigned int size){

      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      while (index < size){
        auto g_ = gravity[index];
        auto manning_ = manning_coeff[index];
        auto miu_ = friction_coeff[index];
        auto h_ = h[index];
        auto hU_ = hU[index];
        auto hC_ = hC[index];
        auto hU_advection_ = hU_advection[index];
        Scalar h_small = 1e-10;
        if (h_ < h_small){
          hU[index] = hU[index]; // +dt*hU_advection_;
        }
        else{
          // ---- Manning step ----
          Scalar C_ = hC_ / h_;
          Scalar rho_mix = rho_solid*C_ + rho_water*(1.0 - C_);
          //updating velocity rather than discharge
          Scalar C_f = g_*manning_*manning_*pow(h_, -4.0 / 3.0);
          Vector2 U_ = hU_ / h_;
          Vector2 acc_ = hU_advection_ / h_;
          Vector2 U_1 = ManningNewton(acc_, C_f, dt, U_);
          hU_ = U_1*h_;
          // ---- Mohr-Coulomb step ----
          Vector2 _friction_force = 0.0;
          U_ = hU_ / h_;
          if (norm(U_) <= 1e-6){
            _friction_force = 0.0;
          }
          else{
            _friction_force = -1.0*(rho_solid - rho_water) / rho_mix*g_*hC_*miu_*U_ / norm(U_);
          }
          ////constrain friction force
          if (dot(_friction_force*dt, _friction_force*dt) >= dot(hU_, hU_)){
            _friction_force = -1.0 / dt*hU_;
          }
          hU[index] = hU_ + _friction_force*dt;
        }
        index += blockDim.x * gridDim.x;
      }
    }

    void cuFrictionManningMCImplicit(Scalar dt, Scalar rho_water, Scalar rho_solid, cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& manning_coeff, cuFvMappedField<Scalar, on_cell>& friction_coeff, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Scalar, on_cell>& hC, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Vector, on_cell>& hU_advection){

      cuFrictionManningMCImplicitKernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> > (
        gravity.data.dev_ptr(), 
        manning_coeff.data.dev_ptr(), 
        friction_coeff.data.dev_ptr(),
        h.data.dev_ptr(), 
        hC.data.dev_ptr(),
        hU.data.dev_ptr(), 
        hU_advection.data.dev_ptr(),
        rho_water,
        rho_solid,
        dt, 
        h.data.size());

    }

  }//end of namespace fv

}//end of namespace GC


