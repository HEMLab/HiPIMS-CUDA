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
 \file riemann.cc
 \brief Source file for riemann solvers 

  \version 0.1
  \author xilin xia
*/

#include "riemann.h"
#include <cmath>
#include <algorithm>
#include "Scalar.h"

namespace GC{
  // namespace for finite volume method
  namespace fv{
    RiemannFluxSWEs hllcRiemannSolverSWEs(const Scalar& g,const ScalarRiemannState& h,const VectorRiemannState& q){
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
        s_L = std::min(u_L - a_L, u_star - a_star);
      }
      
      if(h_R <= h_small){
        s_R = u_L + 2.0*a_L;
      }else{
        s_R = std::max(u_R + a_R, u_star + a_star);
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

    RiemannFluxSWEs hllcRiemannSolverNonUniGravitySWEs(const ScalarRiemannState& g,const ScalarRiemannState& h,const VectorRiemannState& q){
      Scalar g_L = g.L;
      Scalar g_R = g.R;
      Scalar h_L = h.L;
      Scalar h_R = h.R;
      Scalar qx_L = q.L.x;
      Scalar qy_L = q.L.y;
      Scalar qx_R = q.R.x;
      Scalar qy_R = q.R.y;
      Scalar u_L, u_R, v_L, v_R;
      Scalar h_small = 1.0e-10;
      Scalar grav = 9.81;

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
      
      Scalar a_L = sqrt(g_L*h_L);
      Scalar a_R = sqrt(g_L*h_R);
      
      Scalar h_star = pow((a_L + a_R)/2.0 + (u_L - u_R)/4.0, 2.0)/grav;
      Scalar u_star = (u_L + u_R)/2.0 + a_L - a_R;
      Scalar a_star = sqrt(grav*h_star);

      Scalar s_L, s_R;
      if(h_L <= h_small){
        s_L = u_R - 2.0*a_R;
      }else{
        s_L = std::min(u_L - a_L, u_star - a_star);
      }
      
      if(h_R <= h_small){
        s_R = u_L + 2.0*a_L;
      }else{
        s_R = std::max(u_R + a_R, u_star + a_star);
      }

      Scalar s_M = (s_L*h_R*(u_R - s_R) - s_R*h_L*(u_L - s_L))/(h_R*(u_R - s_R) - h_L*(u_L - s_L));
    
      Scalar h_flux_L = qx_L;
      Scalar qx_flux_L = u_L*qx_L + 0.5*g_L*h_L*h_L;
      Scalar qy_flux_L = u_L*qy_L;

      Scalar h_flux_R = qx_R;
      Scalar qx_flux_R = u_R*qx_R + 0.5*g_R*h_R*h_R;
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

  }//end of namespace fv

}//end of namespace GC


