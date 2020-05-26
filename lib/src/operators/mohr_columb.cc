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
\file friction.h
\brief Source file for Mohr Columb friction formula

\version 0.1
\author xilin xia

*/

#include <cmath>
#include "mohr_columb.h"


namespace GC{

  namespace fv{

    Vector MCFriction(Scalar& g, Scalar& miu, Scalar& h, Vector& hU, Vector& z_grad){
      Scalar small_error = 1e-10;

      Vector u,f;
      if(h <= small_error){
        u = 0.0;
      }else{
        u = hU/h;
      }

      Scalar phi = sqrt(dot(u, u) + pow(dot(u, z_grad), 2))/sqrt(1.0 + dot(z_grad, z_grad));
      if (phi <= small_error){
        f = 0.0;
      }
      else{
        f = -1.0*g*h*miu*u / phi;
      }
      return f;
      //if (u.x == 0.0){
      //  f.x = miu*g*h;
      //  if(a.x > 0.0){
      //    f.x = -f.x;
      //  }else if(a.x == 0.0){
      //    f.x = 0.0;
      //  }
      //}else if (u.x > 0.0){
      //  f.x = -miu*g*h*u.x/norm(u);
      //}else{
      //  f.x = miu*g*h*u.x/norm(u);
      //}
      //if (u.y == 0.0){
      //  f.y = miu*g*h;
      //  if(a.y > 0.0){
      //    f.y = -f.y;
      //  }else if(a.y == 0.0){
      //    f.y = 0.0;
      //  }
      //}else if (u.y > 0.0){
      //  f.y = -miu*g*h*u.y/norm(u);
      //}else{
      //  f.y = miu*g*h*u.y/norm(u);
      //}
      //Scalar t_mid, fx_1, fx_2, fy_1, fy_2;
      //if (a.x != -f.x){
      //  t_mid = -hU.x/(a.x + f.x); 
      //  if (0 <= t_mid && t_mid <=dt){
      //    fx_1 = f.x*t_mid;
      //    if(fabs(f.x) > fabs(a.x)){
      //      fx_2 = -a.x*(dt-t_mid);
      //    }else{
      //      if(a.x >= 0.0){
      //        fx_2 = -f.x*(dt-t_mid);
      //      }else{
      //        fx_2 = f.x*(dt-t_mid);
      //      }
      //    }
      //    f.x = (fx_1 + fx_2)/dt;
      //  }
      //}
      //if (a.y != -f.y){
      //  t_mid = -hU.y/(a.y + f.y); 
      //  if (0 <= t_mid && t_mid <=dt){
      //    fy_1 = f.y*t_mid;
      //    if(fabs(f.y) > fabs(a.y)){
      //      fy_2 = -a.y*(dt-t_mid);
      //    }else{
      //      if(a.y >= 0.0){
      //        fy_2 = -f.y*(dt-t_mid);
      //      }else{
      //        fy_2 = f.y*(dt-t_mid);
      //      }
      //    }
      //    f.y = (fy_1 + fy_2)/dt;
      //  }
      //}
      //if (norm(u) == 0){
      //  f = 0.0;
      //}else{
      //  f = -miu*g*h*u/norm(u);
      //}
      //return f;
    }

  }//end of namespace fv

}//end of namespace GC
