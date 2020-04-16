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
  \flie cuda_riemann_solvers.cuh
  \brief Header file for cuda riemann solvers
  \version 0.1
  \author xilin xia
*/

#ifndef CUDA_RIEMANN_SOLVERS_CUH
#define CUDA_RIEMANN_SOLVERS_CUH

 __device__ RiemannFluxSWEs cuLaxFriedrichsRiemannSolverSWEs(const Scalar& g,const ScalarRiemannState& h,const VectorRiemannState& q){
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

      Scalar s_max = fmax(fabs(u_L) + sqrt(g*h_L) ,fabs(u_R) + sqrt(g*h_R));  

      Scalar h_flux = 0.5*(h_L*u_L + h_R*u_R) - 0.5*s_max*(h_R - h_L);
      Scalar qx_flux = 0.5*(h_L*u_L*u_L + 0.5*g*h_L*h_L + h_R*u_R*u_R + 0.5*g*h_R*h_R) - 0.5*s_max*(h_R*u_R - h_L*u_L);
      Scalar qy_flux = 0.5*(h_L*v_L*u_L + h_R*v_R*u_R) - 0.5*s_max*(h_R*v_R - h_L*v_L);

      return RiemannFluxSWEs(h_flux, Vector2(qx_flux, qy_flux));
    };


    __device__ RiemannFluxSWEs cuHLLCRiemannSolverSWEs(const Scalar& g,const ScalarRiemannState& h,const VectorRiemannState& q){
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

    __device__ RiemannFluxSWEs cuHLLCRiemannSolverMSWEs(const Scalar& g, const Scalar& z_f, const ScalarRiemannState& h, const ScalarRiemannState& eta, const VectorRiemannState& u){
      Scalar h_L = h.L;
      Scalar h_R = h.R;
      Scalar eta_L = eta.L;
      Scalar eta_R = eta.R;
      Scalar u_L = u.L.x;
      Scalar v_L = u.L.y;
      Scalar u_R = u.R.x;
      Scalar v_R = u.R.y;
      Scalar qx_L = h_L*u_L;
      Scalar qx_R = h_R*u_R;
      Scalar qy_L = h_L*v_L;
      Scalar qy_R = h_R*v_R;
      Scalar h_small = 1.0e-10;

      if (h_L <= h_small && h_R <= h_small){
        return RiemannFluxSWEs(0.0, Vector2(0.5*g*eta_L*eta_L - g*z_f*eta_L, 0.0));
      }

      Scalar a_L = sqrt(g*h_L);
      Scalar a_R = sqrt(g*h_R);

      Scalar h_star = pow((a_L + a_R) / 2.0 + (u_L - u_R) / 4.0, 2.0) / g;
      Scalar u_star = (u_L + u_R) / 2.0 + a_L - a_R;
      Scalar a_star = sqrt(g*h_star);

      Scalar s_L, s_R;
      if (h_L <= h_small){
        s_L = u_R - 2.0*a_R;
      }
      else{
        s_L = fmin(u_L - a_L, u_star - a_star);
      }

      if (h_R <= h_small){
        s_R = u_L + 2.0*a_L;
      }
      else{
        s_R = fmax(u_R + a_R, u_star + a_star);
      }

      Scalar s_M = (s_L*h_R*(u_R - s_R) - s_R*h_L*(u_L - s_L)) / (h_R*(u_R - s_R) - h_L*(u_L - s_L));

      Scalar h_flux_L = qx_L;
      Scalar qx_flux_L = u_L*qx_L + 0.5*g*eta_L*eta_L - g*z_f*eta_L;
      Scalar qy_flux_L = u_L*qy_L;

      Scalar h_flux_R = qx_R;
      Scalar qx_flux_R = u_R*qx_R + 0.5*g*eta_R*eta_R - g*z_f*eta_R;
      Scalar qy_flux_R = u_R*qy_R;

      Scalar h_flux_M = (s_R*h_flux_L - s_L*h_flux_R + s_L*s_R*(eta_R - eta_L)) / (s_R - s_L);
      Scalar qx_flux_M = (s_R*qx_flux_L - s_L*qx_flux_R + s_L*s_R*(qx_R - qx_L)) / (s_R - s_L);

      Scalar h_flux, qx_flux, qy_flux;
      if (0.0 <= s_L){
        h_flux = h_flux_L;
        qx_flux = qx_flux_L;
        qy_flux = qy_flux_L;
      }
      else if (s_L <= 0.0 && 0.0 <= s_M){
        h_flux = h_flux_M;
        qx_flux = qx_flux_M;
        qy_flux = h_flux_M*v_L;
      }
      else if (s_M <= 0.0 && 0.0 <= s_R){
        h_flux = h_flux_M;
        qx_flux = qx_flux_M;
        qy_flux = h_flux_M*v_R;
      }
      else{
        h_flux = h_flux_R;
        qx_flux = qx_flux_R;
        qy_flux = qy_flux_R;
      }
      return RiemannFluxSWEs(h_flux, Vector2(qx_flux, qy_flux));
    };

__device__ void GeorgeRiemannType(double hL, double hR, double uL, double uR, double drytol, double g, double* _hm, double* _s1m, double* _s2m){

  double hm,u1m,u2m,um,delu,s1m,s2m;
  double h_max,h_min,h0,F_max,F_min,dfdh,F0,slope,gL,gR;

  //Test for Riemann structure

  h_min = fmin(hR,hL);
  h_max = fmax(hR,hL);
  delu=uR-uL;


  if (h_min < drytol){
    hm=0.0;
    um=0.0;
    s1m=uR+uL-2.0*sqrt(g*hR)+2.0*sqrt(g*hL);
    s2m=uR+uL-2.0*sqrt(g*hR)+2.0*sqrt(g*hL);
  }else{
    F_min= delu+2.0*(sqrt(g*h_min)-sqrt(g*h_max));
    F_max= delu + (h_max-h_min)*(sqrt(.5*g*(h_max+h_min)/(h_max*h_min)));
    if (F_min > 0.0){ //2-rarefactions
      hm=(1.0/(16.0*g))*pow(fmax(0.0,-delu+2.0*(sqrt(g*hL)+sqrt(g*hR))),2.0);
      double sign = 0.0;
      if(hm > 0.0){
        sign = 1.0;
      }else if(hm < 0.0){
        sign = -1.0;
      }
      um=sign*(uL+2.0*(sqrt(g*hL)-sqrt(g*hm)));
      s1m=uL+2.0*sqrt(g*hL)-3.0*sqrt(g*hm);
      s2m=uR-2.0*sqrt(g*hR)+3.0*sqrt(g*hm);
    }else if (F_max < 0.0){ //2 shocks
      //root finding using a Newton iteration on sqrt(h)===
      h0=h_max;
      double h1 = h_max;
      int k = 0;
      while(true){
        k++;
        gL=sqrt(0.5*g*(1.0/h0 + 1.0/hL));
        gR=sqrt(0.5*g*(1.0/h0 + 1.0/hR));
        F0=delu+(h0-hL)*gL + (h0-hR)*gR;
        dfdh=gL-g*(h0-hL)/(4.0*(h0*h0)*gL)+gR-g*(h0-hR)/(4.0*(h0*h0)*gR);
        slope=2.0*sqrt(h0)*dfdh;
        h1=pow((sqrt(h0)-F0/slope),2.0);
        if(fabs(h0-h1) <= 0.0001*fabs(h0) || k > 10){
          break;
        }
        h0 = h1;
      }
      hm=h0;
      u1m=uL-(hm-hL)*sqrt((.50*g)*(1.0/hm + 1.0/hL));
      u2m=uR+(hm-hR)*sqrt((.50*g)*(1.0/hm + 1.0/hR));
      um=.5*(u1m+u2m);
      s1m=u1m-sqrt(g*hm);
      s2m=u2m+sqrt(g*hm);
    }else{ //one shock one rarefaction
      h0=h_min;
      double h1 = h_min;
      int k = 0;
      while(true){
        k++;
        F0=delu + 2.0*(sqrt(g*h0)-sqrt(g*h_max)) + (h0-h_min)*sqrt(.5*g*(1.0/h0+1.0/h_min));
        slope=(F_max-F0)/(h_max-h_min);
        h1=h0-F0/slope;
        if(fabs(h0-h1) <= 0.0001*fabs(h0) || k > 10){
          break;
        }
        h0 = h1;
      }
      hm=h0;
      if (hL > hR){
        um=uL+2.0*sqrt(g*hL)-2.0*sqrt(g*hm);
        s1m=uL+2.0*sqrt(g*hL)-3.0*sqrt(g*hm);
        s2m=uL+2.0*sqrt(g*hL)-sqrt(g*hm);
      }else{
        s2m=uR-2.0*sqrt(g*hR)+3.0*sqrt(g*hm);
        s1m=uR-2.0*sqrt(g*hR)+sqrt(g*hm);
        um=uR-2.0*sqrt(g*hR)+2.0*sqrt(g*hm);
      }
    }
  }

  *_hm = hm;
  *_s1m = s1m;
  *_s2m = s2m;
}


__device__ void GeorgeRiemannJCP(double hL,double hR, double huL,double huR,double bL,
    double bR,double uL,double uR,double vL,double vR,
    double phiL,double phiR,double sE1,double sE2,double drytol,double g,double sw[3],double fw[3][3]){

  //local
  double A[3][3];
  double r[3][3];
  double lambda[3];
  double del[3];
  double beta[3];

  double delh,delhu,delphi,delb,delnorm;
  double criticaltol;
  double s1s2bar,s1s2tilde,hbar,hLstar,hRstar;
  double huRstar,huLstar,uRstar,uLstar,hstarHLL;
  double deldelh,deldelphi;
  double s1m,s2m,hm;
  double det1,det2,det3,determinant;

  bool sonic;
  
  //determine del vectors
  delh = hR-hL;
  delhu = huR-huL;
  delphi = phiR-phiL;
  delb = bR-bL;
  delnorm = pow(delh,2) + pow(delphi,2);

  GeorgeRiemannType(hL, hR, uL, uR, drytol,g,&hm,&s1m,&s2m);

  lambda[0] = fmin(sE1,s2m); //Modified Einfeldt speed
  lambda[2] = fmax(sE2,s1m); //Modified Einfeldt speed
  sE1=lambda[0];
  sE2=lambda[2];
  lambda[1] = 0.0;  // ### Fix to avoid uninitialized value in loop on mw -- Correct?? ###

  hstarHLL = fmax((huL-huR+sE2*hR-sE1*hL)/(sE2-sE1),0.0); // middle state in an HLL solve

  for(int mw = 0; mw < 3; mw++){
    r[0][mw]=1.0;
    r[1][mw]=lambda[mw];
    r[2][mw]=pow((lambda[mw]),2);
  }

  lambda[1] = 0.5*(lambda[0]+lambda[2]);
  r[0][1]=0.0;
  r[1][1]=0.0;
  r[2][1]=1.0;

  //determine the steady state wave -------------------
  criticaltol = 1.0e-6;
  deldelh = -delb;
  deldelphi = -g*0.5*(hR+hL)*delb;

  //determine a few quanitites needed for steady state wave if iterated
  hLstar=hL;
  hRstar=hR;
  uLstar=uL;
  uRstar=uR;
  huLstar=uLstar*hLstar;
  huRstar=uRstar*hRstar;


  hbar =  fmax(0.5*(hLstar+hRstar),0.0);
  s1s2bar = 0.25*pow((uLstar+uRstar),2) - g*hbar;
  s1s2tilde= fmax(0.0,uLstar*uRstar) - g*hbar;

  
  //find if sonic problem
  sonic = false;
  if (fabs(s1s2bar) < criticaltol) sonic= true;
  if (s1s2bar*s1s2tilde < criticaltol) sonic=true;
  if (s1s2bar*sE1*sE2 < criticaltol) sonic = true;
  if (fmin(fabs(sE1),fabs(sE2)) < criticaltol) sonic= true;
  if (sE1 < 0.0 && s1m > 0.0) sonic = true;
  if (sE2 > 0.0 && s2m < 0.0) sonic = true;
  if ((uL+sqrt(g*hL))*(uR+sqrt(g*hR)) < 0.0) sonic= true;
  if ((uL-sqrt(g*hL))*(uR-sqrt(g*hR)) < 0.0) sonic= true;

  //find jump in h, deldelh
  if (sonic){
    deldelh =  -delb;
  }else{
    deldelh = delb*g*hbar/s1s2bar;
  }
  //find bounds in case of critical state resonance, or negative states
  if (sE1 < -1.0*criticaltol && sE2 > criticaltol){
    deldelh = fmin(deldelh,hstarHLL*(sE2-sE1)/sE2);
    deldelh = fmax(deldelh,hstarHLL*(sE2-sE1)/sE1);
  }else if (sE1 >= criticaltol){
    deldelh = fmin(deldelh,hstarHLL*(sE2-sE1)/sE1);
    deldelh = fmax(deldelh,-hL);
  }else if (sE2 <= -1.0*criticaltol){
    deldelh = fmin(deldelh,hR);
    deldelh = fmax(deldelh,hstarHLL*(sE2-sE1)/sE2);
  }

  //find jump in phi, deldelphi
  if (sonic){
    deldelphi = -g*hbar*delb;
  }else{
    deldelphi = -delb*g*hbar*s1s2tilde/s1s2bar;
  }

  //find bounds in case of critical state resonance, or negative states
  deldelphi=fmin(deldelphi,g*fmax(-hLstar*delb,-hRstar*delb));
  deldelphi=fmax(deldelphi,g*fmin(-hLstar*delb,-hRstar*delb));

  del[0]=delh-deldelh;
  del[1]=delhu;
  del[2]=delphi-deldelphi;

  //Determine determinant of eigenvector matrix========
  det1=r[0][0]*(r[1][1]*r[2][2]-r[1][2]*r[2][1]);
  det2=r[0][1]*(r[1][0]*r[2][2]-r[1][2]*r[2][0]);
  det3=r[0][2]*(r[1][0]*r[2][1]-r[1][1]*r[2][0]);
  determinant=det1-det2+det3;

  //solve for beta(k) using Cramers Rule=================
  for(int k=0; k < 3; k++){
    for(int mw=0; mw < 3; mw++){
      A[0][mw]=r[0][mw];
      A[1][mw]=r[1][mw];
      A[2][mw]=r[2][mw];
    }
    A[0][k]=del[0];
    A[1][k]=del[1];
    A[2][k]=del[2];
    det1=A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1]);
    det2=A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0]);
    det3=A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);
    beta[k]=(det1-det2+det3)/determinant;
  }

  for(int mw = 0; mw < 3; mw++){
    sw[mw]=lambda[mw];
    fw[0][mw]=beta[mw]*r[1][mw]; //mass flux
    fw[1][mw]=beta[mw]*r[2][mw]; //normal momentum flux
    fw[2][mw]=beta[mw]*r[1][mw]; //shear momentum flux
  }
  //find transverse components (ie huv jumps).
  fw[2][0]=fw[2][0]*vL;
  fw[2][2]=fw[2][2]*vR;
  fw[2][1]= hR*uR*vR - hL*uL*vL - fw[2][0]- fw[2][2];

}




__device__ RiemannFluxSWEs cuGeorgeRiemannSolverSWEs(const Scalar g, const ScalarRiemannState& h,const VectorRiemannState& q, const ScalarRiemannState& b){

  double drytol = 1e-10;

  double hR,hL,huR,huL,uR,uL,hvR,hvL,vR,vL,phiR,phiL;
  double bR,bL,sL,sR,sRoe1,sRoe2,sE1,sE2,uhat,chat;
  double s1m,s2m;
  double hstar,hstartest;

  //Riemann problem variables
  hL = h.L;
  hR = h.R; 
  huL = q.L.x;
  huR = q.R.x;
  bL = b.L;
  bR = b.R;
  hvL= q.L.y; 
  hvR= q.R.y;

  //check for wet/dry boundary
  if(hR > drytol){
    uR=huR/hR;
    vR=hvR/hR;
    phiR = 0.5*g*hR*hR + huR*huR/hR;
  }else{
    hR = 0.0;
    huR = 0.0;
    hvR = 0.0;
    uR = 0.0;
    vR = 0.0;
    phiR = 0.0;
  }

  if(hL > drytol){
    uL=huL/hL;
    vL=hvL/hL;
    phiL = 0.5*g*hL*hL + huL*huL/hL;
  }else{
    hL=0.0;
    huL=0.0;
    hvL=0.0;
    uL=0.0;
    vL=0.0;
    phiL = 0.0;
  }

  double wall[3];

  wall[0] = 1.0;
  wall[1] = 1.0;
  wall[2] = 1.0;

  if (hR < drytol){
    GeorgeRiemannType(hL, hL, uL, -1.0*uL, drytol,g,&hstar,&s1m,&s2m);
    hstartest=fmax(hL,hstar);
    if (hstartest+bL < bR){ //right state should become ghost values that mirror left for wall problem
      bR=hstartest+bL;
      wall[1]=0.0;
      wall[2]=0.0;
      hR=hL;
      huR=-huL;
      bR=bL;
      phiR=phiL;
      uR=-uL;
      vR=vL;
    }else if(hL+bL < bR){
      bR=hL+bL;
    }
  }else if (hL< drytol){ //right surface is lower than left topo
    GeorgeRiemannType(hR, hR, -1.0*uR, uR, drytol,g,&hstar,&s1m,&s2m);
    hstartest=fmax(hR,hstar);
    if (hstartest+bR < bL){  //left state should become ghost values that mirror right
      bL=hstartest+bR;
      wall[0] = 0.0;
      wall[1] = 0.0;
      hL=hR;
      huL=-huR;
      bL=bR;
      phiL=phiR;
      uL=-uR;
      vL=vR;
    }else if (hR+bR < bL){
      bL=hR+bR;
    }
  }


  //determine wave speeds
  sL=uL-sqrt(g*hL); // 1 wave speed of left state
  sR=uR+sqrt(g*hR); // 2 wave speed of right state

  uhat=(sqrt(g*hL)*uL + sqrt(g*hR)*uR)/(sqrt(g*hR)+sqrt(g*hL)); // Roe average
  chat=sqrt(g*0.5*(hR+hL)); // Roe average
  sRoe1=uhat-chat; // Roe wave speed 1 wave
  sRoe2=uhat+chat; // Roe wave speed 2 wave

  sE1 = fmin(sL,sRoe1); // Eindfeldt speed 1 wave
  sE2 = fmax(sR,sRoe2); // Eindfeldt speed 2 wave

  //--------------------end initializing...finally----------

  double fw[3][3];
  double sw[3];  

  GeorgeRiemannJCP(hL,hR, huL,huR,bL,bR,uL,uR,vL,vR,phiL,phiR,sE1,sE2,drytol,g,sw,fw);
  
  for(int mw = 0; mw < 3; mw++){
    sw[mw]=sw[mw]*wall[mw];
    fw[0][mw]=fw[0][mw]*wall[mw]; 
    fw[1][mw]=fw[1][mw]*wall[mw];
    fw[2][mw]=fw[2][mw]*wall[mw];
  }

  double amdq[3] = {0.0};
  double apdq[3] = {0.0};

  for(int mw=0; mw < 3; mw++){
    if (sw[mw] < 0.0){
      amdq[0] = amdq[0] + fw[0][mw];
      amdq[1] = amdq[1] + fw[1][mw];
      amdq[2] = amdq[2] + fw[2][mw];
    }else if (sw[mw] > 0.0){
      apdq[0] = apdq[0] + fw[0][mw];
      apdq[1] = apdq[1] + fw[1][mw];
      apdq[2] = apdq[2] + fw[2][mw];
    }else{
      apdq[0] = apdq[0] + 0.5*fw[0][mw];
      apdq[1] = apdq[1] + 0.5*fw[1][mw];
      apdq[2] = apdq[2] + 0.5*fw[2][mw];
      amdq[0] = amdq[0] + 0.5*fw[0][mw];
      amdq[1] = amdq[1] + 0.5*fw[1][mw];
      amdq[2] = amdq[2] + 0.5*fw[2][mw];
    }
  }

  //double h_flux = (amdq[0] - apdq[0] + huL + huR)*0.5;

  return RiemannFluxSWEs(amdq[0], Vector2(amdq[1], amdq[2]));

}

#endif
