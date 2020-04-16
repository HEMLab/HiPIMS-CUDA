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
\file cuda_limiter.cu
\brief Source file for gradient operator by cuda

\version 0.1
\author xilin xia

*/


#include "cuda_limiter.h"
#include "cuda_kernel_launch_parameters.h"
//#include "cuda_boundary.h"


namespace GC{

  namespace fv{

    __device__ void NonSimplex(unsigned int index, double standard_matrix[][3], unsigned int M, double &L_x, double &L_y){
      unsigned int h = M - 2;
      unsigned int k = M - 1; //indices for the initial working set, e.g. L_x >= 0, L_y >= 0
      unsigned int to_remove = 0;
      unsigned int to_add = 0;
      L_x = 0.0;
      L_y = 0.0;
      Scalar p_1 = 0.0;
      Scalar p_2 = 0.0;
      //avoiding very small numbers
      for(int i = 0; i < 13; i++){
        for(int j = 0; j < 3; j++){
          if(fabs(standard_matrix[i][j]) <= 1e-11){
            standard_matrix[i][j] = 0.0;
          }
        }
      }

      if (fabs(standard_matrix[M][0]) < 1e-11 && fabs(standard_matrix[M][1]) < 1e-11){
        return;
      }
      while (true){
        Scalar det_A = standard_matrix[h][0] * standard_matrix[k][1] - standard_matrix[h][1] * standard_matrix[k][0];
        if (fabs(det_A) < 1e-22){
          p_1 = standard_matrix[k][1] / (standard_matrix[k][0] * standard_matrix[M][1] - standard_matrix[k][1] * standard_matrix[M][0]);
          p_2 = -standard_matrix[k][0] / (standard_matrix[k][0] * standard_matrix[M][1] - standard_matrix[k][1] * standard_matrix[M][0]);
        }
        else{
          //calculate lambda1 and lambda2
          Scalar lambda_1 = (standard_matrix[M][0] * standard_matrix[k][1] - standard_matrix[M][1] * standard_matrix[k][0]) / det_A;
          Scalar lambda_2 = (standard_matrix[M][1] * standard_matrix[h][0] - standard_matrix[M][0] * standard_matrix[h][1]) / det_A;
          Scalar e_1 = 0.0;
          Scalar e_2 = 0.0;
          Scalar min_lambda = lambda_1;
          if (min_lambda > lambda_2){
            min_lambda = lambda_2;
            e_2 = 1.0;
            to_remove = 1;
          }
          else{
            e_1 = 1.0;
            to_remove = 0;
          }
          if (min_lambda >= 0.0){
            //reaches optimal
            return;
          }
          else{
            //calculate p_1 and p_2
            p_1 = (e_1 * standard_matrix[k][1] - e_2 * standard_matrix[h][1]) / det_A;
            p_2 = (e_2 * standard_matrix[h][0] - e_1 * standard_matrix[k][0]) / det_A;
          }
        }
        Scalar min_sigma = 1e15;
        for (unsigned int j = 0; j < M; j++){
          if (j != h && j != k){
            Scalar step_length = standard_matrix[j][0] * p_1 + standard_matrix[j][1] * p_2;
            Scalar new_sigma = 1e15;
            if (step_length < 0.0){
              new_sigma = -(standard_matrix[j][0] * L_x + standard_matrix[j][1] * L_y - standard_matrix[j][2]) / step_length;
            }
            if (min_sigma > new_sigma){
              min_sigma = new_sigma;
              to_add = j;
            }
          }
        }
        //update the variables
        L_x = L_x + min_sigma*p_1;
        L_y = L_y + min_sigma*p_2;
        //modify the working set
        if (to_remove == 0){
          h = to_add;
        }
        else{
          k = to_add;
        }
      }
    }

    //This kernel only supports cell with maximumly 4 faces and 2D problem
    __global__ void cuGradientLimiterKernelNonSimplex(Scalar* phi_data, Scalar* _phi_bound, unsigned int phi_size, Flag* cell_neigbours_dimensions, ShortDualHandle* cell_neigbours, unsigned int cell_neighbours_length, Vector2* cell_centres, Flag* cell_halffacets, unsigned int cell_halffacets_length, Vector2* face_centres, Vector2* face_normals, Vector2* phi_gradient, unsigned int phi_gradient_size){
      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      double standard_matrix[13][3];
      while (index < phi_size){
        Flag cell_neigbours_number = cell_neigbours_dimensions[index];
        unsigned int M = cell_neigbours_number * 2 + 4; //four constraints plus 4 extra ones
        standard_matrix[M - 1][0] = 1.0; //L_x >= 0.0
        standard_matrix[M - 1][1] = 0.0;
        standard_matrix[M - 1][2] = 0.0;
        standard_matrix[M - 2][0] = 0.0; //L_y >= 0.0
        standard_matrix[M - 2][1] = 1.0;
        standard_matrix[M - 2][2] = 0.0;
        standard_matrix[M - 3][0] = -1.0; //-L_x >= -1.0
        standard_matrix[M - 3][1] = 0.0;
        standard_matrix[M - 3][2] = -1.0;
        standard_matrix[M - 4][0] = 0.0; //-L_y <= -1.0
        standard_matrix[M - 4][1] = -1.0;
        standard_matrix[M - 4][2] = -1.0;
        Vector2 _gradient = phi_gradient[index];
        Vector2 centre_this = cell_centres[index];
        Vector2 direction;
        Scalar phi_this = phi_data[index];
        Scalar c_x = 0.0;  //the coefficients for objective function
        Scalar c_y = 0.0;
        for (Flag i = 0; i < cell_neigbours_number; ++i){
          ShortDualHandle neib = cell_neigbours[i*cell_neighbours_length + index];
          Scalar d_phi;
          if (!neib.is_boundary()){
            Flag id_neib = neib.get_global_id();
            Scalar phi_neib = phi_data[id_neib];
            d_phi = phi_neib - phi_this;
            Vector2 centre_neib = cell_centres[id_neib];
            direction = centre_neib - centre_this;
            //Flag id_face = cell_halffacets[i*cell_halffacets_length + index];
            //Vector2 centre_face = face_centres[id_face];
            //direction = centre_face - centre_this;
          }
          else{
            Flag id_face = cell_halffacets[i*cell_halffacets_length + index];
            Vector2 centre_face = face_centres[id_face];
            direction = 2.0*(centre_face - centre_this);
            Vector2 normal = uni(face_normals[id_face]);
            Flag id_boundary = neib.get_global_id();
//            ShortTripleFlag phi_boundary_type = phi_bound[id_boundary];
//            Scalar phi_bound = cuGetBoundary(phi_boundary_type, phi_this, normal);
            Scalar phi_bound = _phi_bound[id_boundary];
            d_phi = phi_bound - phi_this;
          }
          Scalar sign = 0.0;
          if (d_phi > 0){
            sign = -1.0;
          }
          else{
            sign = 1.0;
          }
          Scalar phi_x = _gradient.x*direction.x;
          Scalar phi_y = _gradient.y*direction.y;
          standard_matrix[2 * i][0] = phi_x*sign;
          standard_matrix[2 * i][1] = phi_y*sign;
          standard_matrix[2 * i][2] = d_phi*sign;
          standard_matrix[2 * i + 1][0] = -phi_x*sign;
          standard_matrix[2 * i + 1][1] = -phi_y*sign;
          standard_matrix[2 * i + 1][2] = 0.0;
          c_x -= fabs(phi_x);
          c_y -= fabs(phi_y);
        }
        standard_matrix[M][0] = c_x; //Obj = sum|(D_x*u_xi)|*L_x + sum|(D_y*u_yi)|*L_x
        standard_matrix[M][1] = c_y;
        standard_matrix[M][2] = 0.0;
        Scalar L_x = 0.0;
        Scalar L_y = 0.0;
        __syncthreads();
        NonSimplex(index, standard_matrix, M, L_x, L_y);
        _gradient.x *= L_x;
        _gradient.y *= L_y;
        phi_gradient[index] = _gradient;
        index += blockDim.x * gridDim.x;
      }
    }

    //This kernel only supports cell with maximumly 4 faces and 2D problem
    __global__ void cuGradientLimiterKernelNonSimplex(Vector2* phi_data, Vector2* _phi_bound, unsigned int phi_size, Flag* cell_neigbours_dimensions, ShortDualHandle* cell_neigbours, unsigned int cell_neighbours_length, Vector2* cell_centres, Flag* cell_halffacets, unsigned int cell_halffacets_length, Vector2* face_centres, Vector2* face_normals, Tensor2* phi_gradient, unsigned int phi_gradient_size){
      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      double standard_matrix[13][3];
      while (index < phi_size){
        Flag cell_neigbours_number = cell_neigbours_dimensions[index];
        unsigned int M = cell_neigbours_number * 2 + 4; //four constraints plus 2 extra ones
        Tensor2 _gradient = phi_gradient[index];

        //-----------------------X axis------------------------------
        standard_matrix[M - 1][0] = 1.0; //L_x >= 0.0
        standard_matrix[M - 1][1] = 0.0;
        standard_matrix[M - 1][2] = 0.0;
        standard_matrix[M - 2][0] = 0.0; //L_y >= 0.0
        standard_matrix[M - 2][1] = 1.0;
        standard_matrix[M - 2][2] = 0.0;
        standard_matrix[M - 3][0] = -1.0; //-L_x >= -1.0
        standard_matrix[M - 3][1] = 0.0;
        standard_matrix[M - 3][2] = -1.0;
        standard_matrix[M - 4][0] = 0.0; //-L_y <= -1.0
        standard_matrix[M - 4][1] = -1.0;
        standard_matrix[M - 4][2] = -1.0;
        Vector2 centre_this = cell_centres[index];
        Vector2 direction;
        Vector2 phi_this = phi_data[index];
        Scalar phi_this_x = phi_this.x;
        Scalar c_x = 0.0;
        Scalar c_y = 0.0;
        for (Flag i = 0; i < cell_neigbours_number; ++i){
          ShortDualHandle neib = cell_neigbours[i*cell_neighbours_length + index];
          Scalar d_phi;
          if (!neib.is_boundary()){
            Flag id_neib = neib.get_global_id();
            Scalar phi_neib_x = phi_data[id_neib].x;
            d_phi = phi_neib_x - phi_this_x;
            Vector2 centre_neib = cell_centres[id_neib];
            direction = centre_neib - centre_this;
            //Flag id_face = cell_halffacets[i*cell_halffacets_length + index];
            //Vector2 centre_face = face_centres[id_face];
            //direction = centre_face - centre_this;
          }
          else{
            Flag id_face = cell_halffacets[i*cell_halffacets_length + index];
            Vector2 centre_face = face_centres[id_face];
            direction = 2.0*(centre_face - centre_this);
            Vector2 normal = uni(face_normals[id_face]);
            Flag id_boundary = neib.get_global_id();
//            ShortTripleFlag phi_boundary_type = phi_bound[id_boundary];
//            Vector2 phi_bound = cuGetBoundary(phi_boundary_type, phi_this, normal);
            Vector2 phi_bound = _phi_bound[id_boundary];
            d_phi = phi_bound.x - phi_this.x;
          }
          Scalar sign = 0.0;
          if (d_phi > 0){
            sign = -1.0;
          }
          else{
            sign = 1.0;
          }
          Scalar phi_x = _gradient.xx*direction.x;
          Scalar phi_y = _gradient.xy*direction.y;
          standard_matrix[2 * i][0] = phi_x*sign;
          standard_matrix[2 * i][1] = phi_y*sign;
          standard_matrix[2 * i][2] = d_phi*sign;
          standard_matrix[2 * i + 1][0] = -phi_x*sign;
          standard_matrix[2 * i + 1][1] = -phi_y*sign;
          standard_matrix[2 * i + 1][2] = 0.0;
          c_x -= fabs(phi_x);
          c_y -= fabs(phi_y);
        }
        standard_matrix[M][0] = c_x; //Obj = sum|(D_x*u_xi)|*L_x + sum|(D_y*u_yi)|*L_x
        standard_matrix[M][1] = c_y;
        standard_matrix[M][2] = 0.0;
        Scalar L_x = 0.0;
        Scalar L_y = 0.0;
        __syncthreads();
        NonSimplex(index, standard_matrix, M, L_x, L_y);
        _gradient.xx *= L_x;
        _gradient.xy *= L_y;

        //-----------------------Y axis------------------------------
        standard_matrix[M - 1][0] = 1.0; //L_x >= 0.0
        standard_matrix[M - 1][1] = 0.0;
        standard_matrix[M - 1][2] = 0.0;
        standard_matrix[M - 2][0] = 0.0; //L_y >= 0.0
        standard_matrix[M - 2][1] = 1.0;
        standard_matrix[M - 2][2] = 0.0;
        standard_matrix[M - 3][0] = -1.0; //-L_x >= -1.0
        standard_matrix[M - 3][1] = 0.0;
        standard_matrix[M - 3][2] = -1.0;
        standard_matrix[M - 4][0] = 0.0; //-L_y <= -1.0
        standard_matrix[M - 4][1] = -1.0;
        standard_matrix[M - 4][2] = -1.0;
        Scalar phi_this_y = phi_this.y;
        c_x = 0.0;
        c_y = 0.0;
        for (Flag i = 0; i < cell_neigbours_number; ++i){
          ShortDualHandle neib = cell_neigbours[i*cell_neighbours_length + index];
          Scalar d_phi;
          if (!neib.is_boundary()){
            Flag id_neib = neib.get_global_id();
            Scalar phi_neib_y = phi_data[id_neib].y;
            d_phi = phi_neib_y - phi_this_y;
            Vector2 centre_neib = cell_centres[id_neib];
            direction = centre_neib - centre_this;
            //Flag id_face = cell_halffacets[i*cell_halffacets_length + index];
            //Vector2 centre_face = face_centres[id_face];
            //direction = centre_face - centre_this;
          }
          else{
            Flag id_face = cell_halffacets[i*cell_halffacets_length + index];
            Vector2 centre_face = face_centres[id_face];
            direction = 2.0*(centre_face - centre_this);
            Vector2 normal = uni(face_normals[id_face]);
            Flag id_boundary = neib.get_global_id();
//            ShortTripleFlag phi_boundary_type = phi_bound[id_boundary];
//            Vector2 phi_bound = cuGetBoundary(phi_boundary_type, phi_this, normal);
            Vector2 phi_bound = _phi_bound[id_boundary];
            d_phi = phi_bound.y - phi_this.y;
          }
          Scalar sign = 0.0;
          if (d_phi > 0){
            sign = -1.0;
          }
          else{
            sign = 1.0;
          }
          Scalar phi_x = _gradient.yx*direction.x;
          Scalar phi_y = _gradient.yy*direction.y;
          standard_matrix[2 * i][0] = phi_x*sign;
          standard_matrix[2 * i][1] = phi_y*sign;
          standard_matrix[2 * i][2] = d_phi*sign;
          standard_matrix[2 * i + 1][0] = -phi_x*sign;
          standard_matrix[2 * i + 1][1] = -phi_y*sign;
          standard_matrix[2 * i + 1][2] = 0.0;
          c_x -= fabs(phi_x);
          c_y -= fabs(phi_y);
        }
        standard_matrix[M][0] = c_x; //Obj = sum|(D_x*u_xi)|*L_x + sum|(D_y*u_yi)|*L_x
        standard_matrix[M][1] = c_y;
        standard_matrix[M][2] = 0.0;
        L_x = 0.0;
        L_y = 0.0;
        __syncthreads();
        NonSimplex(index, standard_matrix, M, L_x, L_y);
        _gradient.yx *= L_x;
        _gradient.yy *= L_y;
        phi_gradient[index] = _gradient;
        index += blockDim.x * gridDim.x;
      }
    }


    __device__ void pivot(double standard_matrix[][3], unsigned int M, unsigned int i, unsigned int j){
      
      double a_ij = standard_matrix[i][j];
      for(unsigned int h = 0; h < M; ++h){
        for(unsigned int k = 0; k < 3; ++k){
          if(h != i && k != j){
            standard_matrix[h][k] -= standard_matrix[h][j]*standard_matrix[i][k]/a_ij; 
          }
        }
      }
      for(unsigned int h = 0; h < M; ++h){
        if(h != i){
          standard_matrix[h][j] = -standard_matrix[h][j]/a_ij;
        }
      }
      for(unsigned int k = 0; k < 3; ++k){
        if(j != k){
          standard_matrix[i][k] = standard_matrix[i][k]/a_ij;
        }
      }
      standard_matrix[i][j] = 1/a_ij;

    }

    
    __device__ void simplex(double standard_matrix[][3], unsigned int M, ShortDualFlag x_vector[], ShortDualFlag y_vector[], ShortDualFlag x_pos[], double &L_x, double &L_y){
      
      //avoiding very small numbers
      for(int i = 0; i < 11; i++){
        for(int j = 0; j < 3; j++){
          if(fabs(standard_matrix[i][j]) <= 1e-14){
            standard_matrix[i][j] = 0.0;
          }
        }
      }

      while(true){
        double b_min = standard_matrix[0][2];
        for(unsigned int i = 0; i < M; ++i){
          if(standard_matrix[i][2] < b_min){
            b_min = standard_matrix[i][2];
          }
        }
        if(b_min >= 0){  // b >= 0
          unsigned int k = 0;
          double c_min, c_1, c_2;
          c_1 = standard_matrix[M][0];
          c_2 = standard_matrix[M][1];
          if(c_1 <= c_2){
            c_min = c_1;
            k = 0;
          }else{
            c_min = c_2;
            k = 1;
          }
          if(c_min >= 0){//Termination criteria
            unsigned int id = x_pos[0].getx();
            if(id == 0){
              L_x = 0.0;
            }else{
              L_x = standard_matrix[x_pos[0].gety()][2];
            }
            id = x_pos[1].getx();
            if(id == 0){
              L_y = 0.0;
            }else{
              L_y = standard_matrix[x_pos[1].gety()][2];
            }
            break;
          }else{// try to pivot
            //find smallest positive b_i/a_ij
            unsigned int h = 0;
            double ba_min = 1e15;
            for(unsigned int i = 0; i < M; ++i){
              double a_ik = standard_matrix[i][k];
              double b_i = standard_matrix[i][2];
              if(a_ik > 0){
                if(b_i/a_ik < ba_min){
                  h = i;
                  ba_min = b_i/a_ik;
                }
              } 
            }
            //pivot
            pivot(standard_matrix, M+1, h, k);
            //update vector x y and positions of variable x
            ShortDualFlag var_in_x = x_vector[k];
            ShortDualFlag var_in_y = y_vector[h];
            x_vector[k] = var_in_y;
            y_vector[h] = var_in_x;
            if(var_in_x.getx() == 0){//is a x-variable, update its position
              unsigned int id = var_in_x.gety();
              x_pos[id].setx(1);
              x_pos[id].sety(h);
            }
            if(var_in_y.getx() == 0){//is a x-variable, update its position
              unsigned int id = var_in_y.gety();
              x_pos[id].setx(0);
              x_pos[id].sety(k);
            }
          }
        }else{ //there exists b < 0
          unsigned int h = 0;
          for(unsigned int i = 0; i < M; i++){
            if(standard_matrix[i][2] < 0){
              h = i;
              break;
            }
          }
          unsigned int k = 0;
          for(unsigned int i = 0; i < 2; i++){
            if(standard_matrix[h][i] < 0){
              k = i;
              break;
            }
          }
          double ba_min = standard_matrix[h][2]/standard_matrix[h][k];
          for(unsigned int i = 0; i < M; i++){
            double a_ik = standard_matrix[i][k];
            double b_i = standard_matrix[i][2];
            if(a_ik > 0 && b_i >= 0){
              if(b_i/a_ik < ba_min){
                h = i;
                ba_min = b_i/a_ik;
              }
            }
          }
          //pivot
          pivot(standard_matrix, M+1, h, k);
          //update vector x y and positions of variable x
          ShortDualFlag var_in_x = x_vector[k];
          ShortDualFlag var_in_y = y_vector[h];
          x_vector[k] = var_in_y;
          y_vector[h] = var_in_x;
          if(var_in_x.getx() == 0){//is a x-variable, update its position
            unsigned int id = var_in_x.gety();
            x_pos[id].setx(1);
            x_pos[id].sety(h);
          }
          if(var_in_y.getx() == 0){//is a x-variable, update its position
            unsigned int id = var_in_y.gety();
            x_pos[id].setx(0);
            x_pos[id].sety(k);
          }
        } 
      }
    }

    //This kernel only supports cell with maximumly 4 faces and 2D problem
    __global__ void cuGradientLimiterKernelSimplex(Scalar* phi_data, Scalar* _phi_bound, unsigned int phi_size, Flag* cell_neigbours_dimensions, ShortDualHandle* cell_neigbours, unsigned int cell_neighbours_length, Vector2* cell_centres, Flag* cell_halffacets, unsigned int cell_halffacets_length, Vector2* face_centres, Vector2* face_normals, Vector2* phi_gradient, unsigned int phi_gradient_size){
      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      double standard_matrix[11][3];
      while(index < phi_size){
        Flag cell_neigbours_number = cell_neigbours_dimensions[index];
        unsigned int M = cell_neigbours_number*2 + 2; //four constraints plus 2 extra ones
      //  unsigned int N = 2; //two variables L_x, L_y
        //Short index: 0--x-variable 1--y-variable, Long index: index
        ShortDualFlag x_vector[2] = {ShortDualFlag(0,0), ShortDualFlag(0,1)};
        ShortDualFlag y_vector[10] = {ShortDualFlag(1,0), ShortDualFlag(1,1), ShortDualFlag(1,2), ShortDualFlag(1,3), ShortDualFlag(1,4), 
                                      ShortDualFlag(1,5), ShortDualFlag(1,6), ShortDualFlag(1,7), ShortDualFlag(1,8), ShortDualFlag(1,9)};
        //Short index: 0--x-vector 1--y-vector, Long index: index
        ShortDualFlag x_pos[2] = {ShortDualFlag(0,0), ShortDualFlag(0,1)};
        standard_matrix[M - 2][0] = 1.0; //L_x <= 1.0
        standard_matrix[M - 2][1] = 0.0; 
        standard_matrix[M - 2][2] = 1.0;
        standard_matrix[M - 1][0] = 0.0; //L_y <= 1.0
        standard_matrix[M - 1][1] = 1.0; 
        standard_matrix[M - 1][2] = 1.0;
        Vector2 _gradient = phi_gradient[index];
 //       standard_matrix[M][0] = -fabs(_gradient.x); //Obj = |D_x| * L_x + |D_y| * L_y 
 //       standard_matrix[M][1] = -fabs(_gradient.y); 
 //       standard_matrix[M][2] = 0.0;
        Vector2 centre_this = cell_centres[index];
        Vector2 direction;
        Scalar phi_this = phi_data[index];
        Scalar c_x = 0.0;  //the coefficients for objective function
        Scalar c_y = 0.0;
        for(Flag i = 0; i < cell_neigbours_number; ++i){
          ShortDualHandle neib = cell_neigbours[i*cell_neighbours_length+index];
          Scalar d_phi;
          if(!neib.is_boundary()){
            Flag id_neib = neib.get_global_id();
            Scalar phi_neib = phi_data[id_neib];
            d_phi = phi_neib - phi_this;
            Vector2 centre_neib = cell_centres[id_neib];
            direction = centre_neib - centre_this;
//            Flag id_face = cell_halffacets[i*cell_halffacets_length + index];
//            Vector2 centre_face = face_centres[id_face];
//            direction = centre_face - centre_this;
          }else{
            Flag id_face = cell_halffacets[i*cell_halffacets_length+index];
            Vector2 centre_face = face_centres[id_face];
            direction = 2.0*(centre_face - centre_this);
            Vector2 normal = uni(face_normals[id_face]);
            Flag id_boundary = neib.get_global_id();
//            ShortTripleFlag phi_boundary_type = phi_bound[id_boundary];
//            Scalar phi_bound = cuGetBoundary(phi_boundary_type, phi_this, normal);
//            Scalar phi_bound = _phi_bound[id_boundary];
            d_phi = 0.0; //phi_bound - phi_this;
          }
          Scalar sign = 0.0;
          if(d_phi > 0){
            sign = 1.0;          
          }else{
            sign = -1.0;
          }
          standard_matrix[2 * i][0] = _gradient.x*direction.x*sign;
          standard_matrix[2 * i][1] = _gradient.y*direction.y*sign;
          standard_matrix[2 * i][2] = d_phi*sign;
          standard_matrix[2 * i + 1][0] = -_gradient.x*direction.x*sign;
          standard_matrix[2 * i + 1][1] = -_gradient.y*direction.y*sign;
          standard_matrix[2 * i + 1][2] = 0.0;
          c_x -= fabs(_gradient.x*direction.x);
          c_y -= fabs(_gradient.y*direction.y);
        }
        standard_matrix[M][0] = c_x; //Obj = sum|(D_x*u_xi)|*L_x + sum|(D_y*u_yi)|*L_x
        standard_matrix[M][1] = c_y; 
        standard_matrix[M][2] = 0.0;
        double L_x = 1.0;
        double L_y = 1.0;
        __syncthreads();
        simplex(standard_matrix, M, x_vector, y_vector, x_pos, L_x, L_y);
        _gradient.x *= L_x;
        _gradient.y *= L_y;
        phi_gradient[index] = _gradient;
        index += blockDim.x * gridDim.x;
      }
    }


    //This kernel only supports cell with maximumly 4 faces and 2D problem
    __global__ void cuGradientLimiterKernelSimplex(Vector2* phi_data, Vector2* _phi_bound, unsigned int phi_size, Flag* cell_neigbours_dimensions, ShortDualHandle* cell_neigbours, unsigned int cell_neighbours_length, Vector2* cell_centres, Flag* cell_halffacets, unsigned int cell_halffacets_length, Vector2* face_centres, Vector2* face_normals, Tensor2* phi_gradient, unsigned int phi_gradient_size){
      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      double standard_matrix[11][3];
      while(index < phi_size){
        Flag cell_neigbours_number = cell_neigbours_dimensions[index];
        unsigned int M = cell_neigbours_number*2 + 2; //four constraints plus 2 extra ones
      //  unsigned int N = 2; //two variables L_x, L_y
        Tensor2 _gradient = phi_gradient[index];
        
      //-----------------------X axis------------------------------

        //Short index: 0--x-variable 1--y-variable, Long index: index
        ShortDualFlag x_vector[2] = {ShortDualFlag(0,0), ShortDualFlag(0,1)};
        ShortDualFlag y_vector[10] = {ShortDualFlag(1,0), ShortDualFlag(1,1), ShortDualFlag(1,2), ShortDualFlag(1,3), ShortDualFlag(1,4), 
                                      ShortDualFlag(1,5), ShortDualFlag(1,6), ShortDualFlag(1,7), ShortDualFlag(1,8), ShortDualFlag(1,9)};
        //Short index: 0--x-vector 1--y-vector, Long index: index
        ShortDualFlag x_pos[2] = {ShortDualFlag(0,0), ShortDualFlag(0,1)};
        standard_matrix[M - 2][0] = 1.0; //L_x <= 1.0
        standard_matrix[M - 2][1] = 0.0; 
        standard_matrix[M - 2][2] = 1.0;
        standard_matrix[M - 1][0] = 0.0; //L_y <= 1.0
        standard_matrix[M - 1][1] = 1.0; 
        standard_matrix[M - 1][2] = 1.0;
//        standard_matrix[M][0] = -fabs(_gradient.xx); //Obj = |D_x| * L_x + |D_y| * L_y 
//        standard_matrix[M][1] = -fabs(_gradient.xy); 
//        standard_matrix[M][2] = 0.0;
        Vector2 centre_this = cell_centres[index];
        Vector2 direction;
        Vector2 phi_this = phi_data[index];
        Scalar phi_this_x = phi_this.x;
        Scalar c_x = 0.0;
        Scalar c_y = 0.0;
        for(Flag i = 0; i < cell_neigbours_number; ++i){
          ShortDualHandle neib = cell_neigbours[i*cell_neighbours_length+index];
          Scalar d_phi;
          if(!neib.is_boundary()){
            Flag id_neib = neib.get_global_id();
            Scalar phi_neib_x = phi_data[id_neib].x;
            d_phi = phi_neib_x - phi_this_x;
            Vector2 centre_neib = cell_centres[id_neib];
            direction = centre_neib - centre_this;
//            Flag id_face = cell_halffacets[i*cell_halffacets_length + index];
//            Vector2 centre_face = face_centres[id_face];
//            direction = centre_face - centre_this;
          }else{
            Flag id_face = cell_halffacets[i*cell_halffacets_length+index];
            Vector2 centre_face = face_centres[id_face];
            direction = 2.0*(centre_face - centre_this);
            Vector2 normal = uni(face_normals[id_face]);
            Flag id_boundary = neib.get_global_id();
//            ShortTripleFlag phi_boundary_type = phi_bound[id_boundary];
//            Vector2 phi_bound = cuGetBoundary(phi_boundary_type, phi_this, normal);
            Vector2 phi_bound = _phi_bound[id_boundary];
            d_phi = phi_bound.x - phi_this.x;
          }
          if(d_phi > 0){
            standard_matrix[2*i][0] = _gradient.xx*direction.x;
            standard_matrix[2*i][1] = _gradient.xy*direction.y;
            standard_matrix[2*i][2] = d_phi;
            standard_matrix[2*i+1][0] = -_gradient.xx*direction.x;
            standard_matrix[2*i+1][1] = -_gradient.xy*direction.y;
            standard_matrix[2*i+1][2] = 0.0;
          }else{
            standard_matrix[2*i][0] = -_gradient.xx*direction.x;
            standard_matrix[2*i][1] = -_gradient.xy*direction.y;
            standard_matrix[2*i][2] = -d_phi;
            standard_matrix[2*i+1][0] = _gradient.xx*direction.x;
            standard_matrix[2*i+1][1] = _gradient.xy*direction.y;
            standard_matrix[2*i+1][2] = 0.0;
          }
          c_x -= fabs(_gradient.xx*direction.x);
          c_y -= fabs(_gradient.xy*direction.y);
        }
        standard_matrix[M][0] = c_x; //Obj = sum|(D_x*u_xi)|*L_x + sum|(D_y*u_yi)|*L_x
        standard_matrix[M][1] = c_y; 
        standard_matrix[M][2] = 0.0;
        double L_x = 0.0;
        double L_y = 0.0;
        __syncthreads();
        simplex(standard_matrix, M, x_vector, y_vector, x_pos, L_x, L_y);
        _gradient.xx *= L_x;
        _gradient.xy *= L_y;

      //-----------------------Y axis------------------------------

        //Short index: 0--x-variable 1--y-variable, Long index: index
        //ShortDualFlag x_vector[2] = {ShortDualFlag(0,0), ShortDualFlag(0,1)};
        //ShortDualFlag y_vector[10] = {ShortDualFlag(1,0), ShortDualFlag(1,1), ShortDualFlag(1,2), ShortDualFlag(1,3), ShortDualFlag(1,4), 
        //                              ShortDualFlag(1,5), ShortDualFlag(1,6), ShortDualFlag(1,7), ShortDualFlag(1,8), ShortDualFlag(1,9)};
        //Short index: 0--x-vector 1--y-vector, Long index: index
        //ShortDualFlag x_pos[2] = {ShortDualFlag(0,0), ShortDualFlag(0,1)};
        for(int i = 0; i < 2; i++){
          x_vector[i] = ShortDualFlag(0,i);
          x_pos[i] = ShortDualFlag(0,i);
        }
        for(int i = 0; i < 10; i++){
          y_vector[i] = ShortDualFlag(1,i);
        }
        standard_matrix[M - 2][0] = 1.0; //L_x <= 1.0
        standard_matrix[M - 2][1] = 0.0; 
        standard_matrix[M - 2][2] = 1.0;
        standard_matrix[M - 1][0] = 0.0; //L_y <= 1.0
        standard_matrix[M - 1][1] = 1.0; 
        standard_matrix[M - 1][2] = 1.0;
//        standard_matrix[M][0] = -fabs(_gradient.yx); //Obj = |D_x| * L_x + |D_y| * L_y 
//        standard_matrix[M][1] = -fabs(_gradient.yy); 
//        standard_matrix[M][2] = 0.0;
//        Vector2 centre_this = cell_centres[index];
//        Vector2 direction;
//        Vector2 phi_this = phi_data[index];
        Scalar phi_this_y = phi_this.y;
        c_x = 0.0;
        c_y = 0.0;
        for(Flag i = 0; i < cell_neigbours_number; ++i){
          ShortDualHandle neib = cell_neigbours[i*cell_neighbours_length+index];
          Scalar d_phi;
          if(!neib.is_boundary()){
            Flag id_neib = neib.get_global_id();
            Scalar phi_neib_y = phi_data[id_neib].y;
            d_phi = phi_neib_y - phi_this_y;
            Vector2 centre_neib = cell_centres[id_neib];
            direction = centre_neib - centre_this;
//            Flag id_face = cell_halffacets[i*cell_halffacets_length + index];
//            Vector2 centre_face = face_centres[id_face];
//            direction = centre_face - centre_this;
          }else{
            Flag id_face = cell_halffacets[i*cell_halffacets_length+index];
            Vector2 centre_face = face_centres[id_face];
            direction = 2.0*(centre_face - centre_this);
            Vector2 normal = uni(face_normals[id_face]);
            Flag id_boundary = neib.get_global_id();
//            ShortTripleFlag phi_boundary_type = phi_bound[id_boundary];
//            Vector2 phi_bound = cuGetBoundary(phi_boundary_type, phi_this, normal);
            Vector2 phi_bound = _phi_bound[id_boundary];
            d_phi = phi_bound.y - phi_this.y;
          }
          if(d_phi > 0){
            standard_matrix[2*i][0] = _gradient.yx*direction.x;
            standard_matrix[2*i][1] = _gradient.yy*direction.y;
            standard_matrix[2*i][2] = d_phi;
            standard_matrix[2*i+1][0] = -_gradient.yx*direction.x;
            standard_matrix[2*i+1][1] = -_gradient.yy*direction.y;
            standard_matrix[2*i+1][2] = 0.0;
          }else{
            standard_matrix[2*i][0] = -_gradient.yx*direction.x;
            standard_matrix[2*i][1] = -_gradient.yy*direction.y;
            standard_matrix[2*i][2] = -d_phi;
            standard_matrix[2*i+1][0] = _gradient.yx*direction.x;
            standard_matrix[2*i+1][1] = _gradient.yy*direction.y;
            standard_matrix[2*i+1][2] = 0.0;
          }
          c_x -= fabs(_gradient.yx*direction.x);
          c_y -= fabs(_gradient.yy*direction.y);
        }
        standard_matrix[M][0] = c_x; //Obj = sum|(D_x*u_xi)|*L_x + sum|(D_y*u_yi)|*L_x
        standard_matrix[M][1] = c_y; 
        standard_matrix[M][2] = 0.0;
        L_x = 0.0;
        L_y = 0.0;
        __syncthreads();
        simplex(standard_matrix, M, x_vector, y_vector, x_pos, L_x, L_y);
        _gradient.yx *= L_x;
        _gradient.yy *= L_y;
        phi_gradient[index] = _gradient;
        index += blockDim.x * gridDim.x;
      }
    }


    void cuGradientLimiter(cuFvMappedField<Scalar, on_cell>& phi, cuFvMappedField<Vector2, on_cell>& phi_gradient){

      auto mesh = phi.mesh;

      cuGradientLimiterKernelNonSimplex<<<BLOCKS_PER_GRID, 128>>>(phi.data.dev_ptr(), phi.boundary_value.dev_ptr(), phi.data.size(), mesh->cell_neighbours.dims_dev_ptr(), mesh->cell_neighbours.dev_ptr(), mesh->cell_neighbours.length(), mesh->cell_centre_positions.dev_ptr(), mesh->cell_halffacets.dev_ptr(), mesh->cell_halffacets.length(), mesh->halffacet_centre_positions.dev_ptr(), mesh->halffacet_normal_directions.dev_ptr(), phi_gradient.data.dev_ptr(), phi_gradient.data.size()); 

    }


    void cuGradientLimiter(cuFvMappedField<Vector2, on_cell>& phi, cuFvMappedField<Tensor2, on_cell>& phi_gradient){

      auto mesh = phi.mesh;

      cuGradientLimiterKernelNonSimplex<<<BLOCKS_PER_GRID, 128>>>(phi.data.dev_ptr(), phi.boundary_value.dev_ptr(), phi.data.size(), mesh->cell_neighbours.dims_dev_ptr(), mesh->cell_neighbours.dev_ptr(), mesh->cell_neighbours.length(), mesh->cell_centre_positions.dev_ptr(), mesh->cell_halffacets.dev_ptr(), mesh->cell_halffacets.length(), mesh->halffacet_centre_positions.dev_ptr(), mesh->halffacet_normal_directions.dev_ptr(), phi_gradient.data.dev_ptr(), phi_gradient.data.size()); 

    }

    __device__ Scalar minmod(const Scalar& upwind, const Scalar& downwind){
      Scalar ratio = 0.0;
      if (fabs(downwind) <= 1e-10){
        ratio = 1.0;
      }
      else{
        ratio = upwind / downwind;
      }

      ratio = fmax(0.0, fmin(1.0, ratio));
      return ratio;
    }

    __global__ void cuGradientLimiterCartesianKernel(Scalar* phi_data, Scalar* _phi_bound, unsigned int phi_size, Flag* cell_neigbours_dimensions, ShortDualHandle* cell_neigbours, unsigned int cell_neighbours_length, Vector2* cell_centres, Flag* cell_halffacets, unsigned int cell_halffacets_length, Vector2* face_centres, Vector2* face_normals, Vector2* phi_gradient, unsigned int phi_gradient_size){
      
      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      while (index < phi_size){
        Scalar phi_this = phi_data[index];
        Vector2 centroid_this = cell_centres[index];
        Vector2 gradient_vector = phi_gradient[index];
        Scalar limit_ratio_x = 1.0;
        Scalar limit_ratio_y = 1.0;
        Scalar _limit_ratio_x = 1.0;
        Scalar _limit_ratio_y = 1.0;
        Flag cell_neigbours_number = cell_neigbours_dimensions[index];
        for (Flag i = 0; i < cell_neigbours_number; ++i){
          ShortDualHandle neib = cell_neigbours[i*cell_neighbours_length + index];
          if (i == 1 || i == 3){ //x-axis
            if (!neib.is_boundary()){ //neighbour is not a boundary
              Flag id_neib = neib.get_global_id();
              Scalar phi_neib = phi_data[id_neib];
              Vector2 centroid_neib = cell_centres[id_neib];
              Vector2 direction_vector = centroid_neib - centroid_this;
              Scalar downwind = dot(gradient_vector, direction_vector);
              Scalar upwind = phi_neib - phi_this;
              _limit_ratio_x = minmod(upwind, downwind);
            }
            else{
              Flag id_boundary = neib.get_global_id();
              Flag id_face = cell_halffacets[i*cell_halffacets_length + index];;
              Vector2 normal = uni(face_normals[id_face]);
//              ShortTripleFlag phi_boundary_type = phi_bound[id_boundary];            
//              Scalar phi_bound = cuGetBoundary(phi_boundary_type, phi_this, normal);
              Scalar phi_bound = _phi_bound[id_boundary];
              Vector2 centroid_boundary = face_centres[id_face];
              Vector2 direction_vector = centroid_boundary - centroid_this;
              Scalar downwind = 2.0*dot(gradient_vector, direction_vector);
              Scalar upwind = phi_bound - phi_this;
              _limit_ratio_x = minmod(upwind, downwind);
            }
            limit_ratio_x = fmin(limit_ratio_x, _limit_ratio_x);
          }
          else{
            if (!neib.is_boundary()){ //neighbour is not a boundary
              Flag id_neib = neib.get_global_id();
              Scalar phi_neib = phi_data[id_neib];
              Vector2 centroid_neib = cell_centres[id_neib];
              Vector2 direction_vector = centroid_neib - centroid_this;
              Scalar downwind = dot(gradient_vector, direction_vector);
              Scalar upwind = phi_neib - phi_this;
              _limit_ratio_y = minmod(upwind, downwind);
            }
            else{
              Flag id_boundary = neib.get_global_id();
              Flag id_face = cell_halffacets[i*cell_halffacets_length + index];;
              Vector2 normal = uni(face_normals[id_face]);
//              ShortTripleFlag phi_boundary_type = phi_bound[id_boundary];
//              Scalar phi_bound = cuGetBoundary(phi_boundary_type, phi_this, normal);
              Scalar phi_bound = _phi_bound[id_boundary];
              Vector2 centroid_boundary = face_centres[id_face];
              Vector2 direction_vector = centroid_boundary - centroid_this;
              Scalar downwind = 2.0*dot(gradient_vector, direction_vector);
              Scalar upwind = phi_bound - phi_this;
              _limit_ratio_y = minmod(upwind, downwind);
            }
            limit_ratio_y = fmin(limit_ratio_y, _limit_ratio_y);
          }
        }
        phi_gradient[index] = Vector2(gradient_vector.x*limit_ratio_x, gradient_vector.y*limit_ratio_y);
        index += blockDim.x * gridDim.x;
      }   
    }



    void cuGradientLimiterCartesian(cuFvMappedField<Scalar, on_cell>& phi, cuFvMappedField<Vector2, on_cell>& phi_gradient){

      auto mesh = phi.mesh;

      cuGradientLimiterCartesianKernel <<<BLOCKS_PER_GRID, THREADS_PER_BLOCK >>>(phi.data.dev_ptr(), phi.boundary_value.dev_ptr(), phi.data.size(), mesh->cell_neighbours.dims_dev_ptr(), mesh->cell_neighbours.dev_ptr(), mesh->cell_neighbours.length(), mesh->cell_centre_positions.dev_ptr(), mesh->cell_halffacets.dev_ptr(), mesh->cell_halffacets.length(), mesh->halffacet_centre_positions.dev_ptr(), mesh->halffacet_normal_directions.dev_ptr(), phi_gradient.data.dev_ptr(), phi_gradient.data.size());

    }

    __global__ void cuGradientLimiterCartesianKernel(Vector2* phi_data, Vector2* _phi_bound, unsigned int phi_size, Flag* cell_neigbours_dimensions, ShortDualHandle* cell_neigbours, unsigned int cell_neighbours_length, Vector2* cell_centres, Flag* cell_halffacets, unsigned int cell_halffacets_length, Vector2* face_centres, Vector2* face_normals, Tensor2* phi_gradient, unsigned int phi_gradient_size){

      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      while (index < phi_size){
        Vector2 phi_this = phi_data[index];
        Vector2 centroid_this = cell_centres[index];
        Tensor2 gradient_vector = phi_gradient[index];
        Scalar limit_ratio_xx = 1.0;
        Scalar limit_ratio_xy = 1.0;
        Scalar _limit_ratio_xx = 1.0;
        Scalar _limit_ratio_xy = 1.0;
        Scalar limit_ratio_yx = 1.0;
        Scalar limit_ratio_yy = 1.0;
        Scalar _limit_ratio_yx = 1.0;
        Scalar _limit_ratio_yy = 1.0;
        Flag cell_neigbours_number = cell_neigbours_dimensions[index];
        for (Flag i = 0; i < cell_neigbours_number; ++i){
          ShortDualHandle neib = cell_neigbours[i*cell_neighbours_length + index];
          if (i == 1 || i == 3){ //x-axis
            if (!neib.is_boundary()){ //neighbour is not a boundary
              Flag id_neib = neib.get_global_id();
              Vector2 phi_neib = phi_data[id_neib];
              Vector2 centroid_neib = cell_centres[id_neib];
              Vector2 direction_vector = centroid_neib - centroid_this;
              Vector2 downwind = dot(gradient_vector, direction_vector);
              Vector2 upwind = phi_neib - phi_this;
              _limit_ratio_xx = minmod(upwind.x, downwind.x);
              _limit_ratio_yx = minmod(upwind.y, downwind.y);
            }
            else{
              Flag id_boundary = neib.get_global_id();
              Flag id_face = cell_halffacets[i*cell_halffacets_length + index];;
              Vector2 normal = uni(face_normals[id_face]);
//              ShortTripleFlag phi_boundary_type = phi_bound[id_boundary];
//              Vector2 phi_bound = cuGetBoundary(phi_boundary_type, phi_this, normal);
              Vector2 phi_bound = _phi_bound[id_boundary];
              Vector2 centroid_boundary = face_centres[id_face];
              Vector2 direction_vector = centroid_boundary - centroid_this;
              Vector2 downwind = 2.0*dot(gradient_vector, direction_vector);
              Vector2 upwind = phi_bound - phi_this;
              _limit_ratio_xx = minmod(upwind.x, downwind.x);
              _limit_ratio_yx = minmod(upwind.y, downwind.y);
            }
            limit_ratio_xx = fmin(limit_ratio_xx, _limit_ratio_xx);
            limit_ratio_yx = fmin(limit_ratio_yx, _limit_ratio_yx);
          }
          else{
            if (!neib.is_boundary()){ //neighbour is not a boundary
              Flag id_neib = neib.get_global_id();
              Vector2 phi_neib = phi_data[id_neib];
              Vector2 centroid_neib = cell_centres[id_neib];
              Vector2 direction_vector = centroid_neib - centroid_this;
              Vector2 downwind = dot(gradient_vector, direction_vector);
              Vector2 upwind = phi_neib - phi_this;
              _limit_ratio_xy = minmod(upwind.x, downwind.x);
              _limit_ratio_yy = minmod(upwind.y, downwind.y);
            }
            else{
              Flag id_boundary = neib.get_global_id();
              Flag id_face = cell_halffacets[i*cell_halffacets_length + index];;
              Vector2 normal = uni(face_normals[id_face]);
//              ShortTripleFlag phi_boundary_type = phi_bound[id_boundary];
//              Vector2 phi_bound = cuGetBoundary(phi_boundary_type, phi_this, normal);
              Vector2 phi_bound = _phi_bound[id_boundary];
              Vector2 centroid_boundary = face_centres[id_face];
              Vector2 direction_vector = centroid_boundary - centroid_this;
              Vector2 downwind = 2.0*dot(gradient_vector, direction_vector);
              Vector2 upwind = phi_bound - phi_this;
              _limit_ratio_xy = minmod(upwind.x, downwind.x);
              _limit_ratio_yy = minmod(upwind.y, downwind.y);
            }
            limit_ratio_xy = fmin(limit_ratio_xy, _limit_ratio_xy);
            limit_ratio_yy = fmin(limit_ratio_yy, _limit_ratio_yy);
          }
        }
        phi_gradient[index] = Tensor2(limit_ratio_xx*gradient_vector.xx,
                                      limit_ratio_xy*gradient_vector.xy,
                                      limit_ratio_yx*gradient_vector.yx,
                                      limit_ratio_yy*gradient_vector.yy);
        index += blockDim.x * gridDim.x;
      }
    }

    void cuGradientLimiterCartesian(cuFvMappedField<Vector2, on_cell>& phi, cuFvMappedField<Tensor2, on_cell>& phi_gradient){

      auto mesh = phi.mesh;

      cuGradientLimiterCartesianKernel <<<BLOCKS_PER_GRID, THREADS_PER_BLOCK >>>(phi.data.dev_ptr(), phi.boundary_value.dev_ptr(), phi.data.size(), mesh->cell_neighbours.dims_dev_ptr(), mesh->cell_neighbours.dev_ptr(), mesh->cell_neighbours.length(), mesh->cell_centre_positions.dev_ptr(), mesh->cell_halffacets.dev_ptr(), mesh->cell_halffacets.length(), mesh->halffacet_centre_positions.dev_ptr(), mesh->halffacet_normal_directions.dev_ptr(), phi_gradient.data.dev_ptr(), phi_gradient.data.size());

    }


  }//end of namespace fv-------

}//end of namespace GC---------

