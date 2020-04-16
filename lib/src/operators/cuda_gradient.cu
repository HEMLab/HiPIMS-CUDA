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
\file cuda_gradient.cu
\brief Source file for gradient operator by cuda

\version 0.1
\author xilin xia

*/

#include "cuda_gradient.h"
#include "cuda_kernel_launch_parameters.h"
//#include "cuda_boundary.h"


namespace GC{

  namespace fv{

    __global__ void cuGradientKernel(Scalar* phi_data, Scalar* phi_bound, unsigned int phi_size, Flag* cell_neigbours_dimensions, ShortDualHandle* cell_neigbours, unsigned int cell_neighbours_length, Flag* cell_halffacets, unsigned int cell_halffacets_length, Scalar* cell_volume,  Scalar* face_area, Vector* face_normal, Vector* phi_gradient){

      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x; 
      while(index < phi_size){
        Scalar phi_this = phi_data[index];
        Scalar volume_this = cell_volume[index];
        Flag cell_neigbours_number = cell_neigbours_dimensions[index];
		    Vector _phi_gradient;
        for(Flag i = 0; i < cell_neigbours_number; ++i){
          ShortDualHandle neib = cell_neigbours[i*cell_neighbours_length+index];
          if(!neib.is_boundary()){
            Flag id_neib = neib.get_global_id();
            Flag id_face = cell_halffacets[i*cell_halffacets_length+index];
            Vector normal = uni(face_normal[id_face]);
            Scalar area = face_area[id_face];
            Scalar phi_neib = phi_data[id_neib];
            _phi_gradient += (phi_neib+phi_this)/2.0*area*normal/volume_this;
          }else{
            Flag id_face = cell_halffacets[i*cell_halffacets_length+index];
            Vector normal = uni(face_normal[id_face]);
            Scalar area = face_area[id_face];
            Flag id_boundary = neib.get_global_id();
//            ShortTripleFlag phi_boundary_type = phi_bound[id_boundary];
            Scalar _phi_bound = phi_bound[id_boundary];
            _phi_gradient += (_phi_bound+phi_this)/2.0*area*normal/volume_this;
          }
        }
        phi_gradient[index] = _phi_gradient;
        index += blockDim.x * gridDim.x;
      }
    }


    void cuGradient(cuFvMappedField<Scalar, on_cell>& phi, cuFvMappedField<Vector, on_cell>& phi_gradient){

      auto mesh = phi.mesh;

      cuGradientKernel<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK>>>(phi.data.dev_ptr(), phi.boundary_value.dev_ptr(), phi.data.size(), mesh->cell_neighbours.dims_dev_ptr(), mesh->cell_neighbours.dev_ptr(), mesh->cell_neighbours.length(), mesh->cell_halffacets.dev_ptr(), mesh->cell_halffacets.length(), mesh->cell_volumes.dev_ptr(),  mesh->halffacet_areas.dev_ptr(), mesh->halffacet_normal_directions.dev_ptr(),  phi_gradient.data.dev_ptr());

    }

    __global__ void cuGradientKernel(Vector* phi_data, Vector* phi_bound, unsigned int phi_size, Flag* cell_neigbours_dimensions, ShortDualHandle* cell_neigbours, unsigned int cell_neighbours_length, Flag* cell_halffacets, unsigned int cell_halffacets_length, Scalar* cell_volume,  Scalar* face_area, Vector* face_normal, Tensor* phi_gradient){
      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x; 
      while(index < phi_size){
        Vector phi_this = phi_data[index];
        Scalar volume_this = cell_volume[index];
        Flag cell_neigbours_number = cell_neigbours_dimensions[index];
		    Tensor _phi_gradient;
        for(Flag i = 0; i < cell_neigbours_number; ++i){
          ShortDualHandle neib = cell_neigbours[i*cell_neighbours_length+index];
          if(!neib.is_boundary()){
            Flag id_neib = neib.get_global_id();
            Flag id_face = cell_halffacets[i*cell_halffacets_length+index];
            Vector normal = uni(face_normal[id_face]);
            Scalar area = face_area[id_face];
            Vector phi_neib = phi_data[id_neib];
            _phi_gradient += product((phi_neib+phi_this)/2.0, normal)*area/volume_this;
          }else{
            Flag id_face = cell_halffacets[i*cell_halffacets_length+index];
            Vector normal = uni(face_normal[id_face]);
            Scalar area = face_area[id_face];
            Flag id_boundary = neib.get_global_id();
//            ShortTripleFlag phi_boundary_type = phi_bound[id_boundary];
            Vector _phi_bound = phi_bound[id_boundary];
            _phi_gradient += product((_phi_bound+phi_this)/2.0, normal)*area/volume_this;
          }
        }
        phi_gradient[index] = _phi_gradient;
        index += blockDim.x * gridDim.x;
      }

    }

    void cuGradient(cuFvMappedField<Vector, on_cell>& phi, cuFvMappedField<Tensor, on_cell>& phi_gradient){

      auto mesh = phi.mesh;

      cuGradientKernel<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK>>>(phi.data.dev_ptr(), phi.boundary_value.dev_ptr(), phi.data.size(), mesh->cell_neighbours.dims_dev_ptr(), mesh->cell_neighbours.dev_ptr(), mesh->cell_neighbours.length(), mesh->cell_halffacets.dev_ptr(), mesh->cell_halffacets.length(), mesh->cell_volumes.dev_ptr(),  mesh->halffacet_areas.dev_ptr(), mesh->halffacet_normal_directions.dev_ptr(),  phi_gradient.data.dev_ptr());

    }

    __device__ Scalar minmod_limiter(Scalar a, Scalar b){
        
      if (a*b <= 0.0){
        return 0.0;
      }
      else if(fabs(a) < fabs(b)){
        return a;
      }
      else{
        return b;
      }

    }

    __device__ Scalar vanLeer_limiter(Scalar a, Scalar b){
      if (a*b <= 0.0){
        return 0.0;
      }
      else{
        return 1.0 / (0.5*(1.0 / a + 1.0 / b));
      }
    }

    __global__ void cuLimitedGradientCartesianKernel(Scalar* phi_data, Scalar* phi_bound, unsigned int phi_size, ShortDualHandle* cell_neigbours, unsigned int cell_neighbours_length, Scalar* cell_volume, Vector* phi_gradient){

      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      while (index < phi_size){
        Scalar phi_this = phi_data[index];
        Scalar volume_this = cell_volume[index];
        Scalar cell_size = sqrt(volume_this);
        Scalar grad_x;
        Scalar grad_y;
        Scalar phi_neib[4];
        for (Flag i = 0; i < 4; ++i){
          ShortDualHandle neib = cell_neigbours[i*cell_neighbours_length + index];
          if (!neib.is_boundary()){
            Flag id_neib = neib.get_global_id();
            phi_neib[i] = phi_data[id_neib];
          }
          else{
            Flag id_boundary = neib.get_global_id();
            phi_neib[i] = phi_bound[id_boundary];
          }
        }
        Scalar grad_x_up = (phi_this - phi_neib[3])/cell_size;
        Scalar grad_x_down = (phi_neib[1] - phi_this)/cell_size;
        grad_x = minmod_limiter(grad_x_up, grad_x_down);
        Scalar grad_y_up = (phi_this - phi_neib[0]) / cell_size;
        Scalar grad_y_down = (phi_neib[2] - phi_this) / cell_size;
        grad_y = minmod_limiter(grad_y_up, grad_y_down);
        phi_gradient[index] = Vector2(grad_x, grad_y);
        index += blockDim.x * gridDim.x;
      }
    }


    void cuLimitedGradientCartesian(cuFvMappedField<Scalar, on_cell>& phi, cuFvMappedField<Vector, on_cell>& phi_gradient){

      auto mesh = phi.mesh;

      cuLimitedGradientCartesianKernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(phi.data.dev_ptr(), 
      phi.boundary_value.dev_ptr(), 
      phi.data.size(), 
      mesh->cell_neighbours.dev_ptr(), 
      mesh->cell_neighbours.length(), 
      mesh->cell_volumes.dev_ptr(), 
      phi_gradient.data.dev_ptr());

    }

    __global__ void cuLimitedGradientCartesianKernel(Vector2* phi_data, Vector2* phi_bound, unsigned int phi_size, ShortDualHandle* cell_neigbours, unsigned int cell_neighbours_length, Scalar* cell_volume, Tensor2* phi_gradient){

      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      while (index < phi_size){
        Vector2 phi_this = phi_data[index];
        Scalar volume_this = cell_volume[index];
        Scalar cell_size = sqrt(volume_this);
        Scalar grad_xx;
        Scalar grad_xy;
        Scalar grad_yx;
        Scalar grad_yy;
        Vector2 phi_neib[4];
        for (Flag i = 0; i < 4; ++i){
          ShortDualHandle neib = cell_neigbours[i*cell_neighbours_length + index];
          if (!neib.is_boundary()){
            Flag id_neib = neib.get_global_id();
            phi_neib[i] = phi_data[id_neib];
//            if (index == 31494){
//              printf("Gradient %d phi_x:%e phi_y:%e\n",i, phi_neib[i].x, phi_neib[i].y);
//            }
          }
          else{
            Flag id_boundary = neib.get_global_id();
            phi_neib[i] = phi_bound[id_boundary];
          }
        }
        Scalar grad_xx_up = (phi_this.x - phi_neib[3].x) / cell_size;
        Scalar grad_xx_down = ((phi_neib[1]).x - phi_this.x) / cell_size;
        grad_xx = minmod_limiter(grad_xx_up, grad_xx_down);
        Scalar grad_xy_up = (phi_this.x - phi_neib[0].x) / cell_size;
        Scalar grad_xy_down = (phi_neib[2].x - phi_this.x) / cell_size;
        grad_xy = minmod_limiter(grad_xy_up, grad_xy_down);
        Scalar grad_yx_up = (phi_this.y - phi_neib[3].y) / cell_size;
        Scalar grad_yx_down = ((phi_neib[1]).y - phi_this.y) / cell_size;
        grad_yx = minmod_limiter(grad_yx_up, grad_yx_down);
        Scalar grad_yy_up = (phi_this.y - phi_neib[0].y) / cell_size;
        Scalar grad_yy_down = (phi_neib[2].y - phi_this.y) / cell_size;
        grad_yy = minmod_limiter(grad_yy_up, grad_yy_down);
//        if (index == 31494){
//          printf("Gradient G_xx:%e G_xy:%e G_yx:%e G_yy:%e\n", grad_xx, grad_yx, grad_xy, grad_yy);
//        }
        phi_gradient[index] = Tensor2(grad_xx, grad_xy, grad_yx, grad_yy);
        index += blockDim.x * gridDim.x;
      }
    }

    void cuLimitedGradientCartesian(cuFvMappedField<Vector2, on_cell>& phi, cuFvMappedField<Tensor2, on_cell>& phi_gradient){

      auto mesh = phi.mesh;

      cuLimitedGradientCartesianKernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(phi.data.dev_ptr(),
        phi.boundary_value.dev_ptr(),
        phi.data.size(),
        mesh->cell_neighbours.dev_ptr(),
        mesh->cell_neighbours.length(),
        mesh->cell_volumes.dev_ptr(),
        phi_gradient.data.dev_ptr());

    }

  }//end of namespace fv-----------

}//--end of namespace GC-----------
