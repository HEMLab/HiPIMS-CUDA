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
\file cuda_hessian.cu
\brief Source file for hessian operator by cuda

\version 0.1
\author xilin xia

*/

#include "cuda_hessian.h"
#include "cuda_kernel_launch_parameters.h"

namespace GC{

  namespace fv{

    __global__ void cuHessianCartesian2DKernel(Scalar* phi_data, Scalar* phi_bound, unsigned int phi_size, Flag* cell_neigbours_dimensions, ShortDualHandle* cell_neigbours, unsigned int cell_neighbours_length, Flag* cell_halffacets, unsigned int cell_halffacets_length, Scalar* cell_volume,  Scalar* face_area, Vector* face_normal, Tensor2* phi_hessian){

      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x; 
      while(index < phi_size){
        Scalar phi_this = phi_data[index];
        Scalar volume_this = cell_volume[index];
        Scalar dx = sqrt(volume_this);
        Scalar dy = sqrt(volume_this);
        //Flag cell_neigbours_number = cell_neigbours_dimensions[index];
		    //Tensor2 _phi_hessian;
        Scalar phi_neib[4];
        Scalar phi_corner[4];
        for(Flag i = 0; i < 4; ++i){
          ShortDualHandle neib = cell_neigbours[i*cell_neighbours_length+index];
          if(!neib.is_boundary()){
            Flag id_neib = neib.get_global_id();
            phi_neib[i] = phi_data[id_neib];
          }else{
            Flag id_face = cell_halffacets[i*cell_halffacets_length+index];
            Vector normal = uni(face_normal[id_face]);
            Flag id_boundary = neib.get_global_id();
//            ShortTripleFlag phi_boundary_type = phi_bound[id_boundary];
//            phi_neib[i] = cuGetBoundary(phi_boundary_type, phi_this, normal);
            phi_neib[i] = phi_bound[id_boundary];
          }
        }
        for(Flag i = 0; i < 4; ++i){
          Flag j = (i + 3) % 4;
          ShortDualHandle neib = cell_neigbours[i*cell_neighbours_length+index];
          if(!neib.is_boundary()){
            Flag id_neib = neib.get_global_id();
            ShortDualHandle corner = cell_neigbours[j*cell_neighbours_length+id_neib];
            if(!corner.is_boundary()){
              Flag id_corner = corner.get_global_id();
              phi_corner[i] = phi_data[id_corner];
            }else{
              Flag id_face = cell_halffacets[j*cell_halffacets_length+id_neib];
              Vector normal = uni(face_normal[id_face]);
              Flag id_boundary = corner.get_global_id();
//              ShortTripleFlag phi_boundary_type = phi_bound[id_boundary];
//              phi_corner[i] = cuGetBoundary(phi_boundary_type, phi_neib[i], normal);
              phi_corner[i] = phi_bound[id_boundary];
            }
          }else{
            Flag id_face = cell_halffacets[i*cell_halffacets_length+index];
            Vector normal = uni(face_normal[id_face]);
            Flag id_boundary = neib.get_global_id();
//            ShortTripleFlag phi_boundary_type = phi_bound[id_boundary];
//            phi_corner[i] = cuGetBoundary(phi_boundary_type, phi_this, normal);
            phi_corner[i] = phi_bound[id_boundary];
          }
        }
        Scalar d_phi_xx = (phi_neib[1] - 2.0*phi_this + phi_neib[3])/dx/dx;
        Scalar d_phi_yy = (phi_neib[2] - 2.0*phi_this + phi_neib[0])/dy/dy;
        Scalar d_phi_xy = (phi_corner[2] - phi_corner[1] - phi_corner[3] + phi_corner[0])/4.0/dx/dy;
        phi_hessian[index] = Tensor2(d_phi_xx, d_phi_xy, d_phi_xy, d_phi_yy);
        index += blockDim.x * gridDim.x;
      }
    }



    void cuHessianCartesian2D(cuFvMappedField<Scalar, on_cell>& phi, cuFvMappedField<Tensor2, on_cell>& phi_hessian){

      auto mesh = phi.mesh;

      cuHessianCartesian2DKernel<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK>>>(phi.data.dev_ptr(), phi.boundary_value.dev_ptr(), phi.data.size(), mesh->cell_neighbours.dims_dev_ptr(), mesh->cell_neighbours.dev_ptr(), mesh->cell_neighbours.length(), mesh->cell_halffacets.dev_ptr(), mesh->cell_halffacets.length(), mesh->cell_volumes.dev_ptr(),  mesh->halffacet_areas.dev_ptr(), mesh->halffacet_normal_directions.dev_ptr(),  phi_hessian.data.dev_ptr());

    }


  }

}
