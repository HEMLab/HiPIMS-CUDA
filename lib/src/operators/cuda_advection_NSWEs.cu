// ======================================================================================
// Name                :    High-Performance Integrated Modelling System
// Description         :    This code pack provides a generic framework for developing 
//                          Geophysical CFD software. Legacy name: GeoClasses
// ======================================================================================
// Version             :    1.0.1 
// Author              :    Xilin Xia
// Create Time         :    2014/10/04
// Update Time         :    2020/04/19
// ======================================================================================
// LICENCE: GPLv3 
// ======================================================================================

/*!
  \flie cuda_advection_NSWEs.cu
  \brief Source file for advection of non hydrostatic shallow water equations

  \version 0.1
  \author xilin xia
*/

#include "cuda_advection_NSWEs.h"
#include "cuda_kernel_launch_parameters.h"
//#include "cuda_boundary.h"
#include "riemann.h"


namespace GC{


  namespace fv{

    
   #include "cuda_riemann_solvers.cuh"


    __global__ void cuAdvectionNSWEs2ndKernel(Scalar* gravity, Scalar* _gravity_bound, Scalar* h, Scalar* _h_bound, Scalar* z, Scalar* _z_bound, Vector* u, Vector* _u_bound, Vector* h_gradient, Vector* eta_gradient, Tensor* u_gradient, unsigned int phi_size, Flag* cell_neigbours_dimensions, ShortDualHandle* cell_neigbours, unsigned int cell_neighbours_length, Flag* cell_halffacets, unsigned int cell_halffacets_length, Vector* cell_centres, Scalar* cell_volume,  Scalar* face_area, Vector* face_normal, Vector* face_centres, Scalar* h_advection, Vector* hU_advection){

      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      while(index < phi_size){
        Scalar g_this = gravity[index];
        Scalar h_this = h[index];
        Scalar z_this = z[index];
        Scalar eta_this = h_this + z_this;
        Vector2 u_this = u[index];
        Vector2 grad_h_this = h_gradient[index];
        Vector2 grad_eta_this = eta_gradient[index];
        Tensor2 grad_u_this = u_gradient[index];
        Scalar volume_this =  cell_volume[index];
        Vector2 cell_centroid_this = cell_centres[index];
        Flag cell_neigbours_number = cell_neigbours_dimensions[index];
        Scalar _h_advection(0.0);
        Vector2 _hU_advection(0.0,0.0);
        for(Flag i = 0; i < cell_neigbours_number; ++i){
          Flag id_face = cell_halffacets[i*cell_halffacets_length+index];
          Vector2 normal = uni(face_normal[id_face]);
          Vector2 shear = uni(perpend(normal));
          Scalar area = face_area[id_face];
          Vector2 face_centroid = face_centres[id_face];
          Vector2 direction_this = face_centroid - cell_centroid_this;
          Scalar _h_this = h_this + dot(grad_h_this,direction_this);
          Scalar _eta_this = eta_this + dot(grad_eta_this,direction_this);
          Scalar _z_this = _eta_this - _h_this;
          Vector2 _u_this = u_this + dot(grad_u_this,direction_this);
          ShortDualHandle neib = cell_neigbours[i*cell_neighbours_length+index];
          if(!neib.is_boundary()){
            Flag id_neib = neib.get_global_id();
            Scalar g_neib = gravity[id_neib];
            Scalar h_neib = h[id_neib];
            Scalar z_neib = z[id_neib];
            Scalar eta_neib = h_neib + z_neib;
            Vector2 u_neib = u[id_neib];
            Vector2 grad_h_neib = h_gradient[id_neib];
            Vector2 grad_eta_neib = eta_gradient[id_neib];
            Tensor2 grad_u_neib = u_gradient[id_neib];
            Vector2 cell_centroid_neib = cell_centres[id_neib];
            Vector2 direction_neib = face_centroid - cell_centroid_neib;
            Scalar _h_neib = h_neib + dot(grad_h_neib, direction_neib);
            Scalar _eta_neib = eta_neib + dot(grad_eta_neib, direction_neib);
            Scalar _z_neib = _eta_neib - _h_neib;
            Vector2 _u_neib = u_neib + dot(grad_u_neib, direction_neib);
            Scalar grav = (g_this + g_neib)/2.0;
            Scalar z_face = fmax(_z_this, _z_neib);
            Scalar h_L = fmax(_h_this + _z_this - z_face, 0.0);
            Scalar h_R = fmax(_h_neib + _z_neib - z_face, 0.0);
            Vector2 u_L(dot(_u_this, normal),dot(_u_this, shear));
            Vector2 u_R(dot(_u_neib, normal),dot(_u_neib, shear));
            auto flux = cuHLLCRiemannSolverSWEs(grav, ScalarRiemannState(h_L, h_R), VectorRiemannState(h_L*u_L, h_R*u_R));
            Vector2 _hU_flux = flux.q.x*normal + flux.q.y*shear;
            Vector2 _z_flux = -0.5*grav*h_L*h_L*normal + 0.5*g_this*(_eta_this - z_this)*(_eta_this - z_this)*normal;
            _h_advection += flux.h*area/volume_this;
            _hU_advection += (_hU_flux + _z_flux)*area/volume_this;
          }else{
            Flag id_boundary = neib.get_global_id();
            Scalar g_bound = _gravity_bound[id_boundary];
            Scalar h_bound = _h_bound[id_boundary];
            Scalar z_bound = _z_bound[id_boundary];
            Scalar grav = (g_this + g_bound)/2.0;
            Scalar z_face = fmax(_z_this, z_bound);
            Scalar h_L = fmax(_h_this + _z_this - z_face, 0.0);
            Scalar h_R = fmax(h_bound + z_bound - z_face, 0.0);
            Vector2 u_bound = _u_bound[id_boundary];
            Vector2 u_L(dot(_u_this, normal),dot(_u_this, shear));
            Vector2 u_R(dot(u_bound, normal),dot(u_bound, shear));
            auto flux = cuHLLCRiemannSolverSWEs(grav, ScalarRiemannState(h_L, h_R), VectorRiemannState(h_L*u_L, h_R*u_R));
            Vector2 _hU_flux = flux.q.x*normal + flux.q.y*shear;
            Vector2 _z_flux = -0.5*grav*h_L*h_L*normal + 0.5*g_this*(_eta_this - z_this)*(_eta_this - z_this)*normal;
            _h_advection += flux.h*area/volume_this;
            _hU_advection += (_hU_flux + _z_flux)*area/volume_this;
          }
        }
        h_advection[index] = _h_advection;
        hU_advection[index] = _hU_advection;
        index += blockDim.x * gridDim.x;
      }
    }

    void cuAdvectionNSWEs2nd(cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Scalar, on_cell>& z, cuFvMappedField<Vector, on_cell>& u, cuFvMappedField<Vector, on_cell>& h_gradient, cuFvMappedField<Vector, on_cell>& eta_gradient, cuFvMappedField<Tensor, on_cell>& u_gradient, cuFvMappedField<Scalar, on_cell>& h_advection, cuFvMappedField<Vector, on_cell>& hU_advection){

      auto mesh = h.mesh;

      cuAdvectionNSWEs2ndKernel<<<BLOCKS_PER_GRID, THREADS_PER_BLOCK>>>(gravity.data.dev_ptr(), gravity.boundary_value.dev_ptr(), h.data.dev_ptr(), h.boundary_value.dev_ptr(), z.data.dev_ptr(), z.boundary_value.dev_ptr(), u.data.dev_ptr(), u.boundary_value.dev_ptr(), h_gradient.data.dev_ptr(), eta_gradient.data.dev_ptr(), u_gradient.data.dev_ptr(), h.data.size(), mesh->cell_neighbours.dims_dev_ptr(), mesh->cell_neighbours.dev_ptr(), mesh->cell_neighbours.length(), mesh->cell_halffacets.dev_ptr(), mesh->cell_halffacets.length(), mesh->cell_centre_positions.dev_ptr(), mesh->cell_volumes.dev_ptr(), mesh->halffacet_areas.dev_ptr(), mesh->halffacet_normal_directions.dev_ptr(), mesh->halffacet_centre_positions.dev_ptr(), h_advection.data.dev_ptr(), hU_advection.data.dev_ptr());

    }

    __global__ void cuAdvectionNSWEs2ndCurvKernel(Scalar* gravity, Scalar* gravity_bound, Scalar* centrifugal, Scalar* h, Scalar* h_bound, Scalar* z, Scalar* z_bound, Vector* u, Vector* u_bound, Vector* h_gradient, Vector* eta_gradient, Tensor* u_gradient, unsigned int phi_size, Flag* cell_neigbours_dimensions, ShortDualHandle* cell_neigbours, unsigned int cell_neighbours_length, Flag* cell_halffacets, unsigned int cell_halffacets_length, Vector* cell_centres, Scalar* cell_volume, Scalar* face_area, Vector* face_normal, Vector* face_centres, Scalar* h_advection, Vector* hU_advection){

      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      while (index < phi_size){
        Scalar g_this = gravity[index];
        Scalar centrifugal_this = centrifugal[index];
        Scalar h_this = h[index];
        Scalar z_this = z[index];
        Scalar eta_this = h_this + z_this;
        Vector2 u_this = u[index];
        Vector2 grad_h_this = h_gradient[index];
        Vector2 grad_eta_this = eta_gradient[index];
        Tensor2 grad_u_this = u_gradient[index];
        Scalar volume_this = cell_volume[index];
        Vector2 cell_centroid_this = cell_centres[index];
        Flag cell_neigbours_number = cell_neigbours_dimensions[index];
        Scalar _h_advection(0.0);
        Vector2 _hU_advection(0.0, 0.0);
        for (Flag i = 0; i < cell_neigbours_number; ++i){
          Flag id_face = cell_halffacets[i*cell_halffacets_length + index];
          Vector2 normal = uni(face_normal[id_face]);
          Vector2 shear = uni(perpend(normal));
          Scalar area = face_area[id_face];
          Vector2 face_centroid = face_centres[id_face];
          Vector2 direction_this = face_centroid - cell_centroid_this;
          Scalar _h_this = h_this + dot(grad_h_this, direction_this);
          Scalar _eta_this = eta_this + dot(grad_eta_this, direction_this);
          Scalar _z_this = _eta_this - _h_this;
          Vector2 _u_this = u_this + dot(grad_u_this, direction_this);
          ShortDualHandle neib = cell_neigbours[i*cell_neighbours_length + index];
          if (!neib.is_boundary()){
            Flag id_neib = neib.get_global_id();
            Scalar g_neib = gravity[id_neib];
            Scalar h_neib = h[id_neib];
            Scalar z_neib = z[id_neib];
            Scalar eta_neib = h_neib + z_neib;
            Vector2 u_neib = u[id_neib];
            Vector2 grad_h_neib = h_gradient[id_neib];
            Vector2 grad_eta_neib = eta_gradient[id_neib];
            Tensor2 grad_u_neib = u_gradient[id_neib];
            Vector2 cell_centroid_neib = cell_centres[id_neib];
            Vector2 direction_neib = face_centroid - cell_centroid_neib;
            Scalar _h_neib = h_neib + dot(grad_h_neib, direction_neib);
            Scalar _eta_neib = eta_neib + dot(grad_eta_neib, direction_neib);
            Scalar _z_neib = _eta_neib - _h_neib;
            Vector2 _u_neib = u_neib + dot(grad_u_neib, direction_neib);
            Scalar grav = (g_this + g_neib) / 2.0;
            Scalar z_face = fmax(_z_this, _z_neib);
            Scalar h_L = fmax(_h_this + _z_this - z_face, 0.0);
            Scalar h_R = fmax(_h_neib + _z_neib - z_face, 0.0);
            Vector2 u_L(dot(_u_this, normal), dot(_u_this, shear));
            Vector2 u_R(dot(_u_neib, normal), dot(_u_neib, shear));
            auto flux = cuHLLCRiemannSolverSWEs(grav, ScalarRiemannState(h_L, h_R), VectorRiemannState(h_L*u_L, h_R*u_R));
            Vector2 _hU_flux = flux.q.x*normal + flux.q.y*shear;
            Vector2 _z_flux = -0.5*(grav + centrifugal_this)*h_L*h_L*normal + 0.5*(g_this + centrifugal_this)*(_eta_this - z_this)*(_eta_this - z_this)*normal;
            _h_advection += flux.h*area / volume_this;
            _hU_advection += (_hU_flux + _z_flux)*area / volume_this;
          }
          else{
            Flag id_boundary = neib.get_global_id();
            Scalar g_bound = gravity_bound[id_boundary];
            Scalar _h_bound = h_bound[id_boundary];
            Scalar _z_bound = z_bound[id_boundary];
            Vector2 _u_bound = u_bound[id_boundary];
            Scalar grav = (g_this + g_bound) / 2.0;
            Scalar z_face = fmax(_z_this, _z_bound);
            Scalar h_L = fmax(_h_this + _z_this - z_face, 0.0);
            Scalar h_R = fmax(_h_bound + _z_bound - z_face, 0.0);
            Vector2 u_L(dot(_u_this, normal), dot(_u_this, shear));
            Vector2 u_R(dot(_u_bound, normal), dot(_u_bound, shear));
            auto flux = cuHLLCRiemannSolverSWEs(grav, ScalarRiemannState(h_L, h_R), VectorRiemannState(h_L*u_L, h_R*u_R));
            Vector2 _hU_flux = flux.q.x*normal + flux.q.y*shear;
            Vector2 _z_flux = -0.5*(grav + centrifugal_this)*h_L*h_L*normal + 0.5*(g_this + centrifugal_this)*(_eta_this - z_this)*(_eta_this - z_this)*normal;
            _h_advection += flux.h*area / volume_this;
            _hU_advection += (_hU_flux + _z_flux)*area / volume_this;
          }
        }
        h_advection[index] = _h_advection;
        hU_advection[index] = _hU_advection;
        index += blockDim.x * gridDim.x;
      }
    }

    void cuAdvectionNSWEs2ndCurv(cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& centrifugal, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Scalar, on_cell>& z, cuFvMappedField<Vector, on_cell>& u, cuFvMappedField<Vector, on_cell>& h_gradient, cuFvMappedField<Vector, on_cell>& eta_gradient, cuFvMappedField<Tensor, on_cell>& u_gradient, cuFvMappedField<Scalar, on_cell>& h_advection, cuFvMappedField<Vector, on_cell>& hU_advection){

      auto mesh = h.mesh;

      cuAdvectionNSWEs2ndCurvKernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(gravity.data.dev_ptr(), gravity.boundary_value.dev_ptr(), centrifugal.data.dev_ptr(), h.data.dev_ptr(), h.boundary_value.dev_ptr(), z.data.dev_ptr(), z.boundary_value.dev_ptr(), u.data.dev_ptr(), u.boundary_value.dev_ptr(), h_gradient.data.dev_ptr(), eta_gradient.data.dev_ptr(), u_gradient.data.dev_ptr(), h.data.size(), mesh->cell_neighbours.dims_dev_ptr(), mesh->cell_neighbours.dev_ptr(), mesh->cell_neighbours.length(), mesh->cell_halffacets.dev_ptr(), mesh->cell_halffacets.length(), mesh->cell_centre_positions.dev_ptr(), mesh->cell_volumes.dev_ptr(), mesh->halffacet_areas.dev_ptr(), mesh->halffacet_normal_directions.dev_ptr(), mesh->halffacet_centre_positions.dev_ptr(), h_advection.data.dev_ptr(), hU_advection.data.dev_ptr());


    }

    __global__ void cuAdvectionNSWEsOrderReducerCartesianKernel(Scalar* h, Scalar* _h_bound, Scalar* z, Scalar* _z_bound, Vector* h_gradient, Vector* eta_gradient, Tensor2* u_grad,  unsigned int phi_size, Flag* cell_neigbours_dimensions, ShortDualHandle* cell_neigbours, unsigned int cell_neighbours_length, Flag* cell_halffacets, unsigned int cell_halffacets_length, Vector* cell_centres, Vector* face_normal, Vector* face_centres){

      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      while (index < phi_size){
        Scalar h_this = h[index];
        Scalar z_this = z[index];
        Scalar eta_this = h_this + z_this;
        Vector2 grad_h_this = h_gradient[index];
        Vector2 grad_eta_this = eta_gradient[index];
        Vector2 cell_centroid_this = cell_centres[index];
        Flag cell_neigbours_number = cell_neigbours_dimensions[index];
        Scalar H[4];
        Scalar ratio_x = 1.0;
        Scalar ratio_y = 1.0;
        for (Flag i = 0; i < cell_neigbours_number; ++i){
          Flag id_face = cell_halffacets[i*cell_halffacets_length + index];
          Vector2 normal = uni(face_normal[id_face]);
          Vector2 face_centroid = face_centres[id_face];
          Vector2 direction_this = face_centroid - cell_centroid_this;         
          Scalar _h_this = h_this + dot(grad_h_this, direction_this);
          Scalar _eta_this = eta_this + dot(grad_eta_this, direction_this);
          Scalar _z_this = _eta_this - _h_this;
          ShortDualHandle neib = cell_neigbours[i*cell_neighbours_length + index];
          if (!neib.is_boundary()){
            Flag id_neib = neib.get_global_id();
            Scalar h_neib = h[id_neib];
            Scalar z_neib = z[id_neib];
            Scalar eta_neib = h_neib + z_neib;
            Vector2 grad_h_neib = h_gradient[id_neib];
            Vector2 grad_eta_neib = eta_gradient[id_neib];
            Vector2 cell_centroid_neib = cell_centres[id_neib];
            Vector2 direction_neib = face_centroid - cell_centroid_neib;
            Scalar _h_neib = h_neib + dot(grad_h_neib, direction_neib);
            Scalar _eta_neib = eta_neib + dot(grad_eta_neib, direction_neib);
            Scalar _z_neib = _eta_neib - _h_neib;
            Scalar z_face = fmax(_z_this, _z_neib);
            Scalar h_L = fmax(_h_this + _z_this - z_face, 0.0);
            H[i] = h_L;
          }
          else{
            Flag id_boundary = neib.get_global_id();
            Scalar h_bound = _h_bound[id_boundary];
            Scalar z_bound = _z_bound[id_boundary];
            Scalar z_face = fmax(_z_this, z_bound);
            Scalar h_L = fmax(_h_this + _z_this - z_face, 0.0);
            Scalar h_R = fmax(h_bound + z_bound - z_face, 0.0);
            H[i] = h_L;
          }
        }
        Scalar sum_h_x = H[1] + H[3];
        Scalar sum_h_y = H[0] + H[2];
        if(fabs(h_this) >= 1e-10){
          ratio_x = sum_h_x/h_this;
          ratio_y = sum_h_y/h_this;
        } 
        Scalar grad_eta = dot(grad_eta_this, grad_eta_this);
        if ((ratio_x <= 0.1 || ratio_y <= 0.1 ) && grad_eta >= 1e-10){
          h_gradient[index] = 0.0;
          eta_gradient[index] = 0.0;
          u_grad[index] = 0.0;
        }
        index += blockDim.x * gridDim.x;
      }
    }

    void cuAdvectionNSWEs2ndRobust(cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Scalar, on_cell>& z, cuFvMappedField<Vector, on_cell>& u, cuFvMappedField<Vector, on_cell>& h_gradient, cuFvMappedField<Vector, on_cell>& eta_gradient, cuFvMappedField<Tensor, on_cell>& u_gradient, cuFvMappedField<Scalar, on_cell>& h_advection, cuFvMappedField<Vector, on_cell>& hU_advection){

      auto mesh = h.mesh;

      cuAdvectionNSWEsOrderReducerCartesianKernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(h.data.dev_ptr(), h.boundary_value.dev_ptr(), z.data.dev_ptr(), z.boundary_value.dev_ptr(), h_gradient.data.dev_ptr(), eta_gradient.data.dev_ptr(), u_gradient.data.dev_ptr(), h.data.size(), mesh->cell_neighbours.dims_dev_ptr(), mesh->cell_neighbours.dev_ptr(), mesh->cell_neighbours.length(), mesh->cell_halffacets.dev_ptr(), mesh->cell_halffacets.length(), mesh->cell_centre_positions.dev_ptr(), mesh->halffacet_normal_directions.dev_ptr(), mesh->halffacet_centre_positions.dev_ptr());

      cuAdvectionNSWEs2ndKernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(gravity.data.dev_ptr(), gravity.boundary_value.dev_ptr(), h.data.dev_ptr(), h.boundary_value.dev_ptr(), z.data.dev_ptr(), z.boundary_value.dev_ptr(), u.data.dev_ptr(), u.boundary_value.dev_ptr(), h_gradient.data.dev_ptr(), eta_gradient.data.dev_ptr(), u_gradient.data.dev_ptr(), h.data.size(), mesh->cell_neighbours.dims_dev_ptr(), mesh->cell_neighbours.dev_ptr(), mesh->cell_neighbours.length(), mesh->cell_halffacets.dev_ptr(), mesh->cell_halffacets.length(), mesh->cell_centre_positions.dev_ptr(), mesh->cell_volumes.dev_ptr(), mesh->halffacet_areas.dev_ptr(), mesh->halffacet_normal_directions.dev_ptr(), mesh->halffacet_centre_positions.dev_ptr(), h_advection.data.dev_ptr(), hU_advection.data.dev_ptr());

    }


    void cuAdvectionNSWEs2ndRobustCurv(cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& centrifugal, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Scalar, on_cell>& z, cuFvMappedField<Vector, on_cell>& u, cuFvMappedField<Vector, on_cell>& h_gradient, cuFvMappedField<Vector, on_cell>& eta_gradient, cuFvMappedField<Tensor, on_cell>& u_gradient, cuFvMappedField<Scalar, on_cell>& h_advection, cuFvMappedField<Vector, on_cell>& hU_advection){

      auto mesh = h.mesh;

      cuAdvectionNSWEsOrderReducerCartesianKernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(h.data.dev_ptr(), h.boundary_value.dev_ptr(), z.data.dev_ptr(), z.boundary_value.dev_ptr(), h_gradient.data.dev_ptr(), eta_gradient.data.dev_ptr(), u_gradient.data.dev_ptr(), h.data.size(), mesh->cell_neighbours.dims_dev_ptr(), mesh->cell_neighbours.dev_ptr(), mesh->cell_neighbours.length(), mesh->cell_halffacets.dev_ptr(), mesh->cell_halffacets.length(), mesh->cell_centre_positions.dev_ptr(), mesh->halffacet_normal_directions.dev_ptr(), mesh->halffacet_centre_positions.dev_ptr());

      cuAdvectionNSWEs2ndCurvKernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(gravity.data.dev_ptr(), gravity.boundary_value.dev_ptr(), centrifugal.data.dev_ptr(), h.data.dev_ptr(), h.boundary_value.dev_ptr(), z.data.dev_ptr(), z.boundary_value.dev_ptr(), u.data.dev_ptr(), u.boundary_value.dev_ptr(), h_gradient.data.dev_ptr(), eta_gradient.data.dev_ptr(), u_gradient.data.dev_ptr(), h.data.size(), mesh->cell_neighbours.dims_dev_ptr(), mesh->cell_neighbours.dev_ptr(), mesh->cell_neighbours.length(), mesh->cell_halffacets.dev_ptr(), mesh->cell_halffacets.length(), mesh->cell_centre_positions.dev_ptr(), mesh->cell_volumes.dev_ptr(), mesh->halffacet_areas.dev_ptr(), mesh->halffacet_normal_directions.dev_ptr(), mesh->halffacet_centre_positions.dev_ptr(), h_advection.data.dev_ptr(), hU_advection.data.dev_ptr());

    }


    __global__ void cuAdvectionNSWEs2ndFastKernel1(Scalar* gravity, Scalar* _gravity_bound, Scalar* h, Scalar* _h_bound, Scalar* z, Scalar* _z_bound, Vector* u, Vector* _u_bound, Vector* h_gradient, Vector* eta_gradient, Tensor* u_gradient, unsigned int phi_size, Flag* cell_neigbours_dimensions, ShortDualHandle* cell_neigbours, unsigned int cell_neighbours_length, Flag* cell_halffacets, unsigned int cell_halffacets_length, Vector* cell_centres, Scalar* cell_volume, Scalar* face_area, Vector* face_normal, Vector* face_centres, Scalar* h_flux, Vector* hU_flux){

      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      while (index < phi_size){
        Scalar g_this = gravity[index];
        Scalar h_this = h[index];
        Scalar z_this = z[index];
        Scalar eta_this = h_this + z_this;
        Vector2 u_this = u[index];
        Vector2 grad_h_this = h_gradient[index];
        Vector2 grad_eta_this = eta_gradient[index];
        Tensor2 grad_u_this = u_gradient[index];
        Scalar volume_this = cell_volume[index];
        Vector2 cell_centroid_this = cell_centres[index];
        Flag cell_neigbours_number = cell_neigbours_dimensions[index];
        for (Flag i = 0; i < cell_neigbours_number; ++i){
          Flag id_face = cell_halffacets[i*cell_halffacets_length + index];
          Vector2 normal = uni(face_normal[id_face]);
          Vector2 shear = uni(perpend(normal));
          Scalar area = face_area[id_face];
          Vector2 face_centroid = face_centres[id_face];
          Vector2 direction_this = face_centroid - cell_centroid_this;
          Scalar _h_this = h_this + dot(grad_h_this, direction_this);
          Scalar _eta_this = eta_this + dot(grad_eta_this, direction_this);
          Scalar _z_this = _eta_this - _h_this;
          Vector2 _u_this = u_this + dot(grad_u_this, direction_this);
          ShortDualHandle neib = cell_neigbours[i*cell_neighbours_length + index];
          if (!neib.is_boundary()){
            Flag id_neib = neib.get_global_id();
            Flag neib_faceid = neib.get_local_id();
            Flag id_neib_face = cell_halffacets[neib_faceid*cell_halffacets_length + id_neib];
            if (index < id_neib){
              Scalar g_neib = gravity[id_neib];
              Scalar h_neib = h[id_neib];
              Scalar z_neib = z[id_neib];
              Scalar eta_neib = h_neib + z_neib;
              Vector2 u_neib = u[id_neib];
              Vector2 grad_h_neib = h_gradient[id_neib];
              Vector2 grad_eta_neib = eta_gradient[id_neib];
              Tensor2 grad_u_neib = u_gradient[id_neib];
              Vector2 cell_centroid_neib = cell_centres[id_neib];
              Scalar volume_neib = cell_volume[id_neib];
              Vector2 direction_neib = face_centroid - cell_centroid_neib;
              Scalar _h_neib = h_neib + dot(grad_h_neib, direction_neib);
              Scalar _eta_neib = eta_neib + dot(grad_eta_neib, direction_neib);
              Scalar _z_neib = _eta_neib - _h_neib;
              Vector2 _u_neib = u_neib + dot(grad_u_neib, direction_neib);
              Scalar grav = (g_this + g_neib) / 2.0;
              Scalar z_face = fmax(_z_this, _z_neib);
              Scalar h_L = fmax(_h_this + _z_this - z_face, 0.0);
              Scalar h_R = fmax(_h_neib + _z_neib - z_face, 0.0);
              Vector2 u_L(dot(_u_this, normal), dot(_u_this, shear));
              Vector2 u_R(dot(_u_neib, normal), dot(_u_neib, shear));
              auto flux = cuHLLCRiemannSolverSWEs(grav, ScalarRiemannState(h_L, h_R), VectorRiemannState(h_L*u_L, h_R*u_R));
              Vector2 _hU_flux = flux.q.x*normal + flux.q.y*shear;
              Vector2 _z_flux_this = -0.5*grav*h_L*h_L*normal + 0.5*g_this*(_eta_this - z_this)*(_eta_this - z_this)*normal;
              h_flux[id_face] = flux.h*area / volume_this;
              hU_flux[id_face] = (_hU_flux + _z_flux_this)*area / volume_this;
              Vector2 _z_flux_neib = -0.5*grav*h_R*h_R*normal + 0.5*g_neib*(_eta_neib - z_neib)*(_eta_neib - z_neib)*normal;
              h_flux[id_neib_face] = -1.0*flux.h*area / volume_neib;
              hU_flux[id_neib_face] = -1.0*(_hU_flux + _z_flux_neib)*area / volume_neib;
            }
          }
          else{
            Flag id_boundary = neib.get_global_id();
            Scalar g_bound = _gravity_bound[id_boundary];
            Scalar h_bound = _h_bound[id_boundary];
            Scalar z_bound = _z_bound[id_boundary];
            Scalar grav = (g_this + g_bound) / 2.0;
            Scalar z_face = fmax(_z_this, z_bound);
            Scalar h_L = fmax(_h_this + _z_this - z_face, 0.0);
            Scalar h_R = fmax(h_bound + z_bound - z_face, 0.0);
            Vector2 u_bound = _u_bound[id_boundary];
            Vector2 u_L(dot(_u_this, normal), dot(_u_this, shear));
            Vector2 u_R(dot(u_bound, normal), dot(u_bound, shear));
            auto flux = cuHLLCRiemannSolverSWEs(grav, ScalarRiemannState(h_L, h_R), VectorRiemannState(h_L*u_L, h_R*u_R));
            Vector2 _hU_flux = flux.q.x*normal + flux.q.y*shear;
            Vector2 _z_flux = -0.5*grav*h_L*h_L*normal + 0.5*g_this*(_eta_this - z_this)*(_eta_this - z_this)*normal;
            h_flux[id_face] = flux.h*area / volume_this;
            hU_flux[id_face] = (_hU_flux + _z_flux)*area / volume_this;
          }
        }
        __syncthreads();
        index += blockDim.x * gridDim.x;
      }
    }

    __global__ void cuAdvectionNSWEs2ndFastKernel2(Scalar* h_flux, Vector* hU_flux, Flag* cell_neigbours_dimensions, ShortDualHandle* cell_neigbours, unsigned int cell_neighbours_length, Flag* cell_halffacets, unsigned int cell_halffacets_length, Scalar* h_advection, Vector* hU_advection, unsigned int cell_list_size){
      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      while (index < cell_list_size){
        Flag cell_neigbours_number = cell_neigbours_dimensions[index];
        Scalar _h_advection(0.0);
        Vector _hU_advection(0.0);
        for (Flag i = 0; i < cell_neigbours_number; ++i){
          Flag id_face = cell_halffacets[i*cell_halffacets_length + index];
          _h_advection += h_flux[id_face];
          _hU_advection += hU_flux[id_face];
        }
        h_advection[index] = _h_advection;
        hU_advection[index] = _hU_advection;
        index += blockDim.x * gridDim.x;
      }
    }

    void cuAdvectionNSWEs2ndFast(cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Scalar, on_cell>& z, cuFvMappedField<Vector, on_cell>& u, cuFvMappedField<Vector, on_cell>& h_gradient, cuFvMappedField<Vector, on_cell>& eta_gradient, cuFvMappedField<Tensor, on_cell>& u_gradient, cuFvMappedField<Scalar, on_halffacet>& h_flux, cuFvMappedField<Vector, on_halffacet>& hU_flux, cuFvMappedField<Scalar, on_cell>& h_advection, cuFvMappedField<Vector, on_cell>& hU_advection){

      auto mesh = h.mesh;

      cuAdvectionNSWEs2ndFastKernel1 << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(gravity.data.dev_ptr(), gravity.boundary_value.dev_ptr(), h.data.dev_ptr(), h.boundary_value.dev_ptr(), z.data.dev_ptr(), z.boundary_value.dev_ptr(), u.data.dev_ptr(), u.boundary_value.dev_ptr(), h_gradient.data.dev_ptr(), eta_gradient.data.dev_ptr(), u_gradient.data.dev_ptr(), h.data.size(), mesh->cell_neighbours.dims_dev_ptr(), mesh->cell_neighbours.dev_ptr(), mesh->cell_neighbours.length(), mesh->cell_halffacets.dev_ptr(), mesh->cell_halffacets.length(), mesh->cell_centre_positions.dev_ptr(), mesh->cell_volumes.dev_ptr(), mesh->halffacet_areas.dev_ptr(), mesh->halffacet_normal_directions.dev_ptr(), mesh->halffacet_centre_positions.dev_ptr(), h_flux.data.dev_ptr(), hU_flux.data.dev_ptr());
      cuAdvectionNSWEs2ndFastKernel2 << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(h_flux.data.dev_ptr(), hU_flux.data.dev_ptr(), mesh->cell_neighbours.dims_dev_ptr(), mesh->cell_neighbours.dev_ptr(), mesh->cell_neighbours.length(), mesh->cell_halffacets.dev_ptr(), mesh->cell_halffacets.length(), h_advection.data.dev_ptr(), hU_advection.data.dev_ptr(), h_advection.data.size());

    }

    __global__ void cuAdvectionMSWEsKernel(Scalar* gravity, Scalar* h, Scalar* _h_bound, Scalar* z, Scalar* _z_bound, Vector* hU, Vector* _hU_bound, unsigned int phi_size, Flag* cell_neigbours_dimensions, ShortDualHandle* cell_neigbours, unsigned int cell_neighbours_length, Flag* cell_halffacets, unsigned int cell_halffacets_length, Scalar* cell_volume, Scalar* face_area, Vector* face_normal, Scalar* h_advection, Vector* hU_advection, Scalar* dt_mass){

      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      Scalar h_small = 1e-10;
      while (index < phi_size){
        Scalar g = gravity[index];
        Scalar h_this = h[index];
        Scalar z_this = z[index];
        Scalar eta_this = h_this + z_this;
        Scalar X_this = h_this / eta_this;
        Vector2 hU_this = hU[index];
        Vector2 u_this = 0.0;
        if (h_this < h_small){
          u_this = 0.0;
        }
        else{
          u_this = hU_this / h_this;
        }
        Scalar volume_this = cell_volume[index];
        Flag cell_neigbours_number = cell_neigbours_dimensions[index];
        Scalar _h_advection(0.0);
        Scalar _h_advection_constraint(0.0);
        Vector2 _hU_advection(0.0, 0.0);
        for (Flag i = 0; i < cell_neigbours_number; ++i){
          Flag id_face = cell_halffacets[i*cell_halffacets_length + index];
          Vector2 normal = uni(face_normal[id_face]);
          Vector2 shear = uni(perpend(normal));
          Scalar area = face_area[id_face];
          ShortDualHandle neib = cell_neigbours[i*cell_neighbours_length + index];
          Scalar h_neib = 0.0;
          Scalar z_neib = 0.0;
          Vector2 hU_neib = 0.0;
          if (!neib.is_boundary()){
            Flag id_neib = neib.get_global_id();
            h_neib = h[id_neib];
            z_neib = z[id_neib];
            hU_neib = hU[id_neib];
          }
          else{
            Flag id_boundary = neib.get_global_id();
            h_neib = _h_bound[id_boundary];
            z_neib = _z_bound[id_boundary];
            hU_neib = _hU_bound[id_boundary];
          }
          Scalar eta_neib = h_neib + z_neib;
          Scalar X_neib = h_neib / eta_neib;
          Vector2 u_neib = 0.0;
          if (h_neib < h_small){
            u_neib = 0.0;
          }
          else{
            u_neib = hU_neib / h_neib;
          }
          Vector2 u_L(dot(u_this, normal), dot(u_this, shear));
          Vector2 u_R(dot(u_neib, normal), dot(u_neib, shear));
          auto flux = cuHLLCRiemannSolverMSWEs(g,0.0, ScalarRiemannState(h_this, h_neib), ScalarRiemannState(eta_this, eta_neib), VectorRiemannState(u_L, u_R));
          Scalar eta_f, X_f;
          if (flux.h > 0){
            eta_f = eta_this;
            X_f = X_this;
          }
          else{
            eta_f = eta_neib;
            X_f = X_neib;
          }
          Scalar _h_flux = flux.h*X_f;
          Vector2 _hU_flux = (flux.q.x*normal + flux.q.y*shear)*X_f;
          Vector2 _z_flux = -0.5*g*eta_this*eta_f*(X_f - X_this)*normal;
          _h_advection += _h_flux*area / volume_this;
          _h_advection_constraint += fmax(0.0, _h_flux)*area / volume_this;
          _hU_advection += (_hU_flux + _z_flux)*area / volume_this;
        }
        h_advection[index] = _h_advection;
        hU_advection[index] = _hU_advection;
        if (_h_advection_constraint <= 1e-10){
          dt_mass[index] = 1e35;
        }
        else{
          dt_mass[index] = eta_this / _h_advection_constraint;
        }
        __syncthreads();
        index += blockDim.x * gridDim.x;
      }
    }

    void cuAdvectionMSWEs(cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Scalar, on_cell>& z, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Scalar, on_cell>& h_advection, cuFvMappedField<Vector, on_cell>& hU_advection, cuFvMappedField<Scalar, on_cell>& dt_mass){

      auto mesh = h.mesh;

      cuAdvectionMSWEsKernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(gravity.data.dev_ptr(),
        h.data.dev_ptr(),
        h.boundary_value.dev_ptr(),
        z.data.dev_ptr(),
        z.boundary_value.dev_ptr(),
        hU.data.dev_ptr(),
        hU.boundary_value.dev_ptr(),
        h.data.size(),
        mesh->cell_neighbours.dims_dev_ptr(),
        mesh->cell_neighbours.dev_ptr(),
        mesh->cell_neighbours.length(),
        mesh->cell_halffacets.dev_ptr(),
        mesh->cell_halffacets.length(),
        mesh->cell_volumes.dev_ptr(),
        mesh->halffacet_areas.dev_ptr(),
        mesh->halffacet_normal_directions.dev_ptr(),
        h_advection.data.dev_ptr(),
        hU_advection.data.dev_ptr(),
        dt_mass.data.dev_ptr());
    }

    __global__ void cuAdvectionSWEsCartesianKernel(Scalar* gravity, Scalar* h, Scalar* _h_bound, Scalar* z, Scalar* _z_bound, Vector* hU, Vector* _hU_bound, unsigned int phi_size, ShortDualHandle* cell_neigbours, unsigned int cell_neighbours_length, Scalar* cell_volume, Scalar* h_advection, Vector* hU_advection){

      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      Scalar h_small = 1e-10;
      Vector2 face_normal[4];
      face_normal[0] = Vector2(0.0, -1.0);
      face_normal[1] = Vector2(1.0, 0.0);
      face_normal[2] = Vector2(0.0, 1.0);
      face_normal[3] = Vector2(-1.0, 0.0);
      while (index < phi_size){
        Scalar g = gravity[index];
        Scalar h_this = h[index];
        Scalar z_this = z[index];
        Vector2 hU_this = hU[index];
        Vector2 u_this = 0.0;
        if (h_this < h_small){
          u_this = 0.0;
        }
        else{
          u_this = hU_this / h_this;
        }
        Scalar volume = cell_volume[index];
        Scalar area = sqrt(volume);
        Scalar _h_advection(0.0);
        Scalar _h_advection_constraint(0.0);
        Vector2 _hU_advection(0.0, 0.0);
        for (Flag i = 0; i < 4; ++i){
          ShortDualHandle neib = cell_neigbours[i*cell_neighbours_length + index];
          Scalar h_neib = 0.0;
          Scalar z_neib = 0.0;
          Vector2 hU_neib = 0.0;
          if (!neib.is_boundary()){
            Flag id_neib = neib.get_global_id();
            h_neib = h[id_neib];
            z_neib = z[id_neib];
            hU_neib = hU[id_neib];
          }
          else{
            Flag id_boundary = neib.get_global_id();
            h_neib = _h_bound[id_boundary];
            z_neib = _z_bound[id_boundary];
            hU_neib = _hU_bound[id_boundary];
          }
          if (h_this < h_small && h_neib < h_small){
            continue;
          }
          Vector2 normal = face_normal[i];
          Vector2 shear = uni(perpend(normal));
          Vector2 u_neib = 0.0;
          if (h_neib < h_small){
            u_neib = 0.0;
          }
          else{
            u_neib = hU_neib / h_neib;
          }
          Vector2 u_L(dot(u_this, normal), dot(u_this, shear));
          Vector2 u_R(dot(u_neib, normal), dot(u_neib, shear));
          Scalar z_f = fmax(z_this, z_neib);
          Scalar h_L = fmax(0.0, h_this + z_this - z_f);
          Scalar h_R = fmax(0.0, h_neib + z_neib - z_f);
          auto flux = cuHLLCRiemannSolverSWEs(g, ScalarRiemannState(h_L, h_R), VectorRiemannState(h_L*u_L, h_R*u_R));
          Vector2 _hU_flux = flux.q.x*normal + flux.q.y*shear;
          Vector2 _z_flux = -0.5*g*h_L*h_L*normal;
          _h_advection += flux.h*area / volume;
          _hU_advection += (_hU_flux + _z_flux)*area / volume;
        }
        h_advection[index] = _h_advection;
        hU_advection[index] = _hU_advection;
        __syncthreads();
        index += blockDim.x * gridDim.x;
      }
    }

    void cuAdvectionSWEsCartesian(cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Scalar, on_cell>& z, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Scalar, on_cell>& h_advection, cuFvMappedField<Vector, on_cell>& hU_advection){

      auto mesh = h.mesh;

      cuAdvectionSWEsCartesianKernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(gravity.data.dev_ptr(),
        h.data.dev_ptr(),
        h.boundary_value.dev_ptr(),
        z.data.dev_ptr(),
        z.boundary_value.dev_ptr(),
        hU.data.dev_ptr(),
        hU.boundary_value.dev_ptr(),
        h.data.size(),
        mesh->cell_neighbours.dev_ptr(),
        mesh->cell_neighbours.length(),
        mesh->cell_volumes.dev_ptr(),
        h_advection.data.dev_ptr(),
        hU_advection.data.dev_ptr());
    }


    //Record flux
    __global__ void cuAdvectionSWEsRFCartesianKernel(Scalar* gravity, Scalar* h, Scalar* _h_bound, Scalar* z, Scalar* _z_bound, Vector* hU, Vector* _hU_bound, unsigned int phi_size, ShortDualHandle* cell_neigbours, unsigned int cell_neighbours_length, Scalar* cell_volume, Scalar* h_advection, Vector* hU_advection, Vector* hU_EW, Vector* hU_NS){

      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      Scalar h_small = 1e-10;
      Vector2 face_normal[4];
      face_normal[0] = Vector2(0.0, -1.0);
      face_normal[1] = Vector2(1.0, 0.0);
      face_normal[2] = Vector2(0.0, 1.0);
      face_normal[3] = Vector2(-1.0, 0.0);
      while (index < phi_size){
        Scalar g = gravity[index];
        Scalar h_this = h[index];
        Scalar z_this = z[index];
        Vector2 hU_this = hU[index];
        Vector2 u_this = 0.0;
        if (h_this < h_small){
          u_this = 0.0;
        }
        else{
          u_this = hU_this / h_this;
        }
        Scalar volume = cell_volume[index];
        Scalar area = sqrt(volume);
        Scalar _h_advection(0.0);
        Scalar _h_advection_constraint(0.0);
        Vector2 _hU_advection(0.0, 0.0);
        for (Flag i = 0; i < 4; ++i){
          ShortDualHandle neib = cell_neigbours[i*cell_neighbours_length + index];
          Scalar h_neib = 0.0;
          Scalar z_neib = 0.0;
          Vector2 hU_neib = 0.0;
          if (!neib.is_boundary()){
            Flag id_neib = neib.get_global_id();
            h_neib = h[id_neib];
            z_neib = z[id_neib];
            hU_neib = hU[id_neib];
          }
          else{
            Flag id_boundary = neib.get_global_id();
            h_neib = _h_bound[id_boundary];
            z_neib = _z_bound[id_boundary];
            //z_neib = z_this;
            hU_neib = _hU_bound[id_boundary];
          }
          if (h_this < h_small && h_neib < h_small){
            continue;
          }
          Vector2 normal = face_normal[i];
          Vector2 shear = uni(perpend(normal));
          Vector2 u_neib = 0.0;
          if (h_neib < h_small){
            u_neib = 0.0;
          }
          else{
            u_neib = hU_neib / h_neib;
          }
          Vector2 u_L(dot(u_this, normal), dot(u_this, shear));
          Vector2 u_R(dot(u_neib, normal), dot(u_neib, shear));
          Scalar z_f = fmax(z_this, z_neib);
          Scalar h_L = fmax(0.0, h_this + z_this - z_f);
          Scalar h_R = fmax(0.0, h_neib + z_neib - z_f);
          auto flux = cuHLLCRiemannSolverSWEs(g, ScalarRiemannState(h_L, h_R), VectorRiemannState(h_L*u_L, h_R*u_R));
          Vector2 _hU_flux = flux.q.x*normal + flux.q.y*shear;
          Vector2 _z_flux = -0.5*g*h_L*h_L*normal;
          _h_advection += flux.h*area / volume;
          if(i == 0){
            hU_NS[index].y = flux.h;
          }
          if(i == 2){
            hU_NS[index].x = flux.h;
          }
          if(i == 1){
            hU_EW[index].x = flux.h;
          }
          if(i == 3){
            hU_EW[index].y = flux.h;
          }
          _hU_advection += (_hU_flux + _z_flux)*area / volume;
        }
        h_advection[index] = _h_advection;
        hU_advection[index] = _hU_advection;
        __syncthreads();
        index += blockDim.x * gridDim.x;
      }
    }

    void cuAdvectionSWEsRFCartesian(cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Scalar, on_cell>& z, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Scalar, on_cell>& h_advection, cuFvMappedField<Vector, on_cell>& hU_advection, cuFvMappedField<Vector, on_cell>& hU_EW, cuFvMappedField<Vector, on_cell>& hU_NS){

      auto mesh = h.mesh;

      cuAdvectionSWEsRFCartesianKernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(gravity.data.dev_ptr(),
        h.data.dev_ptr(),
        h.boundary_value.dev_ptr(),
        z.data.dev_ptr(),
        z.boundary_value.dev_ptr(),
        hU.data.dev_ptr(),
        hU.boundary_value.dev_ptr(),
        h.data.size(),
        mesh->cell_neighbours.dev_ptr(),
        mesh->cell_neighbours.length(),
        mesh->cell_volumes.dev_ptr(),
        h_advection.data.dev_ptr(),
        hU_advection.data.dev_ptr(),
        hU_EW.data.dev_ptr(),
        hU_NS.data.dev_ptr());
    }



    __global__ void cuAdvectionSWEsAPCartesianKernel(Scalar* manning, Scalar* gravity, Scalar* h, Scalar* _h_bound, Scalar* z, Scalar* _z_bound, Vector* hU, Vector* _hU_bound, unsigned int phi_size, ShortDualHandle* cell_neigbours, unsigned int cell_neighbours_length, Scalar* cell_volume, Scalar* h_advection, Vector* hU_advection){

      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      Scalar h_small = 1e-10;
      Vector2 face_normal[4];
      face_normal[0] = Vector2(0.0, -1.0);
      face_normal[1] = Vector2(1.0, 0.0);
      face_normal[2] = Vector2(0.0, 1.0);
      face_normal[3] = Vector2(-1.0, 0.0);
      while (index < phi_size){
        Scalar manning_this = manning[index];
        Scalar g = gravity[index];
        Scalar h_this = h[index];
        Scalar z_this = z[index];
        Vector2 hU_this = hU[index];
        Vector2 u_this = 0.0;
        if (h_this < h_small){
          u_this = 0.0;
        }
        else{
          u_this = hU_this / h_this;
        }
        Scalar volume = cell_volume[index];
        Scalar area = sqrt(volume);
        Scalar _h_advection(0.0);
        Scalar _h_advection_constraint(0.0);
        Scalar Tr_this = sqrt(g*h_this)/(2.0*g*pow(manning_this,2)*norm(u_this)*pow(h_this,-4.0/3.0))/area;
        Vector2 _hU_advection(0.0, 0.0);
        for (Flag i = 0; i < 4; ++i){
          ShortDualHandle neib = cell_neigbours[i*cell_neighbours_length + index];
          Scalar h_neib = 0.0;
          Scalar z_neib = 0.0;
          Vector2 hU_neib = 0.0;
          Scalar manning_neib = 0.0;
          if (!neib.is_boundary()){
            Flag id_neib = neib.get_global_id();
            h_neib = h[id_neib];
            z_neib = z[id_neib];
            hU_neib = hU[id_neib];
            manning_neib = manning[id_neib];   
          }
          else{
            Flag id_boundary = neib.get_global_id();
            h_neib = _h_bound[id_boundary];
            z_neib = _z_bound[id_boundary];
            hU_neib = _hU_bound[id_boundary];
            manning_neib = manning[id_boundary];   
          }
          if (h_this < h_small && h_neib < h_small){
            continue;
          }
          Vector2 normal = face_normal[i];
          Vector2 shear = uni(perpend(normal));
          Vector2 u_neib = 0.0;
          if (h_neib < h_small){
            u_neib = 0.0;
          }
          else{
            u_neib = hU_neib / h_neib;
          }
          Scalar Tr_neib = sqrt(g*h_neib)/(2.0*g*pow(manning_neib,2)*norm(u_neib)*pow(h_neib,-4.0/3.0))/area;          
          Scalar Tr = fmin(1.0, (Tr_this + Tr_neib)*0.5);
          Vector2 u_L(dot(u_this, normal), dot(u_this, shear));
          Vector2 u_R(dot(u_neib, normal), dot(u_neib, shear));
          Scalar z_f = fmax(z_this, z_neib);
          Scalar h_L = fmax(0.0, h_this + z_this - z_f);
          Scalar h_R = fmax(0.0, h_neib + z_neib - z_f);
          auto flux = cuHLLCRiemannSolverSWEs(g, ScalarRiemannState(h_L, h_R), VectorRiemannState(h_L*u_L, h_R*u_R));
          auto flux_relaxation = cuLaxFriedrichsRiemannSolverSWEs(g, ScalarRiemannState(h_L, h_R), VectorRiemannState(h_L*u_L, h_R*u_R));          
          Scalar _h_flux = (1.0 - Tr)*flux_relaxation.h + Tr*flux.h;
          Vector2 _hU_flux = (1.0 - Tr)*(flux_relaxation.q.x*normal + flux_relaxation.q.y*shear) + Tr*(flux.q.x*normal + flux.q.y*shear);
          Vector2 _z_flux = -0.5*g*h_L*h_L*normal;
          _h_advection += _h_flux*area / volume;
          _hU_advection += (_hU_flux + _z_flux)*area / volume;
        }
        h_advection[index] = _h_advection;
        hU_advection[index] = _hU_advection;
        __syncthreads();
        index += blockDim.x * gridDim.x;
      }
    }


    void cuAdvectionSWEsAPCartesian(cuFvMappedField<Scalar, on_cell>& manning, cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Scalar, on_cell>& z, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Scalar, on_cell>& h_advection, cuFvMappedField<Vector, on_cell>& hU_advection){

      auto mesh = h.mesh;

      cuAdvectionSWEsAPCartesianKernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(
        manning.data.dev_ptr(),
        gravity.data.dev_ptr(),
        h.data.dev_ptr(),
        h.boundary_value.dev_ptr(),
        z.data.dev_ptr(),
        z.boundary_value.dev_ptr(),
        hU.data.dev_ptr(),
        hU.boundary_value.dev_ptr(),
        h.data.size(),
        mesh->cell_neighbours.dev_ptr(),
        mesh->cell_neighbours.length(),
        mesh->cell_volumes.dev_ptr(),
        h_advection.data.dev_ptr(),
        hU_advection.data.dev_ptr());
    }


    __global__ void cuAdvectionMSWEsCartesianKernel(Scalar* gravity, Scalar* h, Scalar* _h_bound, Scalar* z, Scalar* _z_bound, Vector* z_gradient, Vector* hU, Vector* _hU_bound, unsigned int phi_size, ShortDualHandle* cell_neigbours, unsigned int cell_neighbours_length, Scalar* cell_volume, Scalar* h_advection, Vector* hU_advection){

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
        if (h_this < h_small){
          u_this = 0.0;
        }
        else{
          u_this = hU_this / h_this;
        }
        Scalar volume = cell_volume[index];
        Scalar area = sqrt(volume);
        Scalar _h_advection(0.0);
        Scalar _h_advection_constraint(0.0);
        Vector2 _hU_advection(0.0, 0.0);
        for (Flag i = 0; i < 4; ++i){
          Vector2 normal = face_normal[i];
          Vector2 shear = face_shear[i];
          ShortDualHandle neib = cell_neigbours[i*cell_neighbours_length + index];
          Scalar g_neib = 0.0;
          Scalar h_neib = 0.0;
          Scalar z_neib = 0.0;
          Vector2 hU_neib = 0.0;
          Vector2 _z_gradient_neib = 0.0;
          if (!neib.is_boundary()){
            Flag id_neib = neib.get_global_id();
            g_neib = gravity[id_neib];
            h_neib = h[id_neib];
            z_neib = z[id_neib];
            hU_neib = hU[id_neib];
            _z_gradient_neib = z_gradient[id_neib];
          }
          else{
            Flag id_boundary = neib.get_global_id();
            g_neib = g_this;
            h_neib = _h_bound[id_boundary];
            z_neib = _z_bound[id_boundary];
            hU_neib = _hU_bound[id_boundary];
          }
          if (h_this < h_small && h_neib < h_small){
            continue;
          }
          Scalar eta_neib = z_neib + h_neib;
          Vector2 u_neib = 0.0;
          if (h_neib < h_small){
            u_neib = 0.0;
          }
          else{
            u_neib = hU_neib / h_neib;
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
          auto flux = cuHLLCRiemannSolverSWEs(g, ScalarRiemannState(h_L, h_R), VectorRiemannState(h_L*u_L, h_R*u_R));
          Scalar _h_flux = flux.h;
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
        }
        h_advection[index] = _h_advection;
        hU_advection[index] = _hU_advection;
        __syncthreads();
        index += blockDim.x * gridDim.x;
      }
    }

    void cuAdvectionMSWEsCartesian(cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Scalar, on_cell>& z, cuFvMappedField<Vector, on_cell>& z_gradient, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Scalar, on_cell>& h_advection, cuFvMappedField<Vector, on_cell>& hU_advection){

      auto mesh = h.mesh;

      cuAdvectionMSWEsCartesianKernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(gravity.data.dev_ptr(),
        h.data.dev_ptr(),
        h.boundary_value.dev_ptr(),
        z.data.dev_ptr(),
        z.boundary_value.dev_ptr(),
        z_gradient.data.dev_ptr(),
        hU.data.dev_ptr(),
        hU.boundary_value.dev_ptr(),
        h.data.size(),
        mesh->cell_neighbours.dev_ptr(),
        mesh->cell_neighbours.length(),
        mesh->cell_volumes.dev_ptr(),
        h_advection.data.dev_ptr(),
        hU_advection.data.dev_ptr());
    }

    __global__ void cuAdvectionNSWEsSRMCartesianKernel(Scalar* gravity, Scalar* centrifugal, Scalar* h, Scalar* _h_bound, Scalar* z, Scalar* _z_bound, Vector* z_gradient, Vector* hU, Vector* _hU_bound, unsigned int phi_size, ShortDualHandle* cell_neigbours, unsigned int cell_neighbours_length, Scalar* cell_volume, Scalar* h_advection, Vector* hU_advection){

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
        Scalar centrifugal_this = centrifugal[index];
        Vector2 _z_gradient_this = z_gradient[index];
        Scalar eta_this = h_this + z_this;
        Vector2 hU_this = hU[index];
        Vector2 u_this = 0.0;
        if (h_this < h_small){
          u_this = 0.0;
        }
        else{
          u_this = hU_this / h_this;
        }
        Scalar volume = cell_volume[index];
        Scalar area = sqrt(volume);
        Scalar _h_advection(0.0);
        Scalar _h_advection_constraint(0.0);
        Vector2 _hU_advection(0.0, 0.0);
        for (Flag i = 0; i < 4; ++i){
          Vector2 normal = face_normal[i];
          Vector2 shear = face_shear[i];
          ShortDualHandle neib = cell_neigbours[i*cell_neighbours_length + index];
          Scalar g_neib = 0.0;
          Scalar h_neib = 0.0;
          Scalar z_neib = 0.0;
          Vector2 hU_neib = 0.0;
          Vector2 _z_gradient_neib = 0.0;
          if (!neib.is_boundary()){
            Flag id_neib = neib.get_global_id();
            g_neib = gravity[id_neib];
            h_neib = h[id_neib];
            z_neib = z[id_neib];
            hU_neib = hU[id_neib];
            _z_gradient_neib = z_gradient[id_neib];
          }
          else{
            Flag id_boundary = neib.get_global_id();
            g_neib = g_this;
            h_neib = _h_bound[id_boundary];
            z_neib = _z_bound[id_boundary];
            hU_neib = _hU_bound[id_boundary];
          }
          if (h_this < h_small && h_neib < h_small){
            continue;
          }
          Scalar eta_neib = z_neib + h_neib;
          Vector2 u_neib = 0.0;
          if (h_neib < h_small){
            u_neib = 0.0;
          }
          else{
            u_neib = hU_neib / h_neib;
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
          Scalar deta_neib = fmax(0.0, fmin(-dz, -eta_neib + eta_this));
          Scalar eta_L = eta_this + deta_this;
          Scalar eta_R = eta_neib + deta_neib;
          Scalar h_L = fmax(0.0, eta_L - z_f);
          Scalar h_R = fmax(0.0, eta_R - z_f);
          Vector2 u_L(dot(u_this, normal), dot(u_this, shear));
          Vector2 u_R(dot(u_neib, normal), dot(u_neib, shear));
          Scalar g = 0.5*(g_this + g_neib);
          auto flux = cuHLLCRiemannSolverSWEs(g, ScalarRiemannState(h_L, h_R), VectorRiemannState(h_L*u_L, h_R*u_R));
          Scalar _h_flux = flux.h;
          Vector2 _hU_flux = (flux.q.x*normal + flux.q.y*shear);
          if (h_neib < h_small){
            delta_z = fmax(0.0, z_f - eta_this);
          }
          else{
            delta_z = fmax(0.0, fmin(dz_clip, z_f - eta_this));
          }
          z_f -= delta_z;
          Scalar eta_f = z_f + h_L;
          Vector2 _z_flux = -0.5*(g + centrifugal_this)*h_L*h_L*normal + (g_this + centrifugal_this)*h_this*eta_f*normal;
          _h_advection += _h_flux*area / volume;
          _hU_advection += (_hU_flux + _z_flux)*area / volume;
        }
        h_advection[index] = _h_advection;
        hU_advection[index] = _hU_advection;
        __syncthreads();
        index += blockDim.x * gridDim.x;
      }
    }

    void cuAdvectionNSWEsSRMCartesian(cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& centrifugal, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Scalar, on_cell>& z, cuFvMappedField<Vector, on_cell>& z_gradient, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Scalar, on_cell>& h_advection, cuFvMappedField<Vector, on_cell>& hU_advection){

      auto mesh = h.mesh;

      cuAdvectionNSWEsSRMCartesianKernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(gravity.data.dev_ptr(),
        centrifugal.data.dev_ptr(),
        h.data.dev_ptr(),
        h.boundary_value.dev_ptr(),
        z.data.dev_ptr(),
        z.boundary_value.dev_ptr(),
        z_gradient.data.dev_ptr(),
        hU.data.dev_ptr(),
        hU.boundary_value.dev_ptr(),
        h.data.size(),
        mesh->cell_neighbours.dev_ptr(),
        mesh->cell_neighbours.length(),
        mesh->cell_volumes.dev_ptr(),
        h_advection.data.dev_ptr(),
        hU_advection.data.dev_ptr());
    }
    __global__ void cuAdvectionNSWEs2ndCartesianKernel(Scalar* gravity, Scalar* _gravity_bound, Scalar* h, Scalar* _h_bound, Scalar* z, Scalar* _z_bound, Vector* u, Vector* _u_bound, Vector* h_gradient, Vector* eta_gradient, Tensor* u_gradient, unsigned int phi_size, ShortDualHandle* cell_neigbours, unsigned int cell_neighbours_length, Scalar* cell_volume, Scalar* h_advection, Vector* hU_advection){

      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      while (index < phi_size){
        Scalar g_this = gravity[index];
        Scalar h_this = h[index];
        Scalar z_this = z[index];
        Scalar eta_this = h_this + z_this;
        Vector2 u_this = u[index];
        Vector2 grad_h_this = h_gradient[index];
        Vector2 grad_eta_this = eta_gradient[index];
        Tensor2 grad_u_this = u_gradient[index];
        Scalar volume_this = cell_volume[index];
        Scalar area = sqrt(volume_this);
        Scalar cell_size = sqrt(volume_this);
        Vector2 face_normal[4];
        face_normal[0] = Vector2(0.0, -1.0);
        face_normal[1] = Vector2(1.0, 0.0);
        face_normal[2] = Vector2(0.0, 1.0);
        face_normal[3] = Vector2(-1.0, 0.0);
        Scalar _h_advection(0.0);
        Vector2 _hU_advection(0.0, 0.0);
        for (Flag i = 0; i < 4; ++i){
          Vector2 normal = face_normal[i];
          Vector2 shear = uni(perpend(normal));
          Vector2 direction_this = 0.5*normal*cell_size;
          Scalar _h_this = h_this + dot(grad_h_this, direction_this);
          Scalar _eta_this = eta_this + dot(grad_eta_this, direction_this);
          Scalar _z_this = _eta_this - _h_this;
          Vector2 _u_this = u_this + dot(grad_u_this, direction_this);
          ShortDualHandle neib = cell_neigbours[i*cell_neighbours_length + index];
          if (!neib.is_boundary()){
            Flag id_neib = neib.get_global_id();
            Scalar g_neib = gravity[id_neib];
            Scalar h_neib = h[id_neib];
            Scalar z_neib = z[id_neib];
            Scalar eta_neib = h_neib + z_neib;
            Vector2 u_neib = u[id_neib];
            Vector2 grad_h_neib = h_gradient[id_neib];
            Vector2 grad_eta_neib = eta_gradient[id_neib];
            Tensor2 grad_u_neib = u_gradient[id_neib];
            Vector2 direction_neib = -0.5*normal*cell_size;
            Scalar _h_neib = h_neib + dot(grad_h_neib, direction_neib);
            Scalar _eta_neib = eta_neib + dot(grad_eta_neib, direction_neib);
            Scalar _z_neib = _eta_neib - _h_neib;
            Vector2 _u_neib = u_neib + dot(grad_u_neib, direction_neib);
            Scalar grav = (g_this + g_neib) / 2.0;
            Scalar z_face = fmax(_z_this, _z_neib);
            Scalar h_L = fmax(_h_this + _z_this - z_face, 0.0);
            Scalar h_R = fmax(_h_neib + _z_neib - z_face, 0.0);
            Vector2 u_L(dot(_u_this, normal), dot(_u_this, shear));
            Vector2 u_R(dot(_u_neib, normal), dot(_u_neib, shear));
            auto flux = cuHLLCRiemannSolverSWEs(grav, ScalarRiemannState(h_L, h_R), VectorRiemannState(h_L*u_L, h_R*u_R));
            Vector2 _hU_flux = flux.q.x*normal + flux.q.y*shear;
            Vector2 _z_flux = -0.5*grav*h_L*h_L*normal + 0.5*g_this*(_eta_this - z_this)*(_eta_this - z_this)*normal;
            _h_advection += flux.h*area / volume_this;
            _hU_advection += (_hU_flux + _z_flux)*area / volume_this;
          }
          else{
            Flag id_boundary = neib.get_global_id();
            Scalar g_bound = _gravity_bound[id_boundary];
            Scalar h_bound = _h_bound[id_boundary];
            Scalar z_bound = _z_bound[id_boundary];
            Scalar grav = (g_this + g_bound) / 2.0;
            Scalar z_face = fmax(_z_this, z_bound);
            Scalar h_L = fmax(_h_this + _z_this - z_face, 0.0);
            Scalar h_R = fmax(h_bound + z_bound - z_face, 0.0);
            Vector2 u_bound = _u_bound[id_boundary];
            Vector2 u_L(dot(_u_this, normal), dot(_u_this, shear));
            Vector2 u_R(dot(u_bound, normal), dot(u_bound, shear));
            auto flux = cuHLLCRiemannSolverSWEs(grav, ScalarRiemannState(h_L, h_R), VectorRiemannState(h_L*u_L, h_R*u_R));
            Vector2 _hU_flux = flux.q.x*normal + flux.q.y*shear;
            Vector2 _z_flux = -0.5*grav*h_L*h_L*normal + 0.5*g_this*(_eta_this - z_this)*(_eta_this - z_this)*normal;
            _h_advection += flux.h*area / volume_this;
            _hU_advection += (_hU_flux + _z_flux)*area / volume_this;
          }
        }
        h_advection[index] = _h_advection;
        hU_advection[index] = _hU_advection;
        index += blockDim.x * gridDim.x;
      }
    }

    __global__ void cuAdvectionNSWEsOrderReducerSimCartesianKernel(Scalar* h, Vector* h_gradient, Vector* eta_gradient, Tensor2* u_gradient, unsigned int phi_size, ShortDualHandle* cell_neigbours, unsigned int cell_neighbours_length){

      unsigned int index = blockDim.x * blockIdx.x + threadIdx.x;
      Scalar h_small = 1e-10;
      Scalar eps = 0.05;
      Scalar h_tres = 1e-3;
      Scalar h_neibs[4];
      while (index < phi_size){
        Scalar h_this = h[index];
        if (h_this <= h_small){
          h_gradient[index] = 0.0;
          eta_gradient[index] = 0.0;
          u_gradient[index] = 0.0;
        }
        else{
          for (Flag i = 0; i < 4; ++i){
            ShortDualHandle neib = cell_neigbours[i*cell_neighbours_length + index];
            if (!neib.is_boundary()){
              Flag id_neib = neib.get_global_id();
              Scalar h_neib = h[id_neib];
              h_neibs[i] = h_neib;
              if (h_neib <= eps*h_this && h_neib > h_small){
                h_gradient[index] = 0.0;
                eta_gradient[index] = 0.0;
                u_gradient[index] = 0.0;
                h_gradient[id_neib] = 0.0;
                eta_gradient[id_neib] = 0.0;
                u_gradient[id_neib] = 0.0;
              }
            }
          }
          if ((h_neibs[0] <= h_small && h_neibs[2] <= h_small) || (h_neibs[1] <= h_small && h_neibs[3] <= h_small)){
            h_gradient[index] = 0.0;
            eta_gradient[index] = 0.0;
            u_gradient[index] = 0.0;
          }
        }
        index += blockDim.x * gridDim.x;
      }
    }

    void cuAdvectionNSWEs2ndCartesian(cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Scalar, on_cell>& z, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Vector, on_cell>& h_gradient, cuFvMappedField<Vector, on_cell>& eta_gradient, cuFvMappedField<Tensor, on_cell>& hU_gradient, cuFvMappedField<Scalar, on_cell>& h_advection, cuFvMappedField<Vector, on_cell>& hU_advection){

      auto mesh = h.mesh;

      cuAdvectionNSWEs2ndCartesianKernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(gravity.data.dev_ptr(), 
      gravity.boundary_value.dev_ptr(), 
      h.data.dev_ptr(), 
      h.boundary_value.dev_ptr(), 
      z.data.dev_ptr(), 
      z.boundary_value.dev_ptr(), 
      hU.data.dev_ptr(), 
      hU.boundary_value.dev_ptr(), 
      h_gradient.data.dev_ptr(), 
      eta_gradient.data.dev_ptr(), 
      hU_gradient.data.dev_ptr(), 
      h.data.size(), 
      mesh->cell_neighbours.dev_ptr(), 
      mesh->cell_neighbours.length(), 
      mesh->cell_volumes.dev_ptr(), 
      h_advection.data.dev_ptr(), 
      hU_advection.data.dev_ptr());

    }

    __global__ void cuAdvectionSWEsGeorgeCartesianKernel(Scalar* gravity, Scalar* h, Scalar* _h_bound, Scalar* z, Scalar* _z_bound, Vector* hU, Vector* _hU_bound, unsigned int phi_size, ShortDualHandle* cell_neigbours, unsigned int cell_neighbours_length, Scalar* cell_volume, Scalar* h_advection, Vector* hU_advection){

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
        Vector2 hU_this = hU[index];
        Vector2 u_this = 0.0;
        if (h_this < h_small){
          u_this = 0.0;
        }
        else{
          u_this = hU_this / h_this;
        }
        Scalar volume = cell_volume[index];
        Scalar area = sqrt(volume);
        Scalar _h_advection(0.0);
        Vector2 _hU_advection(0.0, 0.0);
        for (Flag i = 0; i < 4; ++i){
          Vector2 normal = face_normal[i];
          Vector2 shear = face_shear[i];
          ShortDualHandle neib = cell_neigbours[i*cell_neighbours_length + index];
          Scalar g_neib = 0.0;
          Scalar h_neib = 0.0;
          Scalar z_neib = 0.0;
          Vector2 hU_neib = 0.0;
          if (!neib.is_boundary()){
            Flag id_neib = neib.get_global_id();
            g_neib = gravity[id_neib];
            h_neib = h[id_neib];
            z_neib = z[id_neib];
            hU_neib = hU[id_neib];
          }
          else{
            Flag id_boundary = neib.get_global_id();
            g_neib = g_this;
            h_neib = _h_bound[id_boundary];
            z_neib = _z_bound[id_boundary];
            hU_neib = _hU_bound[id_boundary];
          }
          if (h_this < h_small && h_neib < h_small){
            continue;
          }
          Vector2 u_neib = 0.0;
          if (h_neib < h_small){
            u_neib = 0.0;
          }
          else{
            u_neib = hU_neib / h_neib;
          }
          Scalar h_L = h_this;
          Scalar h_R = h_neib;
          Scalar z_L = z_this;
          Scalar z_R = z_neib;
          Vector2 u_L(dot(u_this, normal), dot(u_this, shear));
          Vector2 u_R(dot(u_neib, normal), dot(u_neib, shear));
          Scalar g = 0.5*(g_this + g_neib);
          auto flux = cuGeorgeRiemannSolverSWEs(g, ScalarRiemannState(h_L, h_R), VectorRiemannState(h_L*u_L, h_R*u_R), ScalarRiemannState(z_L, z_R));
          Scalar _h_flux = flux.h;
          Vector2 _hU_flux = (flux.q.x*normal + flux.q.y*shear);
          _h_advection += _h_flux*area / volume;
          _hU_advection += _hU_flux*area / volume;
        }
        h_advection[index] = _h_advection;
        hU_advection[index] = _hU_advection;
        __syncthreads();
        index += blockDim.x * gridDim.x;
      }
    }


    void cuAdvectionSWEsGeorgeCartesian(cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Scalar, on_cell>& z, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Scalar, on_cell>& h_advection, cuFvMappedField<Vector, on_cell>& hU_advection){

      auto mesh = h.mesh;

      cuAdvectionSWEsGeorgeCartesianKernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(gravity.data.dev_ptr(),
        h.data.dev_ptr(),
        h.boundary_value.dev_ptr(),
        z.data.dev_ptr(),
        z.boundary_value.dev_ptr(),
        hU.data.dev_ptr(),
        hU.boundary_value.dev_ptr(),
        h.data.size(),
        mesh->cell_neighbours.dev_ptr(),
        mesh->cell_neighbours.length(),
        mesh->cell_volumes.dev_ptr(),
        h_advection.data.dev_ptr(),
        hU_advection.data.dev_ptr());
    }

__global__ void cuAdvectionSWEsRFGeorgeCartesianKernel(Scalar* gravity, Scalar* h, Scalar* _h_bound, Scalar* z, Scalar* _z_bound, Vector* hU, Vector* _hU_bound, unsigned int phi_size, ShortDualHandle* cell_neigbours, unsigned int cell_neighbours_length, Scalar* cell_volume, Scalar* h_advection, Vector* hU_advection, Vector* hU_EW, Vector* hU_NS){

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
        Vector2 hU_this = hU[index];
        Vector2 u_this = 0.0;
        if (h_this < h_small){
          u_this = 0.0;
        }
        else{
          u_this = hU_this / h_this;
        }
        Scalar volume = cell_volume[index];
        Scalar area = sqrt(volume);
        Scalar _h_advection(0.0);
        Vector2 _hU_advection(0.0, 0.0);
        for (Flag i = 0; i < 4; ++i){
          Vector2 normal = face_normal[i];
          Vector2 shear = face_shear[i];
          ShortDualHandle neib = cell_neigbours[i*cell_neighbours_length + index];
          Scalar g_neib = 0.0;
          Scalar h_neib = 0.0;
          Scalar z_neib = 0.0;
          Vector2 hU_neib = 0.0;
          if (!neib.is_boundary()){
            Flag id_neib = neib.get_global_id();
            g_neib = gravity[id_neib];
            h_neib = h[id_neib];
            z_neib = z[id_neib];
            hU_neib = hU[id_neib];
          }
          else{
            Flag id_boundary = neib.get_global_id();
            g_neib = g_this;
            h_neib = _h_bound[id_boundary];
            z_neib = _z_bound[id_boundary];
            hU_neib = _hU_bound[id_boundary];
          }
          if (h_this < h_small && h_neib < h_small){
            continue;
          }
          Vector2 u_neib = 0.0;
          if (h_neib < h_small){
            u_neib = 0.0;
          }
          else{
            u_neib = hU_neib / h_neib;
          }
          Scalar h_L = h_this;
          Scalar h_R = h_neib;
          Scalar z_L = z_this;
          Scalar z_R = z_neib;
          Vector2 u_L(dot(u_this, normal), dot(u_this, shear));
          Vector2 u_R(dot(u_neib, normal), dot(u_neib, shear));
          Scalar g = 0.5*(g_this + g_neib);
          auto flux = cuGeorgeRiemannSolverSWEs(g, ScalarRiemannState(h_L, h_R), VectorRiemannState(h_L*u_L, h_R*u_R), ScalarRiemannState(z_L, z_R));
          Scalar _h_flux = flux.h;
          Vector2 _hU_flux = (flux.q.x*normal + flux.q.y*shear);
          if(i == 0){
            hU_NS[index].y = flux.h;
          }
          if(i == 2){
            hU_NS[index].x = flux.h;
          }
          if(i == 1){
            hU_EW[index].x = flux.h;
          }
          if(i == 3){
            hU_EW[index].y = flux.h;
          }
          _h_advection += _h_flux*area / volume;
          _hU_advection += _hU_flux*area / volume;
        }
        h_advection[index] = _h_advection;
        hU_advection[index] = _hU_advection;
        __syncthreads();
        index += blockDim.x * gridDim.x;
      }
    }


    void cuAdvectionSWEsRFGeorgeCartesian(cuFvMappedField<Scalar, on_cell>& gravity, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Scalar, on_cell>& z, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Scalar, on_cell>& h_advection, cuFvMappedField<Vector, on_cell>& hU_advection, cuFvMappedField<Vector, on_cell>& hU_EW, cuFvMappedField<Vector, on_cell>& hU_NS){

      auto mesh = h.mesh;

      cuAdvectionSWEsRFGeorgeCartesianKernel << <BLOCKS_PER_GRID, THREADS_PER_BLOCK >> >(gravity.data.dev_ptr(),
        h.data.dev_ptr(),
        h.boundary_value.dev_ptr(),
        z.data.dev_ptr(),
        z.boundary_value.dev_ptr(),
        hU.data.dev_ptr(),
        hU.boundary_value.dev_ptr(),
        h.data.size(),
        mesh->cell_neighbours.dev_ptr(),
        mesh->cell_neighbours.length(),
        mesh->cell_volumes.dev_ptr(),
        h_advection.data.dev_ptr(),
        hU_advection.data.dev_ptr(),
        hU_EW.data.dev_ptr(),
        hU_NS.data.dev_ptr());
    }

  }//end of namespace fv--------


}//--end of namespace GC--------

