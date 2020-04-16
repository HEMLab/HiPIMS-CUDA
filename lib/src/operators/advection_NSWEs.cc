// ======================================================================================
// Name                :    GeoClasses : Generic Geophysical Flow Modelling Framework
// Description         :    This code pack provides a generic framework for developing 
//                          Geophysical CFD software.
// ======================================================================================
// Version             :    0.1 
// Author              :    Xilin Xia (PhD candidate in Newcastle University)
// Create Time         :    2014/10/04
// Update Time         :    2015/10/09
// ======================================================================================
// Copyright @ Xilin Xia 2015 . All rights reserved.
// ======================================================================================

/*!
  \flie advection_NSWEs_2nd.cc
  \brief Source file for flux in non hydrostatic shallow water equations

  \version 0.1
  \author xilin xia
*/

#include "advection_NSWEs.h"
#include <algorithm>
#include "Flag.h"
#include "boundary.h"

namespace GC{

  namespace fv{

    const Scalar small_value = 1e-10;

    void AdvectionNSWEs_2nd(fvScalarFieldOnCell& g, fvScalarFieldOnCell& h, fvVectorFieldOnCell& hU, fvScalarFieldOnCell& z, fvVectorFieldOnCell& grad_h, fvVectorFieldOnCell& grad_surface, fvTensorFieldOnCell& grad_u, fvScalarFieldOnCell& h_flux, 
                        fvVectorFieldOnCell& hU_flux){
      auto g_begin = g.data_begin();
      auto g_boundary_value_begin = g.boundary_value_begin();
      auto h_begin = h.data_begin();
      auto h_boundary_value_begin = h.boundary_value_begin();
      auto hU_begin = hU.data_begin();
      auto hU_boundary_value_begin = hU.boundary_value_begin();
      auto z_begin = z.data_begin();
      auto z_boundary_value_begin = z.boundary_value_begin();
      auto grad_h_begin = grad_h.data_begin();
      auto grad_surface_begin = grad_surface.data_begin();
      auto grad_u_begin = grad_u.data_begin();
      auto h_flux_begin = h_flux.data_begin();
      auto hU_flux_begin = hU_flux.data_begin();
      auto faces_begin = h.mesh.Cell.HalfFacets.begin();
      auto face_normal_begin = h.mesh.HalfFacet.Normal.begin();
      auto face_area_begin = h.mesh.HalfFacet.Area.begin();
      auto face_centroid_begin = h.mesh.HalfFacet.Centroid.begin();
      auto cell_neighbour_begin = h.mesh.Cell.Neighbours.begin();
      auto cell_volume_begin = h.mesh.Cell.Volume.begin();
      auto cell_centroid_begin = h.mesh.Cell.Centroid.begin();

      //initialize the fluxes
#pragma omp parallel for      
      for(auto h_flux_iter = h_flux.data_begin(); h_flux_iter < h_flux.data_end(); ++h_flux_iter){
            *(h_flux_iter) = 0.0;
      }

#pragma omp parallel for
      for(auto hU_flux_iter = hU_flux.data_begin(); hU_flux_iter < hU_flux.data_end(); ++hU_flux_iter){
            *(hU_flux_iter) = Vector2(0.0, 0.0);
      }

      //iterate all the cells
#pragma omp parallel for firstprivate(g_begin, g_boundary_value_begin, h_begin, h_boundary_value_begin, \
                                      hU_begin, hU_boundary_value_begin, z_begin, z_boundary_value_begin, \
                                      grad_h_begin, grad_surface_begin, grad_u_begin, \
                                      h_flux_begin, hU_flux_begin, faces_begin, face_normal_begin, face_area_begin, \
                                      face_centroid_begin, cell_neighbour_begin, cell_volume_begin, cell_centroid_begin)

      for(auto cell_iter = h.mesh.Cell.Neighbours.begin();
               cell_iter < h.mesh.Cell.Neighbours.end();
               ++cell_iter){
        Flag id_this = cell_iter - cell_neighbour_begin;
        Scalar g_this = *(g_begin + id_this);
        Scalar h_this = *(h_begin + id_this);
        Scalar z_this = *(z_begin + id_this);
        Scalar eta_this = h_this + z_this;
        Vector2 hU_this = *(hU_begin + id_this);
        Vector2 u_this;
        if(fabs(h_this) < small_value){
          u_this = Vector2(0.0, 0.0);
        }else{
          u_this = hU_this/h_this;
        }
        Vector2 grad_h_this = *(grad_h_begin + id_this);
        Vector2 grad_eta_this = *(grad_surface_begin + id_this);
        Tensor2 grad_u_this = *(grad_u_begin + id_this);
        Scalar volume_this = *(cell_volume_begin + id_this);
        Vector2 cell_centroid_this = *(cell_centroid_begin + id_this);
        //iterate all the neighbours
        Flag neib_cnt = 0;
        for(auto neib_iter : *cell_iter){
          Flag id_face = (*(faces_begin + id_this))[neib_cnt];
          Vector2 normal = uni(*(face_normal_begin + id_face));
          Vector2 shear = uni(perpend(normal));
          Scalar area = *(face_area_begin + id_face);
          Vector2 face_centroid = *(face_centroid_begin + id_face);
          Vector2 direction_this = face_centroid - cell_centroid_this;
          Scalar _h_this = h_this + dot(grad_h_this,direction_this);
          Scalar _eta_this = eta_this + dot(grad_eta_this,direction_this);
          Scalar _z_this = _eta_this - _h_this;
          Vector2 _u_this = u_this + dot(grad_u_this,direction_this);
          if (!neib_iter.is_boundary()){
            Flag id_neib = neib_iter.get_global_id();
            if (id_this < id_neib){
              Scalar volume_neib = *(cell_volume_begin + id_neib);
              Scalar g_neib = *(g_begin + id_neib);
              Scalar h_neib = *(h_begin + id_neib);
              Scalar z_neib = *(z_begin + id_neib);
              Scalar eta_neib = h_neib + z_neib;
              Vector2 hU_neib = *(hU_begin + id_neib);
              Vector2 u_neib;
              if(fabs(h_neib) < small_value){
                u_neib = Vector2(0.0, 0.0);
              }else{
                u_neib = hU_neib/h_neib;
              }
              Vector2 grad_h_neib = *(grad_h_begin + id_neib);
              Vector2 grad_eta_neib = *(grad_surface_begin + id_neib);
              Tensor2 grad_u_neib = *(grad_u_begin + id_neib);
              Vector2 cell_centroid_neib = *(cell_centroid_begin + id_neib);
              Vector2 direction_neib = face_centroid - cell_centroid_neib;
              Scalar _h_neib = h_neib + dot(grad_h_neib, direction_neib);
              Scalar _eta_neib = eta_neib + dot(grad_eta_neib, direction_neib);
              Scalar _z_neib = _eta_neib - _h_neib;
              Vector2 _u_neib = u_neib + dot(grad_u_neib, direction_neib);
              Scalar grav = (g_this + g_neib)/2.0;
              Scalar z_face = std::max(_z_this, _z_neib);
              Scalar h_L = std::max(_h_this + _z_this - z_face, 0.0);
              Scalar h_R = std::max(_h_neib + _z_neib - z_face, 0.0);
              Vector2 u_L(dot(_u_this, normal),dot(_u_this, shear));
              Vector2 u_R(dot(_u_neib, normal),dot(_u_neib, shear));
              auto flux = hllcRiemannSolverSWEs(grav, ScalarRiemannState(h_L, h_R), VectorRiemannState(h_L*u_L, h_R*u_R));
              Vector2 _hU_flux = flux.q.x*normal + flux.q.y*shear;
              Vector2 _z_flux_this = -0.5*grav*h_L*h_L*normal +0.5*g_this*(_eta_this - z_this)*(_eta_this - z_this)*normal;
              Vector2 _z_flux_neib = -0.5*grav*h_R*h_R*normal +0.5*g_neib*(_eta_neib - z_neib)*(_eta_neib - z_neib)*normal;
#pragma omp atomic
              *(h_flux_begin + id_this) += flux.h*area/volume_this;
#pragma omp atomic
              (*(hU_flux_begin + id_this)).x += (_hU_flux.x + _z_flux_this.x)*area/volume_this;
#pragma omp atomic
              (*(hU_flux_begin + id_this)).y += (_hU_flux.y + _z_flux_this.y)*area/volume_this;
#pragma omp atomic
              *(h_flux_begin + id_neib) -= flux.h*area/volume_neib;
#pragma omp atomic
              (*(hU_flux_begin + id_neib)).x -= (_hU_flux.x + _z_flux_neib.x)*area/volume_neib;
#pragma omp atomic
              (*(hU_flux_begin + id_neib)).y -= (_hU_flux.y + _z_flux_neib.y)*area/volume_neib;
            }
          }else{
            Flag id_face = (*(faces_begin + id_this))[neib_cnt];
            Vector2 normal = uni(*(face_normal_begin + id_face));
            Vector2 shear = uni(perpend(normal));
            Scalar area = *(face_area_begin + id_face);
            Flag id_boundary = neib_iter.get_global_id();
            Scalar g_bound = *(g_boundary_value_begin + id_boundary);;
            Scalar h_bound = *(h_boundary_value_begin + id_boundary);
            Scalar z_bound = *(z_boundary_value_begin + id_boundary);
            Scalar grav = (g_this + g_bound)/2.0;
            Scalar z_face = std::max(_z_this, z_bound);
            Scalar h_L = std::max(_h_this + _z_this - z_face, 0.0);
            Scalar h_R = std::max(h_bound + z_bound - z_face, 0.0);
            Vector2 hU_bound = *(hU_boundary_value_begin + id_boundary);
            Vector2 u_bound;
            if(fabs(h_bound) < small_value){
              u_bound = Vector2(0.0, 0.0);
            }else{
              u_bound = hU_bound/h_bound;
            }
            Vector2 u_L(dot(_u_this, normal),dot(_u_this, shear));
            Vector2 u_R(dot(u_bound, normal),dot(u_bound, shear));
            auto flux = hllcRiemannSolverSWEs(grav, ScalarRiemannState(h_L, h_R), VectorRiemannState(h_L*u_L, h_R*u_R));
            Vector2 _hU_flux = flux.q.x*normal + flux.q.y*shear;
            Vector2 _z_flux_this = -0.5*grav*h_L*h_L*normal +0.5*g_this*(_eta_this - z_this)*(_eta_this - z_this)*normal;
#pragma omp atomic
            *(h_flux_begin + id_this) += flux.h*area/volume_this;
#pragma omp atomic
            (*(hU_flux_begin + id_this)).x += (_hU_flux.x + _z_flux_this.x)*area/volume_this;
#pragma omp atomic
            (*(hU_flux_begin + id_this)).y += (_hU_flux.y + _z_flux_this.y)*area/volume_this;
          }
          neib_cnt++;
        }
        *(hU_flux_begin + id_this) = *(hU_flux_begin + id_this);
      }

    }

  }//--end of namespace fv------

} //--end of namespace GC-------
