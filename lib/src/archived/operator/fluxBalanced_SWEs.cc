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
  \flie flux_SWEs_MCL.cc
  \brief Source file for shallow water equations flux

  \version 0.1
  \author xilin xia
*/

#include "fluxBalanced_SWEs.h"
#include <algorithm>
#include "Flag.h"
#include "boundary.h"


namespace GC{

  namespace fv{

    void FluxBalancedSWEs_1st::operator()(Scalar g,
                                  fvScalarFieldOnCell& h,
                                  fvVectorFieldOnCell& hU,
                                  fvScalarFieldOnCell& z,
                                  fvScalarFieldOnCell& h_flux,
                                  fvVectorFieldOnCell& hU_flux){
          auto h_begin = h.data_begin();
          auto h_boundary_type_begin = h.boundary_type_begin();
          auto hU_begin = hU.data_begin();
          auto hU_boundary_type_begin = hU.boundary_type_begin();
          auto z_begin = z.data_begin();
          auto z_boundary_type_begin = z.boundary_type_begin();
          auto h_flux_begin = h_flux.data_begin();
          auto hU_flux_begin = hU_flux.data_begin();
          auto faces_begin = h.mesh.Cell.HalfFacets.begin();
          auto face_normal_begin = h.mesh.HalfFacet.Normal.begin();
          auto face_area_begin = h.mesh.HalfFacet.Area.begin();
          auto cell_neighbour_begin = h.mesh.Cell.Neighbours.begin();
          auto cell_volume_begin = h.mesh.Cell.Volume.begin();
          //initialize the fluxes

          for(auto h_flux_iter = h_flux_begin; h_flux_iter < h_flux.data_end(); ++h_flux_iter){
            *(h_flux_iter) = 0.0;
          }

          for(auto hU_flux_iter = hU_flux_begin; hU_flux_iter < hU_flux.data_end(); ++hU_flux_iter){
            *(hU_flux_iter) = Vector2(0.0, 0.0);
          }

          //iterate all the cells
          for(auto cell_iter = h.mesh.Cell.Neighbours.begin();
                   cell_iter < h.mesh.Cell.Neighbours.end();
                   ++cell_iter){
            Flag id_this = cell_iter - cell_neighbour_begin;
            Scalar h_this = *(h_begin + id_this);
            Scalar z_this = *(z_begin + id_this);
            Vector2 hU_this = *(hU_begin + id_this);
            Vector2 u_this, u_neib, u_bound;
            if(h_this == 0.0){
              u_this = 0.0;
            }else{
              u_this = hU_this/h_this;
            }
            Scalar volume_this = *(cell_volume_begin + id_this);
            //iterate all the neighbours
            Flag neib_cnt = 0;
            for(auto neib_iter : *cell_iter){              
              if (!neib_iter.is_boundary()){ //neighbour is not a boundary
                Flag id_neib = neib_iter.get_global_id();
                if (id_this < id_neib){ // preventing repetive calculation
                  Flag id_face = (*(faces_begin + id_this))[neib_cnt];
                  Vector2 normal = uni(*(face_normal_begin + id_face));
                  Vector2 shear = uni(perpend(normal)); 
                  Scalar area = *(face_area_begin + id_face);               
                  Scalar h_neib = *(h_begin + id_neib);
                  Scalar z_neib = *(z_begin + id_neib);
                  Scalar z_face = std::max(z_this, z_neib);
                  Scalar h_L = std::max(h_this + z_this - z_face, 0.0);
                  Scalar h_R = std::max(h_neib + z_neib - z_face, 0.0);
                  Scalar volume_neib = *(cell_volume_begin + id_neib);
                  Vector2 hU_neib = *(hU_begin + id_neib);
                  if(h_neib == 0.0){
                    u_neib = 0.0;
                  }else{
                    u_neib = hU_neib/h_neib;
                  }
                  Vector2 u_L(dot(u_this, normal),dot(u_this, shear));
                  Vector2 u_R(dot(u_neib, normal),dot(u_neib, shear));
                  auto flux = RiemannSolver(g, ScalarRiemannState(h_L, h_R), VectorRiemannState(h_L*u_L, h_R*u_R)); 
                  Vector2 _hU_flux = flux.q.x*normal + flux.q.y*shear;
                  Vector2 _z_flux_this = -0.5*g*h_L*h_L*normal;
                  Vector2 _z_flux_neib = -0.5*g*h_R*h_R*normal;
                  *(h_flux_begin + id_this) += flux.h*area/volume_this;
                  (*(hU_flux_begin + id_this)).x += (_hU_flux.x + _z_flux_this.x)*area/volume_this;
                  (*(hU_flux_begin + id_this)).y += (_hU_flux.y + _z_flux_this.y)*area/volume_this;
                  *(h_flux_begin + id_neib) -= flux.h*area/volume_neib;
                  (*(hU_flux_begin + id_neib)).x -= (_hU_flux.x + _z_flux_neib.x)*area/volume_neib;
                  (*(hU_flux_begin + id_neib)).y -= (_hU_flux.y + _z_flux_neib.y)*area/volume_neib;
                }
              }else{
                Flag id_face = (*(faces_begin + id_this))[neib_cnt];
                Vector2 normal = uni(*(face_normal_begin + id_face));
                Vector2 shear = uni(perpend(normal));
                Scalar area = *(face_area_begin + id_face);
                Flag id_boundary = neib_iter.get_global_id();
                auto h_boundary_type = *(h_boundary_type_begin + id_boundary);
                auto z_boundary_type = *(z_boundary_type_begin + id_boundary);
                auto hU_boundary_type = *(hU_boundary_type_begin + id_boundary);
                Scalar h_bound = GetBoundary(h_boundary_type, h_this, normal);
                Scalar z_bound = GetBoundary(z_boundary_type, z_this, normal);
                Scalar z_face = std::max(z_this, z_bound);
                Scalar h_L = std::max(h_this + z_this - z_face, 0.0);
                Scalar h_R = std::max(h_bound + z_bound - z_face, 0.0);
                Vector2 hU_bound = GetBoundary(hU_boundary_type, hU_this, normal);
                if(h_bound == 0.0){
                  u_bound = 0.0;
                }else{
                  u_bound = hU_bound/h_bound;
                }
                Vector2 u_L(dot(u_this, normal),dot(u_this, shear));
                Vector2 u_R(dot(u_bound, normal),dot(u_bound, shear));
                auto flux = RiemannSolver(g, ScalarRiemannState(h_L, h_R), VectorRiemannState(h_L*u_L, h_R*u_R)); 
                Vector2 _hU_flux = flux.q.x*normal + flux.q.y*shear;
                Vector2 _z_flux_this = -0.5*g*h_L*h_L*normal;
                *(h_flux_begin + id_this) += flux.h*area/volume_this;
                (*(hU_flux_begin + id_this)).x += (_hU_flux.x + _z_flux_this.x)*area/volume_this;
                (*(hU_flux_begin + id_this)).y += (_hU_flux.y + _z_flux_this.y)*area/volume_this;
              }
              neib_cnt++;
            }
          }
    }

  }//end of namespace

}//end of namespace
