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
  \flie fluxNonUniGravity_SWEs.cc
  \brief Source file for shallow water equations flux

  \version 0.1
  \author xilin xia
*/

#include "fluxNonUniGravity_SWEs.h"
#include "Flag.h"
#include "boundary.h"


namespace GC{//--namespace GeoX-------------------------------------

  namespace fv{//--namespace fv---------------------------------------
        
        void FluxNonUniGravitySWEs_1st::operator() (fvScalarFieldOnCell& g,
                                       fvScalarFieldOnCell& h, 
                                       fvVectorFieldOnCell& hU, 
                                       fvScalarFieldOnCell& h_flux,
                                       fvVectorFieldOnCell& hU_flux){
          auto g_begin = g.data_begin();
          auto g_boundary_type_begin = g.boundary_type_begin();
          auto h_begin = h.data_begin();
          auto h_boundary_type_begin = h.boundary_type_begin();
          auto hU_begin = hU.data_begin();
          auto hU_boundary_type_begin = hU.boundary_type_begin();
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

//#pragma omp parallel for firstprivate(h_begin, h_boundary_type_begin, hU_begin, hU_boundary_type_begin, h_flux_begin, hU_flux_begin, faces_begin, face_normal_begin, face_area_begin, cell_neighbour_begin, cell_volume_begin)
          //iterate all the cells
          for(auto cell_iter = h.mesh.Cell.Neighbours.begin();
                   cell_iter < h.mesh.Cell.Neighbours.end();
                   ++cell_iter){
            Flag id_this = cell_iter - cell_neighbour_begin;
            Scalar h_this = *(h_begin + id_this);
            Scalar g_this = *(g_begin + id_this);
            Vector2 hU_this = *(hU_begin + id_this);
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
                  Scalar g_neib = *(g_begin + id_neib);
                  Scalar h_neib = *(h_begin + id_neib);
                  Scalar volume_neib = *(cell_volume_begin + id_neib);
                  Vector2 hU_neib = *(hU_begin + id_neib);
                  Vector2 hU_L(dot(hU_this, normal),dot(hU_this, shear));
                  Vector2 hU_R(dot(hU_neib, normal),dot(hU_neib, shear));
                  auto flux = RiemannSolver(ScalarRiemannState(g_this, g_neib), ScalarRiemannState(h_this, h_neib), VectorRiemannState(hU_L, hU_R)); 
                  Vector2 _hU_flux = flux.q.x*normal + flux.q.y*shear;
//#pragma omp atomic
                  *(h_flux_begin + id_this) += flux.h*area/volume_this;
//#pragma omp atomic
                  (*(hU_flux_begin + id_this)).x += _hU_flux.x*area/volume_this;
//#pragma omp atomic
                  (*(hU_flux_begin + id_this)).y += _hU_flux.y*area/volume_this;
//#pragma omp atomic
                  *(h_flux_begin + id_neib) -= flux.h*area/volume_neib;
//#pragma omp atomic
                  (*(hU_flux_begin + id_neib)).x -= _hU_flux.x*area/volume_neib;
//#pragma omp atomic
                  (*(hU_flux_begin + id_neib)).y -= _hU_flux.y*area/volume_neib;
                }
              }else{
                Flag id_face = (*(faces_begin + id_this))[neib_cnt];
                Vector2 normal = uni(*(face_normal_begin + id_face));
                Vector2 shear = uni(perpend(normal));
                Scalar area = *(face_area_begin + id_face);
                Flag id_boundary = neib_iter.get_global_id();
                auto g_boundary_type = *(g_boundary_type_begin + id_boundary);
                auto h_boundary_type = *(h_boundary_type_begin + id_boundary);
                auto hU_boundary_type = *(hU_boundary_type_begin + id_boundary);
                Scalar g_bound = GetBoundary(g_boundary_type, g_this, normal);
                Scalar h_bound = GetBoundary(h_boundary_type, h_this, normal);
                Vector2 hU_bound = GetBoundary(hU_boundary_type, hU_this, normal);
                Vector2 hU_L(dot(hU_this, normal),dot(hU_this, shear));
                Vector2 hU_R(dot(hU_bound, normal),dot(hU_bound, shear));
                auto flux = RiemannSolver(ScalarRiemannState(g_this, g_bound), ScalarRiemannState(h_this, h_bound), VectorRiemannState(hU_L, hU_R)); 
                Vector2 _hU_flux = flux.q.x*normal + flux.q.y*shear;
//#pragma omp atomic
                *(h_flux_begin + id_this) += flux.h*area/volume_this;
//#pragma omp atomic
                (*(hU_flux_begin + id_this)).x += _hU_flux.x*area/volume_this;
//#pragma omp atomic
                (*(hU_flux_begin + id_this)).y += _hU_flux.y*area/volume_this;
              }
              neib_cnt++;
            }
          }
        }
  } //--end namespace fv----------------------------------------------

}//--end namespace GeoX-----------------------------------------------
  


