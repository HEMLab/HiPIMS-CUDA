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
 \file limiter.cc
 \brief Source file for limiter and limit functions 

  \version 0.1
  \author xilin xia
*/

#include "limiter.h"
#include <cmath>
#include <algorithm>
#include "Scalar.h"
#include "boundary.h"


namespace GC{
  //namespace for finite volume method
  namespace fv{

    Scalar minmod(const Scalar& upwind, const Scalar& downwind){
      Scalar ratio = 0.0;
      if (fabs(downwind) <= 1e-10){
        ratio = 1.0;
      }else{
        ratio = upwind/downwind;
      }

      ratio = std::max( 0.0, std::min(1.0, ratio));
      return ratio;
    }

    fvMappedField<Vector, on_cell>& limiterCartesian::operator() (fvMappedField<Scalar, on_cell>& phi, fvMappedField<Vector, on_cell>& gradient){
      vector_buffer.initialize_by_field(gradient);
      auto phi_begin = phi.data_begin();
      auto phi_boundary_value_begin = phi.boundary_value_begin();
      auto phi_gradient_begin = gradient.data_begin();
      auto limited_gradient_begin = vector_buffer.data_begin();
      auto faces_begin = phi.mesh.Cell.HalfFacets.begin();
      auto face_normal_begin = phi.mesh.HalfFacet.Normal.begin();
      auto face_centroid_begin = phi.mesh.HalfFacet.Centroid.begin();
      auto cell_neighbour_begin = phi.mesh.Cell.Neighbours.begin();
      auto cell_centroid_begin = phi.mesh.Cell.Centroid.begin();

      //initialize the limited gradients
      for(auto limited_gradient_iter = limited_gradient_begin; limited_gradient_iter < vector_buffer.data_end(); ++limited_gradient_iter){
        *(limited_gradient_iter) = Vector(0.0);
      }

      //iterate all the cells
      for(auto cell_iter = phi.mesh.Cell.Neighbours.begin();
               cell_iter < phi.mesh.Cell.Neighbours.end();
               ++cell_iter){
        Flag id_this = cell_iter - cell_neighbour_begin;
        Scalar phi_this = *(phi_begin + id_this);
        Vector centroid_this = *(cell_centroid_begin + id_this);
        Vector gradient_vector = *(phi_gradient_begin + id_this);
        Scalar limit_ratio_x = 1.0;
        Scalar limit_ratio_y = 1.0;
        Scalar _limit_ratio_x = 1.0;
        Scalar _limit_ratio_y = 1.0;
        Flag neib_cnt = 0;
        for(auto neib_iter : *cell_iter){
          if (neib_cnt == 1 || neib_cnt == 3){
            if (!neib_iter.is_boundary()){ //neighbour is not a boundary
              Flag id_neib = neib_iter.get_global_id();
              Scalar phi_neib = *(phi_begin + id_neib);
              Vector centroid_neib = *(cell_centroid_begin + id_neib);
              Vector direction_vector = centroid_neib - centroid_this;
              Scalar downwind = dot(gradient_vector, direction_vector);
              Scalar upwind =  phi_neib - phi_this;
              _limit_ratio_x = limiterFunction(upwind, downwind);
            }else{
              Flag id_face = (*(faces_begin + id_this))[neib_cnt];
              Vector normal = uni(*(face_normal_begin + id_face));
              Flag id_boundary = neib_iter.get_global_id();
              Scalar phi_bound = *(phi_boundary_value_begin + id_boundary);
              Vector centroid_boundary = *(face_centroid_begin + id_face);
              Vector direction_vector = centroid_boundary - centroid_this;
              Scalar downwind = 2.0*dot(gradient_vector, direction_vector);
              Scalar upwind = phi_bound - phi_this;
              _limit_ratio_x = limiterFunction(upwind, downwind);
            }
            limit_ratio_x = std::min(limit_ratio_x, _limit_ratio_x);
          }else{
            if (!neib_iter.is_boundary()){ //neighbour is not a boundary
              Flag id_neib = neib_iter.get_global_id();
              Scalar phi_neib = *(phi_begin + id_neib);
              Vector centroid_neib = *(cell_centroid_begin + id_neib);
              Vector direction_vector = centroid_neib - centroid_this;
              Scalar downwind = dot(gradient_vector, direction_vector);
              Scalar upwind = phi_neib - phi_this;
              _limit_ratio_y = limiterFunction(upwind, downwind);
            }else{
              Flag id_face = (*(faces_begin + id_this))[neib_cnt];
              Vector normal = uni(*(face_normal_begin + id_face));
              Flag id_boundary = neib_iter.get_global_id();
              Scalar phi_bound = *(phi_boundary_value_begin + id_boundary);
              Vector centroid_boundary = *(face_centroid_begin + id_face);
              Vector direction_vector = centroid_boundary - centroid_this;
              Scalar downwind = 2.0*dot(gradient_vector, direction_vector);
              Scalar upwind = phi_bound - phi_this;
              _limit_ratio_y = limiterFunction(upwind, downwind);
            }
            limit_ratio_y = std::min(limit_ratio_y, _limit_ratio_y);
          }
          neib_cnt++;    
        }
        *(limited_gradient_begin + id_this) = Vector(limit_ratio_x*gradient_vector.x, limit_ratio_y*gradient_vector.y);
      }
      return vector_buffer;

    }

    fvMappedField<Tensor2, on_cell>& limiterCartesian::operator() (fvMappedField<Vector2, on_cell>& phi, fvMappedField<Tensor2, on_cell>& gradient){
      tensor2_buffer.initialize_by_field(gradient);
      auto phi_begin = phi.data_begin();
      auto phi_boundary_value_begin = phi.boundary_value_begin();
      auto phi_gradient_begin = gradient.data_begin();
      auto faces_begin = phi.mesh.Cell.HalfFacets.begin();
      auto face_normal_begin = phi.mesh.HalfFacet.Normal.begin();
      auto face_centroid_begin = phi.mesh.HalfFacet.Centroid.begin();
      auto limited_gradient_begin = tensor2_buffer.data_begin();
      auto cell_neighbour_begin = phi.mesh.Cell.Neighbours.begin();
      auto cell_centroid_begin = phi.mesh.Cell.Centroid.begin();

      //initialize the limited gradients
      for(auto limited_gradient_iter = limited_gradient_begin; limited_gradient_iter < tensor2_buffer.data_end(); ++limited_gradient_iter){
        *(limited_gradient_iter) = Tensor2(0.0, 0.0, 0.0, 0.0);
      }

      //iterate all the cells
      for(auto cell_iter = phi.mesh.Cell.Neighbours.begin();
               cell_iter < phi.mesh.Cell.Neighbours.end();
               ++cell_iter){
        Flag id_this = cell_iter - cell_neighbour_begin;
        Vector2 phi_this = *(phi_begin + id_this);
        Vector2 centroid_this = *(cell_centroid_begin + id_this);
        Tensor2 gradient_vector = *(phi_gradient_begin + id_this);
        Scalar limit_ratio_xx = 1.0;
        Scalar limit_ratio_xy = 1.0;
        Scalar _limit_ratio_xx = 1.0;
        Scalar _limit_ratio_xy = 1.0;
        Scalar limit_ratio_yx = 1.0;
        Scalar limit_ratio_yy = 1.0;
        Scalar _limit_ratio_yx = 1.0;
        Scalar _limit_ratio_yy = 1.0;

        Flag neib_cnt = 0;
        for(auto neib_iter : *cell_iter){
          if (neib_cnt == 1 || neib_cnt == 3){
            if (!neib_iter.is_boundary()){ //neighbour is not a boundary
              Flag id_neib = neib_iter.get_global_id();
              Vector2 phi_neib = *(phi_begin + id_neib);
              Vector2 centroid_neib = *(cell_centroid_begin + id_neib);
              Vector2 direction_vector = centroid_neib - centroid_this;
              Vector2 downwind = dot(gradient_vector, direction_vector);
              Vector2 upwind = phi_neib - phi_this;
              _limit_ratio_xx = limiterFunction(upwind.x, downwind.x);
              _limit_ratio_yx = limiterFunction(upwind.y, downwind.y);
            }else{
              Flag id_face = (*(faces_begin + id_this))[neib_cnt];
              Vector normal = uni(*(face_normal_begin + id_face));
              Flag id_boundary = neib_iter.get_global_id();
              Vector2 phi_bound = *(phi_boundary_value_begin + id_boundary);
              Vector2 centroid_boundary = *(face_centroid_begin + id_face);
              Vector2 direction_vector = centroid_boundary - centroid_this;
              Vector2 downwind = 2.0*dot(gradient_vector, direction_vector);
              Vector2 upwind = phi_bound - phi_this;
              _limit_ratio_xx = limiterFunction(upwind.x, downwind.x);
              _limit_ratio_yx = limiterFunction(upwind.y, downwind.y);
            }
            limit_ratio_xx = std::min(limit_ratio_xx, _limit_ratio_xx);
            limit_ratio_yx = std::min(limit_ratio_yx, _limit_ratio_yx);
          }else{
            if (!neib_iter.is_boundary()){ //neighbour is not a boundary
              Flag id_neib = neib_iter.get_global_id();
              Vector2 phi_neib = *(phi_begin + id_neib);
              Vector2 centroid_neib = *(cell_centroid_begin + id_neib);
              Vector2 direction_vector = centroid_neib - centroid_this;
              Vector2 downwind = dot(gradient_vector, direction_vector);
              Vector2 upwind = phi_neib - phi_this;
              _limit_ratio_xy = limiterFunction(upwind.x, downwind.x);
              _limit_ratio_yy = limiterFunction(upwind.y, downwind.y);
            }else{
              Flag id_face = (*(faces_begin + id_this))[neib_cnt];
              Vector normal = uni(*(face_normal_begin + id_face));
              Flag id_boundary = neib_iter.get_global_id();
              Vector2 phi_bound = *(phi_boundary_value_begin + id_boundary);
              Vector2 centroid_boundary = *(face_centroid_begin + id_face);
              Vector2 direction_vector = centroid_boundary - centroid_this;
              Vector2 downwind = 2.0*dot(gradient_vector, direction_vector);
              Vector2 upwind = phi_bound - phi_this;
              _limit_ratio_xy = limiterFunction(upwind.x, downwind.x);
              _limit_ratio_yy = limiterFunction(upwind.y, downwind.y);
            }
            limit_ratio_xy = std::min(limit_ratio_xy, _limit_ratio_xy);
            limit_ratio_yy = std::min(limit_ratio_yy, _limit_ratio_yy);
          }
          neib_cnt++;
        }
        *(limited_gradient_begin + id_this) = Tensor2(limit_ratio_xx*gradient_vector.xx,
                                                      limit_ratio_xy*gradient_vector.xy,
                                                      limit_ratio_yx*gradient_vector.yx,
                                                      limit_ratio_yy*gradient_vector.yy);
      }
      return tensor2_buffer;
    }

    void LimiterCartesian(fvMappedField<Scalar, on_cell>& phi, fvMappedField<Vector, on_cell>& gradient){
      auto phi_begin = phi.data_begin();
      auto phi_boundary_value_begin = phi.boundary_value_begin();
      auto phi_gradient_begin = gradient.data_begin();
      auto faces_begin = phi.mesh.Cell.HalfFacets.begin();
      auto face_normal_begin = phi.mesh.HalfFacet.Normal.begin();
      auto face_centroid_begin = phi.mesh.HalfFacet.Centroid.begin();
      auto cell_neighbour_begin = phi.mesh.Cell.Neighbours.begin();
      auto cell_centroid_begin = phi.mesh.Cell.Centroid.begin();

      //iterate all the cells
#pragma omp parallel for firstprivate (phi_begin,\
                                       phi_boundary_value_begin,\
                                       phi_gradient_begin,\
                                       faces_begin,\
                                       face_normal_begin,\
                                       face_centroid_begin,\
                                       cell_neighbour_begin,\
                                       cell_centroid_begin)

      for (auto cell_iter = phi.mesh.Cell.Neighbours.begin();
        cell_iter < phi.mesh.Cell.Neighbours.end();
        ++cell_iter){
        Flag id_this = cell_iter - cell_neighbour_begin;
        Scalar phi_this = *(phi_begin + id_this);
        Vector centroid_this = *(cell_centroid_begin + id_this);
        Vector gradient_vector = *(phi_gradient_begin + id_this);
        Scalar limit_ratio_x = 1.0;
        Scalar limit_ratio_y = 1.0;
        Scalar _limit_ratio_x = 1.0;
        Scalar _limit_ratio_y = 1.0;
        Flag neib_cnt = 0;
        for (auto neib_iter : *cell_iter){
          if (neib_cnt == 1 || neib_cnt == 3){
            if (!neib_iter.is_boundary()){ //neighbour is not a boundary
              Flag id_neib = neib_iter.get_global_id();
              Scalar phi_neib = *(phi_begin + id_neib);
              Vector centroid_neib = *(cell_centroid_begin + id_neib);
              Vector direction_vector = centroid_neib - centroid_this;
              Scalar downwind = dot(gradient_vector, direction_vector);
              Scalar upwind = phi_neib - phi_this;
              _limit_ratio_x = minmod(upwind, downwind);
            }
            else{
              Flag id_face = (*(faces_begin + id_this))[neib_cnt];
              Vector normal = uni(*(face_normal_begin + id_face));
              Flag id_boundary = neib_iter.get_global_id();
              Scalar phi_bound = *(phi_boundary_value_begin + id_boundary);
              Vector centroid_boundary = *(face_centroid_begin + id_face);
              Vector direction_vector = centroid_boundary - centroid_this;
              Scalar downwind = 2.0*dot(gradient_vector, direction_vector);
              Scalar upwind = phi_bound - phi_this;
              _limit_ratio_x = minmod(upwind, downwind);
            }
            limit_ratio_x = std::min(limit_ratio_x, _limit_ratio_x);
          }
          else{
            if (!neib_iter.is_boundary()){ //neighbour is not a boundary
              Flag id_neib = neib_iter.get_global_id();
              Scalar phi_neib = *(phi_begin + id_neib);
              Vector centroid_neib = *(cell_centroid_begin + id_neib);
              Vector direction_vector = centroid_neib - centroid_this;
              Scalar downwind = dot(gradient_vector, direction_vector);
              Scalar upwind = phi_neib - phi_this;
              _limit_ratio_y = minmod(upwind, downwind);
            }
            else{
              Flag id_face = (*(faces_begin + id_this))[neib_cnt];
              Vector normal = uni(*(face_normal_begin + id_face));
              Flag id_boundary = neib_iter.get_global_id();
              Scalar phi_bound = *(phi_boundary_value_begin + id_boundary);
              Vector centroid_boundary = *(face_centroid_begin + id_face);
              Vector direction_vector = centroid_boundary - centroid_this;
              Scalar downwind = 2.0*dot(gradient_vector, direction_vector);
              Scalar upwind = phi_bound - phi_this;
              _limit_ratio_y = minmod(upwind, downwind);
            }
            limit_ratio_y = std::min(limit_ratio_y, _limit_ratio_y);
          }
          neib_cnt++;
        }
        *(phi_gradient_begin + id_this) = Vector(limit_ratio_x*gradient_vector.x, limit_ratio_y*gradient_vector.y);
      }
    }

    void LimiterCartesian(fvMappedField<Vector2, on_cell>& phi, fvMappedField<Tensor2, on_cell>& gradient){
      auto phi_begin = phi.data_begin();
      auto phi_boundary_value_begin = phi.boundary_value_begin();
      auto phi_gradient_begin = gradient.data_begin();
      auto faces_begin = phi.mesh.Cell.HalfFacets.begin();
      auto face_normal_begin = phi.mesh.HalfFacet.Normal.begin();
      auto face_centroid_begin = phi.mesh.HalfFacet.Centroid.begin();
      auto cell_neighbour_begin = phi.mesh.Cell.Neighbours.begin();
      auto cell_centroid_begin = phi.mesh.Cell.Centroid.begin();

      //iterate all the cells
#pragma omp parallel for firstprivate (phi_begin,\
                                       phi_boundary_value_begin,\
                                       phi_gradient_begin,\
                                       faces_begin,\
                                       face_normal_begin,\
                                       face_centroid_begin,\
                                       cell_neighbour_begin,\
                                       cell_centroid_begin)
      for (auto cell_iter = phi.mesh.Cell.Neighbours.begin();
        cell_iter < phi.mesh.Cell.Neighbours.end();
        ++cell_iter){
        Flag id_this = cell_iter - cell_neighbour_begin;
        Vector2 phi_this = *(phi_begin + id_this);
        Vector2 centroid_this = *(cell_centroid_begin + id_this);
        Tensor2 gradient_vector = *(phi_gradient_begin + id_this);
        Scalar limit_ratio_xx = 1.0;
        Scalar limit_ratio_xy = 1.0;
        Scalar _limit_ratio_xx = 1.0;
        Scalar _limit_ratio_xy = 1.0;
        Scalar limit_ratio_yx = 1.0;
        Scalar limit_ratio_yy = 1.0;
        Scalar _limit_ratio_yx = 1.0;
        Scalar _limit_ratio_yy = 1.0;

        Flag neib_cnt = 0;
        for (auto neib_iter : *cell_iter){
          if (neib_cnt == 1 || neib_cnt == 3){
            if (!neib_iter.is_boundary()){ //neighbour is not a boundary
              Flag id_neib = neib_iter.get_global_id();
              Vector2 phi_neib = *(phi_begin + id_neib);
              Vector2 centroid_neib = *(cell_centroid_begin + id_neib);
              Vector2 direction_vector = centroid_neib - centroid_this;
              Vector2 downwind = dot(gradient_vector, direction_vector);
              Vector2 upwind = phi_neib - phi_this;
              _limit_ratio_xx = minmod(upwind.x, downwind.x);
              _limit_ratio_yx = minmod(upwind.y, downwind.y);
            }
            else{
              Flag id_face = (*(faces_begin + id_this))[neib_cnt];
              Vector normal = uni(*(face_normal_begin + id_face));
              Flag id_boundary = neib_iter.get_global_id();
              Vector2 phi_bound = *(phi_boundary_value_begin + id_boundary);
              Vector2 centroid_boundary = *(face_centroid_begin + id_face);
              Vector2 direction_vector = centroid_boundary - centroid_this;
              Vector2 downwind = 2.0*dot(gradient_vector, direction_vector);
              Vector2 upwind = phi_bound - phi_this;
              _limit_ratio_xx = minmod(upwind.x, downwind.x);
              _limit_ratio_yx = minmod(upwind.y, downwind.y);
            }
            limit_ratio_xx = std::min(limit_ratio_xx, _limit_ratio_xx);
            limit_ratio_yx = std::min(limit_ratio_yx, _limit_ratio_yx);
          }
          else{
            if (!neib_iter.is_boundary()){ //neighbour is not a boundary
              Flag id_neib = neib_iter.get_global_id();
              Vector2 phi_neib = *(phi_begin + id_neib);
              Vector2 centroid_neib = *(cell_centroid_begin + id_neib);
              Vector2 direction_vector = centroid_neib - centroid_this;
              Vector2 downwind = dot(gradient_vector, direction_vector);
              Vector2 upwind = phi_neib - phi_this;
              _limit_ratio_xy = minmod(upwind.x, downwind.x);
              _limit_ratio_yy = minmod(upwind.y, downwind.y);
            }
            else{
              Flag id_face = (*(faces_begin + id_this))[neib_cnt];
              Vector normal = uni(*(face_normal_begin + id_face));
              Flag id_boundary = neib_iter.get_global_id();
              Vector2 phi_bound = *(phi_boundary_value_begin + id_boundary);
              Vector2 centroid_boundary = *(face_centroid_begin + id_face);
              Vector2 direction_vector = centroid_boundary - centroid_this;
              Vector2 downwind = 2.0*dot(gradient_vector, direction_vector);
              Vector2 upwind = phi_bound - phi_this;
              _limit_ratio_xy = minmod(upwind.x, downwind.x);
              _limit_ratio_yy = minmod(upwind.y, downwind.y);
            }
            limit_ratio_xy = std::min(limit_ratio_xy, _limit_ratio_xy);
            limit_ratio_yy = std::min(limit_ratio_yy, _limit_ratio_yy);
          }
          neib_cnt++;
        }
        *(phi_gradient_begin + id_this) = Tensor2(limit_ratio_xx*gradient_vector.xx,
          limit_ratio_xy*gradient_vector.xy,
          limit_ratio_yx*gradient_vector.yx,
          limit_ratio_yy*gradient_vector.yy);
      }
    }

  }//--end of namespace fv-------

}//--end of namespace GC---------
