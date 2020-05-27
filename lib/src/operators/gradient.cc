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
\file gradient.cc
\brief Source file for ddt operator

\version 0.1
\author xilin xia

*/

#include "gradient.h"
#include "boundary.h"

namespace GC{

  namespace fv{

    fvMappedField<Vector, on_cell>& gradient::operator()(fvMappedField<Scalar, on_cell>& phi){
      vector_buffer.initialize_by_field(phi);
      auto phi_begin = phi.data_begin();
      auto phi_boundary_value_begin = phi.boundary_value_begin();
      auto phi_gradient_begin = vector_buffer.data_begin();
      auto faces_begin = phi.mesh.Cell.HalfFacets.begin();
      auto face_normal_begin = phi.mesh.HalfFacet.Normal.begin();
      auto face_area_begin = phi.mesh.HalfFacet.Area.begin();
      auto cell_neighbour_begin = phi.mesh.Cell.Neighbours.begin();
      auto cell_volume_begin = phi.mesh.Cell.Volume.begin();

      //initialize the gradients
      for(auto phi_gradient_iter = phi_gradient_begin; phi_gradient_iter < vector_buffer.data_end(); ++phi_gradient_iter){
        *(phi_gradient_iter) = Vector(0.0);
      }

      //iterate all the cells
      for(auto cell_iter = phi.mesh.Cell.Neighbours.begin();
               cell_iter < phi.mesh.Cell.Neighbours.end();
               ++cell_iter){
        Flag id_this = cell_iter - cell_neighbour_begin;
        Scalar phi_this = *(phi_begin + id_this);
        Scalar volume_this = *(cell_volume_begin + id_this);
        Flag neib_cnt = 0;
        for(auto neib_iter : *cell_iter){
          if (!neib_iter.is_boundary()){ //neighbour is not a boundary
            Flag id_neib = neib_iter.get_global_id();
            Flag id_face = (*(faces_begin + id_this))[neib_cnt];
            Vector normal = uni(*(face_normal_begin + id_face));
            Scalar area = *(face_area_begin + id_face);
            Scalar phi_neib = *(phi_begin + id_neib);
            *(phi_gradient_begin + id_this) += (phi_neib+phi_this)/2.0*area*normal/volume_this;
          }else{
            Flag id_face = (*(faces_begin + id_this))[neib_cnt];
            Vector normal = uni(*(face_normal_begin + id_face));
            Scalar area = *(face_area_begin + id_face);
            Flag id_boundary = neib_iter.get_global_id();
            Scalar phi_bound = *(phi_boundary_value_begin + id_boundary);
            *(phi_gradient_begin + id_this) += (phi_bound+phi_this)/2.0*area*normal/volume_this;
          }
          neib_cnt++;
        }
      }
      return vector_buffer;
    }

    fvMappedField<Tensor, on_cell>& gradient::operator()(fvMappedField<Vector, on_cell>& phi){
      tensor_buffer.initialize_by_field(phi);
      auto phi_begin = phi.data_begin();
      auto phi_boundary_type_begin = phi.boundary_type_begin();
      auto phi_gradient_begin = tensor_buffer.data_begin();
      auto faces_begin = phi.mesh.Cell.HalfFacets.begin();
      auto face_normal_begin = phi.mesh.HalfFacet.Normal.begin();
      auto face_area_begin = phi.mesh.HalfFacet.Area.begin();
      auto cell_neighbour_begin = phi.mesh.Cell.Neighbours.begin();
      auto cell_volume_begin = phi.mesh.Cell.Volume.begin();

      //initialize the gradients
      for(auto phi_gradient_iter = phi_gradient_begin; phi_gradient_iter < tensor_buffer.data_end(); ++phi_gradient_iter){
        *(phi_gradient_iter) = Tensor(0.0, 0.0, 0.0, 0.0);
      }

      //iterate all the cells
      for(auto cell_iter = phi.mesh.Cell.Neighbours.begin();
               cell_iter < phi.mesh.Cell.Neighbours.end();
               ++cell_iter){
        Flag id_this = cell_iter - cell_neighbour_begin;
        Vector phi_this = *(phi_begin + id_this);
        Scalar volume_this = *(cell_volume_begin + id_this);
        Flag neib_cnt = 0;
        for(auto neib_iter : *cell_iter){
          if (!neib_iter.is_boundary()){ //neighbour is not a boundary
            Flag id_neib = neib_iter.get_global_id();
            Flag id_face = (*(faces_begin + id_this))[neib_cnt];
            Vector normal = uni(*(face_normal_begin + id_face));
            Scalar area = *(face_area_begin + id_face);
            Vector phi_neib = *(phi_begin + id_neib);
            *(phi_gradient_begin + id_this) += product((phi_neib+phi_this)/2.0, normal)*area/volume_this;
          }else{
            Flag id_face = (*(faces_begin + id_this))[neib_cnt];
            Vector normal = uni(*(face_normal_begin + id_face));
            Scalar area = *(face_area_begin + id_face);
            Flag id_boundary = neib_iter.get_global_id();
            auto phi_boundary_type = *(phi_boundary_type_begin + id_boundary);
            Vector phi_bound = GetBoundary(phi_boundary_type, phi_this, normal);
            *(phi_gradient_begin + id_this) += product((phi_bound+phi_this)/2.0, normal)*area/volume_this;
          }
          neib_cnt++;
        }
      }
      return tensor_buffer;
    }

    void Gradient(fvMappedField<Scalar, on_cell>& phi, fvMappedField<Vector, on_cell>& phi_grad){

      auto phi_begin = phi.data_begin();
      auto phi_boundary_value_begin = phi.boundary_value_begin();
      auto phi_gradient_begin = phi_grad.data_begin();
      auto faces_begin = phi.mesh.Cell.HalfFacets.begin();
      auto face_normal_begin = phi.mesh.HalfFacet.Normal.begin();
      auto face_area_begin = phi.mesh.HalfFacet.Area.begin();
      auto cell_neighbour_begin = phi.mesh.Cell.Neighbours.begin();
      auto cell_volume_begin = phi.mesh.Cell.Volume.begin();

      //initialize the gradients
#pragma omp parallel for
      for (auto phi_gradient_iter = phi_grad.data_begin(); phi_gradient_iter < phi_grad.data_end(); ++phi_gradient_iter){
        *(phi_gradient_iter) = Vector(0.0);
      }

#pragma omp parallel for firstprivate (phi_begin,\
                                       phi_boundary_value_begin,\
                                       phi_gradient_begin,\
                                       faces_begin,\
                                       face_normal_begin,\
                                       face_area_begin,\
                                       cell_neighbour_begin,\
                                       cell_volume_begin)
      for (auto cell_iter = phi.mesh.Cell.Neighbours.begin();
        cell_iter < phi.mesh.Cell.Neighbours.end();
        ++cell_iter){
        Flag id_this = cell_iter - cell_neighbour_begin;
        Scalar phi_this = *(phi_begin + id_this);
        Scalar volume_this = *(cell_volume_begin + id_this);
        Flag neib_cnt = 0;
        for (auto neib_iter : *cell_iter){
          if (!neib_iter.is_boundary()){ //neighbour is not a boundary
            Flag id_neib = neib_iter.get_global_id();
            Flag id_face = (*(faces_begin + id_this))[neib_cnt];
            Vector normal = uni(*(face_normal_begin + id_face));
            Scalar area = *(face_area_begin + id_face);
            Scalar phi_neib = *(phi_begin + id_neib);
            *(phi_gradient_begin + id_this) += (phi_neib + phi_this) / 2.0*area*normal / volume_this;
          }
          else{
            Flag id_face = (*(faces_begin + id_this))[neib_cnt];
            Vector normal = uni(*(face_normal_begin + id_face));
            Scalar area = *(face_area_begin + id_face);
            Flag id_boundary = neib_iter.get_global_id();
            Scalar phi_bound = *(phi_boundary_value_begin + id_boundary);
            *(phi_gradient_begin + id_this) += (phi_bound + phi_this) / 2.0*area*normal / volume_this;
          }
          neib_cnt++;
        }
      }
    }

    void Gradient(fvMappedField<Vector, on_cell>& phi, fvMappedField<Tensor, on_cell>& phi_grad){

      auto phi_begin = phi.data_begin();
      auto phi_boundary_value_begin = phi.boundary_value_begin();
      auto phi_gradient_begin = phi_grad.data_begin();
      auto faces_begin = phi.mesh.Cell.HalfFacets.begin();
      auto face_normal_begin = phi.mesh.HalfFacet.Normal.begin();
      auto face_area_begin = phi.mesh.HalfFacet.Area.begin();
      auto cell_neighbour_begin = phi.mesh.Cell.Neighbours.begin();
      auto cell_volume_begin = phi.mesh.Cell.Volume.begin();

      //initialize the gradients
#pragma omp parallel for
      for (auto phi_gradient_iter = phi_grad.data_begin(); phi_gradient_iter < phi_grad.data_end(); ++phi_gradient_iter){
        *(phi_gradient_iter) = Tensor(0.0, 0.0, 0.0, 0.0);
      }

#pragma omp parallel for firstprivate (phi_begin,\
                                       phi_boundary_value_begin,\
                                       phi_gradient_begin,\
                                       faces_begin,\
                                       face_normal_begin,\
                                       face_area_begin,\
                                       cell_neighbour_begin,\
                                       cell_volume_begin)
      //iterate all the cells
      for (auto cell_iter = phi.mesh.Cell.Neighbours.begin();
        cell_iter < phi.mesh.Cell.Neighbours.end();
        ++cell_iter){
        Flag id_this = cell_iter - cell_neighbour_begin;
        Vector phi_this = *(phi_begin + id_this);
        Scalar volume_this = *(cell_volume_begin + id_this);
        Flag neib_cnt = 0;
        for (auto neib_iter : *cell_iter){
          if (!neib_iter.is_boundary()){ //neighbour is not a boundary
            Flag id_neib = neib_iter.get_global_id();
            Flag id_face = (*(faces_begin + id_this))[neib_cnt];
            Vector normal = uni(*(face_normal_begin + id_face));
            Scalar area = *(face_area_begin + id_face);
            Vector phi_neib = *(phi_begin + id_neib);
            *(phi_gradient_begin + id_this) += product((phi_neib + phi_this) / 2.0, normal)*area / volume_this;
          }
          else{
            Flag id_face = (*(faces_begin + id_this))[neib_cnt];
            Vector normal = uni(*(face_normal_begin + id_face));
            Scalar area = *(face_area_begin + id_face);
            Flag id_boundary = neib_iter.get_global_id();
            Vector phi_bound = *(phi_boundary_value_begin + id_boundary);
            *(phi_gradient_begin + id_this) += product((phi_bound + phi_this) / 2.0, normal)*area / volume_this;
          }
          neib_cnt++;
        }
      }
    }

  } //end of namespace fv
}//end of namespace GC
