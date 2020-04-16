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
 \file mapped_field.h
 \brief Header file for mapped field class 

 \version 0.1
 \author xilin xia
*/ 



#ifndef MAPPED_FIELD_H
#define MAPPED_FIELD_H

#include <vector>
#include <map>
#include <memory>
#include <iostream>
#include <utility>
#include <algorithm>
#include "Flag.h"
#include "mesh_interface.h"
#include "field_reader.h"
#include "boundary_funcs.h"

namespace GC{

  enum COPY_MODES {full, partial};
  enum MAPPING_MODES{on_vertex, on_halffacet, on_cell};
  
  template <typename T_A, typename T_b, MAPPING_MODES C> class diagonalMatrix;

  ///This class implements field data mapped on finite volume mesh 
  template <typename T, MAPPING_MODES C> class fvMappedField{
    public:
      fvMappedField() = default;

      fvMappedField(const fvMeshQueries& _mesh) :mesh(_mesh), current_t(0.0){
          switch (C){
          case on_cell: //on cell
            data.resize(mesh.Cell.size());
            boundary_type.resize(mesh.Boundary.size()); ///to do: boundary type information on either cell or half facet may not be needed
            boundary_value.resize(mesh.Boundary.size());
            break;
          case on_halffacet: //on halffacet
            data.resize(mesh.HalfFacet.size());
            boundary_type.resize(mesh.Boundary.size());
            boundary_value.resize(mesh.Boundary.size());
            break;
          case on_vertex: //on vertex
            data.resize(mesh.Vertex.size());
            boundary_type.resize(mesh.Boundary.size());
            boundary_value.resize(mesh.Boundary.size());
            break;
          }
      }

      fvMappedField(const fvMeshQueries& _mesh, 
        const fieldReader& field_reader) :mesh(_mesh), current_t(0.0){
        switch (C){
        case on_cell: //on cell
          data.resize(mesh.Cell.size());
          boundary_type.resize(mesh.Boundary.size()); ///to do: boundary type information on either cell or half facet may not be needed
          boundary_value.resize(mesh.Boundary.size());
          break;
        case on_halffacet: //on halffacet
          data.resize(mesh.HalfFacet.size());
          boundary_type.resize(mesh.Boundary.size());
          boundary_value.resize(mesh.Boundary.size());
          break;
        case on_vertex: //on vertex
          data.resize(mesh.Vertex.size());
          boundary_type.resize(mesh.Boundary.size());
          boundary_value.resize(mesh.Boundary.size());
          break;
        }
        Flag element_cnt = 0;
        for (auto data_iter : field_reader.data){
          Flag id_element = data_iter.first;
          if (id_element <= data.size()){
            data[id_element] = data_iter.second;
            element_cnt++;
          }
        }
        Flag boundary_cnt = 0;
        for (auto boundary_type_iter : field_reader.boundary_type){
          Flag id_element = boundary_type_iter.first;
          for (auto boundary_face : *(mesh.Cell.Neighbours.begin() + id_element)){
            if (boundary_face.is_boundary()){
              Flag id_boundary = boundary_face.get_global_id();
              boundary_type[id_boundary] = boundary_type_iter.second;
              boundary_cnt++;
            }
          }
        }
        if (element_cnt != data.size() || boundary_cnt != boundary_type.size()){
          std::cout << "Warning: Incomplete initializing." << std::endl;
        }
        update_boundary_values();
      }

      fvMappedField(const fvMeshQueries& _mesh,
        const completeFieldReader& field_reader) :mesh(_mesh), current_t(0.0){
        Flag region_mask_size = field_reader.region_mask.size();
        switch (C){
        case on_cell: //on cell
          data.resize(mesh.Cell.size());
          boundary_type.resize(mesh.Boundary.size()); ///to do: boundary type information on either cell or half facet may not be needed
          boundary_value.resize(mesh.Boundary.size());
          if (region_mask_size > 0){
            region_mask.resize(mesh.Cell.size());
          }
          break;
        case on_halffacet: //on halffacet
          data.resize(mesh.HalfFacet.size());
          boundary_type.resize(mesh.Boundary.size());
          boundary_value.resize(mesh.Boundary.size());
          if (region_mask_size > 0){
            region_mask.resize(mesh.HalfFacet.size());
          }
          break;
        case on_vertex: //on vertex
          data.resize(mesh.Vertex.size());
          boundary_type.resize(mesh.Boundary.size());
          boundary_value.resize(mesh.Boundary.size());
          if (region_mask_size > 0){
            region_mask.resize(mesh.Vertex.size());
          }
          break;
        }

        Flag element_cnt = 0;
        for (auto data_iter : field_reader.data){
          Flag id_element = data_iter.first;
          if (id_element <= data.size()){
            data[id_element] = data_iter.second;
            element_cnt++;
          }
        }
        Flag boundary_cnt = 0;
        for (auto boundary_type_iter : field_reader.boundary_type){
          Flag id_element = boundary_type_iter.first;
          for (auto boundary_face : *(mesh.Cell.Neighbours.begin() + id_element)){
            if (boundary_face.is_boundary()){
              Flag id_boundary = boundary_face.get_global_id();
              boundary_type[id_boundary] = boundary_type_iter.second;
              boundary_cnt++;
            }
          }
        }

        if (element_cnt != data.size() || boundary_cnt != boundary_type.size()){
          std::cout << "Warning: Incomplete initializing." << std::endl;
        }

        element_cnt = 0;
        for (auto region_mask_iter : field_reader.region_mask){
          Flag id_element = region_mask_iter.first;
          if (id_element <= region_mask.size()){
            region_mask[id_element] = region_mask_iter.second;
            element_cnt++;
          }
        }

        Flag boundary_source_cnt = field_reader.time_series.size();
        time_series.resize(boundary_source_cnt);
        boundary_source.resize(boundary_source_cnt);
        for (unsigned int i = 0; i < boundary_source_cnt; ++i){
          Flag source_value_cnt = field_reader.time_series[i].size();
          time_series[i].resize(source_value_cnt);
          boundary_source[i].resize(source_value_cnt);
          for (unsigned int j = 0; j < source_value_cnt; ++j){
            time_series[i][j] = field_reader.time_series[i][j];
            boundary_source[i][j] = field_reader.boundary_source[i][j];
          }
        }
        update_boundary_values();

        Flag data_source_cnt = field_reader.data_time_series.size();
        data_source_current.resize(data_source_cnt);
        data_time_series.resize(data_source_cnt);
        data_source.resize(data_source_cnt);
        for (unsigned int i = 0; i < data_source_cnt; ++i){
          Flag source_value_cnt = field_reader.data_time_series[i].size();
          data_time_series[i].resize(source_value_cnt);
          data_source[i].resize(source_value_cnt);
          for (unsigned int j = 0; j < source_value_cnt; ++j){
            data_time_series[i][j] = field_reader.data_time_series[i][j];
            data_source[i][j] = field_reader.data_source[i][j];
          }
        }
        update_data_values();
      }

      ///copy constructor, it has two modes: full and partial, default value is full
      fvMappedField(const fvMappedField& other, COPY_MODES c_m = full)
        :mesh(other.mesh), current_t(other.current_t){
           if(c_m == full){
             data = other.data;  //only copy data when copy mode is full
             boundary_type = other.boundary_type;
             boundary_value = other.boundary_value;
             time_series = other.time_series;
             boundary_source = other.boundary_source;
             data_time_series = other.data_time_series;
             data_source = other.data_source;
           }else{
             data = std::vector<T>(other.data.size()); //when copying mode is partial, all members are empty
             boundary_type = std::vector<ShortTripleFlag>(other.boundary_type.size()); //when copying mode is partial, all members are empty
             boundary_value = std::vector<T>(other.boundary_value.size()); //when copying mode is partial, all members are empty
           }
         }
      ///move constructor
      fvMappedField(const fvMappedField&& other) :
        mesh(other.mesh),
        current_t(other.current_t),
        data(std::move(other.data)),
        boundary_type(std::move(other.boundary_type)),
        boundary_value(std::move(other.boundary_value)),
        time_series(std::move(other.time_series)),
        boundary_source(std::move(other.boundary_source)),
        data_time_series(std::move(other.data_time_series)),
        data_source(std::move(other.data_source))
        {}
      ///assignment operator
      fvMappedField& operator= (fvMappedField& rhs){
        if(this != &rhs){
          mesh = rhs.mesh;
          current_t = rhs.current_t;
          data = rhs.data;
          boundary_type = rhs.boundary_type;
          boundary_value = rhs.boundary_value;
          time_series = rhs.time_series;
          boundary_source = rhs.boundary_source;
          data_time_series = rhs.data_time_series;
          data_source = rhs.data_source;
        }
        return *this;
      }
      ///move assignment operator
      fvMappedField& operator= (fvMappedField&& rhs){
        if(this != &rhs){
          mesh = rhs.mesh;
          current_t = rhs.current_t;
          data = std::move(rhs.data);
          boundary_type = std::move(rhs.boundary_type);
          boundary_value = std::move(rhs.boundary_value);
          time_series = std::move(rhs.time_series);
          boundary_source = std::move(rhs.boundary_source);
          data_time_series = std::move(rhs.data_time_series);
          data_source = std::move(rhs.data_source);
        }
        return *this;
      }

    public:
      ///iterator for data
      typedef typename std::vector<T>::iterator data_iterator;
      typedef typename std::vector<T>::const_iterator data_const_iterator;
      data_iterator data_begin(){return data.begin();}
      data_iterator data_end(){return data.end();}

      ///iterator for boundary value
      typedef typename std::vector<T>::iterator boundary_value_iterator;
      typedef typename std::vector<T>::const_iterator boundary_value_const_iterator;
      boundary_value_iterator boundary_value_begin(){ return boundary_value.begin(); }
      boundary_value_iterator boundary_value_end(){ return boundary_value.end(); }

      ///iterator for boundary type
      typedef std::vector<ShortTripleFlag>::iterator boundary_type_iterator;
      typedef std::vector<ShortTripleFlag>::const_iterator boundary_type_const_iterator;
      boundary_type_iterator boundary_type_begin(){return boundary_type.begin();}
      boundary_type_iterator boundary_type_end(){return boundary_type.end();}


      ///iterator for time series
      typedef std::vector< std::vector<Scalar> >::iterator time_series_iterator;
      typedef std::vector< std::vector<Scalar> >::const_iterator time_series_const_iterator;
      time_series_iterator time_series_begin(){ return time_series.begin(); }
      time_series_iterator time_series_end(){ return time_series.end(); }

      ///iterator for boundary source
      typedef typename std::vector< std::vector<T> >::iterator boundary_source_iterator;
      typedef typename std::vector< std::vector<T> >::const_iterator boundary_source_const_iterator;
      boundary_source_iterator boundary_source_begin(){ return boundary_source.begin(); }
      boundary_source_iterator boundary_source_end(){ return boundary_source.end(); }

      ///iterator for region mask
      typedef typename std::vector<Flag>::iterator region_mask_iterator;
      typedef typename std::vector<Flag>::const_iterator region_mask_const_iterator;
      region_mask_iterator region_mask_begin(){ return region_mask.begin(); }
      region_mask_iterator region_mask_end(){ return region_mask.end(); }

      ///iterator for data time series
      typedef std::vector< std::vector<Scalar> >::iterator data_time_series_iterator;
      typedef std::vector< std::vector<Scalar> >::const_iterator data_time_series_const_iterator;
      data_time_series_iterator data_time_series_begin(){ return data_time_series.begin(); }
      data_time_series_iterator data_time_series_end(){ return data_time_series.end(); }

      ///iterator for data source
      typedef typename std::vector< std::vector<T> >::iterator data_source_iterator;
      typedef typename std::vector< std::vector<T> >::const_iterator data_source_const_iterator;
      boundary_source_iterator data_source_begin(){ return data_source.begin(); }
      boundary_source_iterator data_source_end(){ return data_source.end(); }

      
      ///initialize from other field, only copy mesh and boundary type
      template <typename T_phi>
      void initialize_by_field(fvMappedField<T_phi, C>& phi){
        mesh = phi.mesh;
        data.resize(mesh.Cell.size());      
        boundary_type.resize(mesh.Boundary.size());
        boundary_value.resize(mesh.Boundary.size());
      }

      ///set by a single value
      void set_single_value(T value){
        for (auto data_iter = data.begin(); data_iter < data.end(); ++ data_iter){
          *data_iter = value;
        }
        for (auto boundary_value_iter = boundary_value.begin(); boundary_value_iter < boundary_value.end(); ++boundary_value_iter){
          *boundary_value_iter = value;
        }
      }

      ///get the boundary value
      T get_boundary(unsigned int index){ return boundary_value[index]; }

      ///update the boundary values
      void update_boundary_values(){
        auto boundary2opposite_begin = mesh.Boundary.Opposite.begin();
        auto cell_halffacets_begin = mesh.Cell.HalfFacets.begin();
        auto halffacet_normal_begin = mesh.HalfFacet.Normal.begin();
        for (auto boundary_type_iter = boundary_type.begin(); boundary_type_iter < boundary_type.end(); ++boundary_type_iter){
          Flag index = boundary_type_iter - boundary_type.begin();
          auto opposite_cell = *(boundary2opposite_begin + index);
          Flag cell_id = opposite_cell.get_global_id();
          Flag halffacet_cnt = opposite_cell.get_local_id();
          Flag halffacet_id = (*(cell_halffacets_begin + cell_id))[halffacet_cnt];
          auto in_value = data[cell_id];
          auto normal = *(halffacet_normal_begin + halffacet_id);
          auto _boundary_type = *boundary_type_iter;
          Flag primary_type = _boundary_type.getx();
          Flag secondary_type = _boundary_type.gety();
          Flag source_id = _boundary_type.getz();
          if (0 == primary_type){ //cannot evolve
            std::cout << "Type is 0, can not be updated!" << std::endl;
          }
          else if (1 == primary_type){ //fixed gradient 

          }
          else if (2 == primary_type){ //zero gradient or wall
            switch(secondary_type){
            case 0:
              boundary_value[index] = ZeroGradientScalar(in_value, normal);
              break;
            case 1:
              boundary_value[index] = ZeroGradientVector(in_value, normal);
              break;
            case 2:
              boundary_value[index] = WallNonSlip(in_value, normal);
              break;
            case 3:
              boundary_value[index] = WallSlip(in_value, normal);
              break;
            } 
          }
          else if (3 == primary_type){ //fixed value
            switch (secondary_type){
            case 0:
              auto time_second_iter = std::lower_bound(time_series[source_id].begin(), time_series[source_id].end(), current_t);
              auto time_first_iter = time_second_iter;
              if (time_first_iter == time_series[source_id].begin()){
                Flag id = time_first_iter - time_series[source_id].begin();
                auto source_value_first = boundary_source[source_id][id];
                auto source_value_second = boundary_source[source_id][id + 1];
                boundary_value[index] = source_value_first;
              }
              else{
                time_first_iter = time_second_iter - 1;
                Flag id = time_first_iter - time_series[source_id].begin();
                auto source_value_first = boundary_source[source_id][id];
                auto source_value_second = boundary_source[source_id][id + 1];
                if (current_t <= *time_second_iter){
                  boundary_value[index] = source_value_first + (source_value_second - source_value_first)*(current_t - *time_first_iter) / (*time_second_iter - *time_first_iter);
                }
                else{
                  boundary_value[index] = source_value_second;
                }
              }
              break;
            }

          }
          else if (4 == primary_type){//calculated

          }
        }

      }

      ///update the boundary values
      void update_data_values(){
        for (auto data_source_iter = data_source_current.begin(); data_source_iter < data_source_current.end(); ++data_source_iter){
          Flag source_id = data_source_iter - data_source_current.begin();
          auto time_second_iter = std::lower_bound(data_time_series[source_id].begin(), data_time_series[source_id].end(), current_t);
          auto time_first_iter = time_second_iter;
          if (time_first_iter == data_time_series[source_id].begin()){
            Flag id = time_first_iter - data_time_series[source_id].begin();
            auto source_value_first = data_source[source_id][id];
            auto source_value_second = data_source[source_id][id + 1];
            *data_source_iter = source_value_first;
          }
          else{
            time_first_iter = time_second_iter - 1;
            Flag id = time_first_iter - data_time_series[source_id].begin();
            auto source_value_first = data_source[source_id][id];
            auto source_value_second = data_source[source_id][id + 1];
            if (current_t <= *time_second_iter){
              *data_source_iter = source_value_first + (source_value_second - source_value_first)*(current_t - *time_first_iter) / (*time_second_iter - *time_first_iter);
            }
            else{
              *data_source_iter = source_value_second;
            }
          }

        }

        for (auto region_mask_iter = region_mask.begin(); region_mask_iter < region_mask.end(); ++region_mask_iter){
          Flag index = region_mask_iter - region_mask.begin();
          Flag source_id = *region_mask_iter;
          data[index] = data_source_current[source_id];
        }
      }

      ///update time
      void update_time(const Scalar& t){
        current_t = t;
      }
    public:  
      fvMeshQueries mesh;
    private:
      std::vector<ShortTripleFlag> boundary_type;  
      std::vector<T> boundary_value;
      std::vector<T> data;
      std::vector< std::vector<Scalar> > time_series;
      std::vector< std::vector<T> > boundary_source;
      std::vector<Flag> region_mask;
      std::vector< std::vector<Scalar> > data_time_series;
      std::vector< std::vector<T> > data_source;
      std::vector<T> data_source_current;
    private:
      Scalar current_t;
    private:
  };

  typedef fvMappedField<Scalar, on_vertex> fvScalarFieldOnVertex;
  typedef fvMappedField<Scalar, on_halffacet> fvScalarFieldOnHalffacet;
  typedef fvMappedField<Scalar, on_cell> fvScalarFieldOnCell;
  typedef fvMappedField<Vector, on_vertex> fvVectorFieldOnVertex;
  typedef fvMappedField<Vector, on_halffacet> fvVectorFieldOnHalffacet;
  typedef fvMappedField<Vector, on_cell> fvVectorFieldOnCell;
  typedef fvMappedField<Tensor, on_vertex> fvTensorFieldOnVertex;
  typedef fvMappedField<Tensor, on_halffacet> fvTensorFieldOnHalffacet;
  typedef fvMappedField<Tensor, on_cell> fvTensorFieldOnCell;

  template <typename T, MAPPING_MODES C>
  fvMappedField<T, C> operator+ (fvMappedField<T, C>& lhs, fvMappedField<T, C>& rhs){
    fvMappedField<T, C> result(lhs);
    auto rhs_begin = rhs.data_begin();
    auto result_begin = result.data_begin();
    for (auto result_iter = result.data_begin(); result_iter < result.data_end(); ++result_iter){
      auto n = result_iter - result_begin;
      *result_iter += *(rhs_begin + n);
    }
    return result;
  }

  template <typename T, MAPPING_MODES C>
  fvMappedField<T, C> operator+ (fvMappedField<T, C>&& lhs, fvMappedField<T, C>& rhs){
    fvMappedField<T, C> result(std::move(lhs));
    auto rhs_begin = rhs.data_begin();
    auto result_begin = result.data_begin();
    for (auto result_iter = result.data_begin(); result_iter < result.data_end(); ++result_iter){
      auto n = result_iter - result_begin;
      *result_iter += *(rhs_begin + n);
    }
    return result;
  }

  template <typename T, MAPPING_MODES C>
  fvMappedField<T, C> operator+ (fvMappedField<T, C>& lhs, fvMappedField<T, C>&& rhs){
    fvMappedField<T, C> result(std::move(rhs));
    auto lhs_begin = lhs.data_begin();
    auto result_begin = result.data_begin();
    for (auto result_iter = result.data_begin(); result_iter < result.data_end(); ++result_iter){
      auto n = result_iter - result_begin;
      *result_iter = *(lhs_begin + n) + (*result_iter);
    }
    return result;
  }

  template <typename T, MAPPING_MODES C>
  fvMappedField<T, C> operator+ (fvMappedField<T, C>&& lhs, fvMappedField<T, C>&& rhs){
    fvMappedField<T, C> result(std::move(lhs));
    auto rhs_begin = rhs.data_begin();
    auto result_begin = result.data_begin();
    for (auto result_iter = result.data_begin(); result_iter < result.data_end(); ++result_iter){
      auto n = result_iter - result_begin;
      *result_iter += *(rhs_begin + n);
    }
    return result;
  }

  template <typename T, MAPPING_MODES C>
  fvMappedField<T, C> operator- (fvMappedField<T, C>& lhs, fvMappedField<T, C>& rhs){
    fvMappedField<T, C> result(lhs);
    auto rhs_begin = rhs.data_begin();
    auto result_begin = result.data_begin();
    for (auto result_iter = result.data_begin(); result_iter < result.data_end(); ++result_iter){
      auto n = result_iter - result_begin;
      *result_iter -= *(rhs_begin + n);
    }
    return result;
  }

  template <typename T, MAPPING_MODES C>
  fvMappedField<T, C> operator- (fvMappedField<T, C>&& lhs, fvMappedField<T, C>& rhs){
    fvMappedField<T, C> result(std::move(lhs));
    auto rhs_begin = rhs.data_begin();
    auto result_begin = result.data_begin();
    for (auto result_iter = result.data_begin(); result_iter < result.data_end(); ++result_iter){
      auto n = result_iter - result_begin;
      *result_iter -= *(rhs_begin + n);
    }
    return result;
  }

  template <typename T, MAPPING_MODES C>
  fvMappedField<T, C> operator- (fvMappedField<T, C>& lhs, fvMappedField<T, C>&& rhs){
    fvMappedField<T, C> result(std::move(rhs));
    auto lhs_begin = lhs.data_begin();
    auto result_begin = result.data_begin();
    for (auto result_iter = result.data_begin(); result_iter < result.data_end(); ++result_iter){
      auto n = result_iter - result_begin;
      *result_iter = *(lhs_begin + n) - (*result_iter);
    }
    return result;
  }

  template <typename T, MAPPING_MODES C>
  fvMappedField<T, C> operator- (fvMappedField<T, C>&& lhs, fvMappedField<T, C>&& rhs){
    fvMappedField<T, C> result(std::move(lhs));
    auto rhs_begin = rhs.data_begin();
    auto result_begin = result.data_begin();
    for (auto result_iter = result.data_begin(); result_iter < result.data_end(); ++result_iter){
      auto n = result_iter - result_begin;
      *result_iter -= *(rhs_begin + n);
    }
    return result;
  }

  template <typename T, MAPPING_MODES C>
  fvMappedField<T, C> operator- (fvMappedField<T, C>& rhs){
    fvMappedField<T, C> result(rhs);
    auto rhs_begin = rhs.data_begin();
    auto result_begin = result.data_begin();
    for (auto result_iter = result.data_begin(); result_iter < result.data_end(); ++result_iter){
      auto n = result_iter - result_begin;
      *result_iter = -(*(rhs_begin + n));
    }
    return result;
  }

  template <typename T, MAPPING_MODES C>
  fvMappedField<T, C> operator- (fvMappedField<T, C>&& rhs){
    fvMappedField<T, C> result(std::move(rhs));
    auto rhs_begin = rhs.data_begin();
    auto result_begin = result.data_begin();
    for (auto result_iter = result.data_begin(); result_iter < result.data_end(); ++result_iter){
      auto n = result_iter - result_begin;
      *result_iter = -(*(rhs_begin + n));
    }
    return result;
  }
  

  template <typename T, MAPPING_MODES C>
  fvMappedField<T, C> operator* (fvMappedField<Scalar, C>& lhs, fvMappedField<T, C>& rhs){
    fvMappedField<T, C> result(rhs);
    auto lhs_begin = lhs.data_begin();
    auto result_begin = result.data_begin();
    for (auto result_iter = result.data_begin(); result_iter < result.data_end(); ++result_iter){
      auto n = result_iter - result_begin;
      *result_iter = *(lhs_begin + n)*(*result_iter);
    }
    return result;
  }

  template <typename T, MAPPING_MODES C>
  fvMappedField<T, C> operator* (fvMappedField<Scalar, C>& lhs, fvMappedField<T, C>&& rhs){
    fvMappedField<T, C> result(std::move(rhs));
    auto lhs_begin = lhs.data_begin();
    auto result_begin = result.data_begin();
    for (auto result_iter = result.data_begin(); result_iter < result.data_end(); ++result_iter){
      auto n = result_iter - result_begin;
      *result_iter = *(lhs_begin + n)*(*result_iter);
    }
    return result;
  }

}

#endif
