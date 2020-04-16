// ======================================================================================
// Name                :    GeoClasses : Generic Geophysical Flow Modelling Framework
// Description         :    This code pack provides a generic framework for developing 
//                          Geophysical CFD software.
// ======================================================================================
// Version             :    0.1 
// Author              :    Xilin Xia (PhD candidate in Newcastle University)
// Create Time         :    2014/10/04
// Update Time         :    2015/10/29
// ======================================================================================
// Copyright @ Xilin Xia 2015 . All rights reserved.
// ======================================================================================

/*!
\flie field_algebra.h
\brief Header file for field algebra functions

\version 0.1
\author xilin xia
*/

#ifndef FIELD_ALGEBRA_H
#define FIELD_ALGEBRA_H

#include "mapped_field.h"

namespace GC{

  namespace fv{

    template <typename T_a, typename T_b, MAPPING_MODES C, typename F>
    void Unary(fvMappedField<T_a, C>& a, fvMappedField<T_b, C>& b, F func){
      auto a_data_begin = a.data_begin();
      auto b_data_begin = b.data_begin();

#pragma omp parallel for firstprivate (a_data_begin,\
                                       b_data_begin)

      for (auto a_data_iter = a.data_begin(); a_data_iter < a.data_end(); ++a_data_iter){
        Flag id = a_data_iter - a_data_begin;
        auto _a = *a_data_iter;
        *(b_data_begin + id) = func(_a);
      }
      auto a_boundary_value_begin = a.boundary_value_begin();
      auto b_boundary_value_begin = b.boundary_value_begin();

#pragma omp parallel for firstprivate (a_boundary_value_begin,\
                                       b_boundary_value_begin)

      for (auto a_boundary_value_iter = a.boundary_value_begin(); a_boundary_value_iter < a.boundary_value_end(); ++a_boundary_value_iter){
        Flag id = a_boundary_value_iter - a_boundary_value_begin;
        auto _a = *a_boundary_value_iter;
        *(b_boundary_value_begin + id) = func(_a);
      }

    }

    template <typename T_a, typename T_b, typename T_c, MAPPING_MODES C, typename F>
    void Binary(fvMappedField<T_a, C>& a, fvMappedField<T_b, C>& b, fvMappedField<T_c, C>& c, F func){
      auto a_data_begin = a.data_begin();
      auto b_data_begin = b.data_begin();
      auto c_data_begin = c.data_begin();

#pragma omp parallel for firstprivate (a_data_begin,\
                                       b_data_begin,\
                                       c_data_begin)

      for (auto a_data_iter = a.data_begin(); a_data_iter < a.data_end(); ++a_data_iter){
        Flag id = a_data_iter - a_data_begin;
        auto _a = *a_data_iter;
        auto _b = *(b_data_begin + id);
        *(c_data_begin + id) = func(_a, _b);
      }
      auto a_boundary_value_begin = a.boundary_value_begin();
      auto b_boundary_value_begin = b.boundary_value_begin();
      auto c_boundary_value_begin = c.boundary_value_begin();

#pragma omp parallel for firstprivate (a_boundary_value_begin,\
                                       b_boundary_value_begin,\
                                       c_boundary_value_begin)

      for (auto a_boundary_value_iter = a.boundary_value_begin(); a_boundary_value_iter < a.boundary_value_end(); ++a_boundary_value_iter){
        Flag id = a_boundary_value_iter - a_boundary_value_begin;
        auto _a = *a_boundary_value_iter;
        auto _b = *(b_boundary_value_begin + id);
        *(c_boundary_value_begin + id) = func(_a, _b);
      }


    }
  
  }//end of namespace fv

}//--end of namespace GC



#endif