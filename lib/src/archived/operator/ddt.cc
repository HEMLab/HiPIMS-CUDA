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
\file ddt.cc
\brief Source file for ddt operator

\version 0.1
\author xilin xia

*/

#include "ddt.h"

namespace GC{

  namespace fv{
    diagonalMatrix<Scalar, Scalar, on_cell>& ddt_optimized::operator() (Scalar dt, fvMappedField<Scalar, on_cell>& phi){
      scalar_matrix_buffer.initialize_by_field(phi);
      auto phi_begin = phi.data_begin();
      auto A_begin = scalar_matrix_buffer.A.begin();
#pragma omp parallel for firstprivate(phi_begin, A_begin)
      for (auto b_iter = scalar_matrix_buffer.b.begin(); b_iter < scalar_matrix_buffer.b.end(); ++b_iter){
        auto id_cell = b_iter - scalar_matrix_buffer.b.begin();
        Scalar A_value = 1.0 / dt;
        Scalar b_value = *(phi_begin + id_cell) * A_value;
        *(A_begin + id_cell) = A_value;
        *b_iter = b_value;
      }
      return scalar_matrix_buffer;
    }

    diagonalMatrix<Scalar, Vector2, on_cell>& ddt_optimized::operator() (Scalar dt, fvMappedField<Vector2, on_cell>& phi){
      vector2_matrix_buffer.initialize_by_field(phi);
      auto phi_begin = phi.data_begin();
      auto A_begin = vector2_matrix_buffer.A.begin();
#pragma omp parallel for firstprivate(phi_begin, A_begin)
      for (auto b_iter = vector2_matrix_buffer.b.begin(); b_iter < vector2_matrix_buffer.b.end(); ++b_iter){
        auto id_cell = b_iter - vector2_matrix_buffer.b.begin();
        Scalar A_value = 1.0 / dt;
        Vector2 b_value = *(phi_begin + id_cell) * A_value;
        *(A_begin + id_cell) = A_value;
        *b_iter = b_value;
      }
      return vector2_matrix_buffer;
    }

    diagonalMatrix<Scalar, Vector3, on_cell>& ddt_optimized::operator() (Scalar dt, fvMappedField<Vector3, on_cell>& phi){
      vector3_matrix_buffer.initialize_by_field(phi);
      auto phi_begin = phi.data_begin();
      auto A_begin = vector3_matrix_buffer.A.begin();
#pragma omp parallel for firstprivate(phi_begin, A_begin)
      for (auto b_iter = vector3_matrix_buffer.b.begin(); b_iter < vector3_matrix_buffer.b.end(); ++b_iter){
        auto id_cell = b_iter - vector3_matrix_buffer.b.begin();
        Scalar A_value = 1.0 / dt;
        Vector3 b_value = *(phi_begin + id_cell) * A_value;
        *(A_begin + id_cell) = A_value;
        *b_iter = b_value;
      }
      return vector3_matrix_buffer;
    }

  }

}
