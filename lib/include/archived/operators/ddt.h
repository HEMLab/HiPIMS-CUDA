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
\file ddt.h
\brief Header file for ddt operator

\version 0.1
\author xilin xia

*/

#ifndef DDT_H
#define DDT_H

//#include <omp.h>
#include "matrix.h"
#include "mapped_field.h"
#include "Scalar.h"
#include "Vector.h"
#include "Flag.h"

namespace GC{

  namespace fv{

    template <class T>
    diagonalMatrix<Scalar, T, on_cell> ddt(Scalar dt, fvMappedField<T, on_cell>& phi){
      diagonalMatrix<Scalar, T, on_cell> result(phi.mesh.Cell.size());
      auto phi_begin = phi.data_begin();
      auto cell_volume_begin = phi.mesh.Cell.Volume.begin();
      auto A_begin = result.A.begin();
//#pragma omp parallel for firstprivate(phi_begin, cell_volume_begin, A_begin)
      for (auto b_iter = result.b.begin(); b_iter < result.b.end(); ++b_iter){
        auto id_cell = b_iter - result.b.begin();
        Scalar A_value = *(cell_volume_begin + id_cell) / dt;
        T b_value = *(phi_begin + id_cell) * A_value;
        *(A_begin + id_cell) = A_value;
        *b_iter = b_value;
      }
      return result;
    }

    template <class T>
    void ddt_buffered(Scalar dt, fvMappedField<T, on_cell>& phi, diagonalMatrix<Scalar, T, on_cell>& equat){
      auto phi_begin = phi.data_begin();
      auto cell_volume_begin = phi.mesh.Cell.Volume.begin();
      auto A_begin = equat.A.begin();
//#pragma omp parallel for firstprivate(phi_begin, cell_volume_begin, A_begin)
      for (auto b_iter = equat.b.begin(); b_iter < equat.b.end(); ++b_iter){
        auto id_cell = b_iter - equat.b.begin();
        Scalar A_value = *(cell_volume_begin + id_cell) / dt;
        T b_value = *(phi_begin + id_cell) * A_value;
        *(A_begin + id_cell) = A_value;
        *b_iter = b_value;
      }
    }

    class ddt_optimized{
    public:
      diagonalMatrix<Scalar, Scalar, on_cell>& operator() (Scalar dt, fvMappedField<Scalar, on_cell>& phi);
      diagonalMatrix<Scalar, Vector2, on_cell>& operator() (Scalar dt, fvMappedField<Vector2, on_cell>& phi);
      diagonalMatrix<Scalar, Vector3, on_cell>& operator() (Scalar dt, fvMappedField<Vector3, on_cell>& phi);
    private:
      diagonalMatrix<Scalar, Scalar, on_cell> scalar_matrix_buffer;
      diagonalMatrix<Scalar, Vector2, on_cell> vector2_matrix_buffer;
      diagonalMatrix<Scalar, Vector3, on_cell> vector3_matrix_buffer;
    };


  }//--end namespace fv


}//--end namespace GC



#endif
