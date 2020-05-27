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
\file time_control.h
\brief Source file for time control class

\version 0.1
\author xilin xia

*/

#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <utility>
//#include <omp.h>
#include "mapped_field.h"

namespace GC{//--namespace GC

  template <typename T, MAPPING_MODES C> class fvMappedField;
  ///This is a template class for diagonal matrix for explicit methods, it can also be
  //interpreted as a space-time field
  template <typename T_A, typename T_b, MAPPING_MODES C>
  class diagonalMatrix{
  public:
    diagonalMatrix() = default;

  public:
    fvMeshQueries mesh;
    std::vector<T_A> A;
    std::vector<T_b> b;
//    std::vector<ShortTripleFlag> boundary; 

  public:
    template <typename T>
    void initialize_by_field(fvMappedField<T, C>& phi){
      mesh = phi.mesh;
      A.resize(mesh.Cell.size());
      b.resize(mesh.Cell.size());      
      //boundary.resize(mesh.Boundary.size());
      //auto phi_boundary_begin = phi.boundary_type_begin();
      //for(auto boundary_iter = boundary.begin(); boundary_iter < boundary.end(); ++boundary_iter){
      //  auto id = boundary_iter - boundary.begin();
      //  *boundary_iter = *(phi_boundary_begin + id);
      //}
    }
  };


//  template <class T_A, class T_b, MAPPING_MODES C>
//  diagonalMatrix<T_A, T_b, C> operator+ (diagonalMatrix<T_A, T_b, C>& lhs, fvMappedField<T_b, C>& rhs){
//    diagonalMatrix<T_A, T_b, C> result(lhs);
//    auto rhs_begin = rhs.data_begin();
//    auto result_begin = result.b.begin();
////#pragma omp parallel for firstprivate(rhs_begin, result_begin)
//    for (auto result_iter = result.b.begin(); result_iter < result.b.end(); ++result_iter){
//      auto n = result_iter - result_begin;
//      *result_iter -= *(rhs_begin + n); //here we are using - rather than +
//    }
//    return result;
//  }
//
//  template <class T_A, class T_b, MAPPING_MODES C>
//  diagonalMatrix<T_A, T_b, C> operator+ (diagonalMatrix<T_A, T_b, C>&& lhs, fvMappedField<T_b, C>& rhs){
//    diagonalMatrix<T_A, T_b, C> result(std::move(lhs));
//    auto rhs_begin = rhs.data_begin();
//    auto result_begin = result.b.begin();
////#pragma omp parallel for firstprivate(rhs_begin, result_begin)
//    for (auto result_iter = result.b.begin(); result_iter < result.b.end(); ++result_iter){
//      auto n = result_iter - result_begin;
//      *result_iter -= *(rhs_begin + n); //here we are using - rather than +
//    }
//    return result;
//  }

  template <class T_A, class T_b, MAPPING_MODES C>
  diagonalMatrix<T_A, T_b, C>& operator== (diagonalMatrix<T_A, T_b, C>& lhs, fvMappedField<T_b, C>& rhs){
    auto rhs_begin = rhs.data_begin();
    auto lhs_begin = lhs.b.begin();
//#pragma omp parallel for firstprivate(rhs_begin, lhs_begin)
    for (auto lhs_iter = lhs.b.begin(); lhs_iter < lhs.b.end(); ++lhs_iter){
      auto n = lhs_iter - lhs_begin;
      *lhs_iter += *(rhs_begin + n); //here we are using +    
    }
    return lhs;
  }

  template <class T_A, class T_b, MAPPING_MODES C>
  diagonalMatrix<T_A, T_b, C>& operator== (diagonalMatrix<T_A, T_b, C>& lhs, fvMappedField<T_b, C>&& rhs){
    auto rhs_begin = rhs.data_begin();
    auto lhs_begin = lhs.b.begin();
//#pragma omp parallel for firstprivate(rhs_begin, lhs_begin)
    for (auto lhs_iter = lhs.b.begin(); lhs_iter < lhs.b.end(); ++lhs_iter){
      auto n = lhs_iter - lhs_begin;
      *lhs_iter += *(rhs_begin + n); //here we are using +    
    }
    return lhs;
  }
  

//  template <class T_A, class T_b, MAPPING_MODES C>
//  diagonalMatrix<T_A, T_b, C> operator- (diagonalMatrix<T_A, T_b, C>& lhs, fvMappedField<T_b, C>& rhs){
//    diagonalMatrix<T_A, T_b, C> result(lhs);
//    auto rhs_begin = rhs.data_begin();
//    auto result_begin = result.b.begin();
////#pragma omp parallel for firstprivate(rhs_begin, result_begin)
//    for (auto result_iter = result.b.begin(); result_iter < result.b.end(); ++result_iter){
//      auto n = result_iter - result_begin;
//      *result_iter += *(rhs_begin + n); //here we are using + rather than -
//    }
//    return result;
//  }
//
//  template <class T_A, class T_b, MAPPING_MODES C>
//  diagonalMatrix<T_A, T_b, C> operator- (diagonalMatrix<T_A, T_b, C>&& lhs, fvMappedField<T_b, C>& rhs){
//    diagonalMatrix<T_A, T_b, C> result(std::move(lhs));
//    auto rhs_begin = rhs.data_begin();
//    auto result_begin = result.b.begin();
////#pragma omp parallel for firstprivate(rhs_begin, result_begin)
//    for (auto result_iter = result.b.begin(); result_iter < result.b.end(); ++result_iter){
//      auto n = result_iter - result_begin;
//      *result_iter += *(rhs_begin + n); //here we are using + rather than -
//    }
//    return result;
//  }
//
//  template <class T_A, class T_b, MAPPING_MODES C>
//  diagonalMatrix<T_A, T_b, C>& operator-= (diagonalMatrix<T_A, T_b, C>& lhs, fvMappedField<T_b, C>& rhs){
//    auto rhs_begin = rhs.data_begin();
//    auto lhs_begin = lhs.b.begin();
////#pragma omp parallel for firstprivate(rhs_begin, lhs_begin)
//    for (auto lhs_iter = lhs.b.begin(); lhs_iter < lhs.b.end(); ++lhs_iter){
//      auto n = lhs_iter - lhs_begin;
//      *lhs_iter += *(rhs_begin + n); //here we are using + rather than -
//    }
//    return lhs;
//  }

  template <class T, MAPPING_MODES C>
  void Solve(diagonalMatrix<Scalar, T, C>& equat, fvMappedField<T, C>& phi){
//#pragma omp parallel for    
    for (auto phi_iter = phi.data_begin(); phi_iter < phi.data_end(); ++phi_iter){
      auto id = phi_iter - phi.data_begin();
      *phi_iter = equat.b[id] / equat.A[id];
      id++;
    }
  }

  template <class T, MAPPING_MODES C>
  void Solve(diagonalMatrix<Scalar, T, C>&& equat, fvMappedField<T, C>& phi){
//#pragma omp parallel for    
    for (auto phi_iter = phi.data_begin(); phi_iter < phi.data_end(); ++phi_iter){
      auto id = phi_iter - phi.data_begin();
      *phi_iter = equat.b[id] / equat.A[id];
      id++;
    }
  }
}//--end namespace GC

#endif
