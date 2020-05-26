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
\file friction.h
\brief Header file for friction operator

\version 0.1
\author xilin xia

*/

#ifndef FRICTION_H
#define FRICTION_H

#include "Flag.h"
#include "Scalar.h"
#include "Vector.h"
#include <functional>
#include "mapped_field.h"
#include "matrix.h"

namespace GC{

  namespace fv{

    ///This is a function class for calculation friction force
    class friction{
      public:
        friction(std::function<Vector(Scalar& g, Scalar& miu, Scalar& h, Vector& hU, Vector& a, Scalar& dt)> _friction_formula);
        fvMappedField<Vector, on_cell>& operator() (fvScalarFieldOnCell& grav, 
                         fvScalarFieldOnCell& friction_coeff,
                         fvScalarFieldOnCell& h, 
                         fvVectorFieldOnCell& hU,
                         diagonalMatrix<Scalar, Vector, on_cell>& eqn);

      private:
        std::function<Vector(Scalar& g, Scalar& miu, Scalar& h, Vector& hU, Vector& a, Scalar& dt)> friction_formula;
        fvMappedField<Vector, on_cell> buffer;
    };

    template<typename F>
    void Friction(fvScalarFieldOnCell& grav, fvScalarFieldOnCell& friction_coeff, fvScalarFieldOnCell& h, fvVectorFieldOnCell& hU, fvVectorFieldOnCell& z_grad, fvVectorFieldOnCell& friction_force, F friction_formula){
      auto g_begin = grav.data_begin();
      auto miu_begin = friction_coeff.data_begin();
      auto h_begin = h.data_begin();
      auto hU_begin = hU.data_begin();
      auto z_grad_begin = z_grad.data_begin();
      auto friction_force_begin = friction_force.data_begin();

#pragma omp parallel for firstprivate (g_begin,\
                                       miu_begin,\
                                       h_begin,\
                                       hU_begin,\
                                       z_grad_begin,\
                                       friction_force_begin,\
                                       friction_formula)

      for (auto friction_iter = friction_force.data_begin(); friction_iter < friction_force.data_end(); ++friction_iter){
        auto id = friction_iter - friction_force_begin;
        auto g = *(g_begin + id);
        auto miu = *(miu_begin + id);
        auto h_ = *(h_begin + id);
        auto hU_ = *(hU_begin + id);
        auto z_grad_ = *(z_grad_begin + id);
        *friction_iter = friction_formula(g, miu, h_, hU_, z_grad_);
      }
    }
    

  }//end of namespace fv

}//end of namespace GC



#endif

