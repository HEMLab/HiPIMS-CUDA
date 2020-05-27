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
\brief Source file for friction operator

\version 0.1
\author xilin xia

*/

#include "friction.h"

namespace GC{

  namespace fv{

    friction::friction(std::function<Vector(Scalar& g, Scalar& miu, Scalar& h, Vector& hU, Vector& a, Scalar& dt)> _friction_formula):friction_formula(_friction_formula){}
    fvMappedField<Vector, on_cell>& friction::operator() (fvScalarFieldOnCell& grav, 
                         fvScalarFieldOnCell& friction_coeff,
                         fvScalarFieldOnCell& h, 
                         fvVectorFieldOnCell& hU,
                         diagonalMatrix<Scalar, Vector, on_cell>& eqn){
      buffer.initialize_by_field(hU);
      auto g_begin = grav.data_begin();
      auto miu_begin = friction_coeff.data_begin();
      auto h_begin = h.data_begin();
      auto hU_begin = hU.data_begin();
      auto A_begin = eqn.A.begin();
      auto b_begin = eqn.b.begin(); ///b in equation Ax = b, it has phi/dt + dphidt
      auto friction_begin = buffer.data_begin();
      for (auto friction_iter = buffer.data_begin(); friction_iter < buffer.data_end(); ++friction_iter){
        auto id = friction_iter - friction_begin;
        auto g = *(g_begin + id);
        auto miu = *(miu_begin + id);
        auto h_ = *(h_begin + id);
        auto hU_ = *(hU_begin + id);
        auto A = *(A_begin + id);
        auto b = *(b_begin + id);
        auto dt = 1.0/A;
        auto a = b - hU_/dt;
        *friction_iter = friction_formula(g, miu, h_, hU_, a, dt);
      }
      return buffer;
    }
    
  }//end of namespace fv

}//end of namespace GC


