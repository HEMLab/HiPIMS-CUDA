// ======================================================================================
// Name                :    GeoClasses : Generic Geophysical Flow Modelling Framework
// Description         :    This code pack provides a generic framework for developing 
//                          Geophysical CFD software.
// ======================================================================================
// Version             :    0.1 
// Author              :    Xilin Xia (PhD candidate in Newcastle University)
// Create Time         :    2014/10/04
// Update Time         :    2015/10/15
// ======================================================================================
// Copyright @ Xilin Xia 2015 . All rights reserved.
// ======================================================================================

/*!
\flie integrators.h
\brief Header file for integrator template functions
\version 0.1
\author xilin xia
*/

#ifndef INTEGRATORS_H
#define INTEGRATORS_H

namespace GC{

  namespace fv{

    template <typename T, MAPPING_MODES C>
    void EulerIntegrator(fvMappedField<T, C>& phi, fvMappedField<T, C>& phi_acc, Scalar delta_t, Scalar current_t){
      auto phi_data_begin = phi.data_begin();
      auto phi_acc_data_begin = phi_acc.data_begin();
      for (auto phi_data_iter = phi.data_begin(); phi_data_iter < phi.data_end(); ++phi_data_iter){
        Flag id = phi_data_iter - phi_data_begin;
        auto _phi = *phi_data_iter;
        auto _phi_acc = *(phi_acc_data_begin + id);
        *phi_data_iter = _phi + delta_t*_phi_acc;
      }
      phi.update_time(current_t);
      phi.update_boundary_values();
    }

  }//end of namespace fv

}//--end of namespace GC


#endif