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
  \flie limiter.h
  \brief Header files for limiter of a gradient to prevent oscilation and limit functions

  \version 0.1
  \author xilin xia
*/

#ifndef LIMITER_H
#define LIMITER_H

#include <functional>
#include "matrix.h"
#include "mapped_field.h"
#include "Scalar.h"
#include "Vector.h"
#include "Tensor.h"
#include "Flag.h"



namespace GC{

  namespace fv{

    ///This is a function class for slope limiter
    class limiterCartesian{
      public:
        limiterCartesian(std::function<Scalar(const Scalar& upwind, const Scalar& downwind)> _limiterFunction):limiterFunction(_limiterFunction){}
        fvMappedField<Vector, on_cell>& operator() (fvMappedField<Scalar, on_cell>& phi, fvMappedField<Vector, on_cell>& gradient);
        fvMappedField<Tensor2, on_cell>& operator() (fvMappedField<Vector2, on_cell>& phi, fvMappedField<Tensor2, on_cell>& gradient);
      private:
        fvMappedField<Vector, on_cell> vector_buffer;
        fvMappedField<Tensor2, on_cell> tensor2_buffer;
        std::function<Scalar(const Scalar& upwind, const Scalar& downwind)> limiterFunction;
    };

    ///This function is the minmod limiter
    
    Scalar minmod(const Scalar& upwind, const Scalar& downwind);

    void LimiterCartesian(fvMappedField<Scalar, on_cell>& phi, fvMappedField<Vector, on_cell>& gradient);
    void LimiterCartesian(fvMappedField<Vector2, on_cell>& phi, fvMappedField<Tensor2, on_cell>& gradient);

  }


}


#endif
