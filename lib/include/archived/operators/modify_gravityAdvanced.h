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
 \file modify_gravityAdvanced.h
 \brief Header file for advanced modify gravity function 

 \version 0.1
 \author xilin xia
*/ 

#ifndef MODIFY_GRAVITY_ADVANCED_H
#define MODIFY_GRAVITY_ADVANCED_H

#include "mapped_field.h"

namespace GC{

  namespace fv{

    class modifyGravityAdvanced{
      public:
        fvScalarFieldOnCell& operator() (fvScalarFieldOnCell& gravity, fvVectorFieldOnCell& z_gradient, fvVectorFieldOnCell& h_gradient);
      private:
        fvScalarFieldOnCell buffer;
    };

  }

}


#endif

