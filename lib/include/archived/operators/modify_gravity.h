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
 \file modify_gravity.h
 \brief Header file for modify gravity function 

 \version 0.1
 \author xilin xia
*/ 




#ifndef MODIFY_GRAVITY_H
#define MODIFY_GRAVITY_H

#include "mapped_field.h"

namespace GC{

  namespace fv{

    class modifyGravity{
      public:
        fvScalarFieldOnCell& operator() (fvScalarFieldOnCell& gravity, fvVectorFieldOnCell& z_gradient);
      private:
        fvScalarFieldOnCell buffer;
    };

  } //end of namespace fv

}//end of namespace GC



#endif
