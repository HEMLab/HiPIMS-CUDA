// ======================================================================================
// Name                :    GeoClasses : Generic Geophysical Flow Modelling Framework
// Description         :    This code pack provides a generic framework for developing 
//                          Geophysical CFD software.
// ======================================================================================
// Version             :    0.1 
// Author              :    Xilin Xia (PhD candidate in Newcastle University)
// Create Time         :    2014/10/04
// Update Time         :    2012/10/29
// ======================================================================================
// Copyright @ Xilin Xia 2014 . All rights reserved.
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
