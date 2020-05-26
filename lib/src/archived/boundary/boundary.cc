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
\flie boundary.cc
\brief Source file for boundary functions

\version 0.1
\author xilin xia
*/

#include "boundary.h"
#include "boundary_func_table.h"

namespace GC{

  Vector3 GetBoundary(const ShortTripleFlag& boundary_type,
                      const Vector3& in_value,
                      const Vector3& normal){
    Flag primary_type = boundary_type.getx();
    Flag secondary_type = boundary_type.gety();
    Flag source_id = boundary_type.getz();
    Vector3 out_value;

    if (0 == primary_type){ //cannot evolve

    }
    else if (1 == primary_type){ //fixed gradient 

    }
    else if (2 == primary_type){ //zero gradient or wall
      auto boundary_func_iter = BoundaryFuncTable::ZeroGradientAndWall.find(secondary_type);
      if (boundary_func_iter != BoundaryFuncTable::ZeroGradientAndWall.end()){
        out_value = boundary_func_iter->second(in_value, normal);
      }
    }
    else if (3 == primary_type){ //mixed 

    }
    else if (4 == primary_type){//calculated

    }

    return out_value;

  }


};