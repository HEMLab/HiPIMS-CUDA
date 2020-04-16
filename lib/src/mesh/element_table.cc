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
#include "element_table.h"

namespace GC{


  const std::map<Flag, ElementProperty> ElementPropertyTable =

  std::map<Flag, ElementProperty>
  {
    {0,
     ElementProperty{ 1, 
     { ElementProperty::Facet{ 0, ElementProperty::Vertices{ 0 } } },
     point_volume, point_normal, point_centre}
    }, //0 - point
    {1,
     ElementProperty{ 2, 
     { ElementProperty::Facet{ 0, ElementProperty::Vertices{ 0 } }, 
       ElementProperty::Facet{ 0, ElementProperty::Vertices{ 1 } } },
     segment_volume, segment_normal, segment_centre}
    }, //1 - segment
    {2,
     ElementProperty{ 3, 
     { ElementProperty::Facet{ 1, ElementProperty::Vertices{ 0, 1 } },
       ElementProperty::Facet{ 1, ElementProperty::Vertices{ 1, 2 } }, 
       ElementProperty::Facet{ 1, ElementProperty::Vertices{ 2, 0 } } },
     triangle_volume, triangle_normal, triangle_centre}},  // 2 - triangle
    {3,
     ElementProperty{ 4, 
     { ElementProperty::Facet{ 1, ElementProperty::Vertices{ 0, 1 } }, 
       ElementProperty::Facet{ 1, ElementProperty::Vertices{ 1, 2 } }, 
       ElementProperty::Facet{ 1, ElementProperty::Vertices{ 2, 3 } }, 
       ElementProperty::Facet{ 1, ElementProperty::Vertices{ 3, 0 } } },
     quadrilateral_volume, quadrilateral_normal, quadilateral_centre}} //3 - quadrilateral
  };

}
