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
 \file mesh_interface.h
 \brief Header file for mesh enquiry interface class
        This file provides an interface between field and mesh

*/

#ifndef MESH_INTERFACE_H //header file protector
#define MESH_INTERFACE_H

#include "mesh_fv_queries.h"

namespace GC{

  typedef fvMeshQueries fvMeshQueries;

}

#endif  //end header file protector