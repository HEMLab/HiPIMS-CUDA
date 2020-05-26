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