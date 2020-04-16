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
\file cuda_kernel_launch_parameters.h
\brief Header file for time control class

\version 0.1
\author xilin xia

*/

#ifndef BLOCKS_PER_GRID
#define BLOCKS_PER_GRID 128
#endif

#ifndef THREADS_PER_BLOCK
#define THREADS_PER_BLOCK 256
#endif
