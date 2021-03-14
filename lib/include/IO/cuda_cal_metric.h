// ======================================================================================
// Name                :    High-Performance Integrated Modelling System
// Description         :    This code pack provides a generic framework for developing 
//                          Geophysical CFD software. Legacy name: GeoClasses
// ======================================================================================
// Version             :    1.0.1 
// Author              :    Xilin Xia
// Create Time         :    2014/10/04
// Update Time         :    2021/01/22
// ======================================================================================
// LICENCE: GPLv3 
// ======================================================================================

/*!
\file cuda_cal_metric.h
\brief Header file for cuda calculating metrics for impact assessment

*/

#ifndef CUDA_CAL_METRIC_H
#define CUDA_CAL_METRIC_H

#include "Flag.h"
#include "Scalar.h"
#include "Vector.h"
#include "cuda_mapped_field.h"

namespace GC{
  
  void cuCalDepthDuration(cuFvMappedField<Scalar, on_cell>& h_old, cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Scalar, on_cell>& t_hGTx, Scalar x, Scalar dt);

}

namespace GC{
  
  void cuCalHazardRating(cuFvMappedField<Scalar, on_cell>& h, cuFvMappedField<Vector, on_cell>& hU, cuFvMappedField<Scalar, on_cell>& hRating);

}

#endif