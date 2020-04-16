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
\file time_control.h
\brief Source file for time control class

\version 0.1
\author xilin xia

*/

#include <algorithm>
#include "time_control.h"

namespace GC{

  bool timeControl::is_end(){
    return current_ >= end_;
  }

  Scalar timeControl::current(){
    return current_;
  }

  Scalar timeControl::dt(){
    return dt_;
  }

  void timeControl::forward(){
    dt_ = std::min(dt_, end_ - current_);
    current_ += dt_;
  }

}