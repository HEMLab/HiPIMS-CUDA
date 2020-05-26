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