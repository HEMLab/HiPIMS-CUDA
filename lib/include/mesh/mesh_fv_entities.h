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
 \file mesh_entities.h
 \brief Header file for basic mesh entities

*/

#ifndef MESH_FV_ENTITIES_H
#define MESH_FV_ENTITIES_H

#include <utility>
#include <vector>
#include "Flag.h"
#include "Scalar.h"
#include "Vector.h"
#include "Tensor.h"

namespace GC{

  class ShortDualHandle : public ShortDualFlag{
    public:
      HEMI_DEV_CALLABLE_INLINE_MEMBER
      ShortDualHandle(unsigned int x = 0, unsigned int y = 0):ShortDualFlag(x,y){}
      
      HEMI_DEV_CALLABLE_INLINE_MEMBER
        bool is_boundary(){
        return getx() == 0;
      }

      HEMI_DEV_CALLABLE_INLINE_MEMBER
      void set_boundary_id(unsigned int id){
        setx(0);
        sety(id);
      }

      HEMI_DEV_CALLABLE_INLINE_MEMBER
      void set_local_id(unsigned int id){
        setx(id + 1);
      }

      HEMI_DEV_CALLABLE_INLINE_MEMBER
      void set_global_id(unsigned int id){
        sety(id);
      }

      HEMI_DEV_CALLABLE_INLINE_MEMBER
      unsigned int get_local_id(){
        return getx() - 1;
      }

      HEMI_DEV_CALLABLE_INLINE_MEMBER
      unsigned int get_global_id(){
        return gety();
      }
  };

  typedef Flag Handle;
  typedef std::vector<Handle> HandleSet;
  typedef std::vector<Handle> MultiHandle;
  typedef std::vector<MultiHandle> MultiHandleSet;
  typedef std::vector<ShortDualHandle> ShortDualHandleSet;
  typedef std::vector<ShortDualHandle> CellNeighbour;
  typedef std::vector<CellNeighbour> CellNeighbourSet;
  typedef std::vector<Handle> CellVertice;
  typedef std::vector<CellVertice> CellVerticeSet;
  typedef std::vector<Scalar> scalarProps;
  typedef std::vector<Vector> vectorProps;

}

#endif
