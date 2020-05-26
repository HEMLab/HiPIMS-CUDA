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
  \file mesh_fv_queries.cc
  \brief Source file for mesh enquiry interface class

*/

#include "mesh_fv_queries.h"
#include "mesh_fv.h"

namespace GC{

  //--assignment operator
  fvMeshQueries& fvMeshQueries::operator= (const fvMeshQueries& rhs){
    if(this != &rhs){
      mesh_ = rhs.mesh_;
      Vertex.my_fvMeshQueries = this;    
      Vertex.Position.my_fvMeshQueries = this;
      Vertex.HalfFacet.my_fvMeshQueries = this;
      Vertex.Cell.my_fvMeshQueries = this;
      HalfFacet.my_fvMeshQueries = this;
      HalfFacet.Centroid.my_fvMeshQueries = this;
      HalfFacet.Normal.my_fvMeshQueries = this;
      HalfFacet.Area.my_fvMeshQueries = this;
      HalfFacet.Cell.my_fvMeshQueries = this;
      Cell.my_fvMeshQueries = this;
      Cell.Volume.my_fvMeshQueries = this;
      Cell.Centroid.my_fvMeshQueries = this;
      Cell.Neighbours.my_fvMeshQueries = this;
      Cell.HalfFacets.my_fvMeshQueries = this;
      Boundary.my_fvMeshQueries = this;
      Boundary.Opposite.my_fvMeshQueries = this;
    }
    return *this;
  }

  //--type name
  const char* fvMeshQueries::type_name(){
    return mesh_->type_name();
  }

  //--Vertex enquiries
  //size
  fvMeshQueries::Vertex::size_type
  fvMeshQueries::Vertex::size(){
    return my_fvMeshQueries->mesh_->vertex_size();
  }

  //Position
  fvMeshQueries::Vertex::Position::iterator
  fvMeshQueries::Vertex::Position::begin(){
    return my_fvMeshQueries->mesh_->vertex_positions_begin();
  }


  fvMeshQueries::Vertex::Position::iterator
  fvMeshQueries::Vertex::Position::end(){
    return my_fvMeshQueries->mesh_->vertex_positions_end();
  }

  //Halffacet;

  fvMeshQueries::Vertex::HalfFacet::iterator
  fvMeshQueries::Vertex::HalfFacet::begin(){
    return my_fvMeshQueries->mesh_->vertex2halffacet_begin();
  }

  fvMeshQueries::Vertex::HalfFacet::iterator
  fvMeshQueries::Vertex::HalfFacet::end(){
    return my_fvMeshQueries->mesh_->vertex2halffacet_end();
  }

  //Cell

  fvMeshQueries::Vertex::Cell::iterator
  fvMeshQueries::Vertex::Cell::begin(){
    return my_fvMeshQueries->mesh_->vertex2cells_begin();
  }

  fvMeshQueries::Vertex::Cell::iterator
  fvMeshQueries::Vertex::Cell::end(){
    return my_fvMeshQueries->mesh_->vertex2cells_end();
  }

  //--Half facet enquiries-----------------------------------
  //size
  fvMeshQueries::HalfFacet::size_type
  fvMeshQueries::HalfFacet::size(){
    return my_fvMeshQueries->mesh_->halffacet_size();
  }

  //Centroid
  fvMeshQueries::HalfFacet::Centroid::iterator
  fvMeshQueries::HalfFacet::Centroid::begin(){
    return my_fvMeshQueries->mesh_->halffacet_centre_positions_begin();
  }

  fvMeshQueries::HalfFacet::Centroid::iterator
  fvMeshQueries::HalfFacet::Centroid::end(){
    return my_fvMeshQueries->mesh_->halffacet_centre_positions_end();
  }

  //Normal
  fvMeshQueries::HalfFacet::Normal::iterator
  fvMeshQueries::HalfFacet::Normal::begin(){
    return my_fvMeshQueries->mesh_->halffacet_normal_directions_begin();
  }

  fvMeshQueries::HalfFacet::Normal::iterator
  fvMeshQueries::HalfFacet::Normal::end(){
    return my_fvMeshQueries->mesh_->halffacet_normal_directions_end();
  }

  //Volume
  fvMeshQueries::HalfFacet::Area::iterator
  fvMeshQueries::HalfFacet::Area::begin(){
    return my_fvMeshQueries->mesh_->halffacet_areas_begin();
  }

  fvMeshQueries::HalfFacet::Area::iterator
  fvMeshQueries::HalfFacet::Area::end(){
    return my_fvMeshQueries->mesh_->halffacet_areas_end();
  }

  //Cell
  fvMeshQueries::HalfFacet::Cell::iterator
  fvMeshQueries::HalfFacet::Cell::begin(){
    return my_fvMeshQueries->mesh_->halffacet2cell_begin();
  }

  fvMeshQueries::HalfFacet::Cell::iterator
  fvMeshQueries::HalfFacet::Cell::end(){
    return my_fvMeshQueries->mesh_->halffacet2cell_end();
  }
  //--End half facet enquiries-------------------------------

  //--Cell enquiries--------------------------------------
  //size
  fvMeshQueries::Cell::size_type
  fvMeshQueries::Cell::size(){
    return my_fvMeshQueries->mesh_->cell_size();
  }

  //Centroid
  fvMeshQueries::Cell::Centroid::iterator
  fvMeshQueries::Cell::Centroid::begin(){
    return my_fvMeshQueries->mesh_->cell_centre_positions_begin();
  }

  fvMeshQueries::Cell::Centroid::iterator
  fvMeshQueries::Cell::Centroid::end(){
    return my_fvMeshQueries->mesh_->cell_centre_positions_end();
  }

  //Volume
  fvMeshQueries::Cell::Volume::iterator
  fvMeshQueries::Cell::Volume::begin(){
    return my_fvMeshQueries->mesh_->cell_volumes_begin();
  }

  fvMeshQueries::Cell::Volume::iterator
  fvMeshQueries::Cell::Volume::end(){
    return my_fvMeshQueries->mesh_->cell_volumes_end();
  }

  //Neighbours
  fvMeshQueries::Cell::Neighbours::iterator
  fvMeshQueries::Cell::Neighbours::begin(){
    return my_fvMeshQueries->mesh_->cell_neighbours_begin();
  }

  fvMeshQueries::Cell::Neighbours::iterator
  fvMeshQueries::Cell::Neighbours::end(){
    return my_fvMeshQueries->mesh_->cell_neighbours_end();
  }

  //Halffacets
  fvMeshQueries::Cell::HalfFacets::iterator
  fvMeshQueries::Cell::HalfFacets::begin(){
    return my_fvMeshQueries->mesh_->cell_halffacets_begin();
  }

  fvMeshQueries::Cell::HalfFacets::iterator
  fvMeshQueries::Cell::HalfFacets::end(){
    return my_fvMeshQueries->mesh_->cell_halffacets_end();
  }

  //--end Cell enquiries----------------------------------

  //--Boundary enquiries----------------------------------
  fvMeshQueries::Boundary::size_type
  fvMeshQueries::Boundary::size(){
    return my_fvMeshQueries->mesh_->boundary_size();
  }

  //Opposite
  fvMeshQueries::Boundary::Opposite::iterator
  fvMeshQueries::Boundary::Opposite::begin(){
    return my_fvMeshQueries->mesh_->boundary2opposite_begin();
  }

  fvMeshQueries::Boundary::Opposite::iterator
  fvMeshQueries::Boundary::Opposite::end(){
    return my_fvMeshQueries->mesh_->boundary2opposite_end();
  }

}

