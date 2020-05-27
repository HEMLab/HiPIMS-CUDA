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
 \file mesh_fv_queries.h
 \brief Header file for mesh enquiry interface class

*/

#ifndef MESH_FV_QURIES_H //header file protector
#define MESH_FV_QURIES_H

#include <memory>
#include "mesh_fv_entities.h"

namespace GC{

  class unstructuredFvMesh;

  ///This class provides an interface for all mesh classes
  class fvMeshQueries{
    private:
      std::shared_ptr<unstructuredFvMesh> mesh_; 
    public:
      fvMeshQueries(std::shared_ptr<unstructuredFvMesh> _mesh = nullptr):mesh_(_mesh), Vertex(this), HalfFacet(this), Cell(this), Boundary(this){}
      fvMeshQueries(const fvMeshQueries& other):mesh_(other.mesh_),Vertex(this), HalfFacet(this), Cell(this), Boundary(this){}
      fvMeshQueries& operator= (const fvMeshQueries& rhs);
    public:
      ///This function returns the type name of the mesh
      const char* type_name();
      ///This class is a memberspace for Vertex, it provides a wrapped interface for all enquires regarding a vertex   
      class Vertex{
        private:
          friend fvMeshQueries;
          fvMeshQueries* my_fvMeshQueries;
          Vertex(fvMeshQueries* _my_fvMeshQueries):
            my_fvMeshQueries(_my_fvMeshQueries),
            Position(my_fvMeshQueries),
            HalfFacet(my_fvMeshQueries),
            Cell(my_fvMeshQueries){}
        public:
          typedef size_t size_type;
          size_type size();
          ///This class is a container for position of this vertex 
          class Position{
            public:
              typedef vectorProps::iterator iterator;
              typedef vectorProps::const_iterator const_iterator;  
              iterator begin();
              iterator end();
            private:
              friend fvMeshQueries;
              Position(fvMeshQueries* _my_fvMeshQueries):my_fvMeshQueries(_my_fvMeshQueries){}
              fvMeshQueries* my_fvMeshQueries;
          }Position; //--end container-----------------------------------

          ///This class is a container for the half facet start from this vertex 
          class HalfFacet{
            public:
              typedef ShortDualHandleSet::iterator iterator;
              typedef ShortDualHandleSet::const_iterator const_iterator;  
              iterator begin();
              iterator end();
            private:
              friend fvMeshQueries;
              HalfFacet(fvMeshQueries* _my_fvMeshQueries):my_fvMeshQueries(_my_fvMeshQueries){}
              fvMeshQueries* my_fvMeshQueries;
          }HalfFacet; //--end container------------------------------------

          ///This class is a container for the cells incident of this vertex 
          class Cell{
            public:
              typedef MultiHandleSet::iterator iterator;
              typedef MultiHandleSet::const_iterator const_iterator;  
              iterator begin();
              iterator end();
            private:
              friend fvMeshQueries;
              Cell(fvMeshQueries* _my_fvMeshQueries):my_fvMeshQueries(_my_fvMeshQueries){}
              fvMeshQueries* my_fvMeshQueries;
          }Cell; //--end container------------------------------------
      }Vertex;//--end container------------------------------------
      
      ///This class is a memberspace for half-facet, it provides a wrapped interface for all enquires regarding half-facet
      class HalfFacet{
        private:
          friend fvMeshQueries;
          fvMeshQueries* my_fvMeshQueries;
          HalfFacet(fvMeshQueries* _my_fvMeshQueries):
            my_fvMeshQueries(_my_fvMeshQueries),
            Centroid(my_fvMeshQueries),
            Normal(my_fvMeshQueries),
            Area(my_fvMeshQueries),
            Cell(my_fvMeshQueries){}
        public:
          typedef size_t size_type;
          size_type size();
          ///This class is a container for the centoids of this halffacet 
          class Centroid{
            public:
              typedef vectorProps::iterator iterator;
              typedef vectorProps::const_iterator const_iterator;  
              iterator begin();
              iterator end();
            private:
              friend fvMeshQueries;
              Centroid(fvMeshQueries* _my_fvMeshQueries):my_fvMeshQueries(_my_fvMeshQueries){}
              fvMeshQueries* my_fvMeshQueries;
          }Centroid; //--end container------------------------------------

          ///This class is a container for the normal directions of this halffacet 
          class Normal{
            public:
              typedef vectorProps::iterator iterator;
              typedef vectorProps::const_iterator const_iterator;  
              iterator begin();
              iterator end();
            private:
              friend fvMeshQueries;
              Normal(fvMeshQueries* _my_fvMeshQueries):my_fvMeshQueries(_my_fvMeshQueries){}
              fvMeshQueries* my_fvMeshQueries;
          }Normal; //--end container------------------------------------

          ///This class is a container for the volume of this halffacet 
          class Area{
            public:
              typedef scalarProps::iterator iterator;
              typedef scalarProps::const_iterator const_iterator;  
              iterator begin();
              iterator end();
            private:
              friend fvMeshQueries;
              Area(fvMeshQueries* _my_fvMeshQueries):my_fvMeshQueries(_my_fvMeshQueries){}
              fvMeshQueries* my_fvMeshQueries;
          }Area; //--end container------------------------------------

          ///This class is a container for the cells this halffacet belongs to 
          class Cell{
            public:
              typedef ShortDualHandleSet::iterator iterator;
              typedef ShortDualHandleSet::const_iterator const_iterator;  
              iterator begin();
              iterator end();
            private:
              friend fvMeshQueries;
              Cell(fvMeshQueries* _my_fvMeshQueries):my_fvMeshQueries(_my_fvMeshQueries){}
              fvMeshQueries* my_fvMeshQueries;
          }Cell; //--end container------------------------------------
      }HalfFacet;//--end memberspace------------------------------------

      ///This class is a memberspace for cell, it provides a wrapped interface for all enquires regarding element
      class Cell{
        private:
          friend fvMeshQueries;
          fvMeshQueries* my_fvMeshQueries;
          Cell(fvMeshQueries* _my_fvMeshQueries):
            my_fvMeshQueries(_my_fvMeshQueries), 
            Centroid(my_fvMeshQueries),
            Volume(my_fvMeshQueries),
            Neighbours(my_fvMeshQueries),
            HalfFacets(my_fvMeshQueries){}
        public:
          typedef size_t size_type;
          size_type size();
          ///This class is a container for the centroid of this cell 
          class Centroid{
            public:
              typedef vectorProps::iterator iterator;
              typedef vectorProps::const_iterator const_iterator;  
              iterator begin();
              iterator end();
            private:
              friend fvMeshQueries;
              Centroid(fvMeshQueries* _my_fvMeshQueries):my_fvMeshQueries(_my_fvMeshQueries){}
              fvMeshQueries* my_fvMeshQueries;
          }Centroid; //--end container------------------------------------

          ///This class is a container for the volume of this cell
          class Volume{
            public:
              typedef scalarProps::iterator iterator;
              typedef scalarProps::const_iterator const_iterator;  
              iterator begin();
              iterator end();
            private:
              friend fvMeshQueries;
              Volume(fvMeshQueries* _my_fvMeshQueries):my_fvMeshQueries(_my_fvMeshQueries){}
              fvMeshQueries* my_fvMeshQueries;
          }Volume; //--end container------------------------------------

          ///This class is a container for cell neighbours
          class Neighbours{
            public:
              typedef CellNeighbourSet::iterator iterator;
              typedef CellNeighbourSet::const_iterator const_iterator;  
              iterator begin();
              iterator end();
            private:
              friend fvMeshQueries;
              Neighbours(fvMeshQueries* _my_fvMeshQueries):my_fvMeshQueries(_my_fvMeshQueries){}
              fvMeshQueries* my_fvMeshQueries;
          }Neighbours; //--end container------------------------------------

          ///This class is a container for half-facets it contains
          class HalfFacets{
            public:
              typedef MultiHandleSet::iterator iterator;
              typedef MultiHandleSet::const_iterator const_iterator;  
              iterator begin();
              iterator end();
            private:
              friend fvMeshQueries;
              HalfFacets(fvMeshQueries* _my_fvMeshQueries):my_fvMeshQueries(_my_fvMeshQueries){}
              fvMeshQueries* my_fvMeshQueries;
          }HalfFacets; //--end container------------------------------------
      }Cell; //--end memberspace----------------------------------------------

      ///This class is a memberspace for boundary, it provides a wrapped interface for all enquires regarding element
      class Boundary{
        private:
          friend fvMeshQueries;
          fvMeshQueries* my_fvMeshQueries;
          Boundary(fvMeshQueries* _my_fvMeshQueries):
            my_fvMeshQueries(_my_fvMeshQueries),
            Opposite(my_fvMeshQueries)
            {}
        public:
          typedef size_t size_type;
          size_type size();
          ///This class is a container for position of this vertex 
          class Opposite{
            public:
              typedef ShortDualHandleSet::iterator iterator;
              typedef ShortDualHandleSet::const_iterator const_iterator;  
              iterator begin();
              iterator end();
            private:
              friend fvMeshQueries;
              Opposite(fvMeshQueries* _my_fvMeshQueries):my_fvMeshQueries(_my_fvMeshQueries){}
              fvMeshQueries* my_fvMeshQueries;
          }Opposite; //--end container-----------------------------------
      }Boundary; //--end memberspace-------------------------------------
  };

}
#endif  //end header file protector
