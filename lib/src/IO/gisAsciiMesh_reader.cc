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
  \file gisAscii_reader.cc
  \brief Source file for arcgis ascii file reader class

*/


#include "gisAsciiMesh_reader.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include "Flag.h"
#include "Scalar.h"
#include "Vector.h"

namespace GC{

  ///constructor
  gisAsciiMeshReader::gisAsciiMeshReader(const char* filename){
    readin(filename);
  }


  ///This function reads element connectivity table from arcgis ascii file
  void gisAsciiMeshReader::readin(const char* filename){
    std::ios_base::sync_with_stdio(false);
    std::ifstream input;
    input.open(filename);
    if(!input){
      std::cout<<"error: unable to open input file: "<<filename<<std::endl;
    }
    std::string line;
    std::string word;
    //read in the head information
    size_t n_cols, n_rows;
    getline(input, line);
    std::istringstream stream;
    stream.str(line);
    stream>>word>>n_cols;
    stream.str(std::string());
    stream.clear();
    getline(input, line);
    stream.str(line);
    stream>>word>>n_rows;
    stream.str(std::string());
    stream.clear();
    Scalar x_ori, y_ori, cell_size, nodata_value;
    getline(input, line);
    stream.str(line);
    stream>>word>>x_ori;
    stream.str(std::string());
    stream.clear();
    getline(input, line);
    stream.str(line);
    stream>>word>>y_ori;
    stream.str(std::string());
    stream.clear();
    getline(input, line);
    stream.str(line);
    stream>>word>>cell_size;
    stream.str(std::string());
    stream.clear();
    getline(input, line);
    stream.str(line);
    stream>>word>>nodata_value;
    stream.str(std::string());
    stream.clear();


    //construct a M*N cell matrix
    std::vector<std::vector<int> > cell_markers;
    cell_markers.resize(n_rows);
    int valid_cell_cnt = 0;
    for(size_t i = 0; i < n_rows; ++i){
      getline(input, line);
      std::istringstream steam;
      stream.str(line);
      for(size_t j = 0; j < n_cols; ++j){
        Scalar z_elevation;
        stream >> z_elevation;
        if (z_elevation > nodata_value){
          cell_markers[n_rows - 1 - i].push_back(valid_cell_cnt); //1 denotes valid cell
          element_types.push_back(3); //3 - quadrilateral elememt
          valid_cell_cnt++;
        }else{
          cell_markers[n_rows - 1 - i].push_back(-1); //0 denotes invalid cell
        }
      }
      stream.str(std::string());
      stream.clear();
    }

    element_vertices.resize(valid_cell_cnt);

    //construct a (M+1)*(N+1) vertex matrix
    std::vector<std::vector<int> > vertex_markers;
    vertex_markers.resize(n_rows + 1);
    for(size_t i = 0; i < n_rows + 1; ++i){
      for(size_t j = 0; j < n_cols + 1; ++j){
        vertex_markers[i].push_back(-1);
      }
    }

    //mark all the vertices
    for(size_t i = 0; i < n_rows; ++i){
      for(size_t j = 0; j < n_cols; ++j){
        if(cell_markers[i][j] >= 0){
          vertex_markers[i][j] = 0;
          vertex_markers[i][j+1] = 0;
          vertex_markers[i+1][j+1] = 0;
          vertex_markers[i+1][j] = 0;
        }
      }
    }

    //build the vertices list
    int vertex_cnt = 0;
    for(size_t i = 0; i < n_rows + 1; ++i){
      for(size_t j = 0; j < n_cols + 1; ++j){
        if(vertex_markers[i][j] >= 0){
          vertex_markers[i][j] = vertex_cnt++;
          vertex_positions.push_back(Vector2(x_ori + double(j)*cell_size, y_ori + double(i)*cell_size));
        }
      }
    }

    //build the element connectivity table  
    Flag cell_cnt = 0;
    for(size_t i = 0; i < n_rows; ++i){
      for(size_t j = 0; j < n_cols; ++j){
        if(cell_markers[i][j] >= 0){
          element_vertices[cell_cnt].push_back(vertex_markers[i][j]);
          element_vertices[cell_cnt].push_back(vertex_markers[i][j+1]);
          element_vertices[cell_cnt].push_back(vertex_markers[i+1][j+1]);
          element_vertices[cell_cnt].push_back(vertex_markers[i+1][j]);
          cell_cnt++;
        }
      }
    }

  }

}//end of namespace GC
