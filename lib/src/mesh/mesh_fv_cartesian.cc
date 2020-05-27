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
\file mesh_fv_cartesian.cc
\brief Source file for finite volume mesh class

*/

#include "mesh_fv_cartesian.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include "Flag.h"
#include "Scalar.h"
#include "Vector.h"

namespace GC{

  CartesianFvMesh::CartesianFvMesh(const char* filename):type("CartesianFvMesh"){
    
    std::ifstream input;
    input.open(filename);
    if (!input){
      std::cout << "error: unable to open input file: " << filename << std::endl;
    }
    std::string line;
    std::string word;
    //read in the head information
    size_t n_cols, n_rows;
    getline(input, line);
    std::istringstream stream;
    stream.str(line);
    stream >> word >> n_cols;
    stream.str(std::string());
    stream.clear();
    getline(input, line);
    stream.str(line);
    stream >> word >> n_rows;
    stream.str(std::string());
    stream.clear();
    Scalar x_ori, y_ori, cell_size, nodata_value;
    getline(input, line);
    stream.str(line);
    stream >> word >> x_ori;
    stream.str(std::string());
    stream.clear();
    getline(input, line);
    stream.str(line);
    stream >> word >> y_ori;
    stream.str(std::string());
    stream.clear();
    getline(input, line);
    stream.str(line);
    stream >> word >> cell_size;
    stream.str(std::string());
    stream.clear();
    getline(input, line);
    stream.str(line);
    stream >> word >> nodata_value;
    stream.str(std::string());
    stream.clear();


    //construct a M*N cell matrix
    std::vector<std::vector<int> > cell_markers;
    cell_markers.resize(n_rows);
    int valid_cell_cnt = 0;
    for (size_t i = 0; i < n_rows; ++i){
      getline(input, line);
      std::istringstream steam;
      stream.str(line);
      for (size_t j = 0; j < n_cols; ++j){
        Scalar z_elevation;
        stream >> z_elevation;
        if (z_elevation > nodata_value){
          cell_markers[i].push_back(valid_cell_cnt); //1 denotes valid cell
          valid_cell_cnt++;
        }
        else{
          cell_markers[i].push_back(-1); //0 denotes invalid cell
        }
      }
      stream.str(std::string());
      stream.clear();
    }

    cell_neighbours.resize(valid_cell_cnt);
    cell_volumes.resize(valid_cell_cnt);
    cell_centre_positions.resize(valid_cell_cnt);
    cell_halffacets.resize(valid_cell_cnt);
    halffacet_normal_directions.resize(4);
    halffacet_areas.resize(4);

    valid_cell_cnt = 0;
    for (size_t i = 0; i < n_rows; ++i){

      for (size_t j = 0; j < n_cols; ++j){
         
        int cell_id = cell_markers[n_rows - 1 - i][j];
        if (cell_id >= 0){
          cell_markers[n_rows - 1 - i][j] = valid_cell_cnt;
          valid_cell_cnt++;
        }

      }

    }

    int boundary_cnt = 0;
    for (size_t i = 0; i < n_rows; ++i){

      for (size_t j = 0; j < n_cols; ++j){        
        int cell_id = cell_markers[i][j];
        if (cell_id >= 0){
          ShortDualHandle cell_handle;
          cell_handle.set_global_id(cell_id);
          //south
          cell_handle.set_local_id(0);
          if (i == n_rows-1){//south boundary
            ShortDualHandle _boundary_handle(0, boundary_cnt);
            cell_neighbours[cell_id].push_back(_boundary_handle);
            boundary2opposite_handles.push_back(cell_handle);
            boundary_cnt++;
          }
          else{
            int neighbour_id = cell_markers[i + 1][j];
            if (neighbour_id < 0){
              ShortDualHandle _boundary_handle(0, boundary_cnt);
              cell_neighbours[cell_id].push_back(_boundary_handle);
              boundary2opposite_handles.push_back(cell_handle);
              boundary_cnt++;
            }
            else{
              ShortDualHandle _neighbour_handle;
              _neighbour_handle.set_local_id(2);
              _neighbour_handle.set_global_id(neighbour_id);
              cell_neighbours[cell_id].push_back(_neighbour_handle);
            }
          }
          //east
          cell_handle.set_local_id(1);
          if (j == n_cols-1){//east boundary
            ShortDualHandle _boundary_handle(0, boundary_cnt);
            cell_neighbours[cell_id].push_back(_boundary_handle);
            boundary2opposite_handles.push_back(cell_handle);
            boundary_cnt++;
          }
          else{
            int neighbour_id = cell_markers[i][j+1];
            if (neighbour_id < 0){
              ShortDualHandle _boundary_handle(0, boundary_cnt);
              cell_neighbours[cell_id].push_back(_boundary_handle);
              boundary2opposite_handles.push_back(cell_handle);
              boundary_cnt++;
            }
            else{
              ShortDualHandle _neighbour_handle;
              _neighbour_handle.set_local_id(3);
              _neighbour_handle.set_global_id(neighbour_id);
              cell_neighbours[cell_id].push_back(_neighbour_handle);
            }
          }
          //north
          cell_handle.set_local_id(2);
          if (i == 0){//north boundary
            ShortDualHandle _boundary_handle(0, boundary_cnt);
            cell_neighbours[cell_id].push_back(_boundary_handle);
            boundary2opposite_handles.push_back(cell_handle);
            boundary_cnt++;
          }
          else{
            int neighbour_id = cell_markers[i - 1][j];
            if (neighbour_id < 0){
              ShortDualHandle _boundary_handle(0, boundary_cnt);
              cell_neighbours[cell_id].push_back(_boundary_handle);
              boundary2opposite_handles.push_back(cell_handle);
              boundary_cnt++;
            }
            else{
              ShortDualHandle _neighbour_handle;
              _neighbour_handle.set_local_id(0);
              _neighbour_handle.set_global_id(neighbour_id);
              cell_neighbours[cell_id].push_back(_neighbour_handle);
            }
          }
          //west
          cell_handle.set_local_id(3);
          if (j == 0){//west boundary
            ShortDualHandle _boundary_handle(0, boundary_cnt);
            cell_neighbours[cell_id].push_back(_boundary_handle);
            boundary2opposite_handles.push_back(cell_handle);
            boundary_cnt++;
          }
          else{
            int neighbour_id = cell_markers[i][j-1];
            if (neighbour_id < 0){
              ShortDualHandle _boundary_handle(0, boundary_cnt);
              cell_neighbours[cell_id].push_back(_boundary_handle);
              boundary2opposite_handles.push_back(cell_handle);
              boundary_cnt++;
            }
            else{
              ShortDualHandle _neighbour_handle;
              _neighbour_handle.set_local_id(1);
              _neighbour_handle.set_global_id(neighbour_id);
              cell_neighbours[cell_id].push_back(_neighbour_handle);
            }
          }
          Scalar x_pos = x_ori + double(j)*cell_size;
          Scalar y_pos = y_ori + double(n_rows -1 - i)*cell_size;
          Vector2 position(x_pos, y_pos);
          cell_centre_positions[cell_id] = position;
          Scalar volume = cell_size*cell_size;
          cell_volumes[cell_id]=volume;
        }
      }

    }

    cell_neighbours.shrink_to_fit();
    cell_centre_positions.shrink_to_fit();
    cell_volumes.shrink_to_fit();

    cell_markers.resize(0);
    cell_markers.shrink_to_fit();

    Vector2 normal[4];
    normal[0] = Vector2(0.0, -1.0);
    normal[1] = Vector2(1.0, 0.0);
    normal[2] = Vector2(0.0, 1.0);
    normal[3] = Vector2(-1.0, 0.0);
    for (int i = 0; i < cell_neighbours.size(); ++i){
      for (int j = 0; j < 4; j++){
        cell_halffacets[i].push_back(j);
      }
    }

    for (int j = 0; j < 4; j++){
      halffacet_normal_directions.push_back(normal[j]);
      halffacet_areas.push_back(cell_size);
    }


    cell_halffacets.shrink_to_fit();
    halffacet_normal_directions.shrink_to_fit();
    halffacet_areas.shrink_to_fit();


  }
}