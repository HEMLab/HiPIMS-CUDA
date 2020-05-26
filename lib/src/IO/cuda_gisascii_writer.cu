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
  \file cuda_gisascii_writer.cc
  \brief Source file for writing output as arcgis ascii files

*/

#include "cuda_gisascii_writer.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>


namespace GC{

  cuGisAsciiWriter::cuGisAsciiWriter(const char* filename){

    std::ifstream input;
    input.open(filename);
    if(!input){
      std::cout<<"error: unable to open input file: "<<filename<<std::endl;
    }
    
    std::string line;
    std::string word;
    //read in the head information
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

    //resize data
    data.resize(n_cols*n_rows);
    
    //initialise the matrix
    int valid_cell_cnt = 0;
    for(int i = 0; i < n_rows; ++i){
      getline(input, line);
      std::istringstream steam;
      stream.str(line);
      for(int j = 0; j < n_cols; ++j){
        Scalar z_elevation;
        stream >> z_elevation;
        if (z_elevation > nodata_value){
          data[i*n_cols + j] = 0.0;          
          valid_cell_cnt++;
        }else{
          data[i*n_cols + j] = nodata_value;          
        }
      }
      stream.str(std::string());
      stream.clear();
    }

    indices.reserve(valid_cell_cnt);
    for(int i = n_rows-1; i >= 0; --i){
      for (int j = 0; j < n_cols; ++j){
        if (data[i*n_cols + j] > nodata_value){
          indices.push_back(i*n_cols + j);
        }
      }
    }

  }
  

  void cuGisAsciiWriter::write(cuFvMappedField<Scalar, on_cell>& phi, const char* filename, Scalar t, const char* directory){

	  phi.data.sync();    
	  auto size = phi.data.size();
	  auto data_array = phi.data.host_ptr();

    for (int i = 0; i < size; ++i){
      int index = indices[i];
      data[index] = data_array[i];
    }

    FILE * output;
    std::ostringstream out_time;
	  out_time << t;
    std::string name = std::string(directory) + std::string(filename) + "_" + out_time.str() + ".asc";
    output = fopen(name.c_str(), "w");

    //write header
    fprintf(output, "ncols %d\n", n_cols);
    fprintf(output, "nrows %d\n", n_rows);
    fprintf(output, "xllcorner %f\n", x_ori);
    fprintf(output, "yllcorner %f\n", y_ori);
    fprintf(output, "cellsize %f\n", cell_size);
    fprintf(output, "NODATA_value %f\n", nodata_value);

    //write data
    for (int i = 0; i < n_rows; ++i){
      for(int j = 0; j < n_cols; ++j){
        fprintf(output, "%f ", data[i*n_cols + j]);
      } 
      fprintf(output, "\n");
    }
    fclose(output);

  }



}
