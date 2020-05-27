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
  \file field_reader.cc
  \brief Source file for basic field file reader class

*/

#include "gisAsciiField_reader.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>

namespace GC{

  gisAsciiFieldReader::gisAsciiFieldReader(const char* filename){
    readin(filename);
  }

  void gisAsciiFieldReader::readin(const char* filename){
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

  }

}// end of namespace GC

