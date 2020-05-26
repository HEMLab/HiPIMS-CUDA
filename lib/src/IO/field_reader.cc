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

#include "field_reader.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>

namespace GC{

  fieldReader::fieldReader(const char* filename){
    readin_field(filename);
  }

  void fieldReader::readin_field(const char* filename){
    std::ios_base::sync_with_stdio(false);
    std::ifstream input;
    input.open(filename);
    if (!input){
      std::cout << "error: unable to open input file: " << filename << std::endl;
      return;
    }
    std::string line;
    std::string word;
    //reading elements
    getline(input, line);
    size_t n_elements;
    getline(input, line);
    std::istringstream stream;
    stream.str(line);
    stream >> n_elements;
    getline(input, line);
    stream.str(std::string());
    stream.clear();
    Vector3 value;
    Flag id_element;
    for (size_t i = 0; i < n_elements; i++){
      getline(input, line);
      stream.str(line);
      value = 0.0;
      stream >> id_element >> value;
      data.insert({ id_element, value });
      stream.str(std::string());
      stream.clear();		//clear the istringstream
    }
    getline(input, line);
    getline(input, line);
    size_t n_boundaries;
    stream.str(line);
    stream >> n_boundaries;
    getline(input, line);
    stream.str(std::string());
    stream.clear();
    Flag id_boundary, primary_type, secondary_type, source_id;
    for (size_t i = 0; i < n_boundaries; i++){
      getline(input, line);
      stream.str(line);
      primary_type = 0;
      secondary_type = 0;
      source_id = 0;
      stream >> id_boundary >> primary_type >> secondary_type >>source_id;
      boundary_type.insert({ id_boundary, ShortTripleFlag(primary_type, secondary_type, source_id) });
      stream.str(std::string());
      stream.clear();		//clear the istringstream
    }
  }

  bool completeFieldReader::readin_region_mask(const char* filename){
    std::ios_base::sync_with_stdio(false);
    std::ifstream input;
    input.open(filename);
    if (!input){
      return false;
    }
    std::string line;
    std::string word;
    //reading elements
    getline(input, line);
    size_t n_elements;
    getline(input, line);
    std::istringstream stream;
    stream.str(line);
    stream >> n_elements;
    getline(input, line);
    stream.str(std::string());
    stream.clear();
    Flag value;
    Flag id_element;
    for (size_t i = 0; i < n_elements; i++){
      getline(input, line);
      stream.str(line);
      value = 0;
      stream >> id_element >> value;
      region_mask.insert({ id_element, value });
      stream.str(std::string());
      stream.clear();		//clear the istringstream
    }
    return true;
  }

  completeFieldReader::completeFieldReader(const char* path, const char* field_name){
    std::string filename = std::string(path) + std::string(field_name) + ".dat";
    readin_field(filename.c_str());
    for (unsigned int i = 0;; ++i){
      std::string filename;
      std::ostringstream file_id;
      file_id << i;
      filename = std::string(path) + std::string(field_name) + "_BC_" + file_id.str() + ".dat";
      bool success = readin_boundary_source(filename.c_str(), i);
      if (!success){
        break;
      }
    }
    filename = std::string(path) + std::string(field_name) + "_mask.dat";
    if (readin_region_mask(filename.c_str())){
      std::string filename = std::string(path) + std::string(field_name) + "_source_all.dat";
      std::ifstream input(filename);
      if(input){
        int num_of_regions;
        std::istringstream stream;
        std::string line;
        Scalar t;
        Scalar value;
        getline(input, line);
        stream.str(line);
        stream >> num_of_regions;
        stream.str(std::string());
        stream.clear();		//clear the istringstream
        data_time_series.resize(num_of_regions);
        data_source.resize(num_of_regions);
        while (!input.eof()){
          getline(input, line);
          stream.str(line);
          stream >> t;
          for (int i = 0; i < num_of_regions; i++){
            stream >> value;
            data_time_series[i].push_back(t);
            data_source[i].push_back(value);
          }
          stream.str(std::string());
          stream.clear();		//clear the istringstream
        }
      }
      else{
        for (unsigned int i = 0;; ++i){
          std::string filename;
          std::ostringstream file_id;
          file_id << i;
          filename = std::string(path) + std::string(field_name) + "_source_" + file_id.str() + ".dat";
          bool success = readin_data_source(filename.c_str(), i);
          if (!success){
            break;
          }
        }
      }
    }
  }

  bool completeFieldReader::readin_boundary_source(const char* filename, const unsigned int cnt){
    std::ios_base::sync_with_stdio(false);
    std::ifstream input;
    input.open(filename);
    if (!input){
      return false;
    }
    else{
      time_series.resize(cnt + 1);
      boundary_source.resize(cnt + 1);
      std::istringstream stream;
      std::string line;
      Scalar t;
      Vector3 value;
      while (!input.eof()){
        getline(input, line);
        stream.str(line);
        stream >> t >> value;
        time_series[cnt].push_back(t);
        boundary_source[cnt].push_back(value);
        stream.str(std::string());
        stream.clear();		//clear the istringstream
      }
      return true;
    }
  }

  bool completeFieldReader::readin_data_source(const char* filename, const unsigned int cnt){
    std::ios_base::sync_with_stdio(false);
    std::ifstream input;
    input.open(filename);
    if (!input){
      return false;
    }
    else{
      data_time_series.resize(cnt + 1);
      data_source.resize(cnt + 1);
      std::istringstream stream;
      std::string line;
      Scalar t;
      Vector3 value;
      while (!input.eof()){
        getline(input, line);
        stream.str(line);
        stream >> t >> value;
        data_time_series[cnt].push_back(t);
        data_source[cnt].push_back(value);
        stream.str(std::string());
        stream.clear();		//clear the istringstream
      }
      return true;
    }
  }

}
