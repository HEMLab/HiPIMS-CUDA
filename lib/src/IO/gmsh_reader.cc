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
  \file gmsh_reader.cc
  \brief Source file for gmsh file reader class

*/


#include "gmsh_reader.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include "Flag.h"
#include "Vector.h"

namespace GC{

  ///constructor
  gmshReader::gmshReader(const char* filename){
    readin(filename);
  }


  //This function reads in element connectivity table from gmsh type file
  void gmshReader::readin(const char* filename){
    std::ifstream input;
    input.open(filename);
    if(!input){
      std::cout<<"error: unable to open input file: "<<filename<<std::endl;
    }
    std::string line;
    std::string word;
    //reading the head
    getline(input, line);
    getline(input, line);
    getline(input, line);
    //reading nodes
    getline(input, line);
    size_t nNodes;
    getline(input, line);
    std::istringstream stream;
    stream.str(line);
    stream>>nNodes;
    stream.str(std::string());
    stream.clear();
    Vector coord;
    for(size_t i = 0; i < nNodes; i++){
      getline(input, line);
      stream.str(line);
      stream>>word>>coord;
      vertex_positions.push_back(coord);
      stream.str(std::string());
      stream.clear();		//clear the istringstream
    }
    getline(input, line);
    //reading elements
    getline(input, line);
    size_t nElements;
    getline(input, line);
    stream.str(line);
    stream>>nElements;
    stream.str(std::string());
    stream.clear();
    element_vertices.resize(nElements);	//allocate memory for all the elements
    Flag eType, vertexID;
    for(size_t i = 0; i <nElements; i++){
      getline(input, line);
      std::istringstream steam;
      stream.str(line);
      stream>>word>>eType;
      size_t nTags;
      stream>>nTags;
      for(size_t j = 0; j < nTags; j++){
        stream>>word;
      }
      element_types.push_back(eType);
      while(stream>>vertexID){
        element_vertices[i].push_back(vertexID - 1); //reduce 1 here
      }
      stream.str(std::string());
      stream.clear();
   //   std::cout<<element_vertices[i][0]<<element_vertices[i][1]<<element_vertices[i][2]<<std::endl;
    }
    getline(input, line);
  }	

}
