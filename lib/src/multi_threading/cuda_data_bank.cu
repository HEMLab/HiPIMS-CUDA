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
  \file cuda_data_bank.cu
  \brief Source file for data bank for syncronizing data between GPUs

  \version 1.0
  \author xilin xia

*/


#include "cuda_data_bank.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include "hemi_error.h"

namespace GC{

  cuDataBankBranch::cuDataBankBranch(unsigned int _size_lower2receive, unsigned int _size_lower2send, unsigned int _size_upper2receive, unsigned int _size_upper2send)
    :size_lower2receive(_size_lower2receive), size_lower2send(_size_lower2send), size_upper2receive(_size_upper2receive), size_upper2send(_size_upper2send){
    lower_status = 0;
    data_lower.resize(size_lower2receive);
    indices_lower2receive.resize(size_lower2receive);
    indices_lower2send.resize(size_lower2send);
    upper_status = 0;
    data_upper.resize(size_upper2receive);
    indices_upper2receive.resize(size_upper2receive);
    indices_upper2send.resize(size_upper2send);
  }

  cuDataBank::cuDataBank(const char* file_name, std::vector<int> device_list){
    std::ifstream input;
    input.open(file_name);
    if (!input){
      std::cout << "error: unable to open input file: " << file_name << std::endl;
      return;
    }
    std::string line;
    std::string word;
    //reading elements
    getline(input, line);
    unsigned int _num_domains;
    getline(input, line);
    std::istringstream stream;
    stream.str(line);
    stream >> _num_domains;
    num_domains = _num_domains;
    stream.str(std::string());
    stream.clear();
    if (device_list.size() < num_domains){
      std::cout << "Fatal error: number of devices is smaller than number of domains" << std::endl;
      return;
    }
    banks.resize(num_domains);
    for (unsigned int i = 0; i < num_domains; i++){
      int dev_id = device_list[i];
      checkCuda(cudaSetDevice(dev_id));
      thrust::host_vector<unsigned int> _indices_lower2receive;
      thrust::host_vector<unsigned int> _indices_lower2send;
      thrust::host_vector<unsigned int> _indices_upper2receive;
      thrust::host_vector<unsigned int> _indices_upper2send;
      getline(input, line);
      //line 1
      getline(input, line);
      stream.str(line);
      unsigned int index;
      while (stream >> index){
        _indices_lower2receive.push_back(index);
      }
      stream.str(std::string());
      stream.clear();
      //line 2
      getline(input, line);
      stream.str(line);
      while (stream >> index){
        _indices_lower2send.push_back(index);
      }
      stream.str(std::string());
      stream.clear();
      //line 3
      getline(input, line);
      stream.str(line);
      while (stream >> index){
        _indices_upper2receive.push_back(index);
      }
      stream.str(std::string());
      stream.clear();
      //line 4
      getline(input, line);
      stream.str(line);
      while (stream >> index){
        _indices_upper2send.push_back(index);
      }
      stream.str(std::string());
      stream.clear();
      unsigned int size_lower2receive = _indices_lower2receive.size();
      unsigned int size_lower2send = _indices_lower2send.size();
      unsigned int size_upper2receive = _indices_upper2receive.size();
      unsigned int size_upper2send = _indices_upper2send.size();
      banks[i] = new cuDataBankBranch(size_lower2receive, size_lower2send, size_upper2receive, size_upper2send);
      (*banks[i]).indices_lower2receive = _indices_lower2receive;
      (*banks[i]).indices_lower2send = _indices_lower2send;
      (*banks[i]).indices_upper2receive = _indices_upper2receive;
      (*banks[i]).indices_upper2send = _indices_upper2send;
    }

  }

}
