// ======================================================================================
// Name                :    GeoClasses : Generic Geophysical Flow Modelling Framework
// Description         :    This code pack provides a generic framework for developing 
//                          Geophysical CFD software.
// ======================================================================================
// Version             :    0.1 
// Author              :    Xilin Xia (PhD candidate in Newcastle University)
// Create Time         :    2014/10/04
// Update Time         :    2016/04/14
// ======================================================================================
// Copyright @ Xilin Xia 2015 . All rights reserved.
// ======================================================================================

/*!
\file cuda_device_query.cu
\brief Header file for device query function

\version 0.1
\author xilin xia
*/

#include "cuda_device_query.h"
#include <cuda_runtime_api.h>
#include <cstdio>
#include <iostream>

namespace GC{

  void deviceQuery(){
    printf("Device Querying...\n");
    int dev_count;
    cudaGetDeviceCount(&dev_count);
    cudaDeviceProp prop;
    std::cout << "---General Information for device---" << std::endl;
    for (int i = 0; i < dev_count; i++){
      cudaGetDeviceProperties(&prop, i);
      std::cout << "Number:   " << i << "   Device Model:   " << prop.name << std::endl;
    }

  choose_device:

    std::cout << "Choosing device, please type in device ID and press ENTER:";
    int chosenID;
    std::cin >> chosenID;
    if (chosenID >= dev_count || chosenID < 0){
      std::cout << "Error: invalid ID, please chose again." << std::endl;
      goto choose_device;
    }
    cudaSetDevice(chosenID);
    std::cout << "Device   " << chosenID << "   is chosen." << std::endl;
  }

}