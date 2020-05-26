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
  \file cuda_data_bank.h
  \brief Header file for data bank for syncronizing data between GPUs
  \version 1.0
  \author xilin xia
*/

#ifndef CUDA_DATA_BANKS_H
#define CUDA_DATA_BANKS_H

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <atomic>
#include "Flag.h"
#include "Scalar.h"
#include "Vector.h"
#include "mapped_field.h"
#include "cuda_mapped_field.h"

namespace GC{

  struct cuDataBankBranch{
  public:
    cuDataBankBranch(unsigned int _size_lower2receive, unsigned int _size_lower2send, unsigned int _size_upper2receive, unsigned int _size_upper2send);
    thrust::device_vector<Vector3> data_lower;
    std::atomic<int> lower_status; // 0 -- empty 1 -- full
    thrust::host_vector<int> indices_lower2receive; //row near boundary
    thrust::host_vector<int> indices_lower2send;
    thrust::device_vector<Vector3> data_upper;
    std::atomic<int> upper_status;
    thrust::host_vector<int> indices_upper2receive; //row near boundary
    thrust::host_vector<int> indices_upper2send;
  private:
    unsigned int size_lower2receive;
    unsigned int size_lower2send;
    unsigned int size_upper2receive;
    unsigned int size_upper2send;
  };

  class cuDataBank{
  public:
    cuDataBank(const char* file_name, std::vector<int> device_list);
    unsigned int domains_size(){return num_domains;};
    typedef cuDataBankBranch* p_cuDataBankBranch;
    std::vector< p_cuDataBankBranch > banks;
  private:
    unsigned int num_domains;
  };

  template <typename T, MAPPING_MODES C>
  void CollectAndSend(cuFvMappedField<T, C>& phi, cuDataBank& bank, unsigned int domain_id, std::vector<int> device_list){
    //Collect
    //lower boundary
    int dev_id_src = device_list[domain_id];
    if (domain_id > 0){ //lowest bank does not have to send lower bank branch
      unsigned int lower_start_id = (*(bank.banks)[domain_id]).indices_lower2send.front();
      unsigned int lower_end_id = (*(bank.banks)[domain_id]).indices_lower2send.back();
      int dev_id_target = device_list[domain_id - 1];
      checkCuda(cudaSetDevice(dev_id_target));
      checkCuda(cudaDeviceEnablePeerAccess(dev_id_src, 0));
      checkCuda(cudaSetDevice(dev_id_src));
      checkCuda(cudaDeviceEnablePeerAccess(dev_id_target, 0));
      void* p_bank = thrust::raw_pointer_cast((*(bank.banks)[domain_id - 1]).data_upper.data()); //target: upper
      while (true){ //keep trying until the bank is empty
        if ((*(bank.banks)[domain_id - 1]).upper_status == 0){ //upper bank of a lower domain is empty
          checkCuda(cudaMemcpy(p_bank, phi.data.dev_ptr() + lower_start_id, (lower_end_id - lower_start_id + 1)*sizeof(T), cudaMemcpyDefault));
          //printf("sending data to lower domain.\n");
          (*(bank.banks)[domain_id - 1]).upper_status = 1; //change status to full
          break;
        }
      }
    }
    if (domain_id < bank.domains_size() - 1){ //up bank dose not have to send upper bank branch
      unsigned int upper_start_id = (*(bank.banks)[domain_id]).indices_upper2send.front();
      unsigned int upper_end_id = (*(bank.banks)[domain_id]).indices_upper2send.back();
      int dev_id_target = device_list[domain_id + 1];
      checkCuda(cudaSetDevice(dev_id_target));
      checkCuda(cudaDeviceEnablePeerAccess(dev_id_src, 0));
      checkCuda(cudaSetDevice(dev_id_src));
      checkCuda(cudaDeviceEnablePeerAccess(dev_id_target, 0));
      void* p_bank = thrust::raw_pointer_cast((*(bank.banks)[domain_id + 1]).data_lower.data()); //target:: lower
      while (true){ //keep trying until the bank is empty
        if ((*(bank.banks)[domain_id + 1]).lower_status == 0){ //lower bank of a upper domain is empty
          checkCuda(cudaMemcpy(p_bank, phi.data.dev_ptr() + upper_start_id, (upper_end_id - upper_start_id + 1)*sizeof(T), cudaMemcpyDefault));
          //printf("sending data to upper domain.\n");
          (*(bank.banks)[domain_id + 1]).lower_status = 1; //change status to full
          break;
        }
      }
    }
  };

  template <typename T, MAPPING_MODES C>
  void ReceiveAndDispatch(cuFvMappedField<T, C>& phi, cuDataBank& bank, unsigned int domain_id, std::vector<int> device_list){
    //Receive and dispatch
    if (domain_id < bank.domains_size() - 1){//receiving from the upper bank
      while (true){ //keep trying until the bank is full
        void* p_bank = thrust::raw_pointer_cast((*(bank.banks)[domain_id]).data_upper.data());
        unsigned int upper_start_id = (*(bank.banks)[domain_id]).indices_upper2receive.front();
        unsigned int upper_end_id = (*(bank.banks)[domain_id]).indices_upper2receive.back();
        if ((*(bank.banks)[domain_id]).upper_status == 1){ //lower bank of a upper domain is full
          checkCuda(cudaMemcpy(phi.data.dev_ptr() + upper_start_id, p_bank, (upper_end_id - upper_start_id + 1)*sizeof(T), cudaMemcpyDefault));
          //printf("receiving data from upper domain.\n");
          (*(bank.banks)[domain_id]).upper_status = 0; //change status to empty
          break;
        }
      }
    }
    if (domain_id > 0){//receiving from the lower bank
      void* p_bank = thrust::raw_pointer_cast((*(bank.banks)[domain_id]).data_lower.data());
      unsigned int lower_start_id = (*(bank.banks)[domain_id]).indices_lower2receive.front();
      unsigned int lower_end_id = (*(bank.banks)[domain_id]).indices_lower2receive.back();
      while (true){ //keep trying until the bank is full
        if ((*(bank.banks)[domain_id]).lower_status == 1){ //lower bank of a upper domain is full
          checkCuda(cudaMemcpy(phi.data.dev_ptr() + lower_start_id, p_bank, (lower_end_id - lower_start_id + 1)*sizeof(T), cudaMemcpyDefault));
          //printf("receiving data from lower domain.\n");
          (*(bank.banks)[domain_id]).lower_status = 0; //change status to empty
          break;
        }
      }
    }
  }

}

#endif