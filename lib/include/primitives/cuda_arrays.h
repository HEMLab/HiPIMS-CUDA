// ======================================================================================
// Name                :    GeoClasses : Generic Geophysical Flow Modelling Framework
// Description         :    This code pack provides a generic framework for developing 
//                          Geophysical CFD software.
// ======================================================================================
// Version             :    0.2 
// Author              :    Xilin Xia (PhD candidate in Newcastle University)
// Create Time         :    2014/10/04
// Update Time         :    2015/10/08
// ======================================================================================
// Copyright @ Xilin Xia 2015 . All rights reserved.
// ======================================================================================

/*!
 \file cuda_arrays.h
 \brief Header file for cuda array class 

*/

#ifndef CUDA_ARRAYS_H
#define CUDA_ARRAYS_H

#include "mesh_fv_entities.h"
#include <cuda_runtime_api.h>
#include <cuda.h>
#include <algorithm>
#include "hemi_error.h"
#include "Flag.h"
#include "Scalar.h"
#include "Vector.h"

namespace GC{


  template <typename T>
  class cuArray{
    public:
      cuArray(unsigned int _size = 0):size_(_size){
        data_host = new T[size_];
        checkCuda(cudaMalloc((void**)&data_dev, size_*sizeof(T)));
        checkCuda(cudaMemset(data_dev, 0, size_*sizeof(T)));
      };
      void resize(unsigned int new_size){
        if(size_ != new_size){
          T *tmp;
          checkCuda(cudaMalloc((void**)&tmp, new_size*sizeof(T)));
          checkCuda(cudaMemset(tmp, 0, new_size*sizeof(T)));
          if(size_ <= new_size){
            checkCuda(cudaMemcpy(tmp, data_dev, size_*sizeof(T), cudaMemcpyDeviceToDevice));
          }else{
            checkCuda(cudaMemcpy(tmp, data_dev, new_size*sizeof(T), cudaMemcpyDeviceToDevice));
          } 
          checkCuda(cudaFree(data_dev));
          data_dev = tmp;
          size_ = new_size;
          delete data_host;
          data_host = new T[size_];
          checkCuda(cudaMemcpy(data_host, data_dev, size_*sizeof(T), cudaMemcpyDeviceToHost));
        }
      }
      void initialize_from_host(typename std::vector<T>::iterator begin, typename std::vector<T>::iterator end){
        unsigned int host_size = end - begin;
        resize(host_size);
        for(auto iter = begin; iter < end; ++iter){
          auto id = iter - begin;
          data_host[id] = *iter;
        }
        checkCuda(cudaMemcpy(data_dev, data_host, host_size*sizeof(T), cudaMemcpyHostToDevice));
      }
      void sync(){
        checkCuda(cudaMemcpy(data_host, data_dev, size_*sizeof(T), cudaMemcpyDeviceToHost));
      }

      void sync_by_host(){
        checkCuda(cudaMemcpy(data_dev, data_host, size_*sizeof(T), cudaMemcpyHostToDevice));
      }
      ///assignment operator
      cuArray& operator= (const cuArray& rhs){
        if(this != &rhs){
          resize(rhs.size());
          checkCuda(cudaMemcpy(data_dev, rhs.dev_ptr(), size_*sizeof(T), cudaMemcpyDeviceToDevice));
          checkCuda(cudaMemcpy(data_host, rhs.host_ptr(), size_*sizeof(T), cudaMemcpyHostToHost));
        }
        return *this;
      }
      ///copy constructor
      cuArray(const cuArray& other){
        resize(other.size());
        checkCuda(cudaMemcpy(data_dev, other.dev_ptr(), size_*sizeof(T), cudaMemcpyDeviceToDevice));
        checkCuda(cudaMemcpy(data_host, other.host_ptr(), size_*sizeof(T), cudaMemcpyHostToHost));
      }
    public:
      T* dev_ptr(){return data_dev;}
      T* host_ptr(){return data_host;}
      T* dev_ptr() const{return data_dev;}
      T* host_ptr() const{return data_host;}
      unsigned int size(){return size_;}
      unsigned int size() const{return size_;}
    private:
      unsigned int size_;
      T *data_dev;
      T *data_host;
    public:
      ~cuArray(){
        delete[] data_host;
        checkCuda(cudaFree(data_dev));
      }
  };

  template <typename T>
  class cu2dArray{
    public:
      cu2dArray(unsigned int _length = 0, unsigned int _width = 0):length_(_length), width_(_width){
        dimension_list_host = new Flag[length_];
        data_host = new T[length_*width_];
        checkCuda(cudaMalloc((void**)&dimension_list_dev, length_*sizeof(Flag)));
        checkCuda(cudaMalloc((void**)&data_dev, length_*width_*sizeof(T)));
        checkCuda(cudaMemset(dimension_list_dev, 0, length_*sizeof(Flag)));
        checkCuda(cudaMemset(data_dev, 0, length_*width_*sizeof(T)));
      };

      void resize(unsigned int new_length, unsigned int new_width){
        T *tmp;
        checkCuda(cudaMalloc((void**)&tmp, new_length*new_width*sizeof(T)));
        checkCuda(cudaMemset(tmp, 0, new_length*new_width*sizeof(T)));
        if(length_ <= new_length){
          if(width_ <= new_width){
            for(unsigned int i = 0; i < width_; ++i){
              checkCuda(cudaMemcpy(tmp+i*new_length, data_dev+i*length_, length_*sizeof(T), cudaMemcpyDeviceToDevice));
            }
          }else{
            for(unsigned int i = 0; i < new_width; ++i){
              checkCuda(cudaMemcpy(tmp+i*new_length, data_dev+i*length_, length_*sizeof(T), cudaMemcpyDeviceToDevice));
            }
          }
        }else{
          if(width_ <= new_width){
            for(unsigned int i = 0; i < width_; ++i){
              checkCuda(cudaMemcpy(tmp+i*new_length, data_dev+i*length_, new_length*sizeof(T), cudaMemcpyDeviceToDevice));
            }
          }else{
            for(unsigned int i = 0; i < new_width; ++i){
              checkCuda(cudaMemcpy(tmp+i*new_length, data_dev+i*length_, new_length*sizeof(T), cudaMemcpyDeviceToDevice));
            }
          }
        }
        if(length_ != new_length){
          Flag *tmp;
          checkCuda(cudaMalloc((void**)&tmp, new_length*sizeof(Flag)));
          checkCuda(cudaMemset(tmp, 0, new_length*sizeof(Flag)));
          if(length_ <= new_length){
            checkCuda(cudaMemcpy(tmp, dimension_list_dev, length_*sizeof(Flag), cudaMemcpyDeviceToDevice));
          }else{
            checkCuda(cudaMemcpy(tmp, dimension_list_dev, new_length*sizeof(Flag), cudaMemcpyDeviceToDevice));
          } 
          checkCuda(cudaFree(dimension_list_dev));
          dimension_list_dev = tmp;
          delete dimension_list_host;
          dimension_list_host = new Flag[new_length];
          checkCuda(cudaMemcpy(dimension_list_host, dimension_list_dev, new_length*sizeof(Flag), cudaMemcpyDeviceToHost));
          for (unsigned int i = 0; i < new_length; ++i){
            unsigned int dim = dimension_list_host[i];
            dimension_list_host[i] = std::min(dim, new_width);
          }
          checkCuda(cudaMemcpy(dimension_list_dev, dimension_list_host, new_length*sizeof(Flag), cudaMemcpyHostToDevice));
        }
        checkCuda(cudaFree(data_dev));
        data_dev = tmp;
        length_ = new_length;
        width_ = new_width;
        delete data_host;
        data_host = new T[length_*width_];
        checkCuda(cudaMemcpy(data_host, data_dev, length_*width_*sizeof(T), cudaMemcpyDeviceToHost));
      }

      void initialize_from_host(typename std::vector< typename std::vector<T> >::iterator begin, typename std::vector< typename std::vector<T> >::iterator end){
        unsigned int host_length = end - begin;
        unsigned int host_width = 0;
        Flag *dimension_list_tmp = new Flag[host_length];

        for(auto iter = begin; iter < end; ++iter){
          auto id = iter - begin;
          unsigned int width_tmp = (*iter).size();
          dimension_list_tmp[id] = width_tmp;
          host_width = std::max(host_width, width_tmp);
        }

        T *data_tmp = new T[host_length*host_width];

        for(auto x_iter = begin; x_iter < end; ++x_iter){
          unsigned int x_id = x_iter - begin;
          unsigned int y_id = 0;
          for(auto y_iter : *x_iter){
            data_tmp[y_id*host_length + x_id] = y_iter;
            y_id++;
          }
        }

        resize(host_length, host_width);
        checkCuda(cudaMemcpy(dimension_list_dev, dimension_list_tmp, host_length*sizeof(Flag), cudaMemcpyHostToDevice));
        checkCuda(cudaMemcpy(data_dev, data_tmp, host_length*host_width*sizeof(T), cudaMemcpyHostToDevice));
        checkCuda(cudaMemcpy(dimension_list_host, dimension_list_tmp, host_length*sizeof(Flag), cudaMemcpyHostToHost));
        checkCuda(cudaMemcpy(data_host, data_tmp, host_length*host_width*sizeof(T), cudaMemcpyHostToHost));

        delete dimension_list_tmp;
        delete data_tmp;
      }

      void sync(){
        checkCuda(cudaMemcpy(dimension_list_host, dimension_list_dev, length_*sizeof(Flag), cudaMemcpyDeviceToHost));
        checkCuda(cudaMemcpy(data_host, data_dev, length_*width_*sizeof(T), cudaMemcpyDeviceToHost));
      }
      ///assignment operator
      cu2dArray& operator= (const cu2dArray& rhs){
        if(this != &rhs){
          resize(rhs.length(), rhs.width());
          checkCuda(cudaMemcpy(dimension_list_dev, rhs.dims_dev_ptr(), length_*sizeof(Flag), cudaMemcpyDeviceToDevice));
          checkCuda(cudaMemcpy(dimension_list_host, rhs.dims_host_ptr(), length_*sizeof(Flag), cudaMemcpyHostToHost));
          checkCuda(cudaMemcpy(data_dev, rhs.dev_ptr(), length_*width_*sizeof(T), cudaMemcpyDeviceToDevice));
          checkCuda(cudaMemcpy(data_host, rhs.host_ptr(), length_*width_*sizeof(T), cudaMemcpyHostToHost));
        }
        return *this;
      }
      ///copy constructor
      cu2dArray(const cu2dArray& other){
        resize(other.length(), other.width());
        checkCuda(cudaMemcpy(dimension_list_dev, other.dims_dev_ptr(), length_*sizeof(Flag), cudaMemcpyDeviceToDevice));
        checkCuda(cudaMemcpy(dimension_list_host, other.dims_host_ptr(), length_*sizeof(Flag), cudaMemcpyHostToHost));
        checkCuda(cudaMemcpy(data_dev, other.dev_ptr(), length_*width_*sizeof(T), cudaMemcpyDeviceToDevice));
        checkCuda(cudaMemcpy(data_host, other.host_ptr(), length_*width_*sizeof(T), cudaMemcpyHostToHost));
      }

    public:
      unsigned int length() {return length_;}
      unsigned int width() {return width_;}
      unsigned int length() const {return length_;}
      unsigned int width() const {return width_;}
    private:
      unsigned int length_;
      unsigned int width_;
    public:
      T* dev_ptr(){return data_dev;}
      T* host_ptr(){return data_host;}
      Flag* dims_dev_ptr(){return dimension_list_dev;}
      Flag* dims_host_ptr(){return dimension_list_host;}
      T* dev_ptr() const{return data_dev;}
      T* host_ptr() const{return data_host;}
      Flag* dims_dev_ptr() const{return dimension_list_dev;}
      Flag* dims_host_ptr() const{return dimension_list_host;}
    private:
      T *data_dev;
      T *data_host;
      Flag *dimension_list_dev;
      Flag *dimension_list_host;
    public:
      ~cu2dArray(){
        delete[] data_host;
        delete[] dimension_list_host;
        checkCuda(cudaFree(data_dev));
        checkCuda(cudaFree(dimension_list_dev));
      }
  };

} //--end of namespace GC------------


#endif
