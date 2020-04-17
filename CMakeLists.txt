cmake_minimum_required(VERSION 2.6)
#project name
project(HiPIMs)

find_package(CUDA)

IF(WIN32)
    set(PLATFORM_CXX_FLAGS /Zi)
    set(CUDA_SEPARABLE_COMPILATION ON)
    set(CUDA_NVCC_FLAGS -arch=sm_35;--expt-extended-lambda)
ENDIF(WIN32)

IF(UNIX)
    set(PLATFORM_CXX_FLAGS "-std=c++11 -O3")
    #set(CUDA_SEPARABLE_COMPILATION ON)
    set(CUDA_NVCC_FLAGS  -arch=sm_35;--expt-extended-lambda)
ENDIF(UNIX)

list(APPEND CMAKE_CXX_FLAGS ${PLATFORM_CXX_FLAGS})

message("CXX  Flags: " ${CMAKE_CXX_FLAGS})
message("NVCC Flags: " ${CUDA_NVCC_FLAGS})


add_subdirectory(lib)
add_subdirectory(apps/UrbanFloodSimulator)
add_subdirectory(apps/MultiGPUsUrbanFloodSimulator)
add_subdirectory(apps/cudaIntegratedShallowFlowSolver)
