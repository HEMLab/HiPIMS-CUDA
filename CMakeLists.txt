cmake_minimum_required(VERSION 3.10)

project(HiPIMs LANGUAGES CXX CUDA)

set(CMAKE_CUDA_STANDARD 14)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_ARCHITECTURES 80)
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --extended-lambda")
set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()

# Enable debug symbols across all platforms
IF(WIN32)
    set(PLATFORM_CXX_FLAGS /Zi)
ENDIF(WIN32)

IF(UNIX)
    set(PLATFORM_CXX_FLAGS "-pthread -std=c++11 -O3")
ENDIF(UNIX)

list(APPEND CMAKE_CXX_FLAGS ${PLATFORM_CXX_FLAGS})

message("CXX  Flags: " ${CMAKE_CXX_FLAGS})
message("NVCC Flags: " ${CUDA_NVCC_FLAGS})

add_subdirectory(lib)
add_subdirectory(apps/cudaFloodSolver)
add_subdirectory(apps/cudaMultiGPUsFloodSolver)
add_subdirectory(apps/cudaIntegratedShallowFlowSolver)
