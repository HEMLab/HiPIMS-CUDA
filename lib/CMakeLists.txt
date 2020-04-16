cmake_minimum_required(VERSION 2.8)

# nvcc flags                                                                                                             

#set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-arch=sm_20)
#set (CUDA_NVCC_FLAGS "-arch=sm_20" CACHE STRING "nvcc flags" FORCE)

file(GLOB_RECURSE gc_srcs "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc")
file(GLOB_RECURSE gc_headers "${CMAKE_CURRENT_SOURCE_DIR}/include/*.h")

set (GC_INCLUDE_DIRS "")
foreach (_header_file ${gc_headers})
	get_filename_component(_dir ${_header_file} DIRECTORY)
	list (APPEND GC_INCLUDE_DIRS ${_dir})
endforeach()
list(REMOVE_DUPLICATES GC_INCLUDE_DIRS)

file(GLOB_RECURSE gc_cuda_srcs "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cu")

#set(GC_INCLUDE_DIRS${GC_INCLUDE_DIRS} PARENT_SCOPE)

if(TARGET lib)

message("lib is already defined")

else()

include_directories(${GC_INCLUDE_DIRS})
add_library(lib STATIC ${gc_srcs})
#cuda_add_library(cuda_lib STATIC ${gc_cuda_srcs})
cuda_add_library(cuda_lib  STATIC ${gc_cuda_srcs})

endif()
