cmake_minimum_required(VERSION 3.10)
project(wave_eq CXX)

include_directories("${PROJECT_BINARY_DIR}/include")
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

if(NOT CMAKE_BUILD_TYPE) # force set build type in cache to debug
  set(CMAKE_BUILD_TYPE "Debug" CACHE STRING
      "Choose the type of build, options are: Debug Release
RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# force set gui setting
SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g -ggdb -D_DEBUG" CACHE STRING "Flags used by the linker during debug builds" FORCE)
SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall" CACHE STRING "Flags used by the linker during release builds" FORCE)

#MKL(BLAS) library
#set(MKL_LD "/opt" CACHE PATH "path to mkl library")
#set(MKL_INCLUDE "/opt" CACHE PATH "path to mkl include")
#include_directories(${MKL_INCLUDE})
#set(CXX_LD ${CXX_LD} ${MKL_LD})
#set(CXX_LIB ${CXX_LIB} mkl_intel_lp64 mkl_tbb_thread mkl_core tbb pthread m dl)

#MPI
find_package(MPI) 
set(METIS_INCLUDE "/opt" CACHE PATH "path to metis include")
set(METIS_LD "/opt" CACHE PATH "path to metis library")
set(PARMETIS_INCLUDE "/opt" CACHE PATH "path to parmetis include")
set(PARMETIS_LD "/opt" CACHE PATH "path to parmetis library")
if(${MPI_FOUND})
      include_directories(${MPI_CXX_INCLUDE_PATH} ${METIS_INCLUDE} ${PARMETIS_INCLUDE})
      set(CXX_LD ${CXX_LD} ${PARMETIS_LD} ${METIS_LD})
      set(CXX_LIB ${CXX_LIB} ${MPI_CXX_LIBRARIES} parmetis metis)
else()
      message(SEND_ERROR  "Cannot find MPI library, please specify manually.")
endif()

#source file list
aux_source_directory(${PROJECT_SOURCE_DIR}/src SRCLIST)

# use C++14
add_definitions(-std=c++14)

#build
LINK_DIRECTORIES(${CXX_LD})
add_executable(wave_eq ${SRCLIST})
target_link_libraries(wave_eq ${CXX_LIB})
