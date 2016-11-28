
cmake_minimum_required (VERSION 2.6)
project (femocs C CXX Fortran)
SET(CMAKE_Fortran_COMPILER gfortran-5)

# specify the name of the library
SET(LIB_NAME "femocs")

# Specify folders where libfemocs and its .mod-file will be written
# Those paths are (not like others) specified in relation to the location of Makefile
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ..)
SET(CMAKE_Fortran_MODULE_DIRECTORY ..)

# Specify the path to deal.II files
SET(DEAL_II_DIR ../dealii)

# Find source files
FILE(GLOB_RECURSE TARGET_SRC "../src/*.cpp" "../src/*.f90" "../heating/source/*")

# Add the directory with headers to search path
INCLUDE_DIRECTORIES(../include ../heating/include)

# specify whether to build in debug or release mode
SET(CMAKE_BUILD_TYPE Release)
#SET(CMAKE_BUILD_TYPE Debug)

# Specify warnings and optimization level
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -w")
SET(CMAKE_CXX_FLAGS_DEBUG "-g2 -O0 -w") #-Wall -Wpedantic")

FIND_PACKAGE(deal.II 8.4 REQUIRED HINTS ${DEAL_II_DIR})
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\nCould not find deal.II from" ${DEAL_II_DIR} "\n")
ENDIF()

FIND_PACKAGE(OpenMP)
IF(OPENMP_FOUND)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()
ADD_LIBRARY (${LIB_NAME} STATIC ${TARGET_SRC})
DEAL_II_SETUP_TARGET (${LIB_NAME})
