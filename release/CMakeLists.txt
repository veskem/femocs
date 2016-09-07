
cmake_minimum_required (VERSION 2.6)
project (femocs C CXX Fortran)
SET(CMAKE_Fortran_COMPILER gfortran-5)

# specify the name of the library
SET(LIB_NAME "femocs")

# Specify folders where libfemocs and its .mod-file will be written
# Those paths are (not like others) specified in relation to the location of Makefile
SET (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ../../lib)
SET (CMAKE_Fortran_MODULE_DIRECTORY ../../include)

# Specify the path to deal.II files
SET(DEAL_II_DIR ../deal.II)

# Find source files
FILE(GLOB_RECURSE TARGET_SRC "../src/*")

# Add the directory with headers to search path
INCLUDE_DIRECTORIES(../include)

# specify whether to build in debug or release mode
SET(CMAKE_BUILD_TYPE Release)
#SET(CMAKE_BUILD_TYPE Debug)

# Specify warnings and optimization level
SET(CMAKE_CXX_FLAGS_RELEASE "-w -O3")
#SET(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wpedantic -O2")
SET(CMAKE_CXX_FLAGS_DEBUG "-g2 -w -O0")

FIND_PACKAGE(deal.II 8.3 REQUIRED HINTS ${DEAL_II_DIR})
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
