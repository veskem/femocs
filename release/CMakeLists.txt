CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

# specify the path to deal.II files
SET(DEAL_II_DIR ../dealii)

# specify the name of the final executable
SET(TARGET "femocs")

# Specify the path for the location of executable
SET (CMAKE_RUNTIME_OUTPUT_DIRECTORY ..)

# path to source files
FILE(GLOB_RECURSE SOURCE "../src/*.cpp")

# Add the directory with headers to search path
INCLUDE_DIRECTORIES(../include ../tetgen)

# path to libraries
LINK_DIRECTORIES(../lib/)

# specify the names of the libraries
SET(LIBRARIES tet)

# specify whether to build in debug or release mode
SET(CMAKE_BUILD_TYPE Release)
#SET(CMAKE_BUILD_TYPE Debug)

# disable or enable all warnings
SET(CMAKE_CXX_FLAGS_RELEASE "-w -O3")
SET(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wpedantic -O2")
# ==============================================

FIND_PACKAGE(deal.II 8.4 REQUIRED
  HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\nCould not locate a (sufficiently recent)"
    "version of deal.II from" ${DEAL_II_DIR} "\n")
ENDIF()

FIND_PACKAGE(OpenMP)
IF(OPENMP_FOUND)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()
ADD_EXECUTABLE(${TARGET} ../Standalone.cpp ${SOURCE})
TARGET_LINK_LIBRARIES(${TARGET} ${LIBRARIES})
DEAL_II_SETUP_TARGET(${TARGET})
