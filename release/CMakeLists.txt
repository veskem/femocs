CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

# path to deal.II files
SET(DEAL_II_DIR ../dealii)

# path for the location of executable
SET (CMAKE_RUNTIME_OUTPUT_DIRECTORY ..)

# name of the final executable
SET(TARGET "femocs")

# path to file with main function
#SET(MAIN "../Test.cpp")

# path to source files
FILE(GLOB_RECURSE SOURCE "../src/*.cpp" "../heating/source/*")

# directories with headers
INCLUDE_DIRECTORIES(../include ../heating/include)

# path to libraries
LINK_DIRECTORIES(../lib)

# names of the libraries
SET(LIBRARIES tet)

# Specify whether to build in debug or release mode
#SET(CMAKE_BUILD_TYPE Release)
#SET(CMAKE_BUILD_TYPE Debug)

# Specify warnings and optimization level
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -w")
SET(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -w") # "-Wall -Wpedantic")
# ==============================================

FIND_PACKAGE(deal.II 8.4 REQUIRED HINTS ${DEAL_II_DIR})
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\nCould not find deal.II from" ${DEAL_II_DIR} "\n")
ENDIF()

FIND_PACKAGE(OpenMP)
IF(OPENMP_FOUND)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()
ADD_EXECUTABLE(${TARGET} ${MAIN} ${SOURCE})
TARGET_LINK_LIBRARIES(${TARGET} ${LIBRARIES})
DEAL_II_SETUP_TARGET(${TARGET})
