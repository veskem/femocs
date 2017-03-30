
cmake_minimum_required (VERSION 2.8.8)
project (femocs C CXX Fortran)
#SET(CMAKE_Fortran_COMPILER gfortran)
SET(CMAKE_CXX_COMPILER gcc)

# specify the name of the library
SET(LIB_NAME "femocs")

# Specify folders where libfemocs and its .mod-file will be written
# Those paths are (not like others) specified in relation to the location of Makefile
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ..)
SET(CMAKE_Fortran_MODULE_DIRECTORY ..)

# Specify the path to deal.II files
SET(DEAL_II_DIR ../dealii)

# Specify the path to CGAL files
SET(CGAL_DIR ../cgal)

# Find source files
FILE(GLOB_RECURSE TARGET_SRC "../src/*.cpp" "../src/*.f90" "../heating/source/*")

# Add the directory with headers to search path
IF (USE_CGAL)
    INCLUDE_DIRECTORIES(../include ../heating/include ../cgal/include)
ELSE()
    INCLUDE_DIRECTORIES(../include ../heating/include)
ENDIF()

# specify whether to build in debug or release mode
#SET(CMAKE_BUILD_TYPE Release)
#SET(CMAKE_BUILD_TYPE Debug)

# ==============================================

FIND_PACKAGE(deal.II 8.4 REQUIRED HINTS ${DEAL_II_DIR})
IF ( NOT ${deal.II_FOUND} )
  MESSAGE(FATAL_ERROR "Deal.II not found from " ${DEAL_II_DIR})
ENDIF()

IF (USE_CGAL)
    FIND_PACKAGE( CGAL HINTS ${CGAL_DIR} )
    IF (CGAL_FOUND)
        MESSAGE("Using CGAL from " ${CGAL_DIR}) 
    ELSE()
        MESSAGE(FATAL_ERROR "CGAL not found") 
    ENDIF()
ELSE()
    MESSAGE("Ignoring CGAL!")
ENDIF()

FIND_PACKAGE( Boost REQUIRED )
IF ( NOT Boost_FOUND )
    MESSAGE(FATAL_ERROR "Boost not found") 
ENDIF()

FIND_PACKAGE(OpenMP)
IF ( NOT OPENMP_FOUND )
    MESSAGE(FATAL_ERROR "OpenMP not found") 
ENDIF()

IF (USE_CGAL)
    INCLUDE( ${CGAL_USE_FILE} )
ENDIF()

# Specify warnings and optimization level
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -w")
SET(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -w") # "-Wall -Wpedantic")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -Dvariable=true")

# ==============================================

DEAL_II_INITIALIZE_CACHED_VARIABLES()
ADD_LIBRARY (${LIB_NAME} STATIC ${TARGET_SRC})
DEAL_II_SETUP_TARGET (${LIB_NAME})
