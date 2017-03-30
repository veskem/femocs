
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)
#SET(CMAKE_CXX_COMPILER gcc)

# path to deal.II files
SET(DEAL_II_DIR ../dealii)

# path to CGAL files
SET(CGAL_DIR ../cgal)

# path for the location of executable
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ..)

# name of the final executable
SET(TARGET "femocs")

# path to file with main function
#SET(MAIN "../Test.cpp")

# path to source files
FILE(GLOB_RECURSE SOURCE "../src/*.cpp" "../heating/source/*")

# directories with headers
IF (USE_CGAL)
    INCLUDE_DIRECTORIES(../include ../heating/include ../dealii/include ../cgal/include)
ELSE()
    INCLUDE_DIRECTORIES(../include ../heating/include ../dealii/include)
ENDIF()

# path to libraries
LINK_DIRECTORIES(../lib)

# names of the libraries
SET(LIBRARIES tet)

## Specify whether to build in debug or release mode
# SET(CMAKE_BUILD_TYPE Release)
# SET(CMAKE_BUILD_TYPE Debug)

# ==============================================

FIND_PACKAGE(deal.II 8.4 REQUIRED HINTS ${DEAL_II_DIR})
IF ( NOT deal.II_FOUND )
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

#SET(warnings "") #"-Wno-overloaded-virtual -Wno-vla -Wno-reorder -Wno-pedantic -Wno-unused-value -Wno-unused-result -Wno-unused-variable -Wno-unused-parameter -Wno-sign-compare")

# Specify warnings and optimization level
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -w")
SET(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wpedantic")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -Dvariable=true")

# ==============================================

DEAL_II_INITIALIZE_CACHED_VARIABLES()
ADD_EXECUTABLE(${TARGET} ${MAIN} ${SOURCE})
TARGET_LINK_LIBRARIES(${TARGET} ${LIBRARIES})
DEAL_II_SETUP_TARGET(${TARGET})
