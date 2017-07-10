
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

# ==============================================
# ======= Variables to be set on cmd-line ======
# ==============================================

## Set compilers
#SET(CMAKE_C_COMPILER gcc)
#SET(CMAKE_CXX_COMPILER c++)

## Specify whether to build in debug or release mode
# SET(CMAKE_BUILD_TYPE Release)
# SET(CMAKE_BUILD_TYPE Debug)

# ==============================================
# ======= Name & path for input & output =======
# ==============================================

# name of the final executable
SET(TARGET "femocs")

# path for the location of executable
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ..)

# path to deal.II files
SET(DEAL_II_DIR ../dealii)

# path to CGAL files
SET(CGAL_DIR ../cgal)

# path to file with main function
#SET(MAIN "../Test.cpp")

# path to source files
FILE(GLOB_RECURSE SOURCE "../src/*.cpp" "../deal-solver/source/*")

# path to the directories with headers
IF (LIBCGAL AND CMAKE_BUILD_TYPE MATCHES Release)
    INCLUDE_DIRECTORIES(../include ../deal-solver/include ../dealii/include ../cgal/include ../GETELEC/modules)
ELSE()
    INCLUDE_DIRECTORIES(../include ../deal-solver/include ../dealii/include ../GETELEC/modules)
ENDIF()

# path to libraries
LINK_DIRECTORIES(../lib ../GETELEC/lib)

# names of the libraries
SET(LIBRARIES tet getelec slatec)

# ==============================================
# =========== Locate dependencies ==============
# ==============================================

FIND_PACKAGE(deal.II 8.4 REQUIRED HINTS ${DEAL_II_DIR})
IF ( NOT deal.II_FOUND )
    MESSAGE(FATAL_ERROR "Deal.II not found from " ${DEAL_II_DIR})
ENDIF()

IF (LIBCGAL AND CMAKE_BUILD_TYPE MATCHES Release)
    FIND_PACKAGE( CGAL HINTS ${CGAL_DIR} )
    IF (CGAL_FOUND)
        MESSAGE("Using CGAL from " ${CGAL_DIR})
        INCLUDE( ${CGAL_USE_FILE} )
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_CGAL=true")
    ELSE()
        MESSAGE(FATAL_ERROR "CGAL not found from " ${CGAL_DIR}) 
    ENDIF()
ELSE()
    MESSAGE(STATUS "Not using CGAL!")
ENDIF()

FIND_PACKAGE( Boost REQUIRED )
IF ( NOT Boost_FOUND )
    MESSAGE(FATAL_ERROR "Boost not found") 
ENDIF()
    
FIND_PACKAGE(OpenMP)
IF ( NOT OPENMP_FOUND )
    MESSAGE(FATAL_ERROR "OpenMP not found") 
ENDIF()

# ==============================================
# ============= Compilation flags ==============
# ==============================================

# Specify warnings and optimization level
# Suppressing individual warnings does not work because they are also controlled by deal.II and CGAL cmake-s

#SET(warnings "-Wno-overloaded-virtual -Wno-type-limits -Wno-reorder -Wno-pedantic -Wno-unused-value -Wno-unused-result -Wno-unused-variable -Wno-unused-parameter -Wno-sign-compare")
SET(warnings "-Wall -Wpedantic")

SET(CMAKE_CXX_FLAGS_RELEASE "-O3 -w")
SET(CMAKE_CXX_FLAGS_DEBUG "-g -O0 ${warnings}")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

# ==============================================
# =============== Start building ===============
# ==============================================

DEAL_II_INITIALIZE_CACHED_VARIABLES()
ADD_EXECUTABLE(${TARGET} ${MAIN} ${SOURCE})
TARGET_LINK_LIBRARIES(${TARGET} ${LIBRARIES})
DEAL_II_SETUP_TARGET(${TARGET})
