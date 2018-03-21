#############################################################################
#           CMake file for testing system dependencies of Femocs
#              In case something is missing, the suggestion 
#              for installing it in Debian/Ubuntu is provided
#############################################################################

cmake_minimum_required (VERSION 2.8.8)

enable_language (Fortran)

# make FindXYZ.cmake files visible for CMake
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})

#############################################################################
# check for versions and availability of various modules
#############################################################################
# c/c++ compiler
IF (NOT CMAKE_COMPILER_IS_GNUCC OR
    (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.9) OR 
    (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 5.9))
  MESSAGE(FATAL_ERROR
          "You are using unsupported compiler with version ${CMAKE_CXX_COMPILER_VERSION}\n"
          "Please install gcc-5 by running\n"
          "sudo apt install gcc-5")
ENDIF()

# fortran compiler
IF (NOT (${CMAKE_Fortran_COMPILER} MATCHES "gfortran"))
  MESSAGE(FATAL_ERROR
          "You are using unsupported Fortran compiler ${CMAKE_Fortran_COMPILER}\n"
          "Please install gfortran by running\n"
          "sudo apt install gfortran") 
ENDIF()

# OpenMP (shared memory parallelization), required by Femocs
FIND_PACKAGE(OpenMP)
IF (NOT OPENMP_FOUND)
  MESSAGE(FATAL_ERROR
          "OpenMP not found. Please re-install gcc") 
ENDIF()

# Boost (linear algebra, rnd nr generation, multithreading etc in c++)
FIND_PACKAGE(Boost)
IF (NOT Boost_FOUND)
  MESSAGE(FATAL_ERROR
          "Boost not found. Please install it by running\n"
          "sudo apt install libboost-all-dev") 
ENDIF()

IF (${machine} MATCHES "ubuntu" OR ${machine} MATCHES "alcyone")
    # Blas, required by Lapack
    FIND_PACKAGE(BLAS)
    IF (NOT BLAS_FOUND)
      MESSAGE(FATAL_ERROR
              "Blas not found. Please install it by running\n"
              "sudo apt install libblas-dev") 
    ENDIF()

    # Lapack (linear algebra package), required in many places
    FIND_PACKAGE(LAPACK)
    IF (NOT LAPACK_FOUND)
      MESSAGE(FATAL_ERROR
              "Lapack not found. Please install it by running\n"
              "sudo apt install liblapack-dev") 
    ENDIF()
ENDIF()

IF (${machine} MATCHES "ubuntu")
    # Tbb (Intel Threading Building Blocks), required by Deal.II
    FIND_PACKAGE(TBB)
    IF (NOT TBB_FOUND)
      MESSAGE(FATAL_ERROR
              "libtbb not found. Please install it by running\n"
              "sudo apt install libtbb-dev") 
    ENDIF()

    # Z, required by Deal.II
    FIND_PACKAGE(Z)
    IF (NOT Z_FOUND)
      MESSAGE(FATAL_ERROR
              "libz not found. Please install it by running\n"
              "sudo apt install libz-dev") 
    ENDIF()

    ## For some reason the check is not working, but hints are still valid
    # netcdf, required by Deal.II
    #SET (NETCDF_CXX "YES")
    #FIND_PACKAGE(NetCDF COMPONENTS CXX)
    #IF (NOT NetCDF_FOUND)
    #  MESSAGE(FATAL_ERROR
    #          "libnetcdf_c++ not found. Please install it by running\n"
    #          "sudo apt install libnetcdf-c++4\n"
    #          "sudo ln -s /usr/lib/*/libnetcdf_c++.so.4 /usr/lib/libnetcdf_c++.so") 
    #ENDIF()

    #  mpfr, required by CGAL
    FIND_PACKAGE(MPFR)
    IF (NOT MPFR_FOUND)
      MESSAGE(FATAL_ERROR
              "libmpfr not found. Please install it by running\n"
              "sudo apt install libmpfr-dev") 
    ENDIF()
ENDIF()

