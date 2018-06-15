# Set  the minimum  required version  of cmake  for a  project.
cmake_minimum_required(VERSION 2.6)
#project(tetgen C CXX)
#set(CMAKE_C_COMPILER gcc-5)
#set(CMAKE_CXX_COMPILER gcc-5)

#Add a library to the project using the specified source files. 
# In Linux/Unix, it will creates the libtet.a
add_library(tet STATIC Tetgen.cxx Predicates.cxx)
