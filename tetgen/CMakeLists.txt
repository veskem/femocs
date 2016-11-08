# Set  the minimum  required version  of cmake  for a  project.
cmake_minimum_required(VERSION 2.6)

#Add a library to the project using the specified source files. 
# In Linux/Unix, it will creates the libtet.a
add_library(tet STATIC Tetgen.cxx Predicates.cxx)
