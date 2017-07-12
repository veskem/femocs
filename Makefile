#############################################################################
#     Makefile for building FEMOCS - Finite Elements on Crystal Surfaces
#                             Mihkel Veske 2016
#
#############################################################################
# Before running make taito, run
#   module load gcc/5.3.0 intelmpi/5.1.3
# Before running make alcyone, run
#   module load PrgEnv-gnu gcc/5.1.0

include release/makefile.defs

all: lib

lib: femocs_lib
femocs_lib:
	make -f release/makefile.lib mode=Release

test_f90: femocs_lib femocs_f90
femocs_f90:
	make -f release/makefile.exec main=${FMAIN} cf=${FCFLAGS} compiler=${F90}

test_c: femocs_lib femocs_c
femocs_c:
	make -f release/makefile.exec main=${CMAIN} cf=${CCFLAGS} compiler=${CC}

solver: femocs_lib femocs_solver
femocs_solver:
	make -f release/makefile.exec main=${SOLVERMAIN} cf=${CXXCFLAGS} compiler=${CXX}

release: femocs_lib femocs_release
femocs_release:
	make  -f release/makefile.exec main=${CXXMAIN} cf=${CXXCFLAGS} compiler=${CXX}

debug: femocs_debug
femocs_debug:
	make -f release/makefile.lib mode=Debug
	make -s -f release/makefile.exec release/femocs.debug main=${CXXMAIN} cf=${CXXCFLAGS} compiler=${CXX}

doc: femocs_doc
femocs_doc:
	make -f release/makefile.doc
	
ubuntu:
	make -s -f release/makefile.cgal release
	make -f release/makefile.install

taito:
	make -s -f release/makefile.cgal taito
	make -f release/makefile.install

alcyone:
	make -f release/makefile.install

clean:
	make -s -f release/makefile.lib clean
	make -s -f release/makefile.exec clean
	make -s -f release/makefile.doc clean

clean-all:
	make -s -f release/makefile.lib clean-all
	make -s -f release/makefile.exec clean-all
	make -s -f release/makefile.install clean-all
	make -s -f release/makefile.cgal clean-all
	make -s -f release/makefile.doc clean-all

help:
	@echo ''
	@echo 'make all        pick default build type for Femocs'
	@echo 'make ubuntu     build in Ubuntu desktop all the external libraries that Femocs needs'
	@echo 'make taito      build in CSC Taito cluster all the external libraries that Femocs needs'
	@echo 'make alcyone    build in Alcyone cluster all the external libraries that Femocs needs'
	@echo 'make lib        build Femocs as static library'
	@echo 'make solver     build Femocs executable from main file in the FEM solver module'
	@echo 'make release    build Femocs executable from c++ main with highest optimization level'
	@echo 'make debug      build Femocs executable from c++ main with debugging features enabled'
	@echo 'make test_f90   build Femocs executable from Fortran main'
	@echo 'make test_c     build Femocs executable from C main'
	@echo 'make doc        generate Femocs documentation in html and pdf format
	@echo 'make clean      delete key files excluding installed libraries to start building from the scratch'
	@echo 'make clean-all  delete all the files and folders produced during the make process'
	@echo ''
