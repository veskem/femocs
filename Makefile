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

lib: lib/libfemocs.a

lib/libfemocs.a: src/* include/* deal-solver/source/* deal-solver/include/*
	make -f release/makefile.lib mode=Release

test_f90: ${MAIN_F90} lib/libfemocs.a
	make -f release/makefile.test_f90

test_c: ${MAIN_C} lib/libfemocs.a
	make -f release/makefile.test_c main_file=${MAIN_C} options=${OPT}

solver: ${MAIN_SOLVER} lib/libfemocs.a
	make -f release/makefile.test_c main_file=${MAIN_SOLVER} options=${OPT}

release: ${MAIN_CPP} lib/libfemocs.a
	make -f release/makefile.test_c main_file=${MAIN_CPP} options=${OPT}

debug: ${MAIN_CPP} src/* include/* deal-solver/source/* deal-solver/include/*
	make -f release/makefile.lib mode=Debug
	make -f release/makefile.test_c main_file=${MAIN_CPP} options=${OPT_DEBUG}

ubuntu:
	make -s -f release/makefile.cgal release
	make -f release/makefile.install

taito:
	make -s -f release/makefile.cgal taito
	make -f release/makefile.install

alcyone:
	make -f release/makefile.install

doc: doc/femocs.pdf
doc/femocs.pdf:
	cd doc; doxygen femocs.doxyfile; cd latex; make; cp refman.pdf ../femocs.pdf

clean:
	make -s -f release/makefile.lib clean
	make -s -f release/makefile.test_f90 clean
	make -s -f release/makefile.test_c clean
	rm -rf doc/html doc/latex *.pdf

clean-all:
	make -s -f release/makefile.lib clean-all
	make -s -f release/makefile.test_f90 clean-all
	make -s -f release/makefile.test_c clean-all
	make -s -f release/makefile.install clean-all
	make -s -f release/makefile.cgal clean-all
	rm -rf doc/html doc/latex *.pdf

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
