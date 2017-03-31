#############################################################################
#     Makefile for building FEMOCS - Finite Elements on Crystal Surfaces
#                             Mihkel Veske 2016
#
#############################################################################
# Before running make taito, run
#   module load gcc/5.3.0 intelmpi/5.1.3

include release/makefile.defs

all: lib

test_f90: lib/libfemocs.a femocs.test_f90
femocs.test_f90: ${MAIN_F90} src/* include/*
	make -f release/makefile.test_f90
	
test_c: lib/libfemocs.a femocs.test_c
femocs.test_c: ${MAIN_C} src/* include/*
	make -f release/makefile.test_c

lib: lib/libfemocs.a 
lib/libfemocs.a: src/* include/* 
	make -f release/makefile.lib mode=Release

release: femocs.release
femocs.release: ${MAIN_CPP} src/* include/* 
	make -f release/makefile.release mode=Release

heating: femocs.heating
femocs.heating: ${MAIN_HEATING} heating/source/* heating/include/* src/* include/*
	make -f release/makefile.heating

debug: 
	make -f release/makefile.debug

ubuntu:
	make -f release/makefile.ubuntu

taito:
	make -f release/makefile.taito

clean:
	make -s -f release/makefile.lib clean
	make -s -f release/makefile.release clean
	make -s -f release/makefile.test_f90 clean
	make -s -f release/makefile.test_c clean
	make -s -f release/makefile.heating clean

clean-all:
	make -s -f release/makefile.lib clean-all
	make -s -f release/makefile.release clean-all
	make -s -f release/makefile.test_f90 clean-all
	make -s -f release/makefile.test_c clean-all
	make -s -f release/makefile.heating clean-all
	make -s -f release/makefile.debug clean-all
	make -s -f release/makefile.ubuntu clean-all
	make -s -f release/makefile.taito clean-all

help:
	@echo ''
	@echo 'make all        pick default build type for Femocs'
	@echo 'make ubuntu     build in Ubuntu desktop all the external libraries that Femocs needs'
	@echo 'make install    build in Taito cluster in CSC all the external libraries that Femocs needs'
	@echo 'make lib        build Femocs as static library'
	@echo 'make heating    build Femocs executable from main file in heating module'
	@echo 'make release    build Femocs executable from c++ main with highest optimization level'
	@echo 'make debug      build Femocs executable from c++ main with debugging features enabled'
	@echo 'make test_f90   build Femocs executable from Fortran main'
	@echo 'make test_c     build Femocs executable from C main'
	@echo 'make clean      delete key files excluding installed libraries to start building from the scratch'
	@echo 'make clean-all  delete all the files and folders produced during the make process'
	@echo ''
