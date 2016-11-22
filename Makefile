#############################################################################
#     Makefile for building FEMOCS - Finite Elements on Crystal Surfaces
#                             Mihkel Veske 2016
#
#############################################################################

include release/makefile.defs

all: fortran

lib: lib/libfemocs.a 
lib/libfemocs.a: src/* include/* 
	make -f release/makefile.lib
		
release: femocs.release
femocs.release: ${MAIN_CPP} src/* include/* 
	make -f release/makefile.release mode=Release
	
debug: femocs.debug
femocs.debug: ${MAIN_CPP} src/* include/* 
	make -f release/makefile.release mode=Debug
	
fortran: femocs.fortran
femocs.fortran: ${MAIN_F90} src/* include/*
	make -f release/makefile.fortran
	
heating: femocs.heating
femocs.heating: ${MAIN_HEATING} heating/source/* heating/include/* src/* include/*
	make -f release/makefile.heating
	
install:
	make -f release/makefile.install
	
clean:
	make -s -f release/makefile.lib clean
	make -s -f release/makefile.release clean
	make -s -f release/makefile.fortran clean
	make -s -f release/makefile.heating clean

clean-all:
	make -s -f release/makefile.lib clean-all
	make -s -f release/makefile.release clean-all
	make -s -f release/makefile.fortran clean-all
	make -s -f release/makefile.heating clean-all
	make -s -f release/makefile.install clean-all

help:
	@echo ''
	@echo 'make all        pick default build type for Femocs'
	@echo 'make lib        build Femocs as static library'
	@echo 'make heating    build Femocs executable together with 3D heating module'
	@echo 'make release    build Femocs executable from c++ main with highest optimization level'
	@echo 'make debug      build Femocs executable from c++ main with debugging features enabled'
	@echo 'make fortran    build Femocs executable from Fortran main'
	@echo 'make install    build all the external libraries that Femocs needs'
	@echo 'make clean      delete key files excluding installed libraries to start building from the scratch'
	@echo 'make clean-all  delete all the files and folders produced during the make process'
	@echo ''
