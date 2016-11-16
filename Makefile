#############################################################################
#   	  Makefile for building Finite Elements on Crystal Surfaces
#                              Mihkel Veske 2016
# 
#############################################################################

include release/makefile.defs

all: debug 

lib: lib/libfemocs.a 
lib/libfemocs.a: src/* include/* 
	make -f release/makefile.lib
		
release: femocs.release
femocs.release: ${MAIN_CPP} src/* include/* 
	make -f release/makefile.release
	
debug: femocs.debug
femocs.debug: ${MAIN_CPP} src/* include/* 
	make -f release/makefile.debug
	
fortran: femocs.fortran
femocs.fortran: ${MAIN_F90} src/* include/*
	make -f release/makefile.fortran
	
install:
	make -f release/makefile.install
	
clean:
	make -s -f release/makefile.lib clean
	make -s -f release/makefile.debug clean
	make -s -f release/makefile.release clean
	make -s -f release/makefile.fortran clean

clean-all:
	make -s -f release/makefile.lib clean-all
	make -s -f release/makefile.debug clean-all
	make -s -f release/makefile.release clean-all
	make -s -f release/makefile.fortran clean-all
	make -s -f release/makefile.install clean-all
help:
	@echo ''
	@echo 'make all        pick default build type for Femocs'
	@echo 'make lib        build Femocs as static library'
	@echo 'make release    build Femocs executable with highest optimization level'
	@echo 'make debug      build Femocs executable with debugging features enabled'
	@echo 'make fortran    build Femocs executable by building Femocs as a static library and using the libray in Fortran code'
	@echo 'make install    build all the external libraries that Femocs needs'
	@echo 'make clean      delete key files excluding installed libraries to start building from the scratch'
	@echo 'make clean-all  delete all the files and folders produced during make process'
	@echo ''
