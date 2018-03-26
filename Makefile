#############################################################################
#     Makefile for building FEMOCS - Finite Elements on Crystal Surfaces
#                             Mihkel Veske 2016
#
#############################################################################
# Before running make taito, run
#   module load gcc/5.3.0 intelmpi/5.1.3
# Before running make alcyone, run
#   module load PrgEnv-gnu gcc/5.1.0

include build/makefile.defs

all: lib

lib: femocs_lib
femocs_lib:
	make -s -f build/makefile.lib lib/libfemocs.a

dlib: femocs_debug_lib
femocs_debug_lib:
	make -s -f build/makefile.lib lib/libfemocs_debug.a

test_f90: femocs_lib femocs_f90
femocs_f90:
	make -f build/makefile.main main=${FMAIN} cf=${FCFLAGS} compiler=${F90}

test_c: femocs_lib femocs_c
femocs_c:
	make -f build/makefile.main main=${CMAIN} cf=${CCFLAGS} compiler=${CC}

test_cpp: femocs_lib femocs_cpp
femocs_cpp:
	make -f build/makefile.main main=${CXXMAIN} cf=${CXXCFLAGS} compiler=${CXX}
	
release: femocs_lib femocs_release
femocs_release:
	make -f build/makefile.main main=${CXXMAIN} cf=${CXXCFLAGS} compiler=${CXX}

debug: femocs_debug_lib femocs_debug
femocs_debug:
	make -f build/makefile.main build/femocs_debug main=${CXXMAIN} cf=${CXXCFLAGS} compiler=${CXX}

#release: femocs_release
#femocs_release:
#	make -s -f build/makefile.exec build/femocs
#
#debug: femocs_debug
#femocs_debug:
#	make -s -f build/makefile.exec build/femocs_debug

doc: femocs_doc
femocs_doc:
	make -f build/makefile.doc

ubuntu:
	mkdir -p share/.build && cd share/.build && rm -rf * && cmake .. -Dmachine=ubuntu
	#make -s -f build/makefile.cgal release
	make -f build/makefile.install

taito:
	#make -s -f build/makefile.cgal taito
	make -f build/makefile.install

alcyone:
	make -f build/makefile.install

clean:
	make -s -f build/makefile.lib clean
	make -s -f build/makefile.exec clean
	make -s -f build/makefile.main clean
	make -s -f build/makefile.doc clean

clean-all:
	@chmod +x ./build/makefile.clean-all
	./build/makefile.clean-all
	@chmod -x ./build/makefile.clean-all

help:
	@echo ''
	@echo 'make all        pick default build type for Femocs'
	@echo ''
	@echo 'make ubuntu     build in Ubuntu desktop all the external libraries that Femocs needs'
	@echo 'make taito      build in CSC Taito cluster all the external libraries that Femocs needs'
	@echo 'make alcyone    build in Alcyone cluster all the external libraries that Femocs needs'
	@echo ''
	@echo 'make lib        build Femocs as static library with maximum optimization level'
	@echo 'make dlib       build Femocs as static library with debugging features enabled'
	@echo 'make release    build Femocs executable from c++ main with highest optimization level'
	@echo 'make debug      build Femocs executable from c++ main with debugging features enabled'
	@echo ''
	@echo 'make test_f90   build Femocs executable from Fortran main'
	@echo 'make test_c     build Femocs executable from C main'
	@echo 'make test_cpp   build Femocs executable from C++ main'
	@echo ''
	@echo 'make doc        generate Femocs documentation in html and pdf format'
	@echo 'make clean      delete key files excluding installed libraries to start building from the scratch'
	@echo 'make clean-all  delete all the files and folders produced during the make process'
	@echo ''
