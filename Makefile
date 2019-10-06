#############################################################################
#     Makefile for building FEMOCS - Finite Elements on Crystal Surfaces
#                             Mihkel Veske 2016
#
#############################################################################
# Before running make install-taito or running Femocs in Taito, run
#   source ./build/load_modules.sh taito
#
# Before running make install-alcyone or running Femocs in Alcyone, run
#   source ./build/load_modules.sh alcyone
#############################################################################

include build/makefile.defs


all: lib

lib: femocs_lib
femocs_lib:
	@make --no-print-directory -f build/makefile.lib build_type=Release

dlib: femocs_dlib
femocs_dlib:
	@make --no-print-directory -f build/makefile.lib build_type=Debug

test_f90: femocs_lib femocs_f90
femocs_f90:
	@make --no-print-directory -f build/makefile.main main=${FMAIN} compiler=${F90}

test_c: femocs_lib femocs_c
femocs_c:
	@make --no-print-directory -f build/makefile.main main=${CMAIN} compiler=${CC}

test_cpp: femocs_lib femocs_release

release: femocs_lib femocs_release
femocs_release:
	@make --no-print-directory -f build/makefile.main main=${CXXMAIN} compiler=${CXX}

debug: femocs_dlib femocs_debug
femocs_debug:
	make -s -f build/makefile.main build/femocs_debug main=${CXXMAIN} compiler=${CXX}

#release: femocs_release
#femocs_release:
#	make -s -f build/makefile.exec build/femocs build_type=Release
#
#debug: femocs_debug
#femocs_debug:
#	make -s -f build/makefile.exec build/femocs build_type=Debug

doc: femocs_doc
femocs_doc:
	make -f build/makefile.doc

install-ubuntu:
	@chmod +x ./build/install.sh
	@./build/install.sh ubuntu

install-taito:
	@chmod +x ./build/install.sh
	@./build/install.sh taito

install-alcyone:
	@chmod +x ./build/install.sh
	@./build/install.sh alcyone

install-cgal:
	@chmod +x ./build/install.sh
	@./build/install.sh cgal

install-no-cgal:
	@chmod +x ./build/install.sh
	@./build/install.sh no-cgal

uninstall-all:
	@chmod +x ./build/uninstall-all.sh
	./build/uninstall-all.sh

clean:
	make -s -f build/makefile.lib clean
	make -s -f build/makefile.exec clean
	make -s -f build/makefile.main clean
	make -s -f build/makefile.doc clean

clean-all:
	make -s -f build/makefile.lib clean-all
	make -s -f build/makefile.exec clean-all
	make -s -f build/makefile.main clean-all
	make -s -f build/makefile.doc clean-all

help:
	@echo ''
	@echo 'make all        pick default build type for Femocs'
	@echo ''
	@echo 'make install-'
	@echo '       ubuntu   build in Ubuntu desktop Femocs mandatory dependencies'
	@echo '       taito    build in CSC Taito cluster all the external libraries that Femocs needs'
	@echo '       alcyone  build in Alcyone cluster all the external libraries that Femocs needs'
	@echo '       cgal     build CGAL and enable its usage in the code'
	@echo '       no-cgal  disable CGAL usage in the code'
	@echo ''
	@echo 'make uninstall-'
	@echo '       all      remove all installation files'
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
