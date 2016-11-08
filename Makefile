
COMPILER = gfortran-5
DEALII_VER = 8.4.1
MAIN = release/Helmod.f90
NPROCS = 4

WARNINGS = -Wall -Wextra
FCFLAGS = -Llib -Ldealii/lib -Ilib -Iinclude -Idealii/include
LDFLAGS = -lfemocs -ltet -ldeal_II -fopenmp -ltbb -lpthread -lumfpack -lblas -llapack -lm -lz -lstdc++ -lnetcdf_c++ -std=c++11

all: release/femocs

release/femocs: lib/build/Makefile lib/libfemocs.mod lib/libtet.a include/Tetgen.h release/femocs.o src/* include/*
	${COMPILER} release/femocs.o ${WARNINGS} ${FCFLAGS} ${LDFLAGS} -o release/femocs	

release/femocs.o: ${MAIN}
	cd lib/build; make -j${NPROCS}
	${COMPILER} -c ${MAIN} ${WARNINGS} ${FCFLAGS} ${LDFLAGS} -o release/femocs.o
	
lib/libfemocs.mod: src/* include/*
	cd lib; mkdir -p build; cd build; rm * -r; cmake ..

lib/build/Makefile: dealii/lib/libdeal_II.a src/* include/*
	cd lib; mkdir -p build; cd build; rm * -r; cmake ..

include/Tetgen.h:
	cp tetgen/Tetgen.h include/Tetgen.h
	
lib/libtet.a:
	cd tetgen; rm -rf build; mkdir build
	cd tetgen/build; cmake ..; make -j${NPROCS}
	cd tetgen; mv build/libtet.a ../lib; rm -rf build
	
dealii/lib/libdeal_II.a: dealii/dealii-${DEALII_VER}.tar.gz
	cd dealii; mv dealii-${DEALII_VER}.tar.gz ..; rm -rf *; mv ../dealii-${DEALII_VER}.tar.gz .
	cd dealii; tar -xvf dealii-${DEALII_VER}.tar.gz; cd dealii-*; mkdir build
	cd dealii/dealii-*/build; cmake -DCMAKE_BUILD_TYPE=Release -DDEAL_II_STATIC_EXECUTABLE=ON -DCMAKE_INSTALL_PREFIX=../.. ..; make -j${NPROCS} install
	cd dealii; rm -r dealii-${DEALII_VER}; rm -r examples

dealii/dealii-${DEALII_VER}.tar.gz:
	cd dealii; wget https://github.com/dealii/dealii/releases/download/v${DEALII_VER}/dealii-${DEALII_VER}.tar.gz

clean:
	rm release/femocs release/femocs.o lib/libfemocs.a lib/libfemocs.mod
