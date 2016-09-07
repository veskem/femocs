
FC = gfortran-5
DEALII = ../deal.II/lib
EXEC = ../test/Helmod.f90

FCFLAGS = -Wall -Wextra -L../lib -L${DEALII} -I../lib -I../include
LDFLAGS = -lfemocs -ltet -ldeal_II -fopenmp -ltbb -lpthread -lumfpack -lblas -llapack -lm -lz -lstdc++ -lnetcdf_c++ -std=c++11

all: femocs

femocs: femocs.o ../src/* ../include/*
	cd build; make; cd ..
	${FC} femocs.o ${FCFLAGS} ${LDFLAGS} -o femocs
	
femocs.o: ${EXEC}
	${FC} -c ${EXEC} ${FCFLAGS} ${LDFLAGS} -o femocs.o
	
clean:
	rm femocs femocs.o
