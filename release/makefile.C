#############################################################################
#     Makefile for building Femocs executable for C main function
#############################################################################

include release/makefile.defs

all: release/femocs

release/femocs: release/femocs.o
	${CC} $< ${OPT} ${FCFLAGS} ${LDFLAGS} -o $@

release/femocs.o: lib/libfemocs.a lib/build/Makefile ${MAIN_C}
	${CC} -c ${MAIN_C} ${OPT} ${CCFLAGS}  -o $@

lib/libfemocs.a: src/* include/* heating/source/* heating/include/*
	make -f release/makefile.lib
    
lib/build/Makefile:
	cd lib; mkdir -p build; cd build; rm * -r; cmake -DCMAKE_Fortran_COMPILER=${F90} ..

clean:
	rm -f release/femocs release/femocs.o include/*.o src/*.o lib/build/Makefile
    
clean-all:
	rm -rf release/femocs release/femocs.o include/*.o src/*.o lib/build
