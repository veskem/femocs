## Overview
FEMOCS - Finite Elements of Crystal Surfaces - is a C++ library for coupling the effects of high electric field into atomistic simulations. For this FEMOCS imports the atomistic data, generates finite element mesh around it, solves the Laplace, heat and continuity equations and converts the calculated electric field and temperature into atomistic forces and velocities, correspondingly. 

FEMOCS can be used both in static and dynamic simulations. For the latter ones FEMOCS checks the difference between current and previous run and decides whether to perform the full calculation or uses the previous result. Such strategy helps to increase the effective

## Instructions to build the FEMOCS
All the build options are displayed with

    $ make help

To clean FEMOCS directory from previous builds type

    $ make clean     # doesn't clean Deal.II and other library installations
    $ make clean-all # performs full clean-up in FEMOCS library

To build FEMOCS as static library type

    $ make lib

To build FEMOCS as Fortran executable that uses FEMOCS library type

    $ make fortran

To build FEMOCS as C++ execuatable in Release mode without building FEMOCS as library type

    $ make release

To build FEMOCS as C++ execuatable in Debug mode without building FEMOCS as library type

    $ make debug

In **release** mode the FEMOCS is fully optimized, in **debug** mode the optimization is turned off and FEMOCS executable can be run in a debugger (GDB or its modifications).

## Testing and running the FEMOCS
The behavior of FEMOCS executable can be changed by changing the main files in **release** directory (*Test.f90* for mode **fortran** and *Test.cpp* for mode **release** and **debug** ).
Another way to control FEMOCS is to modify the configuration script. The sample script that also has PARCAS and HELMOD commands inside is located in **input/md.in**. 

FEMOCS executable that was built in mode **fortran**, **debug** or **release** should be run in FEMOCS main directory by typing

$ ./release/femocs


