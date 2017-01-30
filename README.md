## Overview
FEMOCS - Finite Elements of Crystal Surfaces - is a C++ library for coupling the effects of high electric field into Molecular Dynamics or kintetic Monte Carlo simulations. FEMOCS has the methods to

* import the atomistic data,
* generate finite element mesh around it,
* solve the Laplace, heat and continuity equations,
* convert calculated electric field and temperature into atomistic forces and velocities
* export results to atomistic simulation. 

FEMOCS can be used both in static and dynamic simulations. For the latter ones FEMOCS checks the difference between current and previous run and decides whether to perform the full calculation or to use the previous result. Such strategy helps to increase the computational efficiency of simulations.

## Instructions to build FEMOCS
All the build options are displayed with

    $ make help

FEMOCS directory is cleaned from previous builds by

    $ make clean     # doesn't clean the installations of Deal.II and other libraries
    $ make clean-all # performs the full clean-up in FEMOCS directory

FEMOCS is built as static library by

    $ make lib

For test purposes FEMOCS could also be built as an executable. The main file from directory **release** is chosen with the *make* argument:

    $ make fortran  # main file in Test.f90, full optimization, no debugger
    $ make release  # main file in Test.cpp, full optimization, no debugger
    $ make debug    # main file in Test.cpp, no optimization, debugger enabled

In **release** mode FEMOCS is fully optimized, in **debug** mode the optimization is turned off and FEMOCS executable can be run in a debugger (GDB or its modifications).

## Testing and running FEMOCS
The behavior of FEMOCS executable can be changed by changing the main files in **release** directory.
Another way to control FEMOCS is to modify the configuration script. The sample script that also has PARCAS and HELMOD commands inside is located in **input/md.in**. 

FEMOCS executable that was built in mode **fortran**, **debug** or **release** should be run in FEMOCS main directory by

    $ ./release/femocs
    
## Visualization and debugging
During the run FEMOCS writes couple of files to the **output** if **file_write = true** in *input/md.in*. Those files can be used to estimate the validity of the results. Files with extension *.xyz* and *.movie* contain atomistic data and can be visualized in [OVITO](https://ovito.org/index.php/download). The files ending with *.vtk* contain the geometric data (mesh elements, faces etc) and can be visualized in [ParaView](http://www.paraview.org/download/). Running the code in debug mode allows to follow the progress line-by-line by using [GDB](https://en.wikipedia.org/wiki/GNU_Debugger) or its graphical analogue in [Eclipse](http://www.eclipse.org/cdt/ ) debugger.
