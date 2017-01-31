## Overview
FEMOCS - Finite Elements on Crystal Surfaces - is a C++ library for coupling the effects of high electric field into Molecular Dynamics or kinetic Monte Carlo simulations. FEMOCS has the methods to

* import the atomistic data,
* generate finite element mesh around it,
* solve the Laplace, heat and continuity equations,
* convert calculated electric field and temperature into atomistic forces and velocities,
* export results to atomistic simulation. 

FEMOCS can be used both in static and dynamic simulations. For the latter ones FEMOCS checks the difference between current and previous run and decides whether to perform the full calculation or to use the previous result. Such strategy helps to significantly increase the computational efficiency of simulations.

## Instructions to build FEMOCS
All the build options are displayed with

    $ make help

FEMOCS directory is cleaned from previous builds by

    $ make clean     # doesn't clean the installations of Deal.II and other libraries
    $ make clean-all # performs the full clean-up in FEMOCS directory

Before building FEMOCS as a library or executable, FEMOCS dependencies must be installed. In Ubuntu desktop this can be done with

    $ make ubuntu
    
In CSC Taito cluster first the appropriate modules must be loaded and after that the build can be performed:

    $ module load gcc/5.3.0 intelmpi/5.1.3
    $ make taito

FEMOCS is built as static library by

    $ make lib

For test purposes FEMOCS could also be built as an executable. The options for this are as follows:

    $ make fortran   # main file in release/Test.f90, full optimization, no debugger
    $ make release   # main file in release/Test.cpp, full optimization, no debugger
    $ make debug     # main file in release/Test.cpp, no optimization, debugger enabled

In **release** mode FEMOCS is fully optimized and emulates a static simulation. Running the code in **debug** mode turns optimization off and allows to follow the simulation line-by-line by using [GDB](https://en.wikipedia.org/wiki/GNU_Debugger), [Eclipse debugger](http://www.eclipse.org/cdt/) or their analogues.
    
## Testing and running FEMOCS
The behavior of FEMOCS can be changed by changing the main files in **release** directory or modifing the configuration script. The sample script that also has PARCAS and HELMOD commands inside is located in **input/md.in**. 

FEMOCS executable that was built in **fortran**, **release** or **debug** mode should be run in FEMOCS main directory by

    $ ./release/femocs
    
## Visualization
During the run FEMOCS writes couple of files to the **output** folder if *n_filewrite > 0* in configuration file. Those files can be used to estimate the validity of the results. Files with extension *.xyz* and *.movie* contain atomistic data and can be visualized in [OVITO](https://ovito.org/index.php/download). *.movie* files contain the data from many timesteps, *.xyz* files are the snapshots of a last full run. The files ending with *.vtk* contain the geometric data (mesh elements, faces etc) that can be visualized in [ParaView](http://www.paraview.org/download/). 
