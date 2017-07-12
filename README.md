## Overview
FEMOCS - Finite Elements on Crystal Surfaces - is a C++ library for coupling the effects of high electric field into Molecular Dynamics or kinetic Monte Carlo simulations. FEMOCS has the methods to

* import the atomistic data,
* generate finite element mesh around it,
* solve the Laplace, heat and continuity equations,
* convert calculated electric field and temperature into atomistic forces and velocities,
* export results to atomistic simulation. 

FEMOCS can be used both in static and dynamic simulations. For the latter ones FEMOCS checks the difference between current and previous run and decides whether to perform the full calculation or to use the previous result. Such strategy helps to significantly increase the computational efficiency of simulations.

## Citing
FEMOCS is an open-source and freely available code. The details about its algorithms are described in a preprint in [arXiv](https://arxiv.org/abs/1706.09661).

When publishing results obtained with the help of FEMOCS, please cite

    Veske, M., Kyritsakis, A., Eimre, K., Zadin, V., Aabloo, A. and Djurabekova, F., 2017. Dynamic coupling of a finite element solver to large-scale atomistic simulations. arXiv preprint arXiv:1706.09661.

## Instructions to build FEMOCS
All the build options are displayed with

    $ make help

FEMOCS directory is cleaned from previous builds by

    $ make clean     # doesn't clean the installations of deal.II and other libraries
    $ make clean-all # performs the full clean-up in FEMOCS directory

Before building FEMOCS as a library or executable, FEMOCS dependencies must be installed. In Ubuntu desktop this can be done with

    $ make ubuntu
    
In clusters first the appropriate modules must be loaded and after that the build can be performed:

    $ module load gcc/5.3.0              # CSC Taito
    $ make taito
    $ module load PrgEnv-gnu gcc/5.1.0   # Alcyone
    $ make alcyone

FEMOCS is built as static library by

    $ make lib

For test purposes FEMOCS could also be built as an executable. The options for this are as follows:

    $ make debug     # main file in build/Test.cpp, minimal optimization, debugger & warnings enabled
    $ make release   # main file in build/Test.cpp, full optimization, no debugger & warnings
    $ make solver    # main file in deal-solver/main/main.cc, full optimization, no debugger & warnings
    $ make test_f90  # main file in build/Test.f90, full optimization, no debugger & warnings
    $ make test_c    # main file in build/Test.c, full optimization, no debugger & warnings

In **release** mode FEMOCS is fully optimized and emulates a static simulation. Running the code in **debug** mode minimizes optimization and allows to follow the simulation line-by-line by using [GDB](https://en.wikipedia.org/wiki/GNU_Debugger), [Eclipse debugger](http://www.eclipse.org/cdt/) or their analogues.
    
## Testing and running FEMOCS
The behavior of FEMOCS can be changed by changing the main files in **build** directory or by modifing the configuration script. The sample script that also has PARCAS and HELMOD commands inside is located in **in/md.in**. 

FEMOCS executable that was built in **debug**, **release**, **solver**, **test_f90** or **test_c** mode should be run in FEMOCS main directory by

    $ ./build/femocs        # release, solver, test_f90 or test_c mode
    $ ./build/femocs.debug  # debug mode
    
## Visualization
During the run FEMOCS writes couple of files to the **out** folder if *n_filewrite > 0* in configuration file. Those files can be used to estimate the validity of the results. Files with extension *xyz* and *movie* contain atomistic data and can be visualized in [OVITO](https://ovito.org/index.php/download). *movie* files contain the data from many timesteps, *xyz* files are the snapshots of a last full run. The files ending with *vtk* contain the geometric data (mesh elements, faces etc) that can be visualized in [ParaView](http://www.paraview.org/download/).

## Documentation
To build FEMOCS documentation, first make sure [Doxygen](http://www.stack.nl/~dimitri/doxygen/download.html) is installed in the system. The documentation in *pdf* and *html* format will be generated with command

    $ make doc
