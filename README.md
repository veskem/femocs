## Overview
FEMOCS - Finite Elements on Crystal Surfaces - is a C++ library for coupling the effects of high electric field into Molecular Dynamics or kinetic Monte Carlo simulations. FEMOCS has the methods to

* import the atomistic data,
* generate finite element mesh around it,
* solve the Laplace, heat and continuity equations,
* convert calculated electric field and temperature into atomistic forces and velocities,
* export results to atomistic simulation. 

FEMOCS can be used both in static and dynamic simulations. For the latter ones FEMOCS checks the difference between current and previous run and decides whether to perform the full calculation or to use the previous result. Such strategy helps to significantly increase the computational efficiency of simulations.

## Citing
FEMOCS is an open-source and freely available code. The details about its algorithms are described in an article in [Journal of Computational Physics](https://doi.org/10.1016/j.jcp.2018.04.031).

When publishing results obtained with the help of FEMOCS, please cite

    Veske, M. et al, 2018. Dynamic coupling of a finite element solver to large-scale atomistic simulations. Journal of Computational Physics, 367, pp.279–294.

## Instructions to build FEMOCS
All the build options are displayed with

    $ make help

FEMOCS directory is cleaned from previous builds by

    $ make clean     # deletes executables and object files
    $ make clean-all # performs the full clean-up in FEMOCS directory

Before building FEMOCS as a library or executable, FEMOCS dependencies must be installed. This can be done with

    $ make install-ubuntu   # Ubuntu desktop            
    $ make install-taito    # CSC Taito
    $ make install-alcyone  # Alcyone
    
The files that were build during installation can be removed by
    
    $ make uninstall-all

FEMOCS is built as static library by

    $ make lib       # fully optimized and no debugger flags
    $ make dlib      # debugging options enabled

For test purposes FEMOCS could also be built as an executable. The options for this are as follows:

    $ make debug     # main file in src/main/Main.cpp, minimal optimization, debugger & warnings enabled
    $ make release   # main file in src/main/Main.cpp, full optimization, no debugger & warnings
    $ make test_f90  # main file in src/main/Main.f90, full optimization, no debugger & warnings
    $ make test_cpp  # main file in src/main/Main.cpp, full optimization, no debugger & warnings
    $ make test_c    # main file in src/main/Main.c, full optimization, no debugger & warnings

In **release** mode FEMOCS is fully optimized and emulates a static simulation. Running the code in **debug** mode minimizes optimization and allows to follow the simulation line-by-line by using [GDB](https://en.wikipedia.org/wiki/GNU_Debugger), [Eclipse debugger](http://www.eclipse.org/cdt/) or their analogues.
    
## Testing and running FEMOCS
The behavior of FEMOCS can be changed by changing the main files in **src** directory or by modifing the configuration script. The sample script is located in **in/sample.in**. 

FEMOCS executable that was built in **debug**, **release**, **test_f90**, **test_c** or **test_cpp** mode should be run in FEMOCS main directory by

    $ ./build/femocs        # release, test_f90, test_c and test_cpp mode
    $ ./build/femocs_debug  # debug mode
    
## Visualization
During the run FEMOCS writes couple of files to the **out** folder if *n_filewrite > 0* in configuration file. Those files can be used to estimate the validity of the results. Files with extension *xyz* and *movie* contain atomistic data and can be visualized in [OVITO](https://ovito.org/index.php/download). *movie* files contain the data from many timesteps, *xyz* files are the snapshots of a last full run. The files ending with *vtk* contain the geometric data (mesh elements, faces etc) that can be visualized in [ParaView](http://www.paraview.org/download/).

## Documentation
To build FEMOCS documentation, first make sure [Doxygen](http://www.stack.nl/~dimitri/doxygen/download.html) is installed in the system. The documentation in *pdf* and *html* format will be generated with command

    $ make doc
    
## Setting up Eclipse
For developers it's highly recommended to do FEMOCS development in an IDE (integrated develpment editor) instead of some conventional text editor. However, to take maximum out of IDE, it is necessary to configure it properly before usage. The fresh copy of [Eclipse](https://www.eclipse.org/downloads/eclipse-packages/) underlines c++ code even if it doesn't contain any errors there. It's because fresh copy of Eclipse doesn't know about *stdc++11* features. Below are the instructions to make Eclipse *Oxygen* aware of them.

**Step 1**

    Window > Preferences
    C/C++ > Build > Settings > Tab [Discovery] > CDT GCC Built-in Compiler Settings

Add the *-std=c++11* flag to *Command to get compiler specs*. After the modification it should be similar to

    ${COMMAND} ${FLAGS} -E -P -v -std=c++11 -dD "${INPUTS}"

**Step 2**

Repeat the same procedure in

    Project > Properties
    C/C++ General > Preprocessor include > Tab [Providers] > CDT GCC Built-in Compiler Settings

**Step 3**

    Project > Properties
    C/C++ General > Path and Symbols > Tab [Symbols] > GNU C++

Add the symbol *__cplusplus* with the value *201103L*

**Step 4**

Clean and rebuild both the project (*Project > Clean*, *Project > Rebuild all*) and the index (*Project > C/C++ Index > Rebuild*) as Eclipse tends to cache error messages and show them even though they are gone after changing the settings.
