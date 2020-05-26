## Overview
FEMOCS - Finite Elements on Crystal Surfaces - is a C++ library for coupling the effects of high electric
field into Molecular Dynamics or kinetic Monte Carlo simulations. Current version of FEMOCS can

* import the atomistic data,
* generate finite element mesh around it,
* solve the Poisson's, heat and continuity equations,
* convert calculated electric field and temperature into atomistic forces and velocities,
* export results to atomistic simulation. 

FEMOCS can be used both in static and dynamic simulations. For the latter ones FEMOCS checks the
difference between current and previous run and decides whether to perform the full calculation or to use
the previous result. Such strategy helps to significantly increase the computational efficiency of
simulations.

## Citing
FEMOCS is an open-source and freely available code. The details about its algorithms are published in
[Journal of Computational Physics](https://doi.org/10.1016/j.jcp.2018.04.031) and
[Physical Review E](https://doi.org/10.1103/PhysRevE.101.053307).
When publishing results obtained with the help of FEMOCS, please cite

    Veske, M. et al, 2020. Dynamic coupling between particle-in-cell and atomistic simulations. Physical Review E, 101(5), p.053307.

## Instructions to build FEMOCS
All the build options are displayed with

    $ make help

FEMOCS directory is cleaned from previous builds with

    $ make clean      # deletes executables and object files
    $ make clean-all  # performs the full clean-up in FEMOCS directory

Before building FEMOCS as a library or executable, FEMOCS dependencies must be installed.
This can be done with

    $ make install-ubuntu   # Ubuntu desktop            
    $ make install-taito    # CSC Taito
    $ make install-alcyone  # Alcyone
    
The files that were built during installation can be removed with
    
    $ make uninstall-all

FEMOCS is built as static library with

    $ make lib       # fully optimized and no debugger flags
    $ make dlib      # debugging options enabled

FEMOCS could also be built as an executable. The options for this are as follows:

    $ make debug     # main file in src/main/Main.cpp, minimal optimization, debugger & warnings enabled
    $ make release   # main file in src/main/Main.cpp, full optimization, no debugger & warnings
    $ make test_f90  # main file in src/main/Main.f90, full optimization, no debugger & warnings
    $ make test_cpp  # main file in src/main/Main.cpp, full optimization, no debugger & warnings
    $ make test_c    # main file in src/main/Main.c, full optimization, no debugger & warnings

In **release** mode FEMOCS is fully optimized and emulates a static simulation. Running the code in
**debug** mode minimizes optimization and allows to follow the simulation line-by-line by using
[GDB](https://en.wikipedia.org/wiki/GNU_Debugger), [Eclipse debugger](http://www.eclipse.org/cdt/)
or their analogues.
    
## Notes about linking & compilation
FEMOCS was developed and tested using **GNU** compiler with version **5**. The earlier and later versions
are likely to cause compilation errors. Therefore, before running any FEMOCS installation or compilation
routine, make sure your machine is using the mentioned compiler. One way to ensure that correct compiler
is used is to modify the *.bashrc* file. Another option is to modify compiler variables in 
*build/makefile.defs*. The most reliable option, that is not always available, thou, is to make **GNU 5**
as a default system-wise compiler. Notice that FEMOCS uses *gcc*, *c++* and *gfortran* compilers that all
must meet the correct version criterion.

While linking FEMOCS library it might happen that some of the libraries that are needed by Deal.II are
not found from the system. Often this issue can be solved by modifing the file *share/makefile.femocs*
that was created after running *make install-machine*. In this file find the libraries the linker is
complaining about and remove them.

To obtain the compilation and linker flags for the external code that is using FEMOCS, add
    
    include path_to_femocs/share/makefile.femocs
    
into your makefile. The variables in this file contain all the includes and paths
that are needed to call FEMOCS externally. Those paths are relative to FEMOCS main directory. One way to
change the origin is to use *patsubst* command. For example, if FEMOCS is located inside a directory
*my_project* and you want to obtain the location or FEMOCS libraries with respect to this directory, add
the following line to your makefile

    NEW_PATHS=$(patsubst -L%, -Lmy_project/%, $(FEMOCS_LIBPATH))

## Testing and running FEMOCS
The behavior of FEMOCS can be changed by changing the main files in **src** directory or by modifing the
configuration script. The sample script is located in **in/sample.in**. 

FEMOCS executable that was built in **debug**, **release**, **test_f90**, **test_c** or **test_cpp** mode
should be run in FEMOCS main directory by

    $ ./build/femocs        # release, test_f90, test_c and test_cpp mode
    $ ./build/femocs_debug  # debug mode
    
## Visualization
During the run FEMOCS writes couple of files to the **out** folder if *n_filewrite > 0* in configuration
file. Those files can be used to estimate the validity of the results. Files with extension *xyz* and
*movie* contain atomistic data and can be visualized in [OVITO](https://ovito.org/index.php/download).
*movie* files contain the data from many timesteps, *xyz* files are the snapshots of a last full run.
The files ending with *vtk* contain the geometric data (mesh elements, faces etc) that can be visualized
in [ParaView](http://www.paraview.org/download/).

## Documentation
To build FEMOCS documentation, first make sure [Doxygen](http://www.stack.nl/~dimitri/doxygen/download.html)
is installed in the system. The documentation in *pdf* and *html* format will be generated with command

    $ make doc
    
## Setting up Eclipse
For developers it's highly recommended to do FEMOCS development in an IDE (integrated develpment editor)
instead of some conventional text editor. However, to take maximum out of IDE, it is necessary to
configure it properly before usage. The fresh copy of [Eclipse](https://www.eclipse.org/downloads/eclipse-packages/)
underlines c++ code even if it doesn't contain any errors there. It's because fresh copy of Eclipse does
not know about *stdc++11* features. Below are the instructions to make Eclipse *Oxygen* aware of them.

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

Clean and rebuild both the project (*Project > Clean*, *Project > Rebuild all*) and the index
(*Project > C/C++ Index > Rebuild*) as Eclipse tends to cache error messages and show them even though
they are gone after changing the settings.
