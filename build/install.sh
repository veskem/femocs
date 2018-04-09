#!/bin/bash
# Compile Femocs compilation and linking flags under different machines

## Optimization and warnings
OPT="-O3 -w"
OPT_DBG="-g -Og -Wall -Wpedantic -Wno-unused-local-typedefs"

## Paths to headers
HEADPATH_ALL="-Iinclude -Ilib -Idealii/include -IGETELEC/modules -std=c++11"
HEADPATH_UBUNTU=""
HEADPATH_TAITO="-I/appl/opt/cluster_studio_xe2018/parallel_studio_xe_2018.1.038/compilers_and_libraries_2018/linux/tbb/include"
HEADPATH_ALCYONE="-I/share/intel/composer_xe_2013.1.117/tbb/include"
HEADPATH_CGAL="-DUSE_CGAL=true -Icgal/include"

## Paths to incorporated libraries
LIBPATH_ALL="-Llib -LGETELEC/lib -Ldealii/lib"
LIBPATH_CGAL="-Lcgal/lib/x86_64-linux-gnu -Lcgal/lib64"

## Incorporated libraries
LIB_ALL="-ltet -ldeal_II -lgetelec -lslatec -fopenmp"
LIB_UBUNTU="-ltbb -llapack -lz -lm -lstdc++ -lnetcdf_c++ -lgfortran"
LIB_TAITO="-lz -lm -lstdc++ -lgfortran"
LIB_ALCYONE="-llapack -lz -lm -lstdc++ -lnetcdf_c++ -lgfortran"
LIB_CGAL="-lCGAL"

## Name of Femocs library in release and debug mode
LIB_FEMOCS="-lfemocs"
LIB_FEMOCS_DBG="-lfemocs_debug"


mode=$1

write_initial_flags () {
    rm -rf share/makefile.femocs
    cp build/makefile.empty share/makefile.femocs
    
    sed -i "/^FEMOCS_OPT=/ s|$|$OPT |" share/makefile.femocs
    sed -i "/^FEMOCS_DOPT=/ s|$|$OPT_DBG |" share/makefile.femocs
    
    sed -i "/^FEMOCS_HEADPATH=/ s|$|$HEADPATH_ALL |" share/makefile.femocs
    sed -i "/^FEMOCS_LIBPATH=/ s|$|$LIBPATH_ALL |" share/makefile.femocs

    sed -i "/^FEMOCS_LIB=/ s|$|$LIB_FEMOCS |" share/makefile.femocs
    sed -i "/^FEMOCS_LIB=/ s|$|$LIB_ALL |" share/makefile.femocs
    
    sed -i "/^FEMOCS_DLIB=/ s|$|$LIB_FEMOCS_DBG |" share/makefile.femocs
    sed -i "/^FEMOCS_DLIB=/ s|$|$LIB_ALL |" share/makefile.femocs
}

print_all() { 
    grep "FEMOCS_OPT=" share/makefile.femocs
    grep "FEMOCS_DOPT=" share/makefile.femocs
    grep "FEMOCS_HEADPATH=" share/makefile.femocs
    grep "FEMOCS_LIBPATH=" share/makefile.femocs
    grep "FEMOCS_LIB=" share/makefile.femocs
    grep "FEMOCS_DLIB=" share/makefile.femocs
    grep "MACHINE=" share/makefile.femocs
}

if (test $mode = ubuntu) then
    if [ ! -f share/makefile.femocs ]; then
        echo "Checking Ubuntu dependencies"
        mkdir -p share/.build && rm -rf share/.build/*
        cd share/.build && cmake .. -Dmachine=ubuntu && cd ../..
        echo "All Femocs external dependencies found"
    fi

    echo -e "\nWriting Ubuntu flags"
    write_initial_flags
    
    sed -i "/^FEMOCS_HEADPATH=/ s|=|=$HEADPATH_UBUNTU |" share/makefile.femocs
    sed -i "/^FEMOCS_LIB=/ s|$|$LIB_UBUNTU |" share/makefile.femocs
    sed -i "/^FEMOCS_DLIB=/ s|$|$LIB_UBUNTU |" share/makefile.femocs
    
    sed -i "/MACHINE=/c\MACHINE=ubuntu" share/makefile.femocs
        
    print_all
    
    echo -e "\nInstalling Femocs dependencies"
    make -f build/makefile.install
fi

if (test $mode = taito) then
    echo "Loading Taito modules"
    module load gcc/5.3.0 intelmpi/5.1.3

    echo -e "\nWriting Taito flags"
    write_initial_flags

    sed -i "/^FEMOCS_HEADPATH=/ s|=|=$HEADPATH_TAITO |" share/makefile.femocs
    sed -i "/^FEMOCS_LIB=/ s|$|$LIB_TAITO |" share/makefile.femocs
    sed -i "/^FEMOCS_DLIB=/ s|$|$LIB_TAITO |" share/makefile.femocs
    
    sed -i "/MACHINE=/c\MACHINE=taito" share/makefile.femocs

    print_all
    
    echo -e "\nInstalling Femocs dependencies"
    make -f build/makefile.install
fi

if (test $mode = alcyone) then
    echo "Loading Alcyone modules"
    module load PrgEnv-gnu gcc/5.1.0

    echo -e "\nWriting Alcyone flags"
    write_initial_flags

    sed -i "/^FEMOCS_HEADPATH=/ s|=|=$HEADPATH_ALCYONE |" share/makefile.femocs
    sed -i "/^FEMOCS_LIB=/ s|$|$LIB_ALCYONE |" share/makefile.femocs
    sed -i "/^FEMOCS_DLIB=/ s|$|$LIB_ALCYONE |" share/makefile.femocs

    sed -i "/MACHINE=/c\MACHINE=alcyone" share/makefile.femocs
    
    print_all
    
    echo -e "\nInstalling Femocs dependencies"
    make -f build/makefile.install
fi

if (test $mode = cgal) then
    echo "Adding CGAL flags"

    sed -i "/^FEMOCS_HEADPATH=/ s|=|=$HEADPATH_CGAL |" share/makefile.femocs
    sed -i "/^FEMOCS_LIBPATH=/ s|$|$LIBPATH_CGAL |" share/makefile.femocs
    sed -i "/^FEMOCS_LIB=/ s|=|=$LIB_CGAL |" share/makefile.femocs
    sed -i "/^FEMOCS_DLIB=/ s|=|=$LIB_CGAL |" share/makefile.femocs

    sed -i "/LIBCGAL=/c\LIBCGAL=1" share/makefile.femocs

    print_all
    grep "LIBCGAL=" share/makefile.femocs
    
    machine=$( grep "MACHINE=" share/makefile.femocs | sed "s|MACHINE=||" )       
    echo -e "\nInstalling CGAL in "$machine
    if (test $machine = taito) then
        make -s -f build/makefile.cgal taito
    else
        make -f build/makefile.cgal
    fi
fi

if (test $mode = no-cgal) then
    echo "Removing CGAL flags"

    sed -i -- "s|$HEADPATH_CGAL ||g" share/makefile.femocs
    sed -i -- "s|$LIBPATH_CGAL ||g" share/makefile.femocs
    sed -i -- "s|$LIB_CGAL ||g" share/makefile.femocs

    sed -i "/LIBCGAL=/c\LIBCGAL=" share/makefile.femocs

    print_all
    grep "LIBCGAL=" share/makefile.femocs
fi
