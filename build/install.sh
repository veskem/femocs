#!/bin/bash
# Compile Femocs compilation and linking flags under different machines

source ./build/makefile.defs

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

    sed -i "/^FEMOCS_DLIB=/ s|$|$LIB_FEMOCS |" share/makefile.femocs
    sed -i "/^FEMOCS_DLIB=/ s|$|$DLIB_ALL |" share/makefile.femocs

    sed -i "/^TETGEN_FLAGS=/ s|$|$TETGEN_FLAGS |" share/makefile.femocs
    sed -i "/^DEALII_FLAGS=/ s|$|$DEALII_FLAGS |" share/makefile.femocs
}

print_all() { 
    grep "FEMOCS_OPT=" share/makefile.femocs
    grep "FEMOCS_DOPT=" share/makefile.femocs
    grep "FEMOCS_HEADPATH=" share/makefile.femocs
    grep "FEMOCS_LIBPATH=" share/makefile.femocs
    grep "FEMOCS_LIB=" share/makefile.femocs
    grep "FEMOCS_DLIB=" share/makefile.femocs
    grep "TETGEN_FLAGS=" share/makefile.femocs
    grep "DEALII_FLAGS=" share/makefile.femocs
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
    echo -e "\nWriting Taito flags"
    write_initial_flags

    sed -i "/^FEMOCS_EXT_HEAD=/ s|=|=$HEADPATH_TAITO |" share/makefile.femocs
    sed -i "/^FEMOCS_EXT_LIB=/ s|=|=$LIBPATH_TAITO |" share/makefile.femocs
    sed -i "/^FEMOCS_LIB=/ s|$|$LIB_TAITO |" share/makefile.femocs
    sed -i "/^FEMOCS_DLIB=/ s|$|$LIB_TAITO |" share/makefile.femocs

    sed -i "/MACHINE=/c\MACHINE=taito" share/makefile.femocs

    print_all

    echo -e "\nInstalling Femocs dependencies"
    make -f build/makefile.install
fi

if (test $mode = alcyone) then
    echo -e "\nWriting Alcyone flags"
    write_initial_flags

    sed -i "/^FEMOCS_EXT_HEAD=/ s|=|=$HEADPATH_ALCYONE |" share/makefile.femocs
    sed -i "/^FEMOCS_EXT_LIB=/ s|=|=$LIBPATH_ALCYONE |" share/makefile.femocs
    sed -i "/^FEMOCS_LIB=/ s|$|$LIB_ALCYONE |" share/makefile.femocs
    sed -i "/^FEMOCS_DLIB=/ s|$|$LIB_ALCYONE |" share/makefile.femocs

    sed -i "/MACHINE=/c\MACHINE=alcyone" share/makefile.femocs

    print_all

    echo -e "\nInstalling Femocs dependencies"
    make -f build/makefile.install
fi

if (test $mode = kale) then
    echo -e "\nWriting Kale flags"
    write_initial_flags

    sed -i "/^FEMOCS_EXT_HEAD=/ s|=|=$HEADPATH_KALE |" share/makefile.femocs
    sed -i "/^FEMOCS_EXT_LIB=/ s|=|=$LIBPATH_KALE |" share/makefile.femocs
    sed -i "/^FEMOCS_LIB=/ s|$|$LIB_KALE |" share/makefile.femocs
    sed -i "/^FEMOCS_DLIB=/ s|$|$LIB_KALE |" share/makefile.femocs

    sed -i "/MACHINE=/c\MACHINE=kale" share/makefile.femocs

    print_all

    echo -e "\nInstalling Femocs dependencies"
    make -f build/makefile.install
fi
