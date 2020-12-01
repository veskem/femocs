#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
source $DIR/makefile.defs

mode=$1

if (test $mode = ubuntu) then
    LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$DIR/../dealii/lib:$DIR/../lib
    export LD_LIBRARY_PATH
fi

if (test $mode = taito) then
    echo "Loading Taito modules"

    source $MODULESHOME/init/bash
    module load gcc/5.3.0 intelmpi/5.1.3

    echo -e "\nExporting TBB path"
    LD_LIBRARY_PATH=${LIBPATH_TAITO:2}:${LD_LIBRARY_PATH}
    export LD_LIBRARY_PATH
fi

if (test $mode = alcyone) then
    echo "Loading Alcyone modules"

    source $MODULESHOME/init/bash
    module load PrgEnv-gnu gcc/5.1.0
fi

if (test $mode = kale) then
    echo "Loading Kale modules"

    source $MODULESHOME/init/bash
    module load fgci-common/1.0
    module load boost/1.70.0
    module load tbb/2019_U4-GCCcore-8.2.0
fi

