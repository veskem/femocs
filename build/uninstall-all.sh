#!/bin/bash

read -p "WARNING! uninstall-all will among the other things remove Deal.II installation. Proceed? (y/n) " answer
if echo ${answer} | grep -iq "^y" ; then
    touch share/makefile.femocs
    set -x  # enable verbosity
    make -s -f build/makefile.install clean-all
    rm -rf share/.build
    rm -rf share/makefile.femocs
fi
