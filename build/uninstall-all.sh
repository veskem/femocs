#!/bin/bash

read -p "WARNING! uninstall-all will remove Deal.II & CGAL installation. Proceed? (y/n) " answer
if echo ${answer} | grep -iq "^y" ; then
    set -x  # enable verbosity
    rm -rf share/.build
    rm -rf share/makefile.femocs
    make -s -f build/makefile.install clean-all
    make -s -f build/makefile.cgal clean-all
fi
