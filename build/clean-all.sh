#!/bin/bash

read -p "WARNING! clean-all will also remove Deal.II installation. Proceed? (y/n) " answer
if echo ${answer} | grep -iq "^y" ; then
    set -x  # enable verbosity
    rm -rf share/.build
    make -s -f build/makefile.lib clean-all
    make -s -f build/makefile.exec clean-all
    make -s -f build/makefile.test clean
    make -s -f build/makefile.install clean-all
    make -s -f build/makefile.cgal clean-all
    make -s -f build/makefile.doc clean-all
fi
