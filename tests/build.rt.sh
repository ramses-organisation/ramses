#!/bin/bash

cd bin
cat Makefile.rt |\
    sed 's/NDIM = 1/NDIM = 3/g' |\
    sed 's/NGROUPS = 0/NGROUPS = 3/g' |\
    sed 's/NIONS = 0/NIONS = 3/g' |\
    sed 's/^F90 =.*$/F90 = mpif90 -frecord-marker=4 -ffree-line-length-none -fbacktrace -g -O -fbounds-check -Wuninitialized -Wall/g' |\
    sed 's/^FFLAGS =.*$/FFLAGS = -x f95-cpp-input -ffpe-trap=zero,underflow,overflow,invalid -finit-real=nan  $(DEFINES)/g' > Makefile_new.rt
make -f Makefile_new.rt clean

make -f Makefile_new.rt
