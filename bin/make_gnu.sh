#!/bin/bash

source /work/opt/local/bin/enable-cpe.sh

#rm -rf ~/work/GIZMO_Tables
#cp -r ~/GIZMO_Tables ~/work/

#module swich PrgEnv-cray PrgEnv-intel
#module add gsl/2.4_intel-18.0
#module add fftw3
#module add fftw/2.1.5.9
#module add gsl/115
#module add cray-hdf5/1.8.16

module swich PrgEnv-cray PrgEnv-gnu
#module add gsl/2.4_gnu-7.3
#module add fftw/2.1.5.9
#module add cray-hdf5/1.10.1.1


module list > m.log 2>& 1

make -f Makefile.xc >> m.log 2>& 1
