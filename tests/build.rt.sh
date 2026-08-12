#!/bin/bash

cd bin
make NDIM=$1 RT=1 NIONS=3 NGROUPS=3 MPI=${MPI:-1} OPENMP=${OPENMP:-0} DEBUG=1
