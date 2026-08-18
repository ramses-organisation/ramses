#!/bin/bash

cd bin
make NDIM=$1 MPI=${MPI:-1} OPENMP=${OPENMP:-0} DEBUG=1
