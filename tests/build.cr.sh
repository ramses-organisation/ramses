#!/bin/bash

cd bin
make NDIM=$1 SOLVER=mhd NCR=1 CRFLX=1 PATCH=../patch/cr_tests MPI=1 DEBUG=1
