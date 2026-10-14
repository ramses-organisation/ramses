#!/bin/bash

cd bin
make NDIM=$1 SOLVER=mhd DO_CR=1 NCR_GROUPS=1 MPI=1 DEBUG=1
