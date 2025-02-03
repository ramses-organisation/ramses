#!/bin/bash

cd bin
make NDIM=$1 RT=1 NIONS=3 NGROUPS=3 MPI=1 DEBUG=1
