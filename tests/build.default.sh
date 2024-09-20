#!/bin/bash

cd bin
make NDIM=$1 MPI=1 DEBUG=1
