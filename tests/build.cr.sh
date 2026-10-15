#!/bin/bash

cd bin
make NDIM=$1 SOLVER=mhd DO_CR=1 NCRGRP=1 MPI=1 DEBUG=1
