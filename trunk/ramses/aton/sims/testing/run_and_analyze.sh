#!/bin/bash

mkdir test
cd test

# This takes several minutes.
mpirun -n 1 ../ramses/bin/ramses3d ../sims/testing/stromgren.nml > log

# Get some interesting data.
../ramses/utils/f90/amr2cell -inp output_00007/ -out output_00007/cells.txt
python ../ramses/utils/py/test5spherical.py < output_00007/cells.txt > profile.txt

cd ..

echo ""
echo ""
echo "********************************************************************************"
echo "Now use gnuplot to compare test/profile.txt to sims/testing/expected_profile.txt"
echo "The columns are r, density, xneutral, pressure, temperature, mach, xion"
echo "e.g. to compare the temperature in gnuplot:"
echo "set log y"
echo "plot 'test/profile.txt' using 1:5, 'sims/testing/expected_profile.txt' using 1:5"

