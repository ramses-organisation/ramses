

# RAMSES utils

RAMSES is shipped with a few Fortran programs to extract basic information from RAMSES outputs. 

## Compiling

Some of the programs can be compiled outside of the box. To compile them, customize the `Makefile` to match your configuration (e.g. your compiler), then run

	make

This will generate the binary files listed bellow. To cleanup, run

	make clean

Please also note that some script aren't included in the Makefile. If you find yourself needing one of them, please feel free to issue a Pull Request including the rule to generate it.

## amr2prof

Computes a radial profile out of the hydro files.

## amr2cylprof

Computes a radial profile (in a cylindrical coordinate system) of the hydro file.

## ramses2tipsy

Converts ramses data to gasoline (tipsy) data.

## amr2map/amr2cube

Projects the hydro values onto a fixed 2D-grid/3D-grid.

## part2map/part2cube

Projects the particles values onto a fixed 2D-grid/3D-grid using CIC interpolation. Please note that the CIC is done at on the grid.

## partlist

Extracts the particles in a region of interest and dump them in a single binary file.