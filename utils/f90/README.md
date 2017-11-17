# Utils

Here is a quick reference of what the scripts in the folder do.

## Compiling

To compile, customize the `Makefile` to match your configuration (e.g. your compiler).

## amr2prof

Computes a radial profile out of the hydro files.

## amr2cylprof

Computes a radial profile (in a cylindrical coordinate system) of the hydro file.

## ramses2tipsy

??

## amr2map/amr2cube

Projects the hydro values onto a fixed 2D-grid/3D-grid.

## part2map/part2cube

Projects the particles values onto a fixed 2D-grid/3D-grid using CIC interpolation. Please note that the CIC is done at on the grid.

## partlist

Extracts the particles in a region of interest and dump them in a single binary file.
