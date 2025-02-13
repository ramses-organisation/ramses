# Implementation details

## Keeping track of variables that are evolved on the grid

Grid variables are stored in the array `uold` and `unew`.
Several parameters are used to keep track of the number of different types of variables and their indices in the arrays.

### Overview of the variables

The Euler variables are the density, velocity and pressure.
The number of Euler variables is indicated by `neul` in the code.
We always have density and pressure, but the amount of velocities to keep track of depends on the number of dimensions of the simulation
and whether you use the HYDRO or MHD solver.
For HYDRO, we have `neul = ndim+2` with `ndim` the amount of spatial dimensions of the simulation.
For MHD simulations, we always need to keep track of the three velocities and so `neul = 5`.

When using the MHD solver, we need additional variables to keep track of the magnetic field in the three dimensions.
The magnetic field (on the left cell face) for the three spatial directions is added after the Euler variables.
The number of Euler variables with addition of the magnetic field is indicated by `nhydro` in the code.
This makes it easy to loop over them:
```
do i=1,nhydro
```
We have
- for HYDRO: `nhydro=neul`
- for MHD: `nhydro=neul+3`

Additionally, the magnetic field on the right cell face is added at the very end of the variable array.
This means that when following the evolution of magnetic fields, we need to store 6 additional variables.

Further, one can add more variables by setting the corresponding parameters in the Makefile
- `NPSCAL`: a number of passive scalars
- `NMETALS`: the number of metal species
- `NENER`: a number of non-thermal energies

The total amount of independent variables is indicated by `nvar`. It includes the hydro variables (Euler variables + 3 left magnetic field), metals, passive scalars and non-thermal energies:
```
NVAR := $(NHYDRO)+$(NENER)+$(NPSCAL)+$(NMETALS)
```
In addition, we defined `nvar_all` which in addition also includes the right magentic field.
So we have
- for HYDRO: `nvar_all = nvar`
- for MHD: `nvar_all = nvar+3`

### Indices
To access the Euler variables in uold and unew, use the in indices:
- density: `1`
- pressure: `neul`
- velocities (HYDRO case): `2` up to `1+ndim`
- velocities (MHD case): `2, 3, 4`

To access the magnetic field:
- on the left side of the cells: `neul + 1`, `neul+2`, `neul+3`
- on the right side of the cells: `nvar+1`, `nvar+2`, `nvar+3`


:::{admonition} TODO
:class: attention
indices of nener and passive scalars.
:::
