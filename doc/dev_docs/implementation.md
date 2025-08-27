# Hydrodynamics in RAMSES (implementation details)

## Variables on the grid

### Number of variables `nvar`
The total amount of independent variables stored on the AMR grid is set at compile time by the parameter `nvar`.
It includes
* the Euler variables: density, velocity and pressure
* the magnetic field vector (on the left cell face), when compiled with MHD
* an optional number of non-thermal energies `NENER`
* an optional number of metal species `NMETALS`
* an optional number of passive scalars `NPSCAL`
```
NVAR := $(NHYDRO)+$(NENER)+$(NPSCAL)+$(NMETALS)
```
:::{admonition} Makefile versus source code variables
:class: attention
`NMETALS` and `NPSCAL` are not passed to the source code. For the hydro-solver there is no distinction between these types of variables and both are treated as passive scalars.
:::

### Hydro data structures `uold` and `unew`
In RAMSES, all variables are stored together in the two-dimensional arrays `uold` and `unew`, which contain the value of each variable in each cell of the AMR grid.
They are defined in `hydro_commons` and allocated in `init_hydro`:
```
real(dp),allocatable,dimension(:,:)::uold,unew

allocate(uold(1:ncell,1:nvar))
allocate(unew(1:ncell,1:nvar))
```

### Accessing variables in `uold` and `unew`
Several parameters are used to keep track of the number of different types of variables and their indices in the arrays `uold` and `unew`.

The number of Euler variables is indicated by `neul` in the code.
We always have density and pressure, but the amount of velocities to keep track of depends on the number of dimensions of the simulation
and whether the code is compile with the HYDRO or MHD solver.
For HYDRO, we have `neul = ndim+2` with `ndim` the amount of spatial dimensions of the simulation.
For MHD simulations, we always need to keep track of the three velocities and so `neul = 5`.

:::{admonition} Conservative versus primitive variables
:class: attention
Remark that `uold` and `unew` actually contain the **conservative** Euler variables $\mathbb{U}$: density, momentum and total energy.
:::

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
This means that when following the evolution of magnetic fields, we store 6 additional variables.
For convenience, the code defines the parameter `nvar_all` which, in addition to the independent `nvar` variables, includes the right magentic field.
So we have
- for HYDRO: `nvar_all = nvar`
- for MHD: `nvar_all = nvar+3`

Passive scalars can be used to keep track of metals and star formation recipe variables. These specific scalars are accessed through the indices `imetal`, `idelay`, `ivirial1`, `ivirial2`, `ixion`, `ichem`.
For more info, see the section on star formation.
In the hydro-solver, these are evolved as regular passive scalars.

:::{admonition} Indices
:class: info
To access the Euler variables in `uold` and `unew`, use the indices:
- density: `1`
- total energy or pressure: `neul`
- momentum or velocities (HYDRO case): `2` up to `1+ndim`
- momentum or velocities (MHD case): `2, 3, 4`

To access the magnetic field:
- on the left side of the cells: `neul + 1`, `neul+2`, `neul+3`
- on the right side of the cells: `nvar+1`, `nvar+2`, `nvar+3`

Additional variables are stored after the left magnetic field in the following order:
- NENER: `inener=nhydro+1` up to `nhydro+nener`
- passive scalars: `nhydro+nener+1` up to `nvar`
:::
