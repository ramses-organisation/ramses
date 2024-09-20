

# Runtime Parameters

The Ramses parameter file is based on the Fortran namelist syntax. The Sod test parameter file is shown below, as it should appear if you edit it.
```
$ more tube1d.nml
This is the RAMSES parameter file for Sod''s shock tube test.

&RUN_PARAMS
hydro=.true.
nsubcycle=3*1,2
/

&AMR_PARAMS
levelmin=3
levelmax=10
ngridmax=2000
nexpand=1
boxlen=1.0
/

&BOUNDARY_PARAMS
nboundary=2
ibound_min=-1,+1
ibound_max=-1,+1
bound_type= 1, 1
/

&INIT_PARAMS
nregion=2
region_type(1)='square'
region_type(2)='square'
x_center=0.25,0.75
length_x=0.5,0.5
d_region=1.0,0.125
u_region=0.0,0.0
p_region=1.0,0.1
/

&OUTPUT_PARAMS
noutput=1
tout=0.245
/

&HYDRO_PARAMS
gamma=1.4
courant_factor=0.8
slope_type=2
riemann='hllc'
/

&REFINE_PARAMS
err_grad_d=0.05
err_grad_u=0.05
err_grad_p=0.05
interpol_var=0
interpol_type=2
/

```

This parameter file is organized in namelist blocks. 
Each block starts with &BLOCK_NAME and ends with the character 
"/". Within each block, you can specify parameter values using
standard Fortran namelist syntax. There are currently 11 different 
parameter blocks implemented in RAMSES.

4 parameters blocks are mandatory and must always be present in the 
parameter file. These are `&RUN_PARAMS`, `&AMR_PARAMS`, `&OUTPUT_PARAMS`
and `&INIT_PARAMS`. The 8 other blocks are optional. They must be present in
the file only if they are relevant to the current run. These are  
`&BOUNDARY_PARAMS`, `&HYDRO_PARAMS`, `&PHYSICS_PARAMS`, `&POISSON_PARAMS`,
`&REFINE_PARAMS`,`&CLUMPFIND_PARAMS`, `&SINK_PARAMS` and `&MOVIE_PARAMS`.

The variables you can set or adjust in each namelist block are described
in more details in the following sections.

```{toctree}
:maxdepth: 1
:caption: Parameters:

Global.md
Amr.md
Init.md
Output.md
Boundary_conditions.md
Hydro_Parameters.md
Physics.md
Poisson.md
Refine.md
PHEW.md
unbinding.md
mergertree.md
Sinks.md
Movies.md
TurbulenceDriving.md
Tracers.md
```
