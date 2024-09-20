

# Global parameters

This block, called `&RUN_PARAMS`, contains the run global control 
parameters. These parameters are now briefly described. 
More thorough explanations will be given in dedicated sections of the wiki.


| Variable name, syntax, default value | Fortran type  | Description               |
|:---------------------------- |:------------- |:------------------------- |
| `cosmo=.false.`              | `logical`     | Activate cosmological supercomoving cooordinates and computes the expansion factor |
| `pic=.false.`                |  `logical`    | Activate Particle-In-Cell solver |
| `sink=.false.`               |  `logical`    | Activate sink particles |
| `clumpfind=.false.`          |  `logical`    | Activate the clump finder |
| `tracer=.false.`             |  `logical`    | Use Monte Carlo tracer particles |
| `unbind=.false.`             |  `logical`    | Activate the particle unbinding for clumps |
| `make_mergertree=.false.`    |  `logical`    | Make merger trees |
| `poisson=.false.`            |  `logical`    | Activate Poisson solver for self-gravity |
| `hydro=.false.`              |  `logical`    | Activate hydro or MHD solver. |
| `rt=.false.`                 |  `logical`    | Activate radiative transfer using CPU-based M1 solver. This solver works on the AMR grid. |
| `aton=.false.`               |  `logical`    | Activate radiative transfer using GPU-based M1 solver. This solver works only on unigrid at `levelmin`. |
| `verbose=.false.`            |  `logical`    | Activate verbose mode. |
| `cost_weighting=.true.`      |  `logical`    | Load balancing based on computational cost, not memory. This is rather expensive in term of memory usage. For memory limited runs, using `cost_weighting=.false.` is better. |
| `nrestart=0`                 |  `integer`    | Output file number from which the code loads backup data and resumes the simulation, The default value, zero, is for a fresh start from the beginning (time=0).   |
| `nrestart_quad=0`                 |  `integer`    | Restart with double precision Hilbert keys. Must be equal to `nrestart`. Default value is 0.|
| `nstepmax=1000000`                 |  `integer`    | Maximum number of coarse time step. This can be used to terminate the simulation after a fixed amount of main steps. | 
| `ncontrol=1`                 |  `integer`    | Frequency of screen output for control lines (to stdout), usually redirected into a log file. | 
| `nremap=0`                   |  `integer`    | Frequency of call, in units of coarse time steps, for the load balancing routine, for MPI runs only. The default value, zero, means never. Load balancing is a slow operation, so use as high a value as possible. | 
| `ordering="hilbert"`         |  `character(len=128)`    | Cell ordering method used in the domain decomposition of the grid, for MPI runs only. Possible values are `hilbert`, `planar` and `angular`. |
| `nsubcycle=2,2,2,2,2,2,`     |  `integer array`    | Number of fine level sub-cycling steps within one coarser level time step. Each value in the array corresponds to a given level of refinement, starting from the coarse grid at `levelmin` up to the finest level at `levelmax`. `nsubcycle(1)=1` means that `levelmin` and `levelmin+1` are synchronized. To enforce single time stepping across the whole AMR hierarchy, you need to set `nsubcycle=1,1,1,1,1,1,1,`|
| `static=.false.`            |  `logical`    | Activate full static mode (RT post processing) |
| `static_dm=.false.`            |  `logical`    | Activate static mode for DM particles only (isolated initial conditions relaxation) |
| `static_stars=.false.`            |  `logical`    | Activate static mode for star particles only (isolated initial conditions relaxation) |
| `remap_pscalar=ndim+3,ndim+4,...,nvar`            |  `integer array`    | Mapping for the passive scalars and non-thermal pressures. Value indicates in which variable the scalar from the restart should be loaded. [0 = ignore this scalar in the restart output, -N = do not read but initialise ivar=N to 0, N = read and initialise ivar=N from the restart output]. For example `remap_pscalar=-6,0,7,8,9` translates to: do not read first restart var but initialise ivar=6 to 0, skip second restart var, read ivar=[8,9,10] from the restart snapshot and store it in the current ivar=[7,8,9].|