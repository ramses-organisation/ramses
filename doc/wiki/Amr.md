
# AMR Parameters

This sets of parameters, contained in the namelist block `&AMR_PARAMS`, controls the AMR grid global properties. Parameters specifying the refinement strategy are described [elsewhere](Refine), and the corresponding namelist block `&REFINE_PARAMS` is used only if `levelmax>levelmin`.

 
| Variable name, syntax, default value | Fortran type  | Description       |
|:---------------------------- |:------------- |:------------------------- |
| `levelmin=1`                 |  `integer`    | Minimum level of refinement. This parameter sets the size of the coarse (or base) grid to `nx=2**levelmin`.|
| `levelmax=1`                 |  `integer`    | Maximum level of refinement. If `levelmax=levelmin`, the simulation will be executed on a standad Cartesian grid of linear size `nx=2**levelmin`|
| `ngridmax=1`                 |  `integer`    | Maximum number of grids (or octs) that can be allocated during the run within each MPI process. |
| `ngridtot=1`                 |  `integer`    | Maximum number of grids (or octs) that can be allocated during the run for the entire set of MPI processes. One has `ngridmax=ngridtot/ncpu`.|
| `npartmax=1`                 |  `integer`    | Maximum number of particles of all types that can be allocated during the run within each MPI process. |
| `nparttot=1`                 |  `integer`    | Maximum number of particles of all types that can be allocated during the run for the entire set of MPI processes. One has `npartmax=nparttot/ncpu`.|
| `nsinkmax=1`                 |  `integer`    | Maximum number of sink particles during the run. Sink particles are duplicated over all MPI processes.|
| `nexpand=1`                  |  `integer`    | Number of times the mesh expansion is applied to the refinement map (see mesh smoothing).|
| `boxlen=1.0`                 |  `real`       | Box size in user's units |
| `nlevel_collapse=3`                 |  `integer`       | Number of levels above `levelmin` to follow initial halos collapse with a constant comoving resolution (`cosmo=.true.` only) |
