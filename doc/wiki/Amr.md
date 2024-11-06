
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


The parameters `ngridtot` and `nparttot` specify the maximum memory allocated for AMR grids
and collisionless particles respectively. These numbers should be greater than or equal to the
actual number of AMR grids and particles used during the course of the simulation.
`ngridtot` stands for the total number of AMR grids allocated over all MPI processes. The
`ngridmax` parameter can be used equivalently, but stands for the local number of AMR grids
within each MPI process. Obviously, one has `ngridtot=ngridmax*ncpu`.



:::{warning}
Recall that, in RAMSES, we call “AMR grid” or “oct” a group of $2^{\mathrm{ndim}}$ cells.
If for some reason, during the course of the execution, the maximum allowed
number of grids or particles has been reached, the simulation stops with the
message:
```
    No more free memory
    Increase ngridmax
```
In this case, don’t panic: just increase ngridmax in the Parameter File and
resume the execution, starting from the last backup file.
:::


## Memory management

The `ngridmax` and `npartmax` parameters (or there `*tot` equivalents) control the memory allocated by RAMSES. It is important to know how
much memory in Gb is actually allocated by RAMSES for a given choice of parameters. This
can be approximated by:
- $0.7(\mathtt{ngridmax}/10^6 ) + 0.7(\mathtt{npartmax}/10^7)$ Gbytes per cpu for pure N -body runs,
- $1.4(\mathtt{ngridmax}/10^6 ) + 0.7(\mathtt{npartmax}/10^7)$ Gb per cpu for N -body and hydro runs,
- $1.0(\mathtt{ngridmax}/10^6 )$ Gb per cpu for pure hydro runs.

Because of MPI communications overheads, the actual memory used can be slightly higher.
Note that these numbers are valid for double precision arithmetic. For single precision runs,
using the preprocessor directive `-DNPRE=4`, you can decrease these figures by 40%.

:::{warning}
The memory occupation will be higher when using modules that add new variables, such as MHD, Radiative transfer, Chemistry, Passive scalars, etc.
:::
