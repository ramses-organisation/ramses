

# Boundary Parameters

This set of parameters, contained in the namelist block `&BOUNDARY_PARAMS`, is used to set up boundary conditions on the current simulation. If this namelist block is absent, periodic boundary conditions are assumed. Setting up other types of boundary conditions in RAMSES is quite complex. The default setting, corresponding to a periodic box, should be sufficient in most cases. The strategy to set up boundary conditions is based on using "ghost regions" outside the computational domain, where flow variables are carefully specified in order to mimic the effect of the chosen type of boundary. Note that the order in which boundary regions are specified in the namelist is very important, especially for reflexive or zero gradient boundaries. Specific examples can be found in the namelist/ directory of the package.
 
| Variable name, syntax, default value | Fortran type | Description |
|:---------------------------- |:------------- |:------------------------- |
| `nboundary=1`  | `integer` | Number of ghost regions used to specify the boundary conditions.|
| `bound_type=0,0,0,` | `integer array` | Type of boundary conditions to apply in the corresponding ghost region. Possible values are: `bound_type=0` (periodic), `bound_type=1` (reflexive), `bound_type=2` (outflow, zero gradients), `bound_type=3` (inflow, user specified).|
| `d_bound=0.0` `u_bound=0.0` `v_bound=0.0` `w_bound=0.0` `p_bound=0.0` | `real arrays` | Flow variables in each ghost region (density, velocities and pressure). They are used only for inflow boundary conditions.|
| `ibound_min=0` `jbound_min=0` `kbound_min=0` | `integer arrays` | Coordinates of the lower, left, bottom corner of each boundary region. Each coordinate lies between -1 and +1 in each direction. |
| `ibound_max=0` `jbound_max=0` `kbound_max=0` | `integer arrays` | Likewise, for the upper, right and upper corner of each boundary region. |