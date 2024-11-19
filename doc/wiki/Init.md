
# Initial Conditions Parameters

This sets of parameters, contained in the namelist block `&INIT_PARAMS`. This is used to set up the initial conditions.

| Variable name, syntax, default value | Fortran type | Description |
|:---------------------------- |:------------- |:------------------------- |
| `nregion=1`  | `integer` | Number of independent regions in the computational box used to set up initial flow variables. |
| `region_type=square`  | `10*char` | Geometry defining each region. `square` defines a generalized ellipsoidal shape, while `point` defines a delta function in the flow. |
| `x_center=0.0`  | `real arrays` | X coordinate of the center of each region. |
| `y_center=0.0`  | `real arrays` | Y coordinate of the center of each region. |
| `z_center=0.0`  | `real arrays` | Z coordinate of the center of each region. |
| `length_x=0.0`  | `real arrays` | Size in X direction of each region. |
| `length_y=0.0`  | `real arrays` | Size in Y direction of each region. |
| `length_z=0.0`  | `real arrays` | Size in Z direction of each region. |
| `exp_region=2.0`  | `real arrays` | Exponent defining the norm used to compute distances for the generalized ellipsoid. `exp_region=2` corresponds to a spheroid, `exp_region=1` to a diamond shape, `exp_region>=10` to a perfect square. |
| `d_region=0.0`  | `real arrays` | Density. For `point` regions this is used to define mass. |
| `u_region=0.0`  | `real arrays` | X velocity. For `point` regions this is used to define velocity. |
| `v_region=0.0`  | `real arrays` | Y velocity. For `point` regions this is used to define velocity. |
| `w_region=0.0`  | `real arrays` | Z velocity. For `point` regions this is used to define velocity. |
| `p_region=0.0`  | `real arrays` | Pressure. For `point` regions this is used to define specific pressure. |
| `filetype=ascii`  | `20*char` | Type of initial conditions file for particles. Possible choices are `ascii` or `grafic`. |
| `aexp_ini=10.0`  | `real` | This parameter sets the starting expansion factor for cosmology runs only. Default value is read in the IC file. |
| `multiple=.false.` | `logical` | If `.true.`, each processors reads its own IC file. For parallel runs only. |
| `initfile= ` | `80*char` | Directory where IC files are stored (when relevant).


## Advanced initial conditions

The `condinit` routine in `hydro/condinit.f90` can be modified to set custom initial conditions.
The calling sequence is `call condinit(x,u,dx,ncell)`, where

- `x` is an input array of cell center positions,
- `u` is an output array containing the volume average of the fluid conservative
variables, namely ($\rho$, $\rho u$, $\rho v$, $\rho w$ and $E$), in this exact order.
If more variables are defined, then the user should exploit this routine to define them too.
- `dx` is a single
real value containing the cell size for all the cells and ncell is the number of cells.

This routine
can be used to set the initial conditions directly with Fortran instructions.


## Input files

Another way to define initial conditions in RAMSES is by using input files (`initfile` parameter) in the grafic format.
