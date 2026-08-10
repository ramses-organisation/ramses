
# Initial Conditions Parameters

This sets of parameters, contained in the namelist block `&INIT_PARAMS`. This is used to set up the initial conditions.

| Variable name, syntax, default value | Fortran type | Description |
|:---------------------------- |:------------- |:------------------------- |
| `condinit_kind=region` | `60*char` | Set which kind of initial condition to use. Default is the region-based initialisation |
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
| `A_region=0.0`  | `real arrays` | X magnetic field. |
| `B_region=0.0`  | `real arrays` | Y magnetic field. |
| `C_region=0.0`  | `real arrays` | Z magnetic field. |
| `filetype=ascii`  | `20*char` | Type of initial conditions file for particles. Possible choices are `ascii`, `grafic`, or `dice`. |
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

:::{versionchanged} 2025
There is now the possibility to change initial condition without recompiling each time.
Just write you initial condition as a new routine and add it to the select/case structure
```
  select case (condinit_kind)

  case('region')
     call condinit_default(x, u, dx, nn)

  case('your_new_routine')
     call condinit_your_new_routine(x, u, dx, nn)

  case DEFAULT
     if (myid == 1.and. first_call)  write(*,*) "[condinit] Void or invalid condinit_kind, using default IC"
     call condinit_default(x, u, dx, nn)

  end select
```
:::


## Input files

Another way to define initial conditions in RAMSES is by using input files (`initfile` parameter) in the grafic format.


## DICE initial conditions

(Inspired by the old DICE wiki)

RAMSES offers the possibility to load Gadget1 and Gadget2 initial conditions created with the [DICE code](https://ascl.net/1607.002). The code reads any kind of particles, and transfer them to the RAMSES particle tree. If gas particles are included, the patch will use a NGP scheme to transfer their mass to the corresponding AMR cells. Once the full gas mass has been transferred to the AMR grid, the gas particles are erased from the particle list.

Dedicated parameters are specified through the `&DICE_PARAMS` namelist block:

| Variable name, syntax, default value | Fortran type | Description |
|:---------------------------- |:------------- |:------------------------- |
| `ic_file` | `512*char` | Name of the initial conditions file |
| `ic_nfile=1` | `integer` |  If greater than one, look for files with name matching ic_file//'.n' |
| `ic_ifout=1` | `integer` |  Change ramses output index for restarts |
| `ic_format` | `512*char` |  Format of the initial conditions. Should be 'Gadget1' or 'Gadget2' |
| `ic_center=0.0,0.0,0.0` | 3*`real` |  Shift center parameter. ICs are automatically shifted with `boxlen/2` |
| `ic_scale_pos=1.0` | `real` |  Scaling factor for the position vector  |
| `ic_scale_vel=1.0` | `real` |  Scaling factor for the velocity vector |
| `ic_scale_mass=1.0` | `real` |  Scaling factor for the mass |
| `ic_scale_u=1.0` | `real` |  Scaling factor for the internal energy |
| `ic_scale_age=1.0` | `real` |  Scaling factor for the particles age |
| `ic_scale_metal=1.0` | `real` |  Scaling factor for the metallicity |
| `ic_head_name='HEAD'` | `4*char` |  Name of the Header datablock (Gadget2 format only) |
| `ic_pos_name='POS '` | `4*char` |  Name of the position vector datablock (Gadget2 format only) |
| `ic_vel_name='VEL '` | `4*char` |  Name of the velocity vector datablock (Gadget2 format only) |
| `ic_mass_name='MASS'` | `4*char` |  Name of the mass datablock (Gadget2 format only) |
| `ic_id_name='ID  '` | `4*char` |  Name of the particle identifier datablock (Gadget2 format only) |
| `ic_u_name='U   '` | `4*char` |  Name of the internal energy datablock (Gadget2 format only) |
| `ic_metal_name='Z   '` | `4*char` |  Name of the metallicity datablock (Gadget2 format only) |
| `ic_age_name='AGE '` | `4*char` |  Name of the particle age datablock (Gadget2 format only) |
| `IG_rho=1.0D-6` | `real` |  Intergalactic gas density |
| `IG_T2=1.0D6` | `real` |  Intergalactic gas temperature |
| `IG_metal=1.0D-4` | `real` |  Intergalactic gas metallicity |
| `amr_struct=.false.` | `logical` |  Reproduce the AMR structure of the Gadget2 file resulting from a ramses to gadget conversion |
| `gadget_scale_l=3.085677581282D21` | `real` | Gadget file length unit in cgs |
| `gadget_scale_v=1.0D5` | `` |  Gadget file velocity unit in cgs |
| `gadget_scale_m=1.9891D43` | `` |  Gadget file mass unit in cgs |
| `gadget_scale_t=3.15360e+13` | `` |  Gadget file time unit in cgs |
| `ic_skip_type        = -1` | 6*`integer` |  Skip specific particle type |
| `cosmo_add_gas_index = -1` | 6*`integer` |  Gas particle type index for cosmo runs |
| `ic_mask_ivar = 0` | `integer` | Some mask refinement parameters |
| `ic_mask_min = 1D40` | `real` | Some mask refinement parameters |
| `ic_mask_max = -1D40` | `real` | Some mask refinement parameters |
| `ic_mask_ptype = -1` | `integer` | ??? |
| `ic_mag_const=0.0,0.0,0.0` | 3*`real` |  Background magnetic field value for x,y,z component |
| `ic_mag_center_x=0.0` | 32*`real` |  x-coordinate of the magnetic field symmetry center |
| `ic_mag_center_y=0.0` | 32*`real` |  y-coordinate of the magnetic field symmetry center |
| `ic_mag_center_z=0.0` | 32*`real` |  z-coordinate of the magnetic field symmetry center |
| `ic_mag_axis_x=0.0` | 32*`real` |  Magnetic field normal vector x-component  |
| `ic_mag_axis_y=0.0` | 32*`real` |  Magnetic field normal vector y-component  |
| `ic_mag_axis_z=1.0` | 32*`real` |  Magnetic field normal vector z-component  |
| `ic_mag_scale_R=1.0` | 32*`real` |  Toroidal magnetic field scalelength |
| `ic_mag_scale_H=1.0` | 32*`real` |  Toroidal magnetic field scaleheight |
| `ic_mag_scale_B=0.0` | 32*`real` |  Foreground toroidal magnetic field value |

When using DICE ICs, the code units have to be specified accordingly in the `&UNITS_PARAMS` block. The following parameters are recommended for simplicity:

```
units_density = 0.677025430198932E-22 ! 1e9 Msol/kpc^3
units_time    = 0.470430312423675E+15 ! G=1
units_length  = 0.308567758128200E+22 ! kpc
```

Similarly, the `&INIT_PARAMS` block should contain `filetype='dice'` to indicate that we want to read DICE ICs. Finally, the usual isolated galaxy setup will use isolated boundary conditions:

```
&BOUNDARY_PARAMS
nboundary=6
bound_type= 2, 2,  2,  2,  2,  2
ibound_min=-1, 1, -1, -1, -1, -1
ibound_max=-1, 1,  1,  1,  1,  1
jbound_min= 0, 0, -1,  1, -1, -1
jbound_max= 0, 0, -1,  1,  1,  1
kbound_min= 0, 0,  0,  0, -1,  1
kbound_max= 0, 0,  0,  0, -1,  1
no_inflow=.true.
/
```
