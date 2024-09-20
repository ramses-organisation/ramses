

# Hydro parameters

This namelist is called &HYDRO_PARAMS, and is used to specify runtime parameters for the Godunov solver. These parameters are quite standard in computational fluid dynamics. We briefly
describe them now


| Variable name, syntax, default value | Fortran type  | Description               |
|:---------------------------- |:------------- |:------------------------- |
| `gamma=1.4`           |  `Real`    | Adiabatic exponent for the perfect gas EOS |
| `gamma_rad=1.333`     |  `Real`    | Adiabatic exponent for the non-thermal pressure EOS (if NENER>1) |
| `courant_factor=0.5`  |  `Real`    | CFL number for time step control (less than 1) |
| `smallr=1d-10 `       |  `Real`    | Minimum density to prevent floating exceptions |
| `smallc=1d-10 `       |  `Real`    | Minimum sound speed to prevent floating exceptions |
| `riemann=’llf’`       |`Character LEN=20`| Name of the desired Riemann solver. Possible choices are ’exact’, ’acoustic’, ’llf’, ’hll’ or ’hllc’ for the hydro solver and ’llf’, ’hll’, ’roe’, ’hlld’, ’upwind’ and ’hydro’ for the MHD solver. |
| `riemann2d=’llf’`     |`Characher LEN=20`| Name of the desired 2D Riemann solver for the induction equation (MHD only). Possible choices are ’upwind’, ’llf’, ’roe’, ’hll’, and ’hlld’. |
| `scheme=’muscl’`       |`Character LEN=20`| Name of the desired Godunov integrator. The hydro solver accepts ’muscl’ (MUSCL-HANCOCK scheme) or ’pldme’ (Collela's PLMDE scheme). The MHD solver accepts ’muscl’ or ’induction’. |
| `niter_riemann=10`    |  `Integer` | Maximum number of iterations used in the exact Riemann solver |
| `slope_type=1` |  `Integer`    | Type of slope limiter used in the Godunov scheme for the piecewise linear reconstruction: slope_type=0: First order scheme, slope_type=1: MinMod limiter, slope_type=2: MonCen limiter, slope_type=3: Multi-dimensional MonCen limiter. In 1D runs only, it is also possible to choose:slope_type=4: Superbee limiter, slope_type=5: Ultrabee limiter |
| `slope_mag_type=1` |  `Integer`    | Type of slope limiter used in the Godunov scheme in the MHD solver |
| `pressure_fix=.false.` |  `logical`    | Activate hybrid scheme (conservative or primitive) for high-Mach flows. Useful to prevent negative temperatures. |
| `beta_fix=0d0` |  `Real`    | With pressure_fix=.true., changes the threshold at which the energy is truncated |
| `difmag=0d0` |  `Real`    | Modifies diffusive flux in the PLDME integrator (hydro only). |
| `eta_mag=0d0` |  `Real`    | Modifies dtdiff in PLDME integrator (MHD only, warning, divide by zero error if eta_mag=0d0) |