
# Non-ideal MHD Params

This set of parameters, contained in the namelist block `&NONIDEALMHD_PARAMS`, controls the non-ideal MHD (NIMHD) effects: **ambipolar diffusion** and **Ohmic dissipation**. The Hall effect is not currently implemented.

Non-ideal MHD requires the code to be compiled with `NIMHD=1` in the Makefile, which implies `SOLVER=mhd` and adds the `-DNIMHD` preprocessor flag. If the `&NONIDEALMHD_PARAMS` block is absent from the namelist, or if both `nambipolar` and `nmagdiffu` are `.false.`, the code runs in ideal MHD and none of the NIMHD machinery is active.

The non-ideal terms are added explicitly to the induction equation and to the total-energy flux inside the MHD Godunov update. Because they are treated explicitly, they introduce a diffusive stability constraint on the time step (`dt` scales like `dx^2`), which is reported at each fine step and can be tuned with the CFL-like coefficients below.

For detailed descriptions of the schemes and the resistivity model, see:

* [[1] Masson et al. 2012, Ambipolar diffusion in low-mass star formation](http://arxiv.org/abs/1208.1150)
* [[2] Marchand et al. 2016, Chemical solver to compute molecule and grain abundances and non-ideal MHD resistivities](http://arxiv.org/abs/1607.07658)

## Effect activation

| Variable name, syntax, default value | Fortran type | Description |
|:------------------------------------ |:------------ |:----------- |
| `nambipolar=.false.` | `logical` | Activate ambipolar diffusion. |
| `nmagdiffu=.false.`  | `logical` | Activate Ohmic dissipation (magnetic diffusion). |

## Resistivity model

| Variable name, syntax, default value | Fortran type | Description |
|:------------------------------------ |:------------ |:----------- |
| `resistivity_method=0` | `integer` | How the resistivities are computed. `0`: use a fixed, user-supplied coefficient (`gammaAD` for ambipolar diffusion, `etaMD` for Ohmic dissipation). `1`: analytical model (placeholder, not yet implemented). `2`: density-, temperature- and field-dependent resistivities read from a tabulated chemistry file (see `resistivity_file`). |
| `gammaAD=1.0` | `real` | Ambipolar diffusion coefficient used when `resistivity_method=0`, in CGS. The ambipolar coefficient is `beta = 1/(gammaAD*rho)`; the conversion to code units is handled internally. |
| `etaMD=1.0` | `real` | Fixed Ohmic (magnetic) diffusion coefficient `eta`, used when `resistivity_method=0`, in code units. |
| `resistivity_file="resnh.dat"` | `character(len=80)` | Path to the tabulated resistivity file, read only when `resistivity_method=2`. The file provides Pedersen, Ohm and Hall conductivities as a function of density and temperature. |
| `rho_threshold=1d-10` | `real` | Lower density floor (in code units) applied when evaluating the ambipolar/Ohmic resistivities. Prevents the diffusion coefficients (and hence the diffusion time step) from becoming unphysically large in very low-density regions. When `resistivity_method=2` it is raised to at least the minimum density of the table. |

## Non-ideal heating

The energy released by the non-ideal terms can be deposited either through the total-energy flux (conservative) or as an explicit source term. The two options are mutually exclusive; setting `nimhdheating_in_flux=.false.` automatically switches the source-term treatment on, and vice versa.

| Variable name, syntax, default value | Fortran type | Description |
|:------------------------------------ |:------------ |:----------- |
| `nimhdheating_in_flux=.true.` | `logical` | Add the non-ideal (ambipolar + Ohmic) heating to the total-energy flux. This is the default, conservative treatment. |
| `nimhdheating_source_term=.false.` | `logical` | Add the non-ideal heating as an explicit source term instead of via the flux. Mutually exclusive with `nimhdheating_in_flux`. |

## Time-step control

| Variable name, syntax, default value | Fortran type | Description |
|:------------------------------------ |:------------ |:----------- |
| `coefad=0.1` | `real` | CFL-like safety factor for the ambipolar diffusion time step (`dt_ad = coefad*dx^2/(B^2*beta)`). Smaller values give a more restrictive (safer) time step. |
| `coefohm=0.05` | `real` | CFL-like safety factor for the Ohmic dissipation time step (`dt_ohm = coefohm*dx^2/eta`). Smaller values give a more restrictive (safer) time step. |
| `nminitimestep=.false.` | `logical` | Activate the time-step "flooring" hack. When `.false.` (default, recommended for the test suite), the global time step is reduced to satisfy the full explicit diffusion CFL condition. When `.true.`, the diffusion time step is prevented from falling below a fraction of the ideal-MHD time step (set by `coefalfven`/`coefdtohm`), and the resistivity is rescaled accordingly to keep the update consistent. Use with care: it trades strict accuracy for speed. |
| `coefalfven=1d-10` | `real` | Only used when `nminitimestep=.true.`. Minimum allowed ratio between the ambipolar diffusion time step and the ideal-MHD time step. The default value effectively disables the flooring (no capping). |
| `coefdtohm=1d-10` | `real` | Same as `coefalfven`, but for the Ohmic dissipation time step. |
