# FLD

## Namelist parameters overview

| Variable name, syntax, default value | Fortran type  | Description               |
|:---------------------------- |:------------- |:------------------------- |
| `grey_rad_transfer=.true`           |  `logical`    | Use the grey approximation, meaning there is only 1 radiation bin covering the entire spectral range from 0 to infinity. Only allowed when `NGRP=1`. |
| `freqs_in_Hz=.true.`| `logical`    | If `.true.`, boundary parameter values are provided as frequencies in Hz. Otherwise, they are given as energies in eV. |
| `numin=1.0d5` | `real` |  The lower boundary of the first radiation group. |
| `numax=1.0d19` | `real` |  The upper boundary of the last radiation group. |
| `split_groups_log=.true.` | `logical`    | Generate logaritmically (`.true.`) or linearly (`.false.`) spaced radiation bins between `numin` and `numax`.  |
| `extra_end_group=.false.` | `logical`    | The last radiation group goes from `numax` to infinity, instead of stopping at `numax`.  |
| `read_groups=.false.` | `logical`    | If `.true.`, read radiation group boundaries from a file named 'groups.dat', instead of generating them automatically.  |
| `rosseland_params=1.0` | array of `real`| Rosseland opacity coefficient's parameters. The opacity is calculated as rosseland_params(1)$\times$ density $^{\mathrm{rosseland\_params}(2)}$ $\times$ Temperature $^{\mathrm{rosseland\_params(3)}}$. Output is given in cm $^{-1}$.
| `planck_params=1.0` | array of `real`| Planck opacity coefficient's parameters. Same as rosseland_params |
| `sublimation_kuiper=.false.` | `logical`    |  Mimicks dust sublimation with decreasing d/g ratio, see Kuiper+10 ApJ 
| `dt_control=.false.` | `logical` | If `.true.`, activate the control of the global timestep. It is used in FLD pure diffusion tests, without `hydro`.  |
`dtdiff_params=1d10`  | array of `real` | If `dt_control=.false.`, the global timestep is set as dtnew(ilevel)=dtdiff_params(1)*dtdiff_params(2)**nstep_coarse
|`epsilon_diff=1d-6`   | `real` | Convergence for the iterative method used in the implicit scheme (CG or BiCGstab). By default, two criteria have to be satisified for convergence L $_\infty$ on the radiative energy variation between two iterations and on the ratio of the L $_2$ norm of the residual over the L $_2$ norm of the initial residual.
|`fld_limiter='nolim'`  | `character(LEN=10)` | Choice of the FLD flux limiter: nolim, levermore or minerbo.
|`Tr_floor=10`          | `real` | Background radiation field temperature - WARNING: it affects the pressure_fix in `set_uold` (`if(Tp_loc .lt. 0.3d0*Tr_floor)`). It is also used to compute `scale_E0=a_R*Tr_floor` which renormalises the radiative nergy in the implicit solver (to have x in Ax=b close to values of 1).
|`robin=1.0d0`          | `real`  | Robin type boundary condition between AMR levels for FLD implicit solver (see Commercon et al. 2014 for details). `robin=0.0` corresponds to Von Neumann BC (flux is stored between levels and total radiative energy is conserved) - this is not robust since energy is stored at wrong locations in the grid. `robin=1.0` corresponds to Dirichlet BC (total radiative energy is NOT conserved) - much more robust. 
|`min_optical_depth=1.d-6`  | `real`| Set the minimum optical depth in the cell. It may accelerate convergence of (Bi)CG in optically thin regions. It also ensures a better conditioning of the matrix to invert.
|`Tray_min=0.5d0`       |`real` | Minimum temperature of the radiativve energy. `eray_min=aR*Tray_min=0.5d0**4` is coomputed and is used for radiative energies as` small_r` is used for density.
|`energy_fix=.false.`   | `logical`| Advect internal energy as an extra variable, automaticlally put a t index `nvar`. PUT a +1 to NVAR in the Makefile if you activate `energy_fix`.
| `store_matrix=.true.` | `logical`| Store the matrix coeeficient in the iterative (Bi)CG solver. It worsk well for low number of FLD groups. For large numbers, it allocates to much memroy. 
|`fit_semenov=.false.` | `logical`|   Use a fit of the Semenov et al. (2003) Rosseland anf Planck dust opacities
|

## Setting radiation groups

The FLD implementation can be used either with the grey approximation ([Commercon et al 2011](https://ui.adsabs.harvard.edu/abs/2011A%26A...529A..35C/abstract), [Commercon et al 2014](https://ui.adsabs.harvard.edu/abs/2014A%26A...563A..11C/abstract)) or with the multigroup implementation ([Gonzalez et al 2015](https://ui.adsabs.harvard.edu/abs/2015A%26A...578A..12G/abstract)).

In the grey approximation, we have only one radiation group which covers the entire spectral range from zero to infinity. To use this approximation, compile with `NGRP=1` and set `grey_rad_transfer=.true.` in the namelist. In this case, a Conjugate Gradient algorithm is used in the implicit scheme that updates the radiative energy. 

With the multigroup implementation, a Stabilised Bi-Conjugate Gradient algorithm is used in the implicit scheme that updates the radiative energy. In multigrooup, we can choose as many radiation groups as we want (or as memory constraints allow). The lower and upper limits of each group are determined from the following namelist parameters:
* freqs_in_Hz
* numin and numax
* split_groups_log and extra_end_group
* read_groups

The group boundary parameters can be provided by the user as either frequencies in Hz (`freqs_in_Hz=.true.`) or energies in eV (`freqs_in_Hz=.false.`).

If there is only a single radiation group (`NGRP=1`), the lower and upper frequency (or energy) can be set with the parameters `numin` and `numax`.

If there is more than one group (`NGRP>1`), you can either
* generate group bounds automatically based on the namelist parameters
* provide a list of boundaries in a file named 'groups.dat' (read_groups=.true.)

In the first case, the parameters `split_groups_log` will determine if the groups are splitted logaritmically (`split_groups_log=.true.`) or linearly (`split_groups_log=.false.`). The boundaries will be generated with equal spacing between `numin` and `numax`.

In some cases, you want to specify an additional radiation bin at the end, containing all radiation beyond the group boundary `numax`. This can be done by `setting extra_end_group=.true.`. Enabling this will create NGRP-1 equally spaced radiation groups between `numin` and `numax`, and add an extra group at the end of the list with boundaries `numax` to infinity. Remark that in practice, infinity is hardcoded to be `frequency_upperlimit=1.0d35`.

For example, the set of namelist parameters
TODO
will create the following groups:
TODO

Providing the group boundaries with a file allows for more flexibility. An example of a user provided groups.dat file and corresponding namelist:
TODO


## Opacities

TODO

## Timestep considerations

TODO

## Numerical parameters for the algorithm

## Radiation feedback from sinks