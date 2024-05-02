

# Radiation Hydrodynamics in RAMSES

Radiation hydrodynamics are implemented in RAMSES, as described in those papers:

* [RAMSES-RT: radiation hydrodynamics in the cosmological context](http://arxiv.org/abs/1304.7126)
* [A scheme for radiation pressure and photon diffusion with the M1 closure in RAMSES-RT](http://arxiv.org/abs/1411.6440)
* [A simple model for molecular hydrogen chemistry coupled to radiation hydrodynamics](http://arxiv.org/abs/1802.00445)

## Compiling for RHD runs
An example makefile for an RHD compilation, `Makefile.rt`, is included under `ramses/bin`. To activate radiative transfer, you must compile with the flag `-DRT`, and some RHD specific .f90 files must be included in the compilation (all of which is included in the `Makefile` example).

For a run with only atomic hydrogen (i.e. ionisation species HI and HII), you should use `NIONS=1`. For including molecular hydrogen, add one to `NIONS`, and for including helium ionization, add two to `NIONS`. The number of hydro variables needs to increase accordingly, but this is done automatically in `Makefile.rt`. You also set the number of photon groups in the Makefile with the `NGROUPS` parameter. When running with IR radiation trapping (`rt_isIRtrap=.true.`), you must also compile with a dedicated non-thermal energy variable, by setting `NENER=1` in `Makefile.rt` (`NVAR` is incremented automatically in the `Makefile`). 

## RHD outputs
The radiation field variables are written separately in each output to files named `rt_XXXXX.outYYYYY`, where XXXXX is the output number and YYYYY the cpu number. The naming convention and format of the files is exactly the same as for the hydro variables. For each photon group, there are four cell variables, `c_r*N`, `Fx`, `Fy`, `Fz`, where `c_r` is the reduced light speed, `N` the photon number density, and the rest are the photon number flux in the x, y, and z directions. The factors for converting those to cgs units are stored in a file named `info_rt_XXXXX.txt` in each output directory (`unit_np` and `unit_pf`), along with the reduced light speed and the photon group properties. The non-thermal energy (trapped IR radiation pressure) is stored in runtime next after the thermal pressure, but in the output it comes just before.

## RHD Runtime parameters 
Radiative transfer is activated by setting `rt=.true.` in `&RUN_PARAMS`. For RT post-processing, set also `static_gas=.true.` in the same namelist. Note that this is not sufficient to turn on the advection of radiation -- this is done in the code (with `rt_advect=.true.`) only if sources of radiation (stars, gas, or idealised sources) are detected in the run.

If `rt=.true.`, non-equilibrium thermochemistry of hydrogen (and optionally helium) is used instead of the default equilibrium chemistry in RAMSES. The equilibrium chemistry can also be turned on without RHD, by setting `neq_chem=.true.` in `&PHYSICS_PARAMS` (but you must still compile with `Makefile.rt`, with `NGROUPS=0`).

For RHD runs, there are two additional dedicated namelists:

* [RT_PARAMS](./RHD_params)
* [RT_GROUPS](./RHD_groups)

## Generation of spectral energy distribution (SED) tables for RHD simulations
The (metallicity and age dependent) radiative luminosity of stellar population particles and photon group properties are calculated on-the-fly in RAMSES runs using SED tables. RAMSES-readable tables can be generated from Starburst99 and BC03 formats using a python utility found in `utils/py/sed_utils.py`. The generated files are linked to the RAMSES run with the `sed_dir` parameter in the `&RT_PARAMS` namelist.