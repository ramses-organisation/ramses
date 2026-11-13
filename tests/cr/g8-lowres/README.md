* Test name: `g8-lowres`
* Dimension: `3`
* Solver: `mhd` (`DO_CR=1 NCRGRP=1 NPSCAL=2`)
* Setup: G8-mass (M200 = 1e10 Msun) isolated disc galaxy from the hosted DICE `dwarf_galaxy` ICs (the same as `star/sfvirialtrue-sfmodel0`), with a weak (~1e-3 muG) toroidal magnetic field that sets the CR diffusion direction, star formation, thermal SN feedback with SN CR injection (`fecr=0.1`), CR transport (`Dcr=1e28` cm^2/s, adaptive reduced light speed `cr_varvmax=.true.` with a ~300 km/s floor) and CR cooling. `cr_advect=.false.` initially: CR transport switches on at the first step because the ICs contain star particles. Low resolution for CI: levelmin 6 / levelmax 9 (781 pc to 98 pc in the 50 kpc box). Run to t = 0.5 code units (~7.5 Myr), `t_sne=3` Myr
* Comparison: `check_solution` regression of output_00003 (gas, CR and particle fields) against `g8-lowres-ref.dat`; tolerance 1e-13 on 2 MPI ranks, 1e-4 on selected fields otherwise
* Purpose: Testing SN CR injection, the CR transport auto-switch, and CR transport and cooling in a galaxy setup
* Keywords: CR, galaxy, supernova feedback
* Reference: RAMSES idealised disc galaxy tutorial ICs
