* Test name: `g8-lowres`
* Dimension: `3`
* Solver: `mhd` (`DO_CR=1 NCR_GROUPS=1 NPSCAL=2`)
* Setup: G8-mass (M200 = 1e10 Msun) isolated disc galaxy from the hosted DICE `dwarf_galaxy` ICs (the same as `star/sfvirialtrue-sfmodel0`), with a ~1 muG toroidal magnetic field, star formation, thermal SN feedback with SN CR injection (`fecr=0.1`), CR transport (`Dcr=1e28` cm^2/s, reduced light speed ~300 km/s) and CR cooling. `cr_advect=.false.` initially: CR transport switches on at the first coarse step because the ICs contain star particles. Run to t = 0.5 code units (~7.5 Myr), `t_sne=3` Myr
* Comparison: `check_solution` regression of output_00003 (gas, CR and particle fields) against `g8-lowres-ref.dat`; tolerance 1e-13 on 2 MPI ranks, 1e-4 on selected fields otherwise
* Purpose: Testing SN CR injection, the CR transport auto-switch, and CR transport and cooling in a galaxy setup
* Keywords: CR, galaxy, supernova feedback
* Reference: RAMSES idealised disc galaxy tutorial ICs
