* Test name: `jiang-423`
* Dimension: `2`
* Solver: `mhd` (`CRPHYS=1 NCR_GROUPS=1 PATCH=../patch/cr_tests`)
* Setup: 2D Sedov-type gas blast (central point energy region) sweeping a uniform CR background (E_cr = 1e-10 floor), with CR-gradient refinement (`err_grad_cr=0.05`), adaptive reduced light speed (`cr_varvmax=.true.`) and `cr_nsubcycle=10`
* Comparison: density, E_cr, pressure and AMR level maps of the last output; `check_solution` regression of output_00001 against `jiang-423-ref.dat`
* Purpose: Testing CR transport in a 2D blast with AMR and adaptive reduced light speed
* Keywords: CR, blast, AMR
* Reference: Jiang & Oh (2018), test 4.2.3
