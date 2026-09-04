* Test name: `jiang-414-vel-refine`
* Dimension: `1`
* Solver: `mhd` (`CRPHYS=1 NCR_GROUPS=1 PATCH=../patch/cr_tests`)
* Setup: Same advected, diffusing CR pulse as `jiang-414-vel`, but with CR-gradient refinement (`err_grad_cr`) and `cr_nsubcycle=10`, so the refined region follows the drifting pulse
* Comparison: E_cr, F_cr and AMR level profiles of all outputs; `check_solution` regression of output_00002 against `jiang-414-vel-refine-ref.dat`
* Purpose: Testing CR advection and diffusion with a moving refined region
* Keywords: CR, advection, diffusion, AMR
* Reference: Jiang & Oh (2018), test 4.1.4
