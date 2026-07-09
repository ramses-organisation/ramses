* Test name: `jiang-414-refine`
* Dimension: `1`
* Solver: `mhd` (`CRPHYS=1 NCR_GROUPS=1 PATCH=../patch/cr_tests`)
* Setup: Same diffusing CR pulse as `jiang-414`, but with CR-gradient refinement (`err_grad_cr=0.03`) and `cr_nsubcycle=10`, so the pulse crosses coarse-fine level boundaries
* Comparison: E_cr, F_cr and AMR level profiles of all outputs; `check_solution` regression of output_00002 against `jiang-414-refine-ref.dat`
* Purpose: Testing CR diffusion across AMR level boundaries and `err_grad_cr` refinement
* Keywords: CR, diffusion, AMR
* Reference: Jiang & Oh (2018), test 4.1.4
