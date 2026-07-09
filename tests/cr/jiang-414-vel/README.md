* Test name: `jiang-414-vel`
* Dimension: `1`
* Solver: `mhd` (`CRPHYS=1 NCR_GROUPS=1 PATCH=../patch/cr_tests`)
* Setup: Gaussian CR-energy pulse advected by a frozen background flow (`u_region=1`, `static_gas=.true.`) while diffusing (`Dcr=0.1`) on a uniform periodic grid; the pulse drifts across snapshots
* Comparison: E_cr and F_cr profiles of all outputs; `check_solution` regression of output_00002 against `jiang-414-vel-ref.dat`
* Purpose: Testing combined CR advection and diffusion
* Keywords: CR, advection, diffusion
* Reference: Jiang & Oh (2018), test 4.1.4
