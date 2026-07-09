* Test name: `jiang-414`
* Dimension: `1`
* Solver: `mhd` (`CRPHYS=1 NCR_GROUPS=1 PATCH=../patch/cr_tests`)
* Setup: Static Gaussian CR-energy pulse in a uniform magnetized gas (periodic box), spreading by pure diffusion (`Dcr=0.1`, streaming disabled) on a uniform grid
* Comparison: E_cr and F_cr profiles of all outputs; `check_solution` regression of output_00002 against `jiang-414-ref.dat`
* Purpose: Testing two-moment CR diffusion in 1D
* Keywords: CR, diffusion
* Reference: Jiang & Oh (2018), test 4.1.4
