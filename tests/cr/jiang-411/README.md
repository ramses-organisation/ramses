* Test name: `jiang-411`
* Dimension: `1`
* Solver: `mhd` (`CRPHYS=1 NCR_GROUPS=1 PATCH=../patch/cr_tests`)
* Setup: Gaussian CR-energy pulse in a static uniform gas; streaming-dominated transport with an imposed Alfven speed `cr_v_alfven=1` and negligible explicit diffusion (`Dcr=1e-8`), reduced light speed = 100
* Comparison: E_cr and F_cr profiles of all outputs; `check_solution` regression of output_00002 against `jiang-411-ref.dat`
* Purpose: Testing two-moment CR streaming in 1D
* Keywords: CR, streaming
* Reference: Jiang & Oh (2018), test 4.1.1
