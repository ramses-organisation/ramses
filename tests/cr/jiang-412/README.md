* Test name: `jiang-412`
* Dimension: `2`
* Solver: `mhd` (`CRPHYS=1 NCR_GROUPS=1 PATCH=../patch/cr_tests`)
* Setup: 2D Gaussian CR-energy pulse in a static gas threaded by a uniform field inclined at 45 degrees (`A_region=B_region=1`); streaming diffusion and streaming heating enabled, so the pulse spreads anisotropically along the field
* Comparison: log E_cr maps of all outputs; `check_solution` regression of output_00002 against `jiang-412-ref.dat` (tolerance 3e-6)
* Purpose: Testing anisotropic CR streaming along an inclined magnetic field in 2D
* Keywords: CR, streaming, anisotropic
* Reference: Jiang & Oh (2018), test 4.1.2
