* Test name: `tp-stream-va15`
* Dimension: `1`
* Solver: `mhd` (`CRPHYS=1 NCR_GROUPS=1 PATCH=../patch/cr_tests`)
* Setup: Two-pressure (gas + CR) shock tube as in `tp-nostream` (E_cr = 3 | 1 at x=5), but with CR streaming and streaming heating at an imposed Alfven speed `cr_v_alfven=1.5`; initial and boundary CR fluxes are set to the streaming equilibrium F = 4/3 v_A E_cr
* Comparison: density, velocity, P_gas, P_cr and F_cr profiles; `check_solution` regression of output_00001 against `tp-stream-va15-ref.dat`
* Purpose: Testing the coupled shock tube with faster CR streaming and streaming heating
* Keywords: CR, shock tube, streaming, coupling
