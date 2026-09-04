* Test name: `jiang-421`
* Dimension: `1`
* Solver: `mhd` (`CRPHYS=1 NCR_GROUPS=1 PATCH=../patch/cr_tests`)
* Setup: Sinusoidal CR-energy profile (E_cr = 20 + 10 sin(pi x)) in a periodic box with live gas; streaming (`cr_v_alfven=1`) plus diffusion (`Dcr=2`) with streaming heating, so the CR wave damps while pushing and heating the gas (`cr_nsubcycle=10`)
* Comparison: E_cr, F_cr, density and velocity profiles of all outputs; `check_solution` regression of output_00002 against `jiang-421-ref.dat`
* Purpose: Testing coupled CR-gas evolution with streaming heating
* Keywords: CR, streaming, coupling
* Reference: Jiang & Oh (2018), test 4.2.1
