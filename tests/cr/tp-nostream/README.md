* Test name: `tp-nostream`
* Dimension: `1`
* Solver: `mhd` (`CRPHYS=1 NCR_GROUPS=1 PATCH=../patch/cr_tests`)
* Setup: Two-pressure (gas + CR) Riemann problem: uniform gas (rho=1, P=1, B~0) with a CR-energy jump E_cr = 3 | 1 at x=5 and zero initial flux; the CR pressure gradient drives waves in the live gas. Pure CR diffusion (`Dcr=1/300`, streaming off), Lax-Friedrichs face flux (`cr_hlle=.false.`) and superluminal flux correction; imposed boundaries hold the initial states
* Comparison: density, velocity, P_gas, P_cr and F_cr profiles; `check_solution` regression of output_00001 against `tp-nostream-ref.dat`
* Purpose: Testing gas-CR momentum and energy coupling in a CR-driven shock tube without streaming
* Keywords: CR, shock tube, coupling
