* Test name: `jiang-422`
* Dimension: `1`
* Solver: `mhd` (`CRPHYS=1 NCR_GROUPS=1 PATCH=../patch/cr_tests`)
* Setup: Colliding flows (u = +10 / -10) carrying CRs (E_cr = 50 with flux F = 4/3 u E_cr in each half); the CRs are compressed in the forming shock and push back on the gas through the momentum and energy exchange terms
* Comparison: density, velocity, P_gas, P_cr and F_cr profiles of all outputs; `check_solution` regression of output_00003 against `jiang-422-ref.dat`
* Purpose: Testing gas-CR momentum and energy exchange in a shock
* Keywords: CR, shock, coupling
* Reference: Jiang & Oh (2018), test 4.2.2
