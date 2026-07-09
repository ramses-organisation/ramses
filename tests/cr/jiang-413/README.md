* Test name: `jiang-413`
* Dimension: `1`
* Solver: `mhd` (`CRPHYS=1 NCR_GROUPS=1 PATCH=../patch/cr_tests`)
* Setup: CRs stream through a static gas with a tanh density bump (0.1 -> 1.0 at x=200, width 25); the CR energy is pinned to E_cr=3 at the left reflexive wall and streams into the box, piling up where the Alfven speed drops in the dense slab (bottleneck effect)
* Comparison: density, F_cr and E_cr profiles of all outputs; `check_solution` regression of output_00003 against `jiang-413-ref.dat`
* Purpose: Testing CR streaming across a density (Alfven-speed) jump
* Keywords: CR, streaming, bottleneck
* Reference: Jiang & Oh (2018), test 4.1.3
