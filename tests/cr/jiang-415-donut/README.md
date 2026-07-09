* Test name: `jiang-415-donut`
* Dimension: `2`
* Solver: `mhd` (`CRPHYS=1 NCR_GROUPS=1 PATCH=../patch/cr_tests`)
* Setup: Same magnetic-loop setup as `jiang-415` (CR arc spreading along circular field lines) with first-order gas slopes (`slope_type=1`) and a single output at t=0.26, by which the arc has wrapped into a closed ring
* Comparison: E_cr map of output_00001; `check_solution` regression against `jiang-415-donut-ref.dat` (tolerance 3e-6)
* Purpose: Testing anisotropic CR transport along curved field lines (fully wrapped ring)
* Keywords: CR, anisotropic diffusion, magnetic loop
* Reference: Jiang & Oh (2018), test 4.1.5
