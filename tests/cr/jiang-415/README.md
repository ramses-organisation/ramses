* Test name: `jiang-415`
* Dimension: `2`
* Solver: `mhd` (`CRPHYS=1 NCR_GROUPS=1 PATCH=../patch/cr_tests`)
* Setup: Static uniform gas threaded by a circular magnetic loop (B = curl A); a CR-energy enhancement (E_cr = 12 over a background of 10) on a small arc of the loop spreads along the curved field lines by anisotropic diffusion (`Dcr=0.1`, `Dcr_perp_factor=1e-6`)
* Comparison: E_cr maps of the last outputs; `check_solution` regression of output_00002 against `jiang-415-ref.dat` (tolerance 3e-6)
* Purpose: Testing anisotropic CR transport along curved field lines (hardest 2D transport test)
* Keywords: CR, anisotropic diffusion, magnetic loop
* Reference: Jiang & Oh (2018), test 4.1.5
