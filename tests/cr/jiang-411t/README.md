* Test name: `jiang-411t`
* Dimension: `1`
* Solver: `mhd` (`CRPHYS=1 NCR_GROUPS=1 PATCH=../patch/cr_tests`)
* Setup: Triangular CR-energy profile streaming at an imposed Alfven speed `cr_v_alfven=1` in a static gas; the expanding plateau is driven by time-dependent analytic boundaries (`bound_type=3`, `cr_boundana`)
* Comparison: E_cr and F_cr profiles overlaid on the analytic plateau solution; `check_solution` regression of output_00001 against `jiang-411t-ref.dat`
* Purpose: Testing CR streaming against an exact analytic solution
* Keywords: CR, streaming, analytic
* Reference: Jiang & Oh (2018), test 4.1.1 (triangular variant)
