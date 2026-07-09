* Test name: `jiang-424`
* Dimension: `1`
* Solver: `mhd` (`CRPHYS=1 NCR_GROUPS=1 PATCH=../patch/cr_tests`)
* Setup: CRs stream into a live gas with a tanh density bump (0.1 -> 10 at x=200); the CR energy is pinned to E_cr=3 at the left reflexive wall and its pressure gradient drives an outflow through the dense slab. Streaming + heating enabled, CR-gradient refinement (`err_grad_cr=0.05`), adaptive reduced light speed (`cr_varvmax=.true.`, base reduced light speed = 10) and superluminal flux correction
* Comparison: density, P_th, E_cr, F_cr and AMR level profiles of all outputs; `check_solution` regression of output_00003 against `jiang-424-ref.dat`
* Purpose: Testing a CR-driven outflow with AMR and adaptive reduced light speed
* Keywords: CR, streaming, outflow, AMR
* Reference: Jiang & Oh (2018), test 4.2.4
