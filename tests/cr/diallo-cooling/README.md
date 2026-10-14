* Test name: `diallo-cooling`
* Dimension: `1`
* Solver: `mhd` (`DO_CR=1 NCR_GROUPS=1`)
* Setup: Gaussian CR-energy profile (peak 2.5e-10 erg/cm^3, sigma = 10 pc) in a 1 kpc periodic box of static gas at n_H = 1 cm^-3, with collisional CR cooling (`cr_cooling=.true.`) and negligible transport (`Dcr=1e20` cm^2/s); outputs at 0, 1, 10 and 100 Myr
* Comparison: at every output, E_cr against the analytic decay E_0(x) exp(-lambda_cr n_H t) (max error < 1e-3 of the peak) and the total CR energy against exp(-lambda_cr n_H t) (< 1e-9); `check_solution` regression of output_00004 against `diallo-cooling-ref.dat`
* Purpose: Testing collisional CR cooling, including the n_H unit conversion
* Keywords: CR, cooling
* Reference: Diallo et al. (2026), tests 3.1.1 and 3.1.4 (single-group analogue)
