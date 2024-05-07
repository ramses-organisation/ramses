

# Poisson Parameters

This namelist, `&POISSON_PARAMS`, is used to specify runtime parameters for the Poisson solver. It is used only if `poisson=.true.` or `pic=.true.`


Two different Poisson solvers are available in RAMSES: conjugate gradient (CG) and multigrid (MG). Unlike the CG solver, MG has an initialization overhead cost (at every call of the solver), but is much more efficient on very big levels with few "holes". The multigrid solver is therefore used for all coarse levels.

In addition, MG can be used on refined levels in conjuction with CG. The parameter `cg_levelmin` selects the Poisson solver as follows:

* Coarse levels are solved with MG
* Refined levels with *l* < `cg_levelmin` are solved with MG
* Refined levels with *l* >=  `cg_levelmin` are solved with CG

| Variable name | Fortran type | Default value  | Description      |
|:------------------- |:-------|:----- |:------------------------- |
| `gravity_type`      | `int`  | 0     | Type of gravity force. Possible choices are: `gravity_type=0` self-gravity (Poisson solver); `gravity_type>0` analytical gravity vector; `gravity_type<0` self-gravity plus additional analytical density profile
| `epsilon`           | `real`  | 1e-4  | Stopping criterion for the iterative Poisson solver: residual 2-norm should be lower than `epsilon` times the right hand side 2-norm.
| `gravity_params`    | `real array`  | 0.0, | Parameters used to define the analytical gravity field (routine `gravana.f90`) or the analytical mass density field (routine `rho_ana.f90`).
| `cg_levelmin`       | `integer`  | 999 | Minimum level from which the Conjugate Gradient solver is used in place of the Multigrid solver.
| `cic_levelmax`      |	`integer`  | 0	 | Maximum level for CIC dark matter interpolation (default `cic_levelmax=nlevelmax`)