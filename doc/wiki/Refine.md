

# Refinement parameters #

The block named `&REFINE_PARAMS` contains the parameters related to grid refinement.

| Variable name | Fortran type | Default value  | Description               |
|:------------------- |:-------|:----- |:------------------------- |
| `x_refine`          | `real array` | 0.0   | Geometry-based strategy: center of the refined region at each level of the AMR grid.
| `y_refine`          | `real array` | 0.0   | Geometry-based strategy: center of the refined region at each level of the AMR grid.
| `z_refine`          | `real array` | 0.0   | Geometry-based strategy: center of the refined region at each level of the AMR grid.
| `r_refine`          | `real array` | 1e10   | Geometry-based strategy: radius of the refined region at each level.
| `a_refine`          | `real array` | 1.0   | Geometry-based strategy: ratio Y/X of the refined region at each level.
| `b_refine`          | `real array` | 1.0   | Geometry-based strategy: ratio Z/X of the refined region at each level.
| `exp_refine`        | `real array` | 2.0   | Geometry-based strategy: exponent of the norm.
| `jeans_refine`      | `real array` | -1.0   | Jeans refinement strategy: each level is refined if the cell size exceeds the local Jeans length divided by jeans_refine(ilevel).
| `mass_cut_refine`   | `real` | -1.0   | Mass threshold for particle-based refinement
| `m_refine`          | `real array` | -1.0   | Quasi-Lagrangian strategy: each level is refined if the baryons mass in a cell exceeds `m_refine(ilevel)*mass_sph`, or if the number of dark matter particles exceeds `m_refine(ilevel)`, whatever the mass is.
| `mass_sph`          | `real` | 0.0   | Quasi-Lagrangian strategy: `mass_sph` is used to set a typical mass scale. For cosmo runs, its value is set automatically.
| `err_grad_d`        | `real` | -1.0  | Discontinuity-based strategy: density gradient relative variations above which a cell is refined
| `err_grad_u`        | `real` | -1.0  | Discontinuity-based strategy: velocity gradient relative variations above which a cell is refined
| `err_grad_p`        | `real` | -1.0  | Discontinuity-based strategy: pressure gradient relative variations above which a cell is refined
| `floor_d`           | `real` | 1e-10 | Density floor below which gradients are ignored
| `floor_u`           | `real` | 1e-10 | Velocity floor below which gradients are ignored
| `floor_p`           | `real` | 1e-10 | Pressure floor below which gradients are ignored
| `ivar_refine`       | `int`  | -1    | Variable index for refinement
| `var_cut_refine`    | `real` | -1.0  | Threshold for variable-based refinement
| `interpol_var`      | `int`  | 0     | Variables used to perform interpolation (prolongation) and averaging (restriction). `interpol_type=0`: conservatives; `interpol_type=1`: primitives
| `interpol_type`     | `int`  | 1     | Type of slope limiter used in the interpolation scheme for newly refined cells. `interpol_type=0`: Straight injection (1st order), `interpol_type=1`: MinMod limiter, `interpol_type=2`: MonCen limiter, `interpol_type=3`: unlimited central slope, `interpol_type=4`: type 3 for velocity and type 2 for density and internal energy (if `interpol_var=2`)
| `sink_refine`       | `bool` | `.false.` | Refines grid around sinks to levelmax