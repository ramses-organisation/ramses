

# Refinement parameters #

The block named `&REFINE_PARAMS` contains the parameters related to grid refinement.

| Variable name | Fortran type | Default value  | Description               |
|:------------------- |:-------|:----- |:------------------------- |
| `x_refine`          | `real array` | 0.0   | Geometry-based strategy: center of the refined region at each level of the AMR grid.
| `y_refine`          | `real array` | 0.0   | Geometry-based strategy: center of the refined region at each level of the AMR grid.
| `z_refine`          | `real array` | 0.0   | Geometry-based strategy: center of the refined region at each level of the AMR grid.
| `r_refine`          | `real array` | 1e10   | Geometry-based strategy: **diameter** (yes, not a radius) of the refined region at each level.
| `a_refine`          | `real array` | 1.0   | Geometry-based strategy: ratio Y/X of the refined region at each level.
| `b_refine`          | `real array` | 1.0   | Geometry-based strategy: ratio Z/X of the refined region at each level.
| `exp_refine`        | `real array` | 2.0   | Geometry-based strategy: exponent of the norm.
| `jeans_refine`      | `real array` | -1.0   | Jeans refinement strategy: each level is refined if the cell size exceeds the local Jeans length divided by jeans_refine(ilevel).
| `collapse_jeans_refine` | `bool` | `.false.` | Use modified Jeans refine with prescribed thermodynamics of protostellar collapse
| `collapse_jeans_Tfloor` | `real` | 10.0d0  | Floor temperature (K) for the collapse_jeans_refine criterion
| `collapse_jeans_rho_low` | `real` | 1.0d-8  |Density (g/cc) below which collapse_jeans_refine uses Tfloor, typically at 2nd core formation
| `collapse_jeans_rho_high` | `real` | 1.0d-5   | Density (g/cc) above which collapse_jeans_refine uses the actual gas T, typically at protostar formation
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
| `interpol_var`      | `int`  | 0     | Variables used to perform interpolation (prolongation) and averaging (restriction). `interpol_var=0`: density, momentum and total energy; `interpol_var=1`: density, momentum and internal energy; `interpol_var=2`: density, velocity and internal energy (hydro solver only). See below for how to choose.
| `interpol_type`     | `int`  | 1     | Type of slope limiter used in the interpolation scheme for newly refined cells. `interpol_type=0`: Straight injection (1st order), `interpol_type=1`: MinMod limiter, `interpol_type=2`: MonCen limiter, `interpol_type=3`: unlimited central slope, `interpol_type=4`: type 3 for velocity and type 2 for density and internal energy (if `interpol_var=2`)
| `sink_refine`       | `bool` | `.false.` | Refines grid around sinks to levelmax

## Choosing `interpol_var` ##

Hydro variables have to be passed between AMR levels by two operations, both of which
run at every timestep and not only when the grid actually changes:

- **Prolongation**, from a parent cell to its children. A cell that is refined is split
  into `2**ndim` children, which are given values interpolated from the parent using the
  slope limiter selected by `interpol_type`. The same interpolation fills the buffer
  cells that a fine grid needs where its neighbour is a coarser cell, when filling the hydro stencil. The routine that takes care of this is `interpol_hydro`.
- **Restriction**, from children back to their parent. A split cell is given the average
  of its `2**ndim` children, so that the coarse value stays consistent with the finer
  solution underneath it. In MHD, the magnetic field on the face of an unrefined cell
  touching a refined neighbour is likewise replaced by the average of the fine faces, so
  that div(B)=0 still holds across the boundary. This is handled by the routine `upload_fine`.

`interpol_var` selects which set of variables these transfers act on, and
the choice is not neutral, because the energy of a cell is a sum of contributions:

```
E_total = E_thermal + E_kinetic + E_magnetic
```

Only one of these can be preserved exactly by the transfer.

- `interpol_var=0` transfers the **total** energy, so it is strictly conserved. The
  kinetic and magnetic energies are transferred separately, and since the operation
  changes them, whatever they gain or lose is taken out of the thermal energy. In other
  words the difference is silently converted into heat. This is the intended behaviour:
  it buys exact energy conservation, at the price of a small spurious heating or cooling
  near refinement boundaries.

- `interpol_var=1` transfers the **internal** energy instead, so the thermal energy, and
  hence the pressure, comes out unchanged. The total energy is then no longer conserved
  exactly. This is the better choice when the pressure matters more than strict energy
  conservation, typically when the magnetic or kinetic energy dominates the thermal one
  and the subtraction above would be badly conditioned.

- `interpol_var=2` behaves like `interpol_var=1` but transfers the velocity rather than
  the momentum. It is implemented in the hydro solver only; the MHD solver accepts `0`
  and `1`.

The MHD induction scheme (`scheme='induction'`) **requires** `interpol_var=1`, and the
code will stop otherwise. That scheme evolves only the magnetic field and re-imposes the
density, velocity and total energy analytically at every step, so there is no total
energy to conserve, and letting the magnetic energy leak into the thermal energy would
corrupt the pressure it is supposed to be holding fixed. The check is only applied when
the run actually refines, since `&REFINE_PARAMS` is not read for a uniform grid.
