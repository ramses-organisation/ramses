
# CR Params

This set of parameters, contained in the namelist block `&CR_PARAMS`, controls cosmic-ray transport in RAMSES. Note that the number of CR groups (`NCR_GROUPS`) is a compilation parameter, to be set in the Makefile together with `CRPHYS=1` (see the [documentation page about CR simulations](./CR)). Only `NCR_GROUPS=1` is currently supported; the Makefile rejects other values. Parameters marked `real array` have one entry per CR group.

For detailed descriptions of the scheme, see:

* [[1] A New Numerical Scheme for Cosmic-Ray Transport (Jiang & Oh 2018)](http://arxiv.org/abs/1712.07117)
* [[2] Fitz Axen et al. (2024), for the collisional cooling](https://arxiv.org/abs/2407.17597)
* [[3] Two-moment cosmic ray transport in RAMSES (Rosdahl et al. 2025)](https://arxiv.org/abs/2509.19447)

| Variable name, syntax, default value | Fortran type  | Description       |
|:---------------------------- |:------------- |:------------------------- |
| `cr_advect=.false.`          |  `logical`    | Master switch for CR transport. Nothing CR-related happens unless this is true. |
| `cr_HLLE=.true.`             |  `logical`    | Use the HLLE intercell flux function. If false, a local Lax-Friedrichs flux is used instead. |
| `cr_isotropic_pressure=.true.` |  `logical`  | Closure relation: isotropic CR pressure `P_cr = E_cr/3` ("P1"). If false, an M1-like Eddington tensor is used, with the anisotropic component aligned with the local magnetic field (see §2 in [1]). |
| `cr_flux_correction=.false.` |  `logical`    | Rescale superluminal CR fluxes back to the causality limit after the transport step. |
| `gamma_cr=1.3333`            |  `real array` | CR adiabatic index of each group (relativistic gas: 4/3). |
| **====================** | **====================** | **CR-gas coupling** |
| `cr_feedback=.true.`         |  `logical`    | Master switch for the CR back-reaction on the gas: momentum from the CR pressure gradient and the corresponding energy exchange. Inert when the gas is static (`static_gas=.true.`). |
| `cr_smallr_decouple=1e4`     |  `real`       | CRs are smoothly decoupled from the gas at low densities: source terms are multiplied by `exp(-smallr*cr_smallr_decouple/rho)`. |
| `cr_efloor=1e-30`            |  `real`       | Floor applied to the CR energy density (code units). |
| **====================** | **====================** | **Timestep and reduced light speed** |
| `cr_c_fraction=1.`           |  `real`       | Reduced light-speed fraction: the CR free-streaming speed is `V_m = cr_c_fraction * c`. Reduce it to keep manageable timestep lengths; must be in (0,1] (see [3]). |
| `cr_nsubcycle=1`             |  `integer`    | Maximum number of CR steps during one hydro timestep. The hydro timestep is limited to `cr_nsubcycle` CR Courant steps `dt_cr = courant_factor*dx/(3*V_m)`. |
| `cr_varvmax=.false.`         |  `logical`    | Adaptive reduced light speed: raise `V_m` per level so it stays at least `cr_varvmax_fudge` times the gas signal speed `dx/(3*dt)` (always capped at c; see [3]). |
| `cr_varvmax_fudge=10.`       |  `real`       | Safety factor for the adaptive reduced light speed. |
| `cr_varvmax_vdvs=.false.`    |  `logical`    | With `cr_varvmax`, also keep `V_m` above the streaming speed `gamma_cr*v_A` and the diffusion speed `Dcr/dx`. |
| **====================** | **====================** | **Scattering, diffusion and streaming** |
| `Dcr=1e29`                   |  `real array` | CR diffusion coefficient along the magnetic field, in cm^2/s. Sets the scattering coefficient `sigma_parallel = 1/(3*Dcr)` of the implicit source terms (eq. 10 in [1]). |
| `DCRmax=1e30`                |  `real`       | Upper limit, in cm^2/s, on the effective streaming-diffusion coefficient (which diverges where the CR gradient vanishes). |
| `Dcr_perp_factor=1e-6`       |  `real array` | Suppression factor of the diffusion coefficient perpendicular to the magnetic field: `Dcr_perp = Dcr_perp_factor * Dcr`. |
| `cr_streaming_diffusion=.false.` |  `logical` | Add the streaming contribution to the parallel scattering coefficient, so CRs stream down their pressure gradient along B at the Alfven speed (eq. 18 in [1]). |
| `cr_streaming_heating=.false.` |  `logical`  | Include the streaming velocity in the source-term velocity, which transfers the streaming work `v_A . grad(P_cr)` to the gas as heat. |
| `cr_v_alfven=0.`             |  `real`       | Imposed Alfven speed in code units, used by idealised tests. The default 0 uses the local `B/sqrt(rho)`. |
| `cr_f_taucell=1.`            |  `real`       | Prefactor of the cell scattering optical depth `tau` used to reduce the transport wave speeds in scattering-thick cells (the correction factor on p. 6 of [1]); lower values are more diffusive but more stable. |
| **====================** | **====================** | **Collisional cooling (Fitz Axen et al. 2024)** |
| `cr_cooling=.false.`         |  `logical`    | Enable CR collisional (Coulomb + hadronic) losses: exponential decay of `E_cr` and `F_cr` at the rate `lambda_cr*n_H`, with `lambda_cr = zeta_cr*(1 + 0.22*cr_ne + 0.125*cr_fneut)`. |
| `zeta_cr=7.51e-16`           |  `real`       | Loss-rate coefficient entering `lambda_cr` (cm^3/s). |
| `cr_ne=1e-3`                 |  `real`       | Number of free electrons per hydrogen nucleus. |
| `cr_fneut=0.875`             |  `real`       | Neutral gas fraction. |
| **====================** | **====================** | **CR initialisation regions** |
| `cr_nregion=0`               |  `integer`    | Number of CR initial-condition regions. The region geometry is independent of the gas `&INIT_PARAMS` regions but follows the same conventions. |
| `cr_region_type='square'`    |  `character array` | Geometry of each region: 'square' (generalized ellipsoid, CR state replaced inside) or 'point' (CIC-deposited point value). |
| `cr_reg_x_center=0.` `cr_reg_y_center=0.` `cr_reg_z_center=0.` |  `real arrays` | Coordinates (0 to boxlen) of the center of each region. |
| `cr_reg_length_x=1e10` `cr_reg_length_y=1e10` `cr_reg_length_z=1e10` |  `real arrays` | Size in all directions of each region. |
| `cr_exp_region=2.`           |  `real array` | Exponent of the norm defining the 'square' regions: 2 is a spheroid, 1 a diamond, >=10 a perfect square. |
| `cr_reg_group=1`             |  `integer array` | Target CR group of each region. Currently unused: `cr_region_u` addresses all group slots directly. |
| `cr_region_u=0.`             |  `real matrix` | CR state of each region, in code units: `cr_region_u(k,1)` is the CR energy density and `cr_region_u(k,2:ncrvar)` the CR flux of region k (per group: energy followed by ndim flux components). |
| **====================** | **====================** | **Boundaries** |
| `cr_boundary_u=0.`           |  `real matrix` | CR state imposed on each `bound_type=3` boundary region, applied by the default `cr_boundana`: `cr_boundary_u(b,1)` is the CR energy density and `cr_boundary_u(b,2:ncrvar)` the CR flux on boundary b (code units). |
| **====================** | **====================** | **Refinement** |
| `err_grad_cr=-1.0`           |  `real array` | Relative CR-energy gradient above which a cell is refined, per CR group. Negative values (the default) disable CR refinement. |
| **====================** | **====================** | **Idealised tests** |
| `jiang_test=''`              |  `character(len=32)` | Test-selection string dispatched by the `patch/cr_tests` initial/boundary conditions (Jiang & Oh 2018 suite and two-pressure shock tubes); empty for production runs. |
