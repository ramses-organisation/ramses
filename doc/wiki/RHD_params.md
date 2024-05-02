

This sets of parameters, contained in the namelist block `&RT_PARAMS`, sets radiative properties for RAMSES RHD runs. Note that the number of photon groups (`NGROUPS`) is a compilation parameter, to be set in the Makefile. 

For detailed descriptions of the concepts described here, see the following papers:

* [[1] RAMSES-RT: radiation hydrodynamics in the cosmological context](http://arxiv.org/abs/1304.7126)
* [[2] A scheme for radiation pressure and photon diffusion with the M1 closure in RAMSES-RT](http://arxiv.org/abs/1411.6440)
* [[3] Galaxies that Shine](http://arxiv.org/abs/1501.04632)
* [[4] A simple model for molecular hydrogen chemistry coupled to radiation hydrodynamics](http://arxiv.org/abs/1802.00445)

 
| Variable name, syntax, default value | Fortran type  | Description       |
|:---------------------------- |:------------- |:------------------------- |
| `X=0.76`                     |  `real`       | Hydrogen mass fraction.|
| `Y=0.24`                     |  `real`       | Helium mass fraction.|
| `isH2=.false.`               |  `logical`    | Include molecular hydrogen? See [4]|
| `isHe=.true.`                |  `logical`    | Include helium ionization?|
| `rt_flux_scheme=’glf’`       |  `character(len=10)`    | Intercell flux function for the advection of radiation (see §3.2 in [1]). Use either ’glf’ or ’hll’. Note that only the GLF flux function is compatible with the inclusion of trapped IR radiation (see §2.2 in [2]). |
| `hll_evals_file=”`           |  `character(len=128)`    | Eigenvalues file, necessary only for the HLL intercell flux. Can also be set by environment variable `RAMSES_HLLFILE`. Such a file can be found with the source code (rt/hll_evals.list).|
| `rt_c_fraction=1.`           |  `real`    | Reduced light speed fraction, for keeping a manageable timestep-length (see §4.1 in [1]). The default corresponds to a full light speed.|
| `rt_courant_factor=0.8`      |  `real`    | Courant factor for photon advection between cells.|
| `rt_nsubcycle=1`             |  `integer` | Maximum number of RT-steps during one hydro/gravity/etc timestep.|
| `rt_smooth=.true.`           |  `logical` | Smooth out the operator splitting of photon advection and thermochemistry by incrementally updating the advected quantities in the chemistry. Usually speeds up the calculation. See §4.4 in [1].  |
| `rt_otsa=.true.`             |  `logical` | Assume the on-the-spot approximation, such that straight-to-ground recombinations in H and He do not emit ionising radiation (i.e. emitted photons are assumed to be instantly absorbed in the same cell -- see §3.3.2 in [1]).|
| `rt_is_init_xion=.false.`    |  `logical` | Set to initialize H and He ionization fractions from local photoionisation equilibrium (PIE), using the temperature, density, and radiation in each cell. Note that this is done by default (even if `rt_is_init_xion=.false.`) when starting simulations from scratch -- this parameter should only be used for resetting the ionisation fractions in _restarts_, which can be useful when postprocessing hydro simulations with radiative transfer. |
| `upload_equilibrium_x=.false.` |  `logical` | Set PIE ionization fractions when merging cells, instead of taking children averages. |
| `rt_is_outflow_bound=.false.` |  `logical` | Force outflow boundary for RT on all box sides, regardless of how boundaries are defined for hydrodynamics. By default, boundaries are the same for RT and hydrodynamics. |
| `rt_Tconst=-1`               |  `real`    | Constant temperature, in Kelvin, to assume for all temperature-dependent interaction rates (to run the first Iliev test). The default negative value means the actual cell temperature is used. |
| `rt_output_coolstats=.false.` |  `logical`    | Write thermochemistry statistics to the standard output.|
| `iPEH_group=0` |  `int`    | Photon group used for photoelectric heating (default no group)|
| `heat_unresolved_HII=0` |  `int`    | Subgrid model for unresolved HII regions (1==thermal heating, 2==non-thermal heating with NENER)|
| **====================** | **====================** | **Stellar emission parameters** | 
| `rt_star=.false.`            |  `logical` | Turn on photon emission from stellar particles. If `rt_star=.true.` and `rt=.true.` (in `RUN_PARAMS`), RT turns on when the first stellar particles are created in a simulation.  |
| `sed_dir=”`                  |  `character(len=128)` | Directory containing spectral energy distribution (SED) model for the stellar emission (see Appendix B in [1]). This can also be set by the environment variable `RAMSES_SED_DIR`.  |
| `sedprops_update=-1`         |  `integer` | Frequency (per coarse timestep) of photon group updates according to SED model (see Appendix B in [1]). The default value of -1 means that the update is never done.  |
| `nSEDgroups=NGROUPS`         |  `integer` | Number of photon groups dedicated to stellar emission, and relevant only if `rt_star=.true.`. These are the first photon groups: the last `NGROUPS-nSEDgroups` do not carry direct stellar radiation (are e.g. for a propagated UV background). |
| `SED_isEgy=.false.`          |  `logical` | Energy-conserving stellar emission when using a SED model. The default is photon number conserving within each photon group. With stellar particles of different ages and metallicities, the particle emission cannot be both energy conserving and number conserving at the same time, because the individual particles' spectral shapes differ, and we must use 'average' shapes for the photon groups. The  luminosity of a particle can either be correct (against the SED model) in terms of energy but not exactly by photon number, by taking the energy luminosity, or in terms of photon number but not exactly energy, by taking the photon number luminosity.   |
| `rt_esc_frac=1.`             |  `real` | 'Escape fraction' of photons from stellar particles, essentially just a multiplication factor for the particle emission.  |
| `convert_birth_times=.false.` |  `logical` | Convert birth times of stellar particles from conformal to proper time. By default the birth times are stored as proper in RHD simulations, for faster estimation of the stellar luminosity as a function of age. However, when postprocessing non-RHD simulations, this parameter is needed, since non-RHD simulations store the birth time in conformal units.  |
| **====================** | **====================** | **UV background parameters** | 
| `rt_UVsrc_nHmax=-1.`         |  `real` | Hydrogen density threshold (number per cubic cm) for non-homogeneous UV emission and propagation. Default value corresponds to no UV propagation. |
| `nUVgroups=0`          |  `integer` | Number of photon groups dedicated to the propagated UV background, and relevant (and set to `NGROUPS`) only if `rt_UVsrc_nHmax>0.`. These are the last photon groups: the first `NGROUPS-nUVgroups` do not carry the UV background (are e.g. for stellar radiation). |
| `uv_file=”`                  |  `character(len=128)` | File containing UV model, which is needed for a propagated UV background (`rt_UVsrc_nHmax>0`), or a homogeneous one (`haardt_madau=.true.` in `$PHYSICS_PARAMS`). This can also be set by the environment variable RAMSES_UV_FILE. |
| **====================** | **====================** | **Radiation pressure and IR radiation parameters** | 
| `rt_isIR=.false.`        | `logical`   |  Assume first photon group represents IR photons, in local thermal equilibrium (LTE) with dust (which is assumed to scale linearly in content with the gas metallicity). All other photon groups are 'reprocessed' locally and energy-conservatively into the IR group when interacting with dust via their `kappaAbs` parameter.|
| `rt_isIRtrap=.false.`        | `logical`   |  Apply trapping of IR photons in optically thick gas, in effect correctly modelling the IR propagation when the mean free path becomes unresolved (see §2.4.2 in [2]). With this set to `.true.`, the code must be compiled with a dedicated non-thermal energy variable, with `NENER>0` and a corresponding increment in `NVAR` (e.g. set `NENER=1` and increase `NVAR` by one, and recompile).|
| `is_kIR_T=.false.`        | `logical`   |  Use functions for the dust absorption and scattering opacities, such that they scale with the radiation temperature squared (eq. 79 in [2]), with a normalisation set by the `kappaAbs` and `kappaSc` parameters in `$RT_GROUPS`. The default is to use constant opacities (set by the same parameters).|
| `rt_isoPress=.false.`        | `logical`   |  Use the 'reduced flux approximation' described in [3] (§2.7), where the radiation in each cell is assumed to be fully directional.|
| `rt_pressBoost=1.`           | `real`      |  Multiplication factor to boost or reduce radiation pressure on gas and dust from the default.|
| `rt_vc=.false.`              | `logical`   |  Include relativistic corrections for the Lorentz boost and work done by photons on the gas (see Appendix A and B in [2]).|
| **====================** | **====================** | **RT refinement parameters** | 
| `rt_err_grad_n=-1.0`  `rt_err_grad_xHI=-1.0`   `rt_err_grad_xHII=-1.0` |  `real` | Discontinuity-based strategy: photon density and ionization fraction gradients above which a cell is refined. |
| `rt_floor_n=1d-10`    `rt_floor_xHI=1d-10`    `rt_floor_xHII=1d-10`    |  `real` | Discontinuity-based strategy: photon density and ionization fraction floor below which gradients are ignored. |
| `rt_refine_aexp=-1.0`        |  `real` | Cosmological expansion at which to turn on RT refinement strategies. |
| **====================** | **====================** | **RT source regions** | 
| `rt_nsource=0`               |  `integer` | Number of independent source (photon emission) regions in the computational box. |
| `rt_source_type=’square’`    |  `character(len=128) array` | Geometry defining each source region. ’square’ defines a generalized ellipsoidal shape with photons injected everywhere inside, ’shell’ defines a finite width spherical shell into which photons are injected, and ’point’ defines a point source. |
| `rt_src_x_center=0.0`  `rt_src_y_center=0.0`  `rt_src_z_center=0.0`    |  `real arrays` | Coordinates (0 to boxlen) of the center of each source region. |
| `rt_src_length_x=0.0` `rt_src_length_y=0.0` `rt_src_length_z=0.0`      |  `real arrays` | Sizes in all directions of each source region. If a spherical shell is used, rt_src_length_x and rt_src_length_y represent outer and inner radius. |
| `rt_exp_source=2.0`          |  `real array` | Exponents defining the norm used to compute distances for the generalized ellipsoid. `rt_exp_source=2` corresponds to a spheroid, `rt_exp_source=1` to a diamond shape, `rt_exp_source=10` to a perfect square. |
| `rt_src_group=1`            |  `integer array` | Photon groups into which photons are emitted in each source region (1 .le. rt_src_group .le. M, where M is the number of groups). |
| `rt_n_source=0.0` `rt_u_source=0.0` `rt_v_source=0.0` `rt_w_source=0.0`    |  `real arrays` | Injection rates, in cgs units, into photon densities and fluxes. |
| **====================** | **====================** | **RT initialisation regions** | 
| `rt_nregion=0`               |  `integer` | Number of independent initial radiation regions in the computational box. |
| `rt_region_type=’square’`    |  `character(len=128) array` | Geometry defining each initial radiation region. ’square’ defines a generalized ellipsoidal shape with photons initialised everywhere inside, ’shell’ defines a finite width spherical shell in which photons are initialised, and ’point’ defines a point flash. |
| `rt_reg_x_center=0.0`  `rt_reg_y_center=0.0`  `rt_reg_z_center=0.0`    |  `real arrays` | Coordinates (0 to boxlen) of the center of each initialisation region. |
| `rt_reg_length_x=0.0` `rt_reg_length_y=0.0` `rt_reg_length_z=0.0`      |  `real arrays` | Sizes in all directions of each initialisation region. If a spherical shell is used, rt_src_length_x and rt_src_length_y represent outer and inner radius. |
| `rt_reg_source=2.0`          |  `real array` | Exponents defining the norm used to compute distances for the generalized ellipsoid. `rt_reg_source=2` corresponds to a spheroid, `rt_reg_source=1` to a diamond shape, `rt_reg_source=10` to a perfect square. |
| `rt_reg_group=1`            |  `integer array` | Photon groups for which photons are initialised in each region (1 .le. rt_reg_group .le. M, where M is the number of groups). |
| `rt_n_region=0.0` `rt_u_region=0.0` `rt_v_region=0.0` `rt_w_region=0.0`    |  `real arrays` | Values in each region, in cgs units, for photon densities and fluxes. |