

In the version of RAMSES after RUM 2017 ([PR #284](https://bitbucket.org/rteyssie/ramses/pull-requests/284) and after), new blocks were introduced instead of one large `&PHYSICS_PARAMS`:

# Parameters

## Cooling parameters ##

The block named `&COOLING_PARAMS` contains the parameters related to cooling / basic chemistry

| Variable name | Fortran type | Default value  | Description               |
|:------------------- |:-------|:----- |:------------------------- |
| `cooling`                | `boolean` | `.false.`  | Enable gas atomic & metal cooling
| `metal`                  | `logical` | `.false.`  | Enable metals advection (requires 1 hydro variable)
| `isothermal`             | `logical` | `.false.`  | (deprecated) Enable equation of state for gas (heating and cooling disabled if `.true.`)
| `barotropic_eos`         | `logical` | `.false.`  | Enable barotropic equation of state for gas (heating and cooling disabled if `.true.`). Replaced 'isothermal'
| `barotropic_eos_form`    | `string`  | `legacy`   | Type of barotropic EOS. Options: isothermal, polytrope, double_polytrope, custom, legacy
| `polytrope_rho`          | `real`    | 1.0d50     | sets rho0 in EOS (density normalisation or knee-density), in g/cm3
| `polytrope_index`        | `real`    | 1.0d0      | sets gamma in EOS (polytropic index)
| `T_eos`                  | `real`    | 10         | sets T0 in EOS (isothermal temperature or temperature normalisation), in K
| `mu_gas`                 | `real`    | 1d0        | molecular weight
| `haardt_madau`           | `logical` | `.false.`  | Enable UV background
| `J21`                    | `real` | 0.0  | UV flux at threshold in 10^21 units
| `a_spec`                 | `real` | 1.0  | Slope of the UV spectrum
| `self_shielding`         | `logical` | `.false.`  | Enable self-shielding for densities above 0.01 g.cm^-3
| `z_ave`                  | `real` | 0.0  | Initial average metal abundance
| `z_reion`                | `real` | 8.5  | Reionization redshift (UV background disabled for higher redshifts)
| `ind_rsink`              | `real` | 4.0  | Number of cells defining the radius of the sphere where AGN feedback is active
| `T2max`                  | `real` | HUGE  | Temperature ceiling for gas heating (heating ceiling if `isothermal=.false.`) 
| `neq_chem`               | `logical` | `.false.`  | Enable non-equilibrium chemistry


## Star formation parameters ##

The block named `&SF_PARAMS` contains the parameters related to star formation

| Variable name | Fortran type | Default value  | Description               |
|:------------------- |:-------|:----- |:------------------------- |
| `m_star`                 | `real` | -1.0 | Star particle mass in units of mass_sph
| `n_star`                 | `real` | 0.1  | Star formation density threshold in H/cc
| `T2_star`                | `real` | 0.0  | Typical ISM polytropic temperature (cooling floor if `isothermal=.false.`)
| `g_star`                 | `real` | 1.6  | Typical ISM polytropic index (cooling floor if `isothermal=.false.`)
| `del_star`               | `real` | 2D2  | Star formation density threshold in critical density (`cosmo=.true.` only)
| `eps_star`               | `real` | 0.0  | Star formation efficiency
| `jeans_ncells`           | `real` | -1.0 | Jeans polytropic equation of state (cooling floor if `isothermal=.false.`)
| `sf_virial`              | `logical` | `.false.`  | Enable turbulent star formation prescriptions (requires 1 hydro variable)
| `sf_trelax`              | `real` | 0.0  | Relaxation time for star formation (`cosmo=.false.` only)
| `sf_tdiss`               | `real` | 0.0  | Dissipation timescale for subgrid turbulence in units of turbulent crossing time (`sf_virial=.true.` only)
| `sf_model`               | `integer` | 3  | Turbulent star formation model (`sf_virial=.true.` only)
| `sf_log_properties`    | `logical` | `.false.`  | Output gas properties of cells undergoing a star particle birth or supernova event
| `sf_compressive`         | `logical` | `.false.`  | Store and advect the compressive and solenoidal turbulent field in two different hydro variables (`sf_virial=.true.` only)

## Feedback parameters ##

The block named `&FEEDBACK_PARAMS` contains the parameters related to stellar feedback

| Variable name | Fortran type | Default value  | Description               |
|:------------------- |:-------|:----- |:------------------------- |
| `eta_sn`                 | `real` | 0.0  | Supernova mass fraction
| `eta_ssn`                 | `real` | 0.95| Single supernova ejected mass fraction (`sf_imf=.true.` only)
| `t_sne `                 | `real` | 10.0  | Supernova blast time in Myr
| `delayed_cooling`        | `logical` | `.false.`  | Enable delayed cooling through passive energy scalar advection (requires 1 hydro variable)
| `t_diss`                 | `real` | 20.0  | Dissipation timescale for supernova feedback in Myr (`delayed_cooling=.true.` only)
| `yield`                  | `real` | 0.0  | Supernova metal yield
| `mass_gmc`               | `real` | 0.0  | Stochastic exploding GMC mass in solar mass
| `kappa_IR`               | `real` | 0.0  | IR dust opacity for supernova feedback
| `f_ek`                   | `real` | 0.0  | Supernova kinetic energy fraction (only between 0 and 1)
| `f_w`                    | `real` | 0.0  | Supernova mass loading factor (`f_ek>0` only)
| `rbubble`                | `real` | 0.0  | Supernova superbubble radius in pc (`f_ek>0` only)
| `ndebris`                | `integer` | 1 | Supernova debris particle number (`f_ek>0` only)
| `mass_star_max`| `real` | 120.0 | Maximum mass of a star to undergo a supernova with `eta_ssn` efficiency in solar mass (`sf_imf=.true.` only)
| `mass_sne_min`| `real` | 10.0 | Minimum mass of a star to undergo a supernova blast in solar mass (`sf_imf=.true.` only)
| `ir_feedback`            | `logical` | `.false.`  | Enable IR feedback from accreting sink particles
| `ir_eff`                 | `real` | 0.75  | Efficiency of the IR feedback on sink particles (`ir_feedback=.true.`only)


## Units parameters ##

The block named `&UNITS_PARAMS` contains the parameters related to units that are set "by hand" (`cosmo=.false.` only!).

| Variable name | Fortran type | Default value  | Description               |
|:------------------- |:-------|:----- |:------------------------- |
| `units_density`          | `real` | 1.0   | Density unit in cgs (only for `cosmo=.false.`, requires G=1)
| `units_time`             | `real` | 1.0   | Time unit in cgs (only for `cosmo=.false.`, requires G=1)
| `units_length`           | `real` | 1.0   | Length unit in cgs (only for `cosmo=.false.`, requires G=1)


## GRACKLE parameters ##

Additionally, if the code has been compiled with GRACKLE=1 in the Makefile and properly linked against the hdf5 and grackle libraries, it is possible to define a `&GRACKLE_PARAMS` block to control the grackle parameters.
Please visit [https://grackle.readthedocs.io/en/grackle-3.0/Parameters.html](https://grackle.readthedocs.io/en/grackle-3.0/Parameters.html) for a more extensive description of the grackle parameters.

| Variable name | Fortran type | Default value  | Description               |
|:------------------- |:-------|:----- |:------------------------- |
| `use_grackle`                                    |`integer`|1| Activate Grackle
| `grackle_with_radiative_cooling`                 |`integer`|1| Include radiative cooling
| `grackle_primordial_chemistry`                   |`integer`|0| Primordial chemistry flag
| `grackle_metal_cooling`                          |`integer`|0| Enable metal cooling using Cloudy tables
| `grackle_UVbackground`                           |`integer`|0| Enable UV background
| `grackle_cmb_temperature_floor`                  |`integer`|1| Enable effective CMB temperature floor
| `grackle_h2_on_dust`                             |`integer`|0| Enable H2 formation on dust grains
| `grackle_photoelectric_heating`                  |`integer`|0| Enable a spatially uniform heating term approximating photo-electric heating
| `grackle_use_volumetric_heating_rate`            |`integer`|0| Signal that an array of volumetric heating rates is being provided
| `grackle_use_specific_heating_rate`              |`integer`|0| Signal that an array of specific heating rates is being provided
| `grackle_three_body_rate`                        |`integer`|0| Flag to control which three-body H2 formation rate is used
| `grackle_cie_cooling`                            |`integer`|0| Enable H2 collision-induced emission cooling
| `grackle_h2_optical_depth_approximation`         |`integer`|0| Enable H2 cooling attenuation
| `grackle_ih2co`                                  |`integer`|1|
| `grackle_ipiht`                                  |`integer`|1|
| `grackle_NumberOfTemperatureBins`                |`integer`|600|
| `grackle_CaseBRecombination`                     |`integer`|0|
| `grackle_Compton_xray_heating`                   |`integer`|0| Flag to enable Compton heating from an X-ray background
| `grackle_LWbackground_sawtooth_suppression`      |`integer`|0| Flag to enable suppression of Lyman-Werner flux due to Lyman-series absorption
| `grackle_NumberOfDustTemperatureBins`            |`integer`|250|
| `grackle_use_radiative_transfer`                 |`integer`|0| Signal that arrays of ionization and heating rates from radiative transfer solutions are being provided
| `grackle_radiative_transfer_coupled_rate_solver` |`integer`|0| Flag that must be enabled to couple the passed radiative transfer fields to the chemistry solver
| `grackle_radiative_transfer_intermediate_step`   |`integer`|0| Flag to enable intermediate stepping in applying radiative transfer fields to chemistry solver
| `grackle_radiative_transfer_hydrogen_only`       |`integer`|0| Flag to only use hydrogen ionization and heating rates from the radiative transfer solutions
| `grackle_self_shielding_method`                  |`integer`|0| Switch to enable approximate self-shielding from the UV background
| `grackle_Gamma`                                  |`real`|5./3.| The ratio of specific heats for an ideal gas
| `grackle_photoelectric_heating_rate`             |`real`|8.5D-26| The heating rate in units of erg cm-3 s-1
| `grackle_HydrogenFractionByMass`                 |`real`|0.76| 
| `grackle_DeuteriumToHydrogenRatio`               |`real`|2.0*3.4e-5|
| `grackle_SolarMetalFractionByMass`               |`real`|0.01295|
| `grackle_TemperatureStart`                       |`real`|1.0|
| `grackle_TemperatureEnd`                         |`real`|1.0D9|
| `grackle_DustTemperatureStart`                   |`real`|1.0|
| `grackle_DustTemperatureEnd`                     |`real`|1500.|
| `grackle_LWbackground_intensity`                 |`real`|0.0| Intensity of a constant Lyman-Werner H2 photo-dissociating radiation field in units of 10-21 erg s-1 cm-2 Hz-1 sr-1
| `grackle_UVbackground_redshift_on`               |`real`|7.0|
| `grackle_UVbackground_redshift_off`              |`real`|0.0|
| `grackle_UVbackground_redshift_fullon`           |`real`|6.0|
| `grackle_UVbackground_redshift_drop`             |`real`|0.0|
| `grackle_cloudy_electron_fraction_factor`        |`real`|9.153959D-3|
| `grackle_data_file`                              |`character`|| Path to the data file containing the metal cooling and UV background HDF5 tables

---

## Physics parameters (LEGACY ONLY) ##

The block named `&PHYSICS_PARAMS` contains the parameters related to physical models.

| Variable name | Fortran type | Default value  | Description               |
|:------------------- |:-------|:----- |:------------------------- |
| `units_density`          | `real` | 1.0   | Density unit in cgs (only for `cosmo=.false.`, requires G=1)
| `units_time`             | `real` | 1.0   | Time unit in cgs (only for `cosmo=.false.`, requires G=1)
| `units_length`           | `real` | 1.0   | Length unit in cgs (only for `cosmo=.false.`, requires G=1)
| `cooling`                | `boolean` | `.false.`  | Enable gas atomic & metal cooling
| `isothermal`             | `logical` | `.false.`  | Enable equation of state for gas (heating and cooling disabled if `.true.`)
| `metal`                  | `logical` | `.false.`  | Enable metals advection (requires 1 hydro variable)
| `haardt_madau`           | `logical` | `.false.`  | Enable UV background
| `self_shielding`         | `logical` | `.false.`  | Enable self-shielding for densities above 0.01 g.cm^-3
| `smbh`                   | `logical` | `.false.`  | Enable super massive black holes for sink particles
| `agn `                   | `logical` | `.false.`  | Enable AGN feedback from super massive black holes
| `neq_chem`               | `logical` | `.false.`  | Enable non-equilibrium chemistry
| `ir_feedback`            | `logical` | `.false.`  | Enable IR feedback from accreting sink particles
| `sf_virial`              | `logical` | `.false.`  | Enable turbulent star formation prescriptions (requires 1 hydro variable)
| `sf_compressive`         | `logical` | `.false.`  | Store and advect the compressive and solenoidal turbulent field in two different hydro variables (`sf_virial=.true.` only)
| `sf_log_properties`    | `logical` | `.false.`  | Output gas properties of cells undergoing a star particle birth or supernova event
| `delayed_cooling`        | `logical` | `.false.`  | Enable delayed cooling through passive energy scalar advection (requires 1 hydro variable)
| `sf_imf`                 | `logical` | `.false.` | Model Initial Mass Function during thermal feedback events
| `T2_star`                | `real` | 0.0  | Typical ISM polytropic temperature (cooling floor if `isothermal=.false.`)
| `g_star`                 | `real` | 1.6  | Typical ISM polytropic index (cooling floor if `isothermal=.false.`)
| `jeans_ncells`           | `real` | -1.0 | Jeans polytropic equation of state (cooling floor if `isothermal=.false.`)
| `T2max`                  | `real` | 1D50  | Temperature ceiling for gas heating (heating ceiling if `isothermal=.false.`) 
| `n_star`                 | `real` | 0.1  | Star formation density threshold in H/cc
| `m_star`                 | `real` | -1.0 | Star particle mass in units of mass_sph
| `del_star`               | `real` | 2D2  | Star formation density threshold in critical density (`cosmo=.true.` only)
| `eps_star`               | `real` | 0.0  | Star formation efficiency
| `t_star`                 | `real` | 0.0  | Star formation time scale in Gyr (only used if `eps_star>0`)
| `sf_trelax`              | `real` | 0.0  | Relaxation time for star formation (`cosmo=.false.` only)
| `sf_tdiss`               | `real` | 0.0  | Dissipation timescale for subgrid turbulence in units of turbulent crossing time (`sf_virial=.true.` only)
| `sf_model`               | `integer` | 3  | Turbulent star formation model (`sf_virial=.true.` only)
| `eta_sn`                 | `real` | 0.0  | Supernova mass fraction
| `eta_ssn`                 | `real` | 0.95| Single supernova ejected mass fraction (`sf_imf=.true.` only)
| `t_sne `                 | `real` | 10.0  | Supernova blast time in Myr
| `t_diss`                 | `real` | 20.0  | Dissipation timescale for supernova feedback in Myr (`delayed_cooling=.true.` only)
| `yield`                  | `real` | 0.0  | Supernova metal yield
| `mass_gmc`               | `real` | 0.0  | Stochastic exploding GMC mass in solar mass
| `kappa_IR`               | `real` | 0.0  | IR dust opacity for supernova feedback
| `f_ek`                   | `real` | 0.0  | Supernova kinetic energy fraction (only between 0 and 1)
| `f_w`                    | `real` | 0.0  | Supernova mass loading factor (`f_ek>0` only)
| `rbubble`                | `real` | 0.0  | Supernova superbubble radius in pc (`f_ek>0` only)
| `ndebris`                | `integer` | 1 | Supernova debris particle number (`f_ek>0` only)
| `mass_star_max`| `real` | 120.0 | Maximum mass of a star to undergo a supernova with `eta_ssn` efficiency in solar mass (`sf_imf=.true.` only)
| `mass_sne_min`| `real` | 10.0 | Minimum mass of a star to undergo a supernova blast in solar mass (`sf_imf=.true.` only)
| `J21`                    | `real` | 0.0  | UV flux at threshold in 10^21 units
| `a_spec`                 | `real` | 1.0  | Slope of the UV spectrum
| `z_ave`                  | `real` | 0.0  | Initial average metal abundance
| `z_reion`                | `real` | 8.5  | Reionization redshift (UV background disabled for higher redshifts)
| `ind_rsink`              | `real` | 4.0  | Number of cells defining the radius of the sphere where AGN feedback is active
| `ir_eff`                 | `real` | 0.75  | Efficiency of the IR feedback on sink particles (`ir_feedback=.true.`only)
