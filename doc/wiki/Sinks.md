

# Sink particles #

The block named `&SINK_PARAMS` contains the parameters related to the sink particle implementation, which can be used for

* star formation [(most recently Bleuler & Teyssier 2015)](http://arxiv.org/pdf/1409.6528v1)
* supermassive black hole evolution and AGN feedback [(most recently Biernacki, Teyssier and Bleuler 2017)](https://arxiv.org/abs/1701.05190)

## Overview of parameters ##

| Variable name         | Fortran type | Default value | Description               |
|:----------------------|:------------ |:------------- |:------------------------- |
| `smbh`                | `boolean`    | `.false.`     | Controls if sink behaves as a star or SMBH
| `agn`                 | `boolean`    | `.false.`     | Controls if SMBH sink produces feedback
| `create_sinks`        | `boolean`    | `.false.`     | Specifies if sinks are formed with clump finder
| `nsinkmax`            | `integer`    | `2000`        | Maximum number of sinks allowed to form
| `mass_sink_direct_force` | `float`   | `-1`          | Mass in Msun above which sinks are treated with direct N-body solver
| `ir_cloud`            | `integer`    | `4`           | Radius of cloud region in unit of grid spacing
| `ir_cloud_massive`    | `integer`    | `3`           | Radius of massive cloud region in unit of grid spacing
| `sink_soft`           | `integer`    | `2`           | Gravitational softening length in dx at levelmax for "direct force" sinks
| `n_sink`              | `float`      | `.false.`     | Sink (as a star) formation density threshold in H/cc
| `rho_sink`            | `float`      | `.false.`     | Sink (as a star) formation density threshold in g/cc
| `d_sink`              | `float`      | `.false.`     | Sink (as a star) formation density threshold in code units
| `clump_core`          | `boolean`    | `.false.`     | Trims the clump (for sinks as stars only)
| `mass_sink_seed`      | `float`      | `0.0`         | Dynamical mass of sink seed in Msun
| `mass_smbh_seed`      | `float`      | `0.0`         | Accretion mass of sink seed in Msun, if 0, then dynamical and accretion masses are equivalent
| `mass_halo_AGN`       | `float`      | `1e10`        | Mass of a halo in which AGN sinks are seeded
| `mass_clump_AGN`      | `float`      | `1e10`        | Mass of a clump in which AGN sinks are seeded
| `accretion_scheme`    | `string`     | `none`        | Accretion scheme: none, bondi
| `eddington_limit`     | `boolean`    | `.false.`     | Controls if accretion rate should be limited by Eddington rate
| `acc_sink_boost`      | `float`      | `1`           | Value of boost factor in Bondi velocity (-1 to depend on density)
| `boost_threshold_density` | `float`  | `0.1`         | Threshold density for boosting, typically the same as n_star; in H/cc
| `bondi_use_vrel`      | `boolean`    | `.false.`     | Excludes relative velocity between gas and sink from the accretion calculations
| `mass_merger_vel_check_AGN` | `float` | `-1`         | Mass above which the check for two sinks' binding energy is applied, in Msun
| `merging_timescale `  | `float`      | `-1`          | Time during which sinks are considered for merging (non-SMBH only) 
| `verbose_AGN`         | `boolean`    | `.false.`     | Controls verbosity of AGN in the log file
| `AGN_fbk_mode_switch_threshold` | `float` | `0.1`    | Controls the AGN feedback switching
| `AGN_fbk_frac_ener`   | `float`      | `1.0`         | Controls what fraction of energy goes into thermal feedback
| `T2_min`              | `float`      | `1e7`         | Minimum temperature to which AGN blast should heat the gas in K
| `T2_AGN`              | `float`      | `1e12`        | Temperature of AGN blasts in K - feedback efficiency
| `AGN_fbk_frac_mom`    | `float`      | `1.0`         | Controls what fraction of energy goes into kinetic feedback
| `epsilon_kin`         | `float`      | `1.0`         | Kinetic feedback coupling efficiency
| `kin_mass_loading`    | `float`      | `100.`        | Kinetic feedback mass loading
| `cone_opening`        | `float`      | `180.`        | Opening angle of the cone through which momentum feedback proceeds



## Example set of parameters for cosmological simulations with AGN feedback ##

The example listed below by no means can fit everyone's needs, but can serve as a minimum starting example.

```{code-block} fortran
#!fortran
&SINK_PARAMS
! General switches
smbh=.true.                 ! turns sinks into SMBH
agn=.true.                  ! enables AGN feedback
create_sinks=.true.         ! enables formation of new sink particles
mass_sink_direct_force=1.0  ! minimum mass of sink to treat it with direct solver, in M_sun

! Seeding masses
mass_sink_seed=1.0d6        ! dynamical mass of sink particle
mass_smbh_seed=1.0d6        ! accretion mass of sink particle
mass_halo_AGN=5.d10         ! minimum mass of PHEW halo in which a sink is seeded
mass_clump_AGN=1.d9         ! minimum mass of PHEW clump in which a sink is seeded

! Accretion
accretion_scheme='bondi'    ! selects Bondi accretion as accretion mode
eddington_limit=.true.      ! enables Eddington limit on accretion
acc_sink_boost=-1           ! boosts accretion according to Booth&Schaye 2009
boost_threshold_density=0.1 ! threshold density for boosting, typically the same as n_star; in H/cc
bondi_use_vrel=.false.      ! excludes relative velocity between gas and sink from the accretion calculations 

! Merging
mass_merger_vel_check=1e8   ! sum of sinks' masses for which velocities are checked upon merging to determine if the system is bound

! Feedback
T2_min=0                    ! if feedback can heat the gas to this temperature then deposit it; here deposits at every fine step
T2_AGN=0.15d12              ! thermal feedback efficiency
AGN_fbk_frac_ener=1.0       ! controls what fraction of energy goes into thermal feedback, here 100%
AGN_fbk_frac_mom=0.0        ! controls what fraction of energy goes into momentum feedback, here 0%
AGN_fbk_mode_switch_threshold=1d-2 ! switches between thermal (above) and momentum modes for this ratio of Bondi-to-Eddington
epsilon_kin=1.              ! momentum feedback efficiency
kin_mass_loading=100.       ! mass loading factor of momentum feedback
cone_opening=90.            ! opening angle of the cone in which momentum feedback is deposited (180. means full sphere)
/
```

Notes:

* if `agn=.false.` all feedback settings are disregarded
* in order to seed the sink *both* `mass_halo_AGN` and `mass_clump_AGN` have to be satisfied (as well as condition of only one sink per halo and minimum density of the halo - at least star forming)
* currently the only other accretion mode besides `bondi` is `none`; please feel free to add more
* to not boost accretion, set `acc_sink_boost=0.0`
* it is *highly advised* to use `T2_min=0.0` in order to have all feedback modes on the finest timestep
* if `chi_switch=0.0`, then the initial values of `AGN_fbk_frac_*` will be used throughout the simulation