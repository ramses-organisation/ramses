

# Turbulence driving #

The block named `&TURB_PARAMS` contains the parameters related to the turbulence driving. Originally implemented by Andrew Mcleod

## Overview of parameters ##

| Variable name         | Fortran type | Default value | Description               |
|:----------------------|:------------ |:------------- |:------------------------- |
| `turb`                | `boolean`    | `.false.`     |  Turn on or off driving
| `turb_seed`           | `integer`    | `-1`          |  Random number generator seed. -1 = random
| `turb_type`           | `integer`    | `1`           |  How the driving changes over time. 1=driven evolving, 3=decaying
| `instant_turb`        | `boolean`    | `.true.`      |  Generate initial turbulence before start
| `comp_frac`           | `float`      | `0.3333`      |  The weight of compressive over solenoidal modes
| `turb_T`              | `float`      | `1`           |  Turbulent velocity auto-correlation time in code units.
| `turb_Ndt`            | `integer`    | `100`         |  Number of timesteps per auto-correlation time
| `turb_rms`            | `float`      | `1`           |  Root-mean-square  turbulent  forcing  in  code  units.
| `turb_min_rho`        | `float`      | `1d-50`       |  Minimum density for turbulence. Not forcing is added onto cellswith a density less than this value.
| `forcing_power_spectrum`  | `string`     | `parabolic`     | Power spectrum type of the forcing, which describes the relative strength of individual modes. Options are: power_law, parabolic, konstandin