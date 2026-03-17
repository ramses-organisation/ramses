# FLD

## Namelist parameters overview

| Variable name, syntax, default value | Fortran type  | Description               |
|:---------------------------- |:------------- |:------------------------- |
| `grey_rad_transfer=.true`           |  `logical`    | Use the grey approximation, meaning there is only 1 radiation bin covering the entire spectral range from 0 to infinity. Only allowed when `NGRP=1`. |
| `freqs_in_Hz=.true.`| `logical`    | If `.true.`, boundary parameter values are provided as frequencies in Hz. Otherwise, they are given as energies in eV. |
| `numin=1.0d5` | `real` |  The lower boundary of the first radiation group. |
| `numax=1.0d19` | `real` |  The upper boundary of the last radiation group. |
| `split_groups_log=.true.` | `logical`    | Generate logaritmically (`.true.`) or linearly (`.false.`) spaced radiation bins between `numin` and `numax`.  |
| `extra_end_group=.false.` | `logical`    | The last radiation group goes from `numax` to infinity, instead of stopping at `numax`.  |
| `read_groups=.false.` | `logical`    | If `.true.`, read radiation group boundaries from a file named 'groups.dat', instead of generating them automatically.  |
| `rosseland_params=1.0` | array of `real`| Rosseland opacity coefficient's parameters |
| `planck_params=1.0` | array of `real`| Planck opacity coefficient's parameters |
| `sublimation_kuiper=.false.` | `logical`    |  Mimicks dust sublimation with decreasing d/g ratio, see Kuiper+10 ApJ |

## Setting radiation groups

The FLD implementation can be used either with the grey approximation ([Commercon et al 2011](https://ui.adsabs.harvard.edu/abs/2011A%26A...529A..35C/abstract)) or with the multigroup implementation ([Gonzalez et al 2015](https://ui.adsabs.harvard.edu/abs/2015A%26A...578A..12G/abstract)).

In the grey approximation, we have only one radiation group which covers the entire spectral range from zero to infinity. To use this approximation, compile with `NGRP=1` and set `grey_rad_transfer=.true.` in the namelist.

With the multigroup implementation, we can choose as many radiation groups as we want (or as memory constraints allow). The lower and upper limits of each group are determined from the following namelist parameters:
* freqs_in_Hz
* numin and numax
* split_groups_log and extra_end_group
* read_groups

The group boundary parameters can be provided by the user as either frequencies in Hz (`freqs_in_Hz=.true.`) or energies in eV (`freqs_in_Hz=.false.`).

If there is only a single radiation group (`NGRP=1`), the lower and upper frequency (or energy) can be set with the parameters `numin` and `numax`.

If there is more than one group (`NGRP>1`), you can either
* generate group bounds automatically based on the namelist parameters
* provide a list of boundaries in a file named 'groups.dat' (read_groups=.true.)

In the first case, the parameters `split_groups_log` will determine if the groups are splitted logaritmically (`split_groups_log=.true.`) or linearly (`split_groups_log=.false.`). The boundaries will be generated with equal spacing between `numin` and `numax`.

In some cases, you want to specify an additional radiation bin at the end, containing all radiation beyond the group boundary `numax`. This can be done by `setting extra_end_group=.true.`. Enabling this will create NGRP-1 equally spaced radiation groups between `numin` and `numax`, and add an extra group at the end of the list with boundaries `numax` to infinity. Remark that in practice, infinity is hardcoded to be `frequency_upperlimit=1.0d35`.

For example, the set of namelist parameters
TODO
will create the following groups:
TODO

Providing the group boundaries with a file allows for more flexibility. An example of a user provided groups.dat file and corresponding namelist:
TODO


## Opacities

TODO

## Timestep considerations

TODO

## Numerical parameters for the algorithm

## Radiation feedback from sinks