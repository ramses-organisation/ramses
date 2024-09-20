

# Monte Carlo tracer particles

The block named `&TRACER_PARAMS` contains the parameters related to the Monte Carlo tracer particles, see
[Cadiou et al 2018](https://arxiv.org/abs/1810.11401)


## Overview of parameters

| Variable name         | Fortran type | Default value | Description               |
|:----------------------|:------------ |:------------- |:------------------------- |
| `MC_tracer`           | `boolean`    | `.false.`     | Activate MC tracers
| `tracer_feed`         | `string`     | `none`        | Filename to read the tracer from
| `tracer_feed_fmt`     | `string`     | `ascii`       | Format of the input (ascii,binary or inplace)
| `tracer_mass`         | `float`      | `-1.0`        | Mass of the tracers, used for outputs and seed
| `tracer_first_balance_levelmin`  | `integer`    | `-1`           | Set to >0 to add more weight on level finer than this
| `tracer_first_balance_part_per_cell`    | `integer`    | `0`           | Typical initial number of parts per cell