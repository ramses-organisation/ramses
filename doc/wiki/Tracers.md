

# Tracer particles

The block named `&TRACER_PARAMS` contains the parameters related to tracer particles.
Tracers particles are used to follow the gas in a Lagrangian way.
There are three different implementation of tracers in Ramses.

- The classical tracers ([Pichon et al 2011](https://ui.adsabs.harvard.edu/abs/2011MNRAS.418.2493P), [Dubois et al 2012](https://ui.adsabs.harvard.edu/abs/2012MNRAS.423.3616D)).
The velocity of the tracers is computed from a CIC-interpolation of the velocity of the cells.
The trajectory is continuous but this kind of tracers is known to eventually decorrelate with the gas flow ([Price & Federrath 2010](https://ui.adsabs.harvard.edu/abs/2010MNRAS.406.1659P), [Genel 2013](https://ui.adsabs.harvard.edu/abs/2013MNRAS.435.1426G)).
They are activated with the `tracer` parameter in the `RUN_PARAMS` block.


- The Monte Carlo tracer particles ([Cadiou et al 2018](https://arxiv.org/abs/1810.11401)).
Each tracer is attached to a cell. Tracers jumps randomly from a cell to another with a probably matching the fraction of the mass flux between the two cell. If there is enough tracer particles per cell, the flow of tracers statistically matches the flow of the gas. The trajectory is discontinuous.
They are activated with the `tracer` parameter in the `&RUN_PARAMS` block, and the `MC_tracers` in the `&TRACERS_PARAMS` bloc.


- The Itô tracer particles ([Moseley et al 2026](https://ui.adsabs.harvard.edu/abs/2026MNRAS.551g1536M)) can be seen as a mix the two methods above. The velocity of the tracers (the drift) is computed from a CIC-interpolation of the fluxes of the cells. Additionally to this drift, each tracer undergo a random kick at each step, so that the diffusion of the gas is reproduced. As a result, the tracers have a continuous trajectory, and the flow of the gas is reproduced almost as well as with the discrete Monte-Carlo tracers. The scheme implemented is RAMSES is the "CIC Itô-2 tracers".
They are activated with the `tracer` parameter in the `&RUN_PARAMS` block, and the `MC_tracers` and `ito_tracers` in the `&TRACERS_PARAMS` bloc.



| Tracers       |  Type         | Trajectory    | Reproduction of the flow  | Activate with  |
|:--------------|:--------------|:------------- |:--------------------------|----------------|
| Classical     | deterministic | Continuous    | Poor                      |`tracer`
| Monte-Carlo   | stochastic    | Discrete      | Very good                 |`tracer` and `MC_tracers`
| Itô           | stochastic    | Continuous    | Good                      |`tracer` and `MC_tracer` and `Itô tracer`




## Overview of parameters of `&TRACERS_PARAMS`

| Variable name         | Fortran type | Default value | Description               |
|:----------------------|:------------ |:------------- |:------------------------- |
| `MC_tracer`           | `boolean`    | `.false.`     | Activate MC tracers
| `ito_tracer`           | `boolean`    | `.false.`    | Activate Itô tracers (requires `MC_tracer`)
| `tracer_feed`         | `string`     | `none`        | Filename to read the tracer from
| `tracer_feed_fmt`     | `string`     | `ascii`       | Format of the input (ascii - for all kinds-, binary or inplace - for discrete Monte-Carlo only)
| `tracer_mass`         | `float`      | `-1.0`        | Mass of the tracers, used for outputs, seed and ascii input
| `tracer_first_balance_levelmin`  | `integer`    | `-1`           | Set to >0 to add more weight on level finer than this
| `tracer_first_balance_part_per_cell`    | `integer`    | `0`           | Typical initial number of parts per cell
