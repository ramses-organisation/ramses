

# Turbulence driving #

The block named `&TURB_PARAMS` contains the parameters related to the turbulence driving. Originally implemented by Andrew Mcleod.

## General principle[^ref]

The model used for turbulent driving in Ramses
is a generalization of the Ornstein-Uhlenbeck process.
The force is computed in Fourier space and then applied to the gas. The evolution
of the Fourier modes $\vec{\hat{f}}$ of the force was obtained via the resolution
of the following stochastic differential equation:

$$
    \mathrm{d}\vec{\hat{f}}(\vec{k}, t) = - \vec{\hat{f}}(\vec{k}, t)\dfrac{\mathrm{d}t}{T}
    + F_0(\vec{k})\vec{P_\chi}\left(\vec{k}\right) \mathrm{d}\vec{W}_t.
$$

In this equation, $\mathrm{d}t$ is the timestep for integration and $T$ is
the autocorrelation timescale. The perturbation $\mathrm{d}\vec{W_t}$ is a small vector
randomly chosen following the Wiener process. The power spectrum $F_0$ is a way to select the relevant mode.

:::{admonition} Example
:class: tip
A parabolic power spectrum between $k=1$ and $k=3$:

$$
    F_0(\vec{k}) =
    \begin{cases}
    1 - \left(\dfrac{\vec{k}}{2\pi} - 2\right)^2\text{ if }
    1 < \dfrac{\vert k \vert}{2\pi} < 3 \\
    0 \text{ if not.}
    \end{cases}
$$
:::

The projection operator $\vec{P_\chi}$ is a weighted sum of the components of
the Helmholtz decomposition of compressive versus solenoidal modes:

$$
    \vec{P_\chi}(\vec{k}) =  (1 - \chi) \vec{P}^{\perp}(\vec{k}) +
    \chi \vec{P}^{\parallel}(\vec{k}) \;,
$$

with $\vec{P}^{\perp}$ and $\vec{P}^{\parallel}$ the projection operators
respectively perpendicular and parallel to $\vec{k}$
and $\chi$ the compressive
driving fraction. This compressive driving fraction
applies only to the driving and is different from the compressive ratio measured
in the velocity field.
The forcing field $\vec{f}(\vec{x}, t)$ is then computed from the Fourier
transform:

$$
\vec{f}(\vec{x}, t) = g(\chi) f_{\mathrm{rms}}  \int\vec{\hat{f}}(\vec{k}, t)
e^{i\vec{k}\cdot x} d^3\vec{k}\;.
$$

The parameter $f_{\mathrm{rms}}$ is directly linked to the power injected by
the turbulent force into the simulation. The $g(\chi)$ factor is an empirical
correction so that the resulting time-averaged root-mean-square of the power of
the Fourier modes is equal to $f_{\mathrm{rms}}$, independently of the compressive
fraction $\chi$.



## Overview of parameters ##

| Variable name         | Fortran type | Default value | Notation | Description |
|:----------------------|:------------ |:------------- |:---------| :------------------------ |
| `turb`                | `boolean`    | `.false.`     |          | Turn on or off driving
| `turb_seed`           | `integer`    | `-1`          |          | Random number generator seed. -1 = random
| `turb_type`           | `integer`    | `1`           |          | How the driving changes over time. 1=driven evolving, 3=decaying
| `instant_turb`        | `boolean`    | `.true.`      |          | Generate initial turbulence before start
| `comp_frac`           | `float`      | `0.3333`      |  $\chi$  | The weight of compressive over solenoidal modes
| `turb_T`              | `float`      | `1`           |  $T$     | Turbulent velocity auto-correlation time in code units.
| `turb_Ndt`            | `integer`    | `100`         |  $T/dt$  | Number of timesteps per auto-correlation time |
| `turb_rms`            | `float`      | `1`           |  $f_\mathrm{rms}$ |Root-mean-square  turbulent  forcing  in  code  units. |
| `turb_min_rho`        | `float`      | `1d-50`       |          | Minimum density for turbulence. Not forcing is added onto cellswith a density less than this value.
| `forcing_power_spectrum`  | `string`     | `parabolic` | $F_0$   | Power spectrum type of the forcing, which describes the relative strength of individual modes. Options are: power_law, parabolic, konstandin

[^ref]: adapted from Brucy et al. 2024.
