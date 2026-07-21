

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
| `turb_type`           | `integer`    | `1`           |          | How the driving changes over time. 1=driven evolving, 2=driven fixed, 3=decaying. See the "Turbulence types" section below
| `instant_turb`        | `boolean`    | `.true.`      |          | Evolve the field to saturation before the run starts. Recommended for `turb_type=1`. See the "Turbulence types" section below
| `comp_frac`           | `float`      | `0.3333`      |  $\chi$  | The weight of compressive over solenoidal modes
| `turb_T`              | `float`      | `1`           |  $T$     | Turbulent velocity auto-correlation time in code units.
| `turb_Ndt`            | `integer`    | `100`         |  $T/dt$  | Number of timesteps per auto-correlation time |
| `turb_rms`            | `float`      | `1`           |  $f_\mathrm{rms}$ |Root-mean-square  turbulent  forcing  in  code  units. |
| `turb_exact_rms`      | `boolean`    | `.false.`     |          | Always use a forcing rms of exactly $\sqrt{N_{dim}}f_\mathrm{rms}$ instead of letting it follow the random draw. Only applicable for `turb_type=1`. See the "Turbulence types" section below
| `turb_min_rho`        | `float`      | `1d-50`       |          | Minimum density for turbulence. Not forcing is added onto cellswith a density less than this value.
| `forcing_power_spectrum`  | `string`     | `parabolic` | $F_0$   | Power spectrum type of the forcing, which describes the relative strength of individual modes. Options are: power_law, parabolic, konstandin


## Turbulence types ##

The `turb_type` parameter selects one of three fairly different behaviours. The
same turbulent field is generated in all three cases; what changes is how it is
used afterwards.

### `turb_type=1` — evolving driving ###

The default. The field evolves through the Ornstein-Uhlenbeck process described
above, being advanced every $dt = T/N_{dt}$ and linearly interpolated in between,
and is applied as an acceleration at every timestep. The timestep is capped at
$T/N_{dt}$ so that no step skips over a field update.

The rms of the field fluctuates by a few percent about
$\sqrt{N_{dim}}\,f_{\mathrm{rms}}$ from one update to the next. This is intended: $f_{\mathrm{rms}}$ sets the *time-averaged* amplitude, and the
scatter averages out over a run. This behavior can be altered by setting `turb_exact_rms=.true.`. This will force the rms to always be exactly $\sqrt{N_{dim}}\,f_{\mathrm{rms}}$ at every step, so the driving
pattern keeps decorrelating while the injected amplitude stays constant.

An evolving acceleration field grown from rest starts at only $\sqrt{1-e^{-2dt/T}}\approx14\%$
of its saturated amplitude and needs about one autocorrelation time $T$ to build
up. The `instant_turb` option (on by default) evolves the field for several $T$
before the run starts, so that it begins already saturated. Leaving it on is
strongly recommended for `turb_type=1`, especially for runs shorter than or
comparable to $T$, which would otherwise spend most of their length in the
spin-up transient. It has no useful effect for the other two types.

### `turb_type=2` — fixed driving ###

One acceleration field is generated at startup and then held, unchanged, for the whole run.
It is applied as an acceleration at every timestep, exactly as for type 1.

Because the field is frozen, the scatter that averages out for type 1 would
instead become a systematic amplitude bias for the whole simulation. To avoid this,
`turb_exact_rms=.true.` is always set to true for this type, so the field is
renormalised at startup so that its rms is exactly $\sqrt{N_{dim}}\,f_{\mathrm{rms}}$.

Taken together, `turb_type` and `turb_exact_rms` are two independent axes — how
the driving *pattern* behaves, and whether its *amplitude* is pinned:

| | `turb_exact_rms=.false.` | `turb_exact_rms=.true.` |
|:---|:---|:---|
| `turb_type=1` | Pattern and amplitude both fluctuate (the raw Ornstein-Uhlenbeck process) | Pattern decorrelates at constant injected amplitude |
| `turb_type=2` | Frozen pattern, amplitude left to the draw (illegal) | Frozen pattern at exactly the requested amplitude |

### `turb_type=3` — decaying turbulence ###

Despite living in the driving module, this is **not** a driven run. The
turbulent field is applied exactly once, by `init_flow_fine`, as the initial
velocity field, and is then zeroed for the rest of the run. Every forcing and
timestep-limiting path is switched off for this type, so the timestep is set by
the ordinary CFL condition rather than being capped at $T/N_{dt}$.

This has an important consequence for `turb_rms`. For types 1 and 2 the field
is an **acceleration**, and the velocity a run reaches is a saturation value set
by driving against dissipation, typically far below $f_{\mathrm{rms}}$. For type
3 the field is applied with an effective timestep of exactly one code time unit,
which turns it into a **velocity**: the box starts at
$|v|_{\mathrm{rms}} = \sqrt{N_{dim}}\,f_{\mathrm{rms}}$ immediately. As for type
2, the field is renormalised at startup, so this initial velocity dispersion,
and hence the initial Mach number, is exact rather than left to the draw.

So `turb_rms` carries over badly between the two. Reusing a driven value for a
decaying run typically gives an initial Mach number one to two orders of
magnitude too high, and a correspondingly tiny CFL timestep. For a decaying run,
choose `turb_rms` from the initial velocity dispersion you want, divided by
$\sqrt{N_{dim}}$.

Note this differs from the more usual way of setting up decaying turbulence, in
which a driven run is evolved to a steady state and the driving is then switched
off. Here the initial field is a single draw of the process rather than a
saturated turbulent state, so it has not developed the density structure or the
velocity correlations of fully developed turbulence. For that reason
`turb_type=3` is best thought of as a quick way to seed a velocity field, not as
a substitute for genuinely decaying turbulence; for the latter, use the recipe
below.

### Decaying turbulence from a developed field ###

To study freely decaying turbulence it is usually better to develop the
turbulence self-consistently first and only then let it decay, rather than to
start from the single draw of `turb_type=3`:

1. Run with `turb_type=1` (keep `instant_turb` on) until the turbulence is fully
   developed — as a rule of thumb a few autocorrelation times $T$, long enough
   for the density structure and velocity correlations to build up. Write an
   output at the point you want the decay to begin.
2. Restart from that output with turbulence switched off, i.e. `turb=.false.`
   in the namelist.

On restart with `turb=.false.` no forcing is applied and the initial velocity
kick is not re-applied (`init_flow_fine` only applies it on a fresh start,
`nrestart==0`), so the developed field simply decays under its own dissipation
from the state in the restart file. This gives a physically meaningful decaying
run, which the `turb_type=3` initial-velocity field does not.


[^ref]: adapted from Brucy et al. 2024.
