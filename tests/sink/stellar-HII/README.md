* Test name: `stellar-HII`
* Dimension: `3`
* Solver: `hydro` + `RT` (`NGROUPS=3`, `NIONS=3`)
* Comparison: gas, RT and movie variables, plus all sink and stellar object columns of `output_00002`
* Purpose: Testing sink merging, stellar object spawning, and ionising radiation from sinks
* Keywords: sink dynamics, sink merging, stellar objects, RT, HII region, PIC

## Setup

250 pc box, uniform 10 H/cc at T = 10 K, `levelmin=6`, `levelmax=8`,
`nlevelmax_sink=7`. Six sinks are read from `ic_sink`; `create_sinks=.false.`,
so no sink ever forms from the gas and there is no formation threshold to make
the run non-reproducible.

| id | mass (Msun) | position | velocity |
|----|-------------|----------|----------|
| 6  | 4079        | (125,125,125) — box centre | at rest |
| 1  | 4079        | (200,200,125) | at rest |
| 4  | 40.79       | (50,200,125)  | at rest |
| 5  | 40.79       | (200,50,125)  | at rest |
| 2  | 20.395      | (55,50,125)   | `vx = -119` |
| 3  | 20.395      | (45,50,125)   | `vx = +119` |

Sinks 2 and 3 are aimed at each other and **merge** during the run, so
`output_00002` contains five sinks (ids 1, 2, 4, 5, 6). The merged sink keeps
the lower id, which is why the namelist notes `sink_id != isink`.

`stellar=.true.` with `stellar_msink_th=300`, so the two massive sinks spawn
**two stellar objects, on sinks 1 and 6**. `mstellarini=50,50,50,50,50,50`
fixes the first six stellar masses at 50 Msun (2039.2 in code units), which
means `sample_powerlaw` still draws from the IMF but the result is discarded.
That matters: the draw uses a per-rank RNG substream and so depends on the
number of MPI ranks. **Do not remove `mstellarini` without addressing that.**

With `rt_sink=.true.`, those two stellar objects drive the HII regions.

## The exact symmetry — the most useful property of this test

Every sink starts at `z = 125` with `vz = 0`, the gas is uniform, and the only
non-zero initial velocities are in-plane. The whole configuration is therefore
**exactly reflection-symmetric about the plane z = 125**.

Since `l = sum m (r x v)`, and under `z -> -z` both `r_z` and `v_z` flip sign
while the in-plane components do not, the `x` and `y` components of `r x v` are
odd and must cancel. So for every sink, at all times:

* `vz = 0`
* `lx = ly = 0`
* `lz` is the only component the physics permits

**Any non-zero `vz`, `lx` or `ly` is purely numerical.** That makes this test a
null test for anything that can break spatial symmetry, and it is how the NGP
deposition bug in `sink_RT_feedback` was found — see the header comment of
`pm/sink_rt_feedback.f90`.

## What the solution should look like

Exact values are in `stellar-HII-ref.dat`. Orders of magnitude at 2 MPI ranks,
from `output_00002/sink_00002.csv`:

* `vz` ~ 2e-8 for **all** sinks. This is the accumulated round-off floor, not
  physics. If it reaches 1e-3 or above, something is breaking the symmetry.
* `lz` is the physical component. Sink 2 dominates at ~0.066 — it ploughs
  through the gas at `vx = -119` before merging. Sinks 4 and 5 sit at ~1e-5
  with opposite signs, sinks 1 and 6 at ~1e-6 to 1e-5.
* `lx`, `ly` ~ 1e-10 for sinks 2, 4 and 5.
* `lx = -ly` at ~2e-4 (sink 1) and ~4e-6 (sink 6). This is a **residual,
  unexplained symmetry violation**, five orders of magnitude smaller than the
  NGP bug that was fixed, and still present with `rt_sink=.false.`. Not yet
  investigated.
* Two stellar objects, both 50 Msun, both born at the same time.

## Triage recipe when this test fails

1. Diff `output_00002/sink_00002.csv` between the two configurations first.
   `vz`, `lx` and `ly` are symmetry-forbidden, so any change there is a
   numerical artefact and not a physics disagreement.
2. Sinks 2, 4 and 5 carry no stellar object and receive no radiation injection.
   They are a **control group**: if they move too, the problem is not in the
   sink RT path.
3. `rt_sink=.false.` removes the sink radiation entirely and isolates it from
   the rest of the sink machinery.

## Known issues

**Does not pass at 4 MPI ranks** (the reference is generated at 2); the
namelist has carried a `TODO: make work on 4 cpus` since the test was added.
Every failing variable is a `sink_*` one, and all of them are either
symmetry-forbidden or noise-dominated: `sink_lx/ly/lz`, `sink_vz`,
`sink_vz_gas`. The physical quantities are bit-identical between 2 and 4 ranks
(`msink`, `dmfsink`, `acc_rate`, `del_mass`, `x`, `y`, `z`, `tform`,
`rho_gas`, `cs**2`, `etherm`), and `sink_vx/vy` and `sink_vx_gas/vy_gas` agree
to ~1e-10 and ~1e-7. So the 4-rank failure is the default 3e-13 tolerance being
applied to quantities that are legitimately zero, not a physics disagreement.

The inflated tolerances on `scalar_*`, `velocity_*` and `photon_flux_*` are
**inherited from RT**, not caused by the sinks: the whole RT family in this
suite runs at ~3e-6 (`stromgren2d`, `stromgren2d-He`, `shadow2d`,
`crossing-beams`), for reasons nobody has investigated.
