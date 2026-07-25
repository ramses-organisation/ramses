# Cosmic-ray (CR) test family

Tests for the two-moment CR transport module (`cr/`). The `jiang-4xx` tests
reproduce the benchmark suite of Jiang & Oh (2018) — the test number is the
paper section number. The `tp-*` tests are two-pressure (gas + CR) shock tubes
following the CR-modified Riemann problems of Thomas, Pfrommer & Pakmor (2021).

See [`doc/wiki/CR.md`](../../doc/wiki/CR.md) for the scheme and
[`doc/wiki/CR_params.md`](../../doc/wiki/CR_params.md) for the `&CR_PARAMS`
reference. The module is currently limited to a single CR group.

## Building and running

```bash
bash tests/build.cr.sh 3          # NDIM=3; SOLVER=mhd CRPHYS=1 NCR_GROUPS=1 MPI=1 DEBUG=1
cd tests && ./run_test_suite.sh -t cr        # add -p 2 to run under MPI
```

CR transport requires `SOLVER=mhd`: the closure and the streaming/scattering
source terms need the magnetic field.

## Tests

| Test | Dim | What it tests | Reference |
|---|---|---|---|
| `jiang-411` | 1D | Streaming of a Gaussian CR pulse in static gas (imposed `cr_v_alfven=1`) | J&O 4.1.1 |
| `jiang-411t` | 1D | Streaming of a triangular profile against the exact analytic solution (time-dependent analytic boundaries) | J&O 4.1.1 |
| `jiang-412` | 2D | Anisotropic streaming along a 45°-inclined B field, with streaming heating | J&O 4.1.2 |
| `jiang-413` | 1D | Streaming across a density bump: CR bottleneck where the Alfvén speed drops | J&O 4.1.3 |
| `jiang-414` | 1D | Pure two-moment diffusion of a Gaussian pulse (`Dcr=0.1`) | J&O 4.1.4 |
| `jiang-414-refine` | 1D | Same, across AMR coarse–fine boundaries with `err_grad_cr` refinement | J&O 4.1.4 |
| `jiang-414-vel` | 1D | Combined advection + diffusion in a frozen background flow, uniform grid | J&O 4.1.4 |
| `jiang-414-vel-refine` | 1D | Combined advection + diffusion in a frozen background flow, with AMR | J&O 4.1.4 |
| `jiang-415` | 2D | Anisotropic diffusion along curved field lines (magnetic loop; hardest transport test) | J&O 4.1.5 |
| `jiang-415-donut` | 2D | Same loop, evolved until the CR arc wraps into a closed ring | J&O 4.1.5 |
| `jiang-421` | 1D | Coupled CR–gas evolution: sinusoidal CR wave damping while heating the gas | J&O 4.2.1 |
| `jiang-422` | 1D | Gas–CR momentum/energy exchange in colliding flows forming a shock | J&O 4.2.2 |
| `jiang-423` | 2D | Sedov-type blast sweeping a CR background, with AMR + adaptive reduced light speed | J&O 4.2.3 |
| `jiang-424` | 1D | CR-driven outflow through a dense slab, with AMR + adaptive reduced light speed | J&O 4.2.4 |
| `tp-nostream` | 1D | Two-pressure (gas + CR) Riemann problem, pure diffusion | TPP 2021 |
| `tp-stream-va075` | 1D | Two-pressure shock tube with streaming + heating at `v_A=0.75` | TPP 2021 |
| `tp-stream-va15` | 1D | Two-pressure shock tube with streaming + heating at `v_A=1.5` | TPP 2021 |

AMR CR tests must set `cr_nsubcycle>1` (jiang-414-refine, jiang-414-vel-refine,
jiang-423 and jiang-424 use 10).

## Regression checks

Each test's `plot-<name>.py` produces a comparison figure and then calls
`visu_ramses.check_solution` against the committed `<name>-ref.dat`. The whole
family uses `tolerance={"all": 1e-12}` on the per-field cell sums. That is
looser than the framework default of `3e-13` because the sums are not
bit-reproducible across MPI domain decompositions, and tighter than the `3e-6`
used by the `rt/` family.

## References

- Jiang, Y.-F., & Oh, S. P. 2018, ApJ, 854, 5 — [arXiv:1712.07117](https://arxiv.org/abs/1712.07117)
- Thomas, T., Pfrommer, C., & Pakmor, R. 2021, MNRAS, 503, 2242 — [arXiv:2010.11960](https://arxiv.org/abs/2010.11960)
