* Test name: `mergertree`
* Dimension: `3`
* Solver: `hydro` (dark matter only, `hydro=.false.`)
* Purpose: Testing the merger tree, and with it the clump finder (PHEW)
  and the particle unbinding
* Keywords: cosmology, pic, poisson, amr, clumpfind, unbinding, mergertree

Setup
-----

A dark-matter-only cosmological run in a 0.84 h-1 Mpc box, from a=0.0066 to
a=0.1355 (z=150 to z=6.4), with `levelmin=6`, `levelmax=8`.

The initial conditions are the same 128^3 grafic set the `poisson/cosmo` test
uses, but we degrade them to 64^3 with
`utils/f90/degrade_grafic.f90` to make the test cheaper.


Expected results
----------------

Outputs are written at six specified expansion factors.
The clump finder reports 40, 71,
98, 120, 126 and 170 halos in the log. The `halo_*.txt*` files hold more than
that - 79, 111, 155, 195, 221 and 260, which is what `halo_number` sums to -
because they list every main halo above `mass_threshold` without applying the
relevance cut that the `clump_*.txt*` files use.

The trees contain 663 entries in total: 395 clumps with a direct progenitor, 224
newly formed, 3 with a progenitor from an older non-adjacent snapshot, and 41
progenitors that merged into a descendant. Six snapshots is the minimum that
produces the non-adjacent case at all; with three it never occurs.

The structural invariants in `plot-mergertree.py` (`tree_bad_*`,
`tree_reused_progenitor`) must be exactly zero, whatever the task count. They
check that

* `desc_mass` equals `desc_npart` times the particle mass,
* every descendant position is inside the box,
* every direct progenitor was an alive clump in the previous snapshot,
* every negative progenitor was an alive clump in an older, non-adjacent one,
* no progenitor is claimed by more than one descendant.

The reference solution was produced on 2 MPI tasks.

Known issues
---------------

`make_mock_galaxies` is off. The stellar masses it produces cannot be checked at
all: the SHAM relation applies a random scatter and `calc_stellar_mass_params`
in `pm/clfind_commons.f90` seeds that generator from `system_clock()`, so they
differ in every run.

`stellar_mass` in `pm/clfind_commons.f90` overflows here, because `exp(10**(-logMM1))` is unrepresentable
for any halo well below M1 - all of them at these redshifts. That is harmless in
IEEE arithmetic (the term it feeds goes to zero) but fatal under the
`-ffpe-trap=overflow` that `DEBUG=1` compiles in.

Even with checks in place for the previous issue, a
`DEBUG=1` build still cannot run this test, for an unrelated reason: it dies
with SIGFPE in `bordercheck` at `pm/unbinding.f90:1464`, comparing a position
against `peak_pos`. That happens with `make_mock_galaxies` on or off, so it is
in the unbinding rather than the merger tree, and is most likely a NaN left by
the `-finit-real=nan` that `DEBUG=1` also sets.


Note on `-fwhole-program`
-------------------------

This test does not pass with a build that uses `-fwhole-program`. That flag lets
gcc assume no symbol is referenced outside the compiled sources, which is false
for an MPI program: `MPI_IN_PLACE` stops comparing equal to the sentinel the MPI
library expects, and every reduction using it silently drops the root task's own
contribution. `write_progenitor_data` then records only one task's progenitor
counts in `progenitorcount.dat`, and the next snapshot reads past the end of its
buffers and aborts in the heap allocator.

`bin/Makefile` therefore uses `-flto` without `-fwhole-program`. Plain `-flto`
was verified not to have the problem.
