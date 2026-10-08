import glob

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import visu_ramses
from matplotlib.colors import LogNorm

# Columns of the ASCII output written by the clump finder and the merger tree
MTREE_COLS = ["clump", "progenitor", "prog_outputnr", "desc_mass", "desc_npart",
              "desc_x", "desc_y", "desc_z", "desc_vx", "desc_vy", "desc_vz"]
HALO_COLS = ["index", "ncell", "peak_x", "peak_y", "peak_z", "rho_max", "mass"]


def read_table(nout, prefix, columns):
    """Concatenate the per-MPI-task ASCII files of one output into a dict."""
    rows = []
    pattern = "output_%05d/%s_%05d.txt*" % (nout, prefix, nout)
    for fname in sorted(glob.glob(pattern)):
        with open(fname) as fh:
            for line in fh:
                fields = line.split()
                # Skip the header line each task writes.
                if len(fields) != len(columns) or not fields[0].lstrip("-").isdigit():
                    continue
                rows.append([float(f) for f in fields])
    table = np.array(rows) if rows else np.zeros((0, len(columns)))
    return {name: table[:, i] for i, name in enumerate(columns)}


def output_numbers():
    return sorted(int(d.split("_")[1]) for d in glob.glob("output_?????"))


def expansion_factor(nout):
    with open("output_%05d/info_%05d.txt" % (nout, nout)) as fh:
        for line in fh:
            if line.startswith("aexp"):
                return float(line.split("=")[1])
    return np.nan


def main_branch(trees, nout, clump):
    """Walk a clump back in time along its main progenitor branch."""
    branch = []
    while nout in trees and clump in trees[nout]:
        row = trees[nout][clump]
        branch.append((nout, row["desc_mass"]))
        progenitor = int(row["progenitor"])
        if progenitor == 0:
            break
        # A negative progenitor lives in an older, non-adjacent snapshot.  The
        # snapshot number has to decrease, otherwise a corrupted tree would send
        # us round in circles.
        if int(row["prog_outputnr"]) >= nout:
            break
        clump, nout = abs(progenitor), int(row["prog_outputnr"])
    return branch


nouts = output_numbers()
aexp = {n: expansion_factor(n) for n in nouts}
halos = {n: read_table(n, "halo", HALO_COLS) for n in nouts}
mtree = {n: read_table(n, "mergertree", MTREE_COLS) for n in nouts}
last = nouts[-1]

# Index the trees by descendant clump id so branches can be followed.
trees = {}
for n in nouts:
    trees[n] = {}
    for i, clump in enumerate(mtree[n]["clump"]):
        if clump > 0:
            trees[n][int(clump)] = {k: mtree[n][k][i] for k in MTREE_COLS}

# Snapshots that actually contain clumps: the first output is written at astart,
# long before anything has collapsed.
with_halos = [n for n in nouts if len(halos[n]["mass"]) > 0]

fig = plt.figure(figsize=(15, 8))
axes = fig.subplots(nrows=2, ncols=3)

# ----------------------------------------------------------------------------
# 1. Projected mass distribution with the halos found at the last output
# ----------------------------------------------------------------------------
data = visu_ramses.load_snapshot(last, read_hydro=False)
xp = data["particle"]["position_x"]
yp = data["particle"]["position_y"]
mp = data["particle"]["mass"]

ax = axes[0][0]
ax.hist2d(xp, yp, weights=mp, bins=128, range=[[0, 1], [0, 1]],
          norm=LogNorm(vmin=8e-6, vmax=8e-4), cmap="bone", edgecolor="face")
# Area proportional to the halo radius so the small ones stay visible.
ax.scatter(halos[last]["peak_x"], halos[last]["peak_y"],
           s=300 * halos[last]["mass"] ** (1.0 / 3.0), facecolor="none",
           edgecolor="orangered", linewidth=0.8)
ax.set_xlim([0, 1])
ax.set_ylim([0, 1])
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title("particles + halos at a=%.3f" % aexp[last])

# ----------------------------------------------------------------------------
# 2. Cumulative halo mass function
# ----------------------------------------------------------------------------
ax = axes[0][1]
for n in with_halos:
    m = np.sort(halos[n]["mass"])[::-1]
    ax.plot(m, np.arange(1, len(m) + 1), label="a=%.3f" % aexp[n])
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("halo mass [code units]")
ax.set_ylabel("N(>M)")
ax.set_title("cumulative halo mass function")
ax.legend(fontsize=7)

# ----------------------------------------------------------------------------
# 3. Growth of the halo population
# ----------------------------------------------------------------------------
ax = axes[0][2]
a = [aexp[n] for n in with_halos]
ax.plot(a, [len(halos[n]["mass"]) for n in with_halos], "o-", color="C0")
ax.set_xlabel("a")
ax.set_ylabel("number of halos", color="C0")
ax.tick_params(axis="y", labelcolor="C0")
twin = ax.twinx()
twin.plot(a, [halos[n]["mass"].sum() for n in with_halos], "s--", color="C3")
twin.set_ylabel("mass in halos [code units]", color="C3")
twin.tick_params(axis="y", labelcolor="C3")
ax.set_title("halo population vs time")

# ----------------------------------------------------------------------------
# 4. Mass growth along the main branch of the most massive halos
# ----------------------------------------------------------------------------
ax = axes[1][0]
massive = sorted(trees[last], key=lambda c: -trees[last][c]["desc_mass"])[:8]
for clump in massive:
    branch = main_branch(trees, last, clump)
    if len(branch) < 2:
        continue
    ax.plot([aexp[n] for n, _ in branch], [m for _, m in branch], "o-", markersize=3)
ax.set_yscale("log")
ax.set_xlabel("a")
ax.set_ylabel("clump mass [code units]")
ax.set_title("main branch of the 8 most massive clumps")

# ----------------------------------------------------------------------------
# 5. Census of the links in the trees
# ----------------------------------------------------------------------------
ax = axes[1][1]
census = {"direct": [], "new": [], "non-adjacent": [], "merged": []}
for n in with_halos:
    clump = mtree[n]["clump"].astype(int)
    prog = mtree[n]["progenitor"].astype(int)
    census["direct"].append(np.count_nonzero((clump > 0) & (prog > 0)))
    census["new"].append(np.count_nonzero((clump > 0) & (prog == 0)))
    census["non-adjacent"].append(np.count_nonzero((clump > 0) & (prog < 0)))
    census["merged"].append(np.count_nonzero(clump < 0))
bottom = np.zeros(len(with_halos))
for label, counts in census.items():
    counts = np.array(counts, dtype=float)
    ax.bar(a, counts, bottom=bottom, width=0.006, label=label)
    bottom += counts
ax.set_xlabel("a")
ax.set_ylabel("number of tree entries")
ax.set_title("merger tree links")
ax.legend(fontsize=7)

# ----------------------------------------------------------------------------
# 6. How far a clump travels between two snapshots
# ----------------------------------------------------------------------------
ax = axes[1][2]
for n in with_halos[1:]:
    steps = []
    for clump, row in trees[n].items():
        progenitor = int(row["progenitor"])
        nprev = int(row["prog_outputnr"])
        if progenitor <= 0 or nprev not in trees:
            continue
        prev = trees[nprev].get(progenitor)
        if prev is None:
            continue
        d = np.array([row["desc_" + a] - prev["desc_" + a] for a in "xyz"])
        d -= np.round(d)                      # periodic box
        steps.append(np.sqrt((d ** 2).sum()))
    if steps:
        ax.hist(steps, bins=np.logspace(-3.5, -0.5, 25), histtype="step",
                label="a=%.3f" % aexp[n])
ax.set_xscale("log")
ax.set_xlabel("progenitor to descendant distance [box units]")
ax.set_ylabel("number of links")
ax.set_title("main branch displacement per snapshot")
ax.legend(fontsize=7)

fig.tight_layout()
fig.savefig("mergertree.pdf", bbox_inches="tight")

# ----------------------------------------------------------------------------
# Structural invariants of the trees
#
# These do not depend on the exact halo catalogue, so they must be exactly zero
# whatever the number of MPI tasks or the rounding of the Poisson solve.  Any
# non-zero value makes check_solution report an infinite error, i.e. the test
# fails, which is what we want for a broken tree rather than a slightly
# different one.
# ----------------------------------------------------------------------------
# All the dark matter particles have the same mass (one per cell of the ICs).
particle_mass = np.median(mp)

alive = {n: set(trees[n]) for n in nouts}
bad_mass = bad_position = bad_direct = bad_nonadjacent = bad_reused = 0
for n in with_halos:
    t = mtree[n]
    npart = t["desc_npart"]
    keep = npart > 0
    bad_mass += np.count_nonzero(
        np.abs(t["desc_mass"][keep] - npart[keep] * particle_mass) >
        1e-6 * t["desc_mass"][keep])
    for axis in "xyz":
        p = t["desc_" + axis]
        bad_position += np.count_nonzero((p < 0.0) | (p > 1.0))

    clump = t["clump"].astype(int)
    prog = t["progenitor"].astype(int)
    prog_out = t["prog_outputnr"].astype(int)

    # A direct progenitor was an alive clump in the previous snapshot.
    direct = (clump > 0) & (prog > 0)
    for p, o in zip(prog[direct], prog_out[direct]):
        if o in alive and p not in alive[o]:
            bad_direct += 1

    # A negative progenitor is a "past merged progenitor": it was alive in
    # prog_outputnr, which has to be an older, non-adjacent snapshot.
    old = (clump > 0) & (prog < 0)
    for p, o in zip(-prog[old], prog_out[old]):
        if o >= n - 1 or (o in alive and p not in alive[o]):
            bad_nonadjacent += 1

    # A progenitor has exactly one descendant, so no progenitor may be claimed
    # twice - neither as a main progenitor (clump > 0) nor as one that merged
    # into the descendant (clump < 0).
    claimed = [(abs(p), o) for p, o in zip(prog, prog_out) if p != 0]
    bad_reused += len(claimed) - len(set(claimed))

# ----------------------------------------------------------------------------
# Reference solution
# ----------------------------------------------------------------------------
to_check = dict(data["particle"])
to_check["time"] = data["data"]["time"]

# Halo catalogue, over every snapshot that has one.
for key in ["mass", "ncell", "peak_x", "peak_y", "peak_z"]:
    to_check["halo_" + key] = np.concatenate([halos[n][key] for n in with_halos])
to_check["halo_number"] = np.array([len(halos[n]["mass"]) for n in with_halos],
                                   dtype=float)

# Merger tree entries, over every snapshot that has one.
for key in ["desc_mass", "desc_npart", "desc_x", "desc_y", "desc_z",
            "desc_vx", "desc_vy", "desc_vz"]:
    to_check["tree_" + key] = np.concatenate([mtree[n][key] for n in with_halos])
for label in ["direct", "new", "non-adjacent", "merged"]:
    to_check["tree_links_" + label.replace("-", "_")] = np.array(census[label],
                                                                dtype=float)

to_check["tree_bad_mass"] = np.array([bad_mass], dtype=float)
to_check["tree_bad_position"] = np.array([bad_position], dtype=float)
to_check["tree_bad_direct_progenitor"] = np.array([bad_direct], dtype=float)
to_check["tree_bad_nonadjacent_progenitor"] = np.array([bad_nonadjacent], dtype=float)
to_check["tree_reused_progenitor"] = np.array([bad_reused], dtype=float)

# The catalogue was measured to be bit-identical on 1, 2 and 4 MPI tasks, but it
# is a discrete quantity: the clump finder compares densities against fixed
# thresholds, so on another compiler a marginal clump can fall on the other side
# of the relevance cut and shift these sums by ~1/1000.  The catalogue and
# everything derived from it therefore gets a 1% tolerance, while the particle
# distribution - a smooth function of the initial conditions - keeps the default
# one.  See the README.
#tolerance = {key: 1.0e-2 for key in to_check
#             if key.startswith(("halo_", "tree_"))}
# The invariants above must stay exactly zero.
#for key in list(tolerance):
#    if key.startswith("tree_bad") or key == "tree_reused_progenitor":
#        del tolerance[key]

visu_ramses.check_solution(to_check, "mergertree", #tolerance=tolerance,
                           overwrite=False)
