import matplotlib as mpl
mpl.use("Agg")
import os, sys, glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "visu"))
import visu_ramses

# ---- CONFIG (per test) ----
title="jiang-412"; cmap, scale, vmin, vmax = "RdGy_r", "log", 1e-6, 1.0
# ---------------------------

snaps = sorted(int(os.path.basename(d).split("_")[-1])
               for d in glob.glob("output_?????") if os.path.isdir(d))
use = snaps[-3:] if len(snaps) > 3 else snaps

def read_cr_mode(si):
    return os.path.exists("output_%05d/cr_file_descriptor.txt" % si)
mode = "separate cr dump" if read_cr_mode(snaps[0]) else "legacy (CR in hydro file)"
nrm = LogNorm(vmin=vmin, vmax=vmax) if scale == "log" else Normalize(vmin=vmin, vmax=vmax)

fig, axes = plt.subplots(1, len(use), figsize=(4.2 * len(use), 3.6), squeeze=False)
axes = axes[0]
for ax, si in zip(axes, use):
    data = visu_ramses.load_snapshot(si, read_cr=read_cr_mode(si))
    d = data["data"]
    t = float(data["info"]["time"])
    x, y, cr = d["x"], d["y"], d["CRegy_01"]
    xu, yu = np.unique(x), np.unique(y)                       # uniform grid (levelmin==levelmax)
    img = np.full((len(yu), len(xu)), np.nan)
    img[np.searchsorted(yu, y), np.searchsorted(xu, x)] = cr
    im = ax.imshow(img, origin="lower", cmap=cmap, norm=nrm,
                   extent=[xu.min(), xu.max(), yu.min(), yu.max()], aspect="equal")
    ax.set_title("t=%.3g" % t)
    fig.colorbar(im, ax=ax, shrink=0.85)
fig.suptitle("%s  [%s]" % (title, mode))
fig.tight_layout()
fig.savefig(title + ".pdf", bbox_inches="tight")
fig.savefig(title + ".png", dpi=110, bbox_inches="tight")
print("wrote %s.pdf  (mode: %s, snapshots: %s)" % (title, mode, use))

# regression check (mirror RT): final-snapshot sums vs committed <title>-ref.dat
visu_ramses.check_solution(data["data"], title, tolerance={"all": 3e-6},
                           overwrite=(os.environ.get("CR_REF_OVERWRITE") == "1"))
