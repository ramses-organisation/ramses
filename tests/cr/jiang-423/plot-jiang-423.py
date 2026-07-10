import matplotlib as mpl
mpl.use("Agg")
import os, glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.interpolate import griddata
import visu_ramses

title = "jiang-423"
# (var, cmap, scale): mirrors the jiang_423 notebook's 4-panel last-snapshot map
quantities = [("density", "viridis", "linear"),
              ("CRegy_01", "RdGy_r", "log"),
              ("pressure", "inferno", "linear"),
              ("level",    "RdGy_r", "linear")]

snaps = sorted(int(os.path.basename(d).split("_")[-1])
               for d in glob.glob("output_?????") if os.path.isdir(d))
si = snaps[-1]

data = visu_ramses.load_snapshot(si, read_cr=True)
d = data["data"]
t = float(data["info"]["time"])
x, y = d["x"], d["y"]
n = 2 ** 7                                                    # finest level (levelmax=7) → 128
xi = np.linspace(x.min(), x.max(), n)
yi = np.linspace(y.min(), y.max(), n)
XI, YI = np.meshgrid(xi, yi)
ext = [x.min(), x.max(), y.min(), y.max()]

fig, ax = plt.subplots(1, 4, figsize=(16, 3.6))
for a, (var, cmap, scale) in zip(ax, quantities):
    vals = d[var]
    img = griddata((x, y), vals, (XI, YI), method="nearest")  # nearest = block rendering, like osyris.map
    if scale == "log":
        pos = vals[vals > 0]
        nrm = LogNorm(vmin=max(pos.min(), vals.max() * 1e-6), vmax=vals.max())
        im = a.imshow(img, origin="lower", cmap=cmap, norm=nrm, extent=ext, aspect="equal")
    else:
        im = a.imshow(img, origin="lower", cmap=cmap, extent=ext, aspect="equal")
    a.set_title(var)
    fig.colorbar(im, ax=a, shrink=0.8)
fig.suptitle("%s  t=%.3g  (snapshot %d)" % (title, t, si))
fig.tight_layout()
fig.savefig(title + ".pdf", bbox_inches="tight")

visu_ramses.check_solution(data["data"], title, tolerance={"all": 1e-8}, overwrite=False)
