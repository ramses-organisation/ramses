import matplotlib as mpl
mpl.use("Agg")
import os, glob
import numpy as np
import matplotlib.pyplot as plt
import visu_ramses

# ---- CONFIG (per test) ----
title="jiang-421"; center=True;
panels=[("CRegy_01", r"$E_{\rm CR}$", False), ("CRflx_01_x", r"$F_{\rm CR}$", False), ("density", r"$\rho$", False), ("velocity_x", r"$v$", False)]
# ---------------------------

snaps = sorted(int(os.path.basename(d).split("_")[-1])
               for d in glob.glob("output_?????") if os.path.isdir(d))

def q(d, key):
    if key == "P_CR":       return d["CRegy_01"] / 3.0    # CR pressure = (gamma_cr-1) e_c, gamma_cr=4/3
    if key == "position_x": return d["x"]
    return d[key]                                          # density, velocity_x, pressure, CRegy_01, CRflx_01_x, level

fig, axes = plt.subplots(len(panels), 1, figsize=(5.2, 1.9 * len(panels)), sharex=True)
if len(panels) == 1:
    axes = [axes]
cols = plt.cm.viridis(np.linspace(0.0, 0.82, len(snaps)))
for i, si in enumerate(snaps):
    data = visu_ramses.load_snapshot(si, read_cr=True)
    d = data["data"]
    t = float(data["info"]["time"])
    x = q(d, "position_x"); srt = np.argsort(x)
    shift = 0.5 * (x.min() + x.max()) if center else 0.0
    for ax, (key, lab, ylog) in zip(axes, panels):
        ax.plot(x[srt] - shift, q(d, key)[srt], color=cols[i], lw=1.3, label="t=%.3g" % t)
        ax.set_ylabel(lab)
        if ylog:
            ax.set_yscale("log")

axes[0].legend(fontsize=8, framealpha=0.5); axes[-1].set_xlabel("x")
axes[0].set_title(title)
fig.tight_layout(); plt.subplots_adjust(hspace=0.08)
fig.savefig(title + ".pdf", bbox_inches="tight")

# regression check (mirror RT): final-snapshot sums vs committed <title>-ref.dat
visu_ramses.check_solution(data["data"], title, tolerance={"all": 1e-8}, overwrite=False)
