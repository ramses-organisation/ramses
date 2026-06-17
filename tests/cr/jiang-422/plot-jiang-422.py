import matplotlib as mpl
mpl.use("Agg")
import os, sys, glob
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "visu"))
import visu_ramses

# ---- CONFIG (per test) ----
title="jiang-422"; center=False; analytic=None
panels=[("density", r"$\rho$", False), ("velocity_x", r"$v$", False), ("pressure", r"$P_{\rm gas}$", False), ("P_CR", r"$P_{\rm CR}$", False), ("CRflx_01_x", r"$F_{\rm CR}$", False)]
# ---------------------------

snaps = sorted(int(os.path.basename(d).split("_")[-1])
               for d in glob.glob("output_?????") if os.path.isdir(d))

def read_cr_mode(si):
    return os.path.exists("output_%05d/cr_file_descriptor.txt" % si)
mode = "separate cr dump" if read_cr_mode(snaps[0]) else "legacy (CR in hydro file)"

def q(d, key):
    if key == "P_CR":       return d["CRegy_01"] / 3.0    # CR pressure = (gamma_cr-1) e_c, gamma_cr=4/3
    if key == "position_x": return d["x"]
    return d[key]                                          # density, velocity_x, pressure, CRegy_01, CRflx_01_x, level

fig, axes = plt.subplots(len(panels), 1, figsize=(5.2, 1.9 * len(panels)), sharex=True)
if len(panels) == 1:
    axes = [axes]
cols = plt.cm.viridis(np.linspace(0.0, 0.82, len(snaps)))
for i, si in enumerate(snaps):
    data = visu_ramses.load_snapshot(si, read_cr=read_cr_mode(si))
    d = data["data"]
    t = float(data["info"]["time"])
    x = q(d, "position_x"); srt = np.argsort(x)
    shift = 0.5 * (x.min() + x.max()) if center else 0.0
    for ax, (key, lab, ylog) in zip(axes, panels):
        ax.plot(x[srt] - shift, q(d, key)[srt], color=cols[i], lw=1.3, label="t=%.3g" % t)
        ax.set_ylabel(lab)
        if ylog:
            ax.set_yscale("log")
    if analytic == "triangular" and t > 0:   # 411t: panels[0]=CRegy_01, panels[1]=CRflx_01_x
        xa = np.linspace(x.min() - shift, x.max() - shift, 600)
        Ea = 2.0 + 4.0/3.0*t - np.abs(xa)
        xm = np.sqrt(16.0/9.0*t*(3.0+t)); Ea[np.abs(xa) < xm] = 2.0 + 4.0/3.0*t - xm
        axes[0].plot(xa, Ea, "k--", lw=0.9, label="analytic")
        Fa = np.zeros_like(xa); Fa[xa >= xm] = 4.0/3.0*Ea[xa >= xm]; Fa[xa <= -xm] = -4.0/3.0*Ea[xa <= -xm]
        axes[1].plot(xa, Fa, "k--", lw=0.9)
axes[0].legend(fontsize=8, framealpha=0.5); axes[-1].set_xlabel("x")
axes[0].set_title("%s  [%s]" % (title, mode))
fig.tight_layout(); plt.subplots_adjust(hspace=0.08)
fig.savefig(title + ".pdf", bbox_inches="tight")
fig.savefig(title + ".png", dpi=110, bbox_inches="tight")
print("wrote %s.pdf  (mode: %s, %d snapshots)" % (title, mode, len(snaps)))

# regression check (mirror RT): final-snapshot sums vs committed <title>-ref.dat
visu_ramses.check_solution(data["data"], title, tolerance={"all": 1e-8},
                           overwrite=(os.environ.get("CR_REF_OVERWRITE") == "1"))
