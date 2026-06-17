import matplotlib as mpl
mpl.use("Agg")
import os, sys, glob
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "visu"))
import visu_ramses

# Mirrors the jiang_414 notebook (osyris -> visu_ramses): 1D CR-pulse advection,
# CRegy_01 (e_c) and CRflx_01_x (F_c) profiles vs x, all snapshots overlaid.
title  = "jiang-414"
panels = [("CRegy_01", r"$e_c$"), ("CRflx_01_x", r"$F_c$")]
center = True   # shift x so the box is centred on 0

snaps = sorted(int(os.path.basename(d).split("_")[-1])
               for d in glob.glob("output_?????") if os.path.isdir(d))

# Auto-detect output mode: separate cr_NNNNN.out dump vs CR embedded in hydro file
def read_cr_mode(si):
    return os.path.exists("output_%05d/cr_file_descriptor.txt" % si)
mode = "separate cr dump" if read_cr_mode(snaps[0]) else "legacy (CR in hydro file)"

fig, axes = plt.subplots(len(panels), 1, figsize=(5.2, 1.9 * len(panels)), sharex=True)
cols = plt.cm.viridis(np.linspace(0.0, 0.82, len(snaps)))
for i, si in enumerate(snaps):
    data = visu_ramses.load_snapshot(si, read_cr=read_cr_mode(si))
    d = data["data"]
    t = data["info"]["time"]
    x = d["x"]
    srt = np.argsort(x)
    shift = 0.5 * (x.min() + x.max()) if center else 0.0
    for ax, (key, lab) in zip(axes, panels):
        ax.plot(x[srt] - shift, d[key][srt], color=cols[i], lw=1.3, label="t=%.3g" % t)
        ax.set_ylabel(lab)
axes[0].legend(fontsize=8, framealpha=0.5)
axes[-1].set_xlabel("x")
axes[0].set_title("%s  [%s]" % (title, mode))
fig.tight_layout()
plt.subplots_adjust(hspace=0.08)
fig.savefig(title + ".pdf", bbox_inches="tight")
fig.savefig(title + ".png", dpi=110, bbox_inches="tight")
print("wrote %s.pdf  (output mode: %s, %d snapshots)" % (title, mode, len(snaps)))

# regression check (mirror RT): final-snapshot sums vs committed <title>-ref.dat
visu_ramses.check_solution(data["data"], title, tolerance={"all": 1e-8},
                           overwrite=(os.environ.get("CR_REF_OVERWRITE") == "1"))
