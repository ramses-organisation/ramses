import sys
from pathlib import Path
import matplotlib as mpl
import numpy as np

mpl.use("Agg")
import matplotlib.pyplot as plt

sys.path.append(str(Path(__file__).resolve().parent.parent.parent / ".osyris/src"))
import osyris

fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(8, 12))
ax2 = ax[0].twinx()

data = osyris.RamsesDataset(2).load()["mesh"]

osyris.plot(data["position_x"], data["density"], ax=ax[0], marker="o", ls="None")
osyris.plot(data["position_x"], data["level"], color="black", ls="dotted", ax=ax2)
osyris.plot(data["position_x"], data["velocity_x"], ax=ax[1], marker="o", ls="None")
osyris.plot(data["position_x"], data["pressure"], ax=ax[2], marker="o", ls="None")

# Read analytical solution
data_ana = np.loadtxt("sod-tube-ana.dat")
x_ana = data_ana[:, 1]
rho_ana = data_ana[:, 3]
u_ana = data_ana[:, 2]
p_ana = data_ana[:, 4]
ax[0].plot(x_ana, rho_ana, color="C1")
ax[1].plot(x_ana, u_ana, color="C1")
ax[2].plot(x_ana, p_ana, color="C1")

fig.subplots_adjust(wspace=0.3)
fig.savefig("sod-tube.pdf", bbox_inches="tight")

# Check results against reference solution
# visu_ramses.check_solution(data["data"], "sod-tube")
