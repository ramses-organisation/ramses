import matplotlib as mpl

mpl.use("Agg")
import numpy as np
import matplotlib.pyplot as plt
import visu_ramses

fig = plt.figure(figsize=(6, 3))
ax1 = plt.subplot(111)
ax5 = ax1.twinx()

# Load initial RAMSES output
data = visu_ramses.load_snapshot(1)
order = data["data"]["x"].argsort()
x_sim = data["data"]["x"][order]
rho_sim = data["data"]["density"][order]
# Initial density
ax1.plot(x_sim, rho_sim, color="red")

# Load final RAMSES output
data = visu_ramses.load_snapshot(2)
order = data["data"]["x"].argsort()
x_sim = data["data"]["x"][order]
amrlev = data["data"]["level"][order]
rho_sim = data["data"]["density"][order]

# final density
ax1.plot(x_sim, rho_sim, 'o', color="black", markerfacecolor='none', markersize=3)
# AMR level
ax5.plot(x_sim, amrlev, color="black", ls="dotted")

ax1.set_xlabel("Distance (cm)")
ax1.set_ylabel("Density (g/cm3)")
ax5.set_ylabel("AMR Level")

fig.subplots_adjust(wspace=0.3)
fig.savefig("advect1d.pdf", bbox_inches="tight")

# correct pressure
data["data"]["pressure"] = np.where(np.abs(data["data"]["pressure"]-1)<1e-10, 1, data["data"]["pressure"])


# Check results against reference solution
visu_ramses.check_solution(data["data"], "advect1d", tolerance={"density":3e-12}, overwrite=False)
