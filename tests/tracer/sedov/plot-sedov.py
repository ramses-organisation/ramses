import matplotlib as mpl

mpl.use("Agg")
import numpy as np
import matplotlib.pyplot as plt
import visu_ramses
from scipy.interpolate import griddata

fig = plt.figure(constrained_layout=True, figsize=(12, 6))
axes = fig.subplots(nrows=2, ncols=3)

# Load RAMSES output
data = visu_ramses.load_snapshot(2)
x = data["data"]["x"]
y = data["data"]["y"]
dx = data["data"]["dx"]
rho = data["data"]["density"]

xp = data["particle"]["position_x"]
yp = data["particle"]["position_y"]
mp = data["particle"]["mass"]

# Check that the tracer total mass agrees with the gas mass within
# 1/sqrt(Ntracer) (~1 sigma for a Poisson distribution)
Mgas = (data["data"]["density"] * data["data"]["dx"]**2).sum()
Mtracer = data["particle"]["mass"].sum()
Ntracer = data["particle"]["mass"].size

rtol = 1/np.sqrt(Ntracer)

np.testing.assert_allclose(Mgas, Mtracer, rtol=rtol)

# Normalize the tracer mass so that their sum equals the box mass
mp = mp * (rho * dx**2).sum() / mp.sum()

xmin = np.amin(x - 0.5 * dx)
xmax = np.amax(x + 0.5 * dx)
ymin = np.amin(y - 0.5 * dx)
ymax = np.amax(y + 0.5 * dx)

nx = 32
dpx = (xmax - xmin) / float(nx)
dpy = (ymax - ymin) / float(nx)
xpx = np.linspace(xmin + 0.5 * dpx, xmax - 0.5 * dpx, nx)
ypx = np.linspace(ymin + 0.5 * dpy, ymax - 0.5 * dpy, nx)
grid_x, grid_y = np.meshgrid(xpx, ypx)
points = np.transpose([x, y])
z1 = griddata(points, rho, (grid_x, grid_y), method="nearest")

xbins = np.linspace(xmin, xmax, nx + 1)
ybins = np.linspace(ymin, ymax, nx + 1)
xpoints = np.transpose([xp, yp])
z2, _1, _2 = np.histogram2d(xp, yp, bins=(xbins, ybins), weights=mp)
z2 *= z1.sum() / z2.sum()
Npart, _1, _2 = np.histogram2d(xp, yp, bins=(xbins, ybins))

im1 = axes[0, 0].imshow(
    z1, origin="lower", aspect="equal", extent=[xmin, xmax, ymin, ymax]
)
im2 = axes[0, 1].imshow(
    z2, origin="lower", aspect="equal", extent=[xmin, xmax, ymin, ymax],
)

arg = (z2 - z1) / z1 * np.sqrt(Npart)
im3 = axes[0, 2].imshow(
    arg,
    origin="lower",
    aspect="equal",
    extent=[xmin, xmax, ymin, ymax],
    cmap="RdBu",
    vmin=-3,
    vmax=3,
)

im4 = axes[1, 0].imshow(
    z1 / z1.T,
    origin="lower",
    aspect="equal",
    extent=[xmin, xmax, ymin, ymax],
    cmap="RdBu",
    vmin=0.9,
    vmax=1.1,
)

im5 = axes[1, 1].imshow(
    z2 / z2.T,
    origin="lower",
    aspect="equal",
    extent=[xmin, xmax, ymin, ymax],
    cmap="RdBu",
    vmin=0.5,
    vmax=1.5,
)

ymean_ref = z1.mean(axis=0)
ymean = z2.mean(axis=0)
yerr = (z2/np.sqrt(Npart)).mean(axis=0)
axes[1, 2].plot(xbins[1:], ymean-ymean_ref, c="C1")
axes[1, 2].fill_between(xbins[1:], ymean-ymean_ref-yerr/2, ymean-ymean_ref+yerr/2, label="Tracer (x)", color="C1", alpha=0.5)

ymean = z2.mean(axis=1)
yerr = (z2/np.sqrt(Npart)).mean(axis=1)
axes[1, 2].plot(xbins[1:], ymean-ymean_ref, c="C2")
axes[1, 2].fill_between(xbins[1:], ymean-ymean_ref-yerr/2, ymean-ymean_ref+yerr/2, label="Tracer (y)", color="C2", alpha=0.5)
axes[1, 2].axhline(0, ls='--', lw=1, c='black')
axes[1, 2].legend(loc="upper left", fontsize=8)

plt.colorbar(im1, ax=axes[0, 0], label="Density")
plt.colorbar(im2, ax=axes[0, 1], label="Tracer density")
plt.colorbar(im3, ax=axes[0, 2], label=r"$\sqrt{N_\mathrm{tracer}}(\rho_\mathrm{tracer}-\rho_\mathrm{gas})/ {\rho_\mathrm{gas}}$")
plt.colorbar(im4, ax=axes[1, 0], label="Difference symmetric")
plt.colorbar(im5, ax=axes[1, 1], label="Difference symmetric")

for ax in axes.flatten():
    ax.set_xlabel("x")
    ax.set_ylabel("y")
axes[1, 2].set_xlabel("x")
axes[1, 2].set_ylabel(r"Dens. diff. ${\rho_\mathrm{gas}-\rho_\mathrm{tracer}}$")
axes[1, 2].yaxis.tick_right()

fig.savefig("sedov.pdf", bbox_inches="tight")

# Check results against reference solution
dt = {
    k: v for k, v in data["particle"].items()
    if k in ("family", "levelp", "mass", "position_x", "position_y", "tag")
}

# Use a relative tolerance within ± 2 Poisson noise
rtol = {
    key: 2/np.sqrt(data["particle"]["identity"].size)
    for key in dt
}
visu_ramses.check_solution(dt, 'sedov', tolerance=rtol)
