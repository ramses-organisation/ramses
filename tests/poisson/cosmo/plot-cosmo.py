import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import visu_ramses
from matplotlib.colors import LogNorm

fig = plt.figure(figsize=(12, 3.75))
axes = fig.subplots(nrows=1, ncols=3)

# Load RAMSES output
data = visu_ramses.load_snapshot(2,read_hydro=False)
xp = data["particle"]["position_x"]
yp = data["particle"]["position_y"]
zp = data["particle"]["position_z"]
mp = data["particle"]["mass"]

im = axes[0].hist2d(xp,yp,weights=mp,bins=128,range=[[0, 1], [0, 1]],norm=LogNorm(vmin=8e-6,vmax=8e-4),cmap='bone',edgecolor='face')
im = axes[1].hist2d(xp,zp,weights=mp,bins=128,range=[[0, 1], [0, 1]],norm=LogNorm(vmin=8e-6,vmax=8e-4),cmap='bone',edgecolor='face')
im = axes[2].hist2d(yp,zp,weights=mp,bins=128,range=[[0, 1], [0, 1]],norm=LogNorm(vmin=8e-6,vmax=8e-4),cmap='bone',edgecolor='face')
#plt.colorbar(im[3], ax=axes[2])
for ax in axes:
    ax.axis('equal')
    ax.set_xlim([0,1])
    ax.set_ylim([0,1])
axes[0].set_xlabel('x')
axes[0].set_ylabel('y')
axes[1].set_xlabel('x')
axes[1].set_ylabel('z')
axes[2].set_xlabel('y')
axes[2].set_ylabel('z')

fig.savefig("cosmo.pdf", bbox_inches="tight")

to_check = data["particle"]
to_check['time'] = data["data"]["time"]

visu_ramses.check_solution(to_check, 'cosmo')
