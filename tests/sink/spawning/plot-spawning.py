# REMARK: it is normal the reference value for density is negative. This is the sum of the log of the cell densities.

import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import visu_ramses
from scipy.interpolate import griddata

fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(12, 8))

# Load RAMSES output
data = visu_ramses.load_snapshot(2)
x      = data["data"]["x"]
y      = data["data"]["y"]
z      = data["data"]["z"]
dx     = data["data"]["dx"]
rho    = data["data"]["density"]

xmin = np.amin(x-0.5*dx)
xmax = np.amax(x+0.5*dx)
ymin = np.amin(y-0.5*dx)
ymax = np.amax(y+0.5*dx)
zmin = np.amin(z-0.5*dx)
zmax = np.amax(z+0.5*dx)

nx  = 2**6
dpx = (xmax-xmin)/float(nx)
dpy = (ymax-ymin)/float(nx)
dpz = (zmax-zmin)/float(nx)
xpx = np.linspace(xmin+0.5*dpx,xmax-0.5*dpx,nx)
ypx = np.linspace(ymin+0.5*dpy,ymax-0.5*dpy,nx)
zpx = np.linspace(zmin+0.5*dpz,zmax-0.5*dpz,nx)
grid_x, grid_y, grid_z = np.meshgrid(xpx,ypx,zpx)
points = np.transpose([x,y,z])
z1 = griddata(points,rho,(grid_x,grid_y, grid_z),method='nearest')


rho_proj2 = np.sum(z1, axis=0) #proj along y-axis
rho_proj1 = np.sum(z1, axis=1) #proj along x-axis
rho_proj3 = np.sum(z1, axis=2) #proj along z-axis

im1 = ax[0].imshow(rho_proj1, origin="lower", aspect='equal', extent=[xmin, xmax, ymin, ymax])
im2 = ax[1].imshow(rho_proj2, origin="lower", aspect='equal', extent=[xmin, xmax, ymin, ymax])
im3 = ax[2].imshow(rho_proj3, origin="lower", aspect='equal', extent=[xmin, xmax, ymin, ymax])

cb = []
cb.append(plt.colorbar(im1, ax=ax[0], label='Density'))
cb.append(plt.colorbar(im2, ax=ax[1], label='Density'))
cb.append(plt.colorbar(im3, ax=ax[2], label='Density'))

#ax[0].scatter(data["sinks"]["y"], data["sinks"]["z"], color='white', marker='o')
#ax[1].scatter(data["sinks"]["x"], data["sinks"]["z"], color='white', marker='o')
#ax[2].scatter(data["sinks"]["x"], data["sinks"]["y"], color='white', marker='o')

ax[0].set_xlabel('y')
ax[0].set_ylabel('z')
ax[1].set_xlabel('x')
ax[1].set_ylabel('z')
ax[2].set_xlabel('x')
ax[2].set_ylabel('y')

for c in cb:
    c.ax.yaxis.set_label_coords(-1.1, 0.5)

fig.savefig('spawning.pdf',bbox_inches='tight')

# Check results against reference solution
for key in data["sinks"].keys():
    data["data"]["sink_"+key] = data["sinks"][key]
visu_ramses.check_solution(data["data"],'spawning', threshold=1e-30, overwrite=False)
