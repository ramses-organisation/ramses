import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import visu_ramses
from scipy.interpolate import griddata
from mpl_toolkits.axes_grid1 import make_axes_locatable

fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(8, 5))

# Load RAMSES output
data = visu_ramses.load_snapshot(2)
x      = data["data"]["x"]
y      = data["data"]["y"]
dx     = data["data"]["dx"]
rho    = data["data"]["density"]
p      = data["data"]["pressure"]
scalar = data["data"]["scalar_00"]

xmin = np.amin(x-0.5*dx)
xmax = np.amax(x+0.5*dx)
ymin = np.amin(y-0.5*dx)
ymax = np.amax(y+0.5*dx)

nx  = 2**9
dpx = (xmax-xmin)/float(nx)
dpy = (ymax-ymin)/float(nx)
xpx = np.linspace(xmin+0.5*dpx,xmax-0.5*dpx,nx)
ypx = np.linspace(ymin+0.5*dpy,ymax-0.5*dpy,nx)
grid_x, grid_y = np.meshgrid(xpx,ypx)
points = np.transpose([x,y])
z1 = griddata(points,rho,(grid_x,grid_y),method='nearest')
z2 = griddata(points,scalar,(grid_x,grid_y),method='nearest')

# Density
imD = ax[0].imshow(z1, origin="lower", aspect='equal',extent=[xmin, xmax, ymin, ymax], norm=mpl.colors.LogNorm())
imS = ax[1].imshow(z2, origin="lower", aspect='equal',extent=[xmin, xmax, ymin, ymax], cmap='RdBu')

ax[0].set_xticks([])
ax[1].set_xticks([])
ax[0].set_yticks([])
ax[1].set_yticks([])

divider = make_axes_locatable(ax[0])
cax = divider.append_axes("right", size="5%", pad=0.05)
cbD =plt.colorbar(imD, cax=cax, label='density')
cbD.ax.yaxis.set_label_coords(2, 0.5)

divider = make_axes_locatable(ax[1])
cax = divider.append_axes("right", size="5%", pad=0.05)
plt.colorbar(imS, cax=cax, label='passive scalar')


fig.savefig('mixing-scalar.pdf',bbox_inches='tight')

# Check results against reference solution
visu_ramses.check_solution(data["data"],'mixing-scalar', tolerance={"all":1e-10})
