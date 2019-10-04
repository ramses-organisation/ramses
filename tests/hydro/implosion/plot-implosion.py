import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import visu_ramses
from scipy.interpolate import griddata

fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(12, 6))

# Load RAMSES output
data = visu_ramses.load_snapshot(2)
x      = data["data"]["x"]
y      = data["data"]["y"]
dx     = data["data"]["dx"]
rho    = data["data"]["density"]
p      = data["data"]["pressure"]
lev    = data["data"]["level"]
u      = np.sqrt(data["data"]["velocity_x"]**2 + data["data"]["velocity_y"]**2)

xmin = np.amin(x-0.5*dx)
xmax = np.amax(x+0.5*dx)
ymin = np.amin(y-0.5*dx)
ymax = np.amax(y+0.5*dx)

nx  = 256
dpx = (xmax-xmin)/float(nx)
dpy = (ymax-ymin)/float(nx)
xpx = np.linspace(xmin+0.5*dpx,xmax-0.5*dpx,nx)
ypx = np.linspace(ymin+0.5*dpy,ymax-0.5*dpy,nx)
grid_x, grid_y = np.meshgrid(xpx,ypx)
points = np.transpose([x,y])
z1 = griddata(points,rho,(grid_x,grid_y),method='nearest')
z2 = griddata(points,u  ,(grid_x,grid_y),method='nearest')
z3 = griddata(points,p  ,(grid_x,grid_y),method='nearest')

im1 = ax[0, 0].imshow(z1, origin="lower", aspect='equal', extent=[xmin, xmax, ymin, ymax])
im2 = ax[0, 1].imshow(z2, origin="lower", aspect='equal', extent=[xmin, xmax, ymin, ymax], cmap='magma')
im3 = ax[0, 2].imshow(z3, origin="lower", aspect='equal', extent=[xmin, xmax, ymin, ymax], cmap='jet')
im4 = ax[1, 0].imshow(z1/z1.T, origin="lower", aspect='equal', extent=[xmin, xmax, ymin, ymax], cmap='RdBu', vmin=0.9, vmax=1.1)
im5 = ax[1, 1].imshow(z2/z2.T, origin="lower", aspect='equal', extent=[xmin, xmax, ymin, ymax], cmap='RdBu', vmin=0.9, vmax=1.1)
im6 = ax[1, 2].imshow(z3/z3.T, origin="lower", aspect='equal', extent=[xmin, xmax, ymin, ymax], cmap='RdBu', vmin=0.9, vmax=1.1)

cb = []
cb.append(plt.colorbar(im1, ax=ax[0, 0], label='Density'))
cb.append(plt.colorbar(im2, ax=ax[0, 1], label='Symmetry difference'))
cb.append(plt.colorbar(im3, ax=ax[0, 2], label='Velocity'))
cb.append(plt.colorbar(im4, ax=ax[1, 0], label='Symmetry difference'))
cb.append(plt.colorbar(im5, ax=ax[1, 1], label='Pressure'))
cb.append(plt.colorbar(im6, ax=ax[1, 2], label='Symmetry difference'))

for c in cb:
    c.ax.yaxis.set_label_coords(-1.1, 0.5)

for a in ax.flatten():
    a.set_xlabel('x')
    a.set_ylabel('y')

fig.savefig('implosion.pdf',bbox_inches='tight')

# Check results against reference solution
visu_ramses.check_solution(data["data"],'implosion')
