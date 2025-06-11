import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import visu_ramses
from scipy.interpolate import griddata

fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(12, 8))

# Load RAMSES output
data = visu_ramses.load_snapshot(2)
x      = data["data"]["x"]
y      = data["data"]["y"]
z      = data["data"]["z"]
dx     = data["data"]["dx"]
rho    = data["data"]["density"]
p      = data["data"]["pressure"]
l      = data["data"]["level"]

xmin = np.amin(x-0.5*dx)
xmax = np.amax(x+0.5*dx)
ymin = np.amin(y-0.5*dx)
ymax = np.amax(y+0.5*dx)
zmin = np.amin(z-0.5*dx)
zmax = np.amax(z+0.5*dx)

nx  = 2**7
dpx = (xmax-xmin)/float(nx)
dpy = (ymax-ymin)/float(nx)
dpz = (zmax-zmin)/float(nx)
xpx = np.linspace(xmin+0.5*dpx,xmax-0.5*dpx,nx)
ypx = np.linspace(ymin+0.5*dpy,ymax-0.5*dpy,nx)
zpx = np.linspace(zmin+0.5*dpz,zmax-0.5*dpz,nx)
grid_x, grid_y, grid_z = np.meshgrid(xpx,ypx,zpx)
points = np.transpose([x,y,z])
z1 = griddata(points,rho,(grid_x,grid_y, grid_z),method='nearest')
z5 = griddata(points,p,(grid_x,grid_y, grid_z),method='nearest')
zl = griddata(points,l,(grid_x,grid_y, grid_z),method='nearest')

rho_proj2 = np.sum(z1, axis=0) #proj along y-axis
rho_proj1 = np.sum(z1, axis=1) #proj along x-axis
rho_proj3 = np.sum(z1, axis=2) #proj along z-axis
p_proj2 = np.sum(z5, axis=0) #proj along y-axis
p_proj1 = np.sum(z5, axis=1) #proj along x-axis
p_proj3 = np.sum(z5, axis=2) #proj along z-axis
l_proj2 = np.max(zl, axis=0) #proj along y-axis
l_proj1 = np.max(zl, axis=1) #proj along x-axis
l_proj3 = np.max(zl, axis=2) #proj along z-axis


im1 = ax[0,0].imshow(rho_proj1, origin="lower", aspect='equal', extent=[xmin, xmax, ymin, ymax])
im2 = ax[0,1].imshow(rho_proj2, origin="lower", aspect='equal', extent=[xmin, xmax, ymin, ymax])
im3 = ax[0,2].imshow(rho_proj3, origin="lower", aspect='equal', extent=[xmin, xmax, ymin, ymax])

im4 = ax[1,0].imshow(rho_proj1/rho_proj1.T, origin="lower", aspect='equal', extent=[xmin, xmax, ymin, ymax], cmap='RdBu', vmin=0.99, vmax=1.01)
im5 = ax[1,1].imshow(rho_proj2/rho_proj2.T, origin="lower", aspect='equal', extent=[xmin, xmax, ymin, ymax], cmap='RdBu', vmin=0.99, vmax=1.01)
im6 = ax[1,2].imshow(rho_proj3/rho_proj3.T, origin="lower", aspect='equal', extent=[xmin, xmax, ymin, ymax], cmap='RdBu', vmin=0.99, vmax=1.01)

im7 = ax[2,0].imshow(l_proj1, origin="lower", aspect='equal', extent=[xmin, xmax, ymin, ymax])
im8 = ax[2,1].imshow(l_proj2, origin="lower", aspect='equal', extent=[xmin, xmax, ymin, ymax])
im9 = ax[2,2].imshow(l_proj3, origin="lower", aspect='equal', extent=[xmin, xmax, ymin, ymax])

cb = []
cb.append(plt.colorbar(im1, ax=ax[0,0], label='Density'))
cb.append(plt.colorbar(im2, ax=ax[0,1], label='Density'))
cb.append(plt.colorbar(im3, ax=ax[0,2], label='Density'))
cb.append(plt.colorbar(im4, ax=ax[1,0], label='Symmetry difference'))
cb.append(plt.colorbar(im5, ax=ax[1,1], label='Symmetry difference'))
cb.append(plt.colorbar(im6, ax=ax[1,2], label='Symmetry difference'))
cb.append(plt.colorbar(im7, ax=ax[2,0], label='Level'))
cb.append(plt.colorbar(im8, ax=ax[2,1], label='Level'))
cb.append(plt.colorbar(im9, ax=ax[2,2], label='Level'))

for i in [0,1,2]:
    ax[i,0].set_xlabel('y')
    ax[i,0].set_ylabel('z')
    ax[i,1].set_xlabel('x')
    ax[i,1].set_ylabel('z')
    ax[i,2].set_xlabel('x')
    ax[i,2].set_ylabel('y')

for c in cb:
    c.ax.yaxis.set_label_coords(-1.1, 0.5)

fig.savefig('sedov3d.pdf',bbox_inches='tight')

# Check results against reference solution
visu_ramses.check_solution(data["data"],'sedov3d', overwrite=False)
