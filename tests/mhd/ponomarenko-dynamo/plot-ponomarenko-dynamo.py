import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import visu_ramses
from scipy.interpolate import griddata

fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(12, 4),constrained_layout=True)

# Load RAMSES output
data = visu_ramses.load_snapshot(2)
x      = data["data"]["x"]
y      = data["data"]["y"]
z      = data["data"]["z"]
dx     = data["data"]["dx"]
b2     = 0.25*(data["data"]["B_x_left"]+data["data"]["B_x_right"])**2 + \
         0.25*(data["data"]["B_y_left"]+data["data"]["B_y_right"])**2 + \
         0.25*(data["data"]["B_z_left"]+data["data"]["B_z_right"])**2

xmin = np.amin(x-0.5*dx)
xmax = np.amax(x+0.5*dx)
ymin = np.amin(y-0.5*dx)
ymax = np.amax(y+0.5*dx)
zmin = np.amin(z-0.5*dx)
zmax = np.amax(z+0.5*dx)

nx  = 64 #128
dpx = (xmax-xmin)/float(nx)
dpy = (ymax-ymin)/float(nx)
dpz = (zmax-zmin)/float(nx)
xpx = np.linspace(xmin+0.5*dpx,xmax-0.5*dpx,nx)
ypx = np.linspace(ymin+0.5*dpy,ymax-0.5*dpy,nx)
zpx = np.linspace(zmin+0.5*dpz,zmax-0.5*dpz,nx)
grid_x, grid_y, grid_z = np.meshgrid(xpx,ypx,zpx)
points = np.transpose([x,y,z])
z1 = griddata(points,b2,(grid_x,grid_y, grid_z),method='nearest')

b2_proj2 = np.sum(z1, axis=0) #proj along y-axis
b2_proj1 = np.sum(z1, axis=1) #proj along x-axis
b2_proj3 = np.sum(z1, axis=2) #proj along z-axis

im1 = ax[0].imshow(np.log10(b2_proj1), origin="lower", aspect='equal', extent=[xmin, xmax, ymin, ymax],vmin=6.5,vmax=11.5)
im2 = ax[1].imshow(np.log10(b2_proj2), origin="lower", aspect='equal', extent=[xmin, xmax, ymin, ymax],vmin=6.5,vmax=11.5)
im3 = ax[2].imshow(np.log10(b2_proj3), origin="lower", aspect='equal', extent=[xmin, xmax, ymin, ymax],vmin=6.5,vmax=11.5)

cb = []
cb.append(plt.colorbar(im1, ax=ax[0], label='log10(B^2)'))
cb.append(plt.colorbar(im2, ax=ax[1], label='log10(B^2)'))
cb.append(plt.colorbar(im3, ax=ax[2], label='log10(B^2)'))

ax[0].set_xlabel('y')
ax[0].set_ylabel('z')
ax[1].set_xlabel('x')
ax[1].set_ylabel('z')
ax[2].set_xlabel('x')
ax[2].set_ylabel('y')

fig.savefig('ponomarenko-dynamo.pdf',bbox_inches='tight')

# Check results against reference solution
visu_ramses.check_solution(data["data"], 'ponomarenko-dynamo') #, tolerance={"all":3.0e-06},overwrite=False)
