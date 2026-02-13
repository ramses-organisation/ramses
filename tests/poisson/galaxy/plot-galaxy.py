import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import visu_ramses
from matplotlib.colors import LogNorm
import numpy as np
from scipy.interpolate import griddata


# Execute the map part of the test

fig = plt.figure(figsize=(6, 4))
axes = fig.subplots(nrows=2, ncols=2)

# Load RAMSES output
data = visu_ramses.load_snapshot(2)
# particles
print(data['data'].keys())
print(data['particle'].keys())
print(set(data['particle']['family']))
print(set(data['particle']['tag']))
xp = data["particle"]["position_x"]
yp = data["particle"]["position_y"]
zp = data["particle"]["position_z"]
mp = data["particle"]["mass"]

# gas
x_all  = data["data"]["x"]
y_all  = data["data"]["y"]
z_all  = data["data"]["z"]
dx_all = data["data"]["dx"]
rho_all= data["data"]["density"]

filt = (
    (x_all > 150) & (x_all < 250) &
    (y_all > 150) & (y_all < 250) &
    (z_all > 150) & (z_all < 250)
)

x = x_all[filt]
y = y_all[filt]
z = z_all[filt]
dx = dx_all[filt]
rho = rho_all[filt]

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
rho_proj3 = np.sum(z1, axis=2)
rho_proj2 = np.transpose(np.sum(z1, axis=1))

# DM particles
im = axes[0,0].hist2d(xp,yp,weights=mp,bins=256,range=[[150, 250], [150, 250]],norm=LogNorm())
im = axes[1,0].hist2d(xp,zp,weights=mp,bins=256,range=[[150, 250], [150, 250]],norm=LogNorm())
#plt.colorbar(im, ax=axes[1])
# Stars

# gas
im2 = axes[0,1].imshow(rho_proj3, origin="lower", aspect='equal', extent=[xmin, xmax, ymin, ymax])
im3 = axes[1,1].imshow(rho_proj2, origin="lower", aspect='equal', extent=[xmin, xmax, ymin, ymax])
# refinement level


#axes[0].scatter(xp,yp,s=1,alpha=0.005)
#axes[1].scatter(xp,zp,s=1,alpha=0.005)


for ax in axes.flatten():
    ax.axis('equal')

fig.savefig("galaxy.pdf", bbox_inches="tight")

to_check = {}
for key in data["data"].keys():
    to_check[key] = data["data"][key]
for key in data["particle"].keys():
    to_check['particle '+key] = data["particle"][key]

# Check results against reference solution
visu_ramses.check_solution(to_check,'galaxy', overwrite=False)
