"""Plot results of sink accretion test"""
import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot as plt
from matplotlib import ticker
import visu_ramses
import numpy as np
from scipy.interpolate import griddata

out_end = 2

# Fundamental constants
G = 6.67259e-8 #cm^3 g^-1 s^-2             # gravitational constant
YR = 3.1556926e7 #s                        # 1 year
MSUN = 1.989e33 #g                         # solar mass
MH = 1.6737236e-24 #g                      # hydrogen mass
KB = 1.38064852e-16 #cm^2 g s^-2 K^-1      # Boltzman constant
PC = 3.0857e18 #cm                         # 1 parsec
AU = 1.49597871e13 #cm                     # 1 astronomical unit

# code units
unit_d=1.66e-24
unit_t=3.004683525921981e15
unit_l=3.08567758128200e+18
unit_v=unit_l/unit_t

data = visu_ramses.load_snapshot(out_end)
#for key in data["stellars"].keys():
#    data["data"]["stellar_"+key] = data["stellars"][key]

x   = data["data"]["x"]
y   = data["data"]["y"]
z   = data["data"]["z"]
dx  = data["data"]["dx"]
rho = data["data"]["density"]
p   = data["data"]["pressure"] * unit_d * unit_l**2 / unit_t**2
cs2 = p/(rho * unit_d)
temperature = cs2 * 2.37 * MH /KB
ps1 = data["data"]["scalar_00"]
ps2 = data["data"]["scalar_01"]
ps3 = data["data"]["scalar_02"]

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
z2 = griddata(points,temperature,(grid_x,grid_y, grid_z),method='nearest')
z3 = griddata(points,ps1,(grid_x,grid_y, grid_z),method='nearest')
z4 = griddata(points,ps2,(grid_x,grid_y, grid_z),method='nearest')
z5 = griddata(points,ps3,(grid_x,grid_y, grid_z),method='nearest')

fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(10, 10))

im1 = ax[0,0].imshow(z1[:,:,int(2**7/2.)], origin="lower", aspect='equal', extent=[xmin, xmax, ymin, ymax])
im2 = ax[0,1].imshow(np.log10(z2[:,:,int(2**7/2.)]), origin="lower", aspect='equal', extent=[xmin, xmax, ymin, ymax])
im3 = ax[1,0].imshow(np.log10(z3[:,:,int(2**7/2.)]), origin="lower", aspect='equal', extent=[xmin, xmax, ymin, ymax])
im4 = ax[1,1].imshow(np.log10(z4[:,:,int(2**7/2.)]), origin="lower", aspect='equal', extent=[xmin, xmax, ymin, ymax])
im5 = ax[1,2].imshow(np.log10(z5[:,:,int(2**7/2.)]), origin="lower", aspect='equal', extent=[xmin, xmax, ymin, ymax])

plt.colorbar(im1, ax=ax[0,0], label='Density [H/cc]')
plt.colorbar(im2, ax=ax[0,1], label='log(temperature [K])')
plt.colorbar(im3, ax=ax[1,0], label='Ions 1')
plt.colorbar(im4, ax=ax[1,1], label='Ions 2')
plt.colorbar(im5, ax=ax[1,2], label='Ions 3')

for axis in [ax[0,0],ax[0,1],ax[1,0],ax[1,1],ax[1,2]]:
    axis.set_axis_off()
    axis.scatter(data["sinks"]['x'],data["sinks"]['y'], s = 20, marker='x', color='red')

fig.savefig('stellar-HII.pdf', bbox_inches="tight")

# Why is this so inaccurate on multiple cores?
red_tol = 1.0e-7
tolerance={"scalar_00":red_tol, "scalar_01":red_tol, "velocity_x":red_tol, "velocity_y":red_tol, "velocity_z":red_tol}
visu_ramses.check_solution(data["data"],'stellar-HII',tolerance=tolerance,overwrite=False)
