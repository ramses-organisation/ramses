# REMARK: it is normal the reference value for density is negative. This is the sum of the log of the cell densities.
#
# Decaying turbulence (turb_type=3). Despite living in the driving module this
# is not a driven run at all: init_flow_fine applies the initial field once as a
# velocity kick, turb_check_time then zeroes it, and every forcing path in
# amr_step/courant_fine/newdt_fine is switched off by 'turb_type /= 3'.
#
# So the bottom panel shows the decay of the kinetic energy, which is the part
# that actually evolves. The turbulent rms is expected to be identically zero at
# every step; a non-zero value would mean forcing is leaking into a run that is
# supposed to be freely decaying, so it is reported rather than plotted.

import matplotlib as mpl
mpl.use('Agg')
import os
import numpy as np
import matplotlib.pyplot as plt
import visu_ramses
import turb_log
from scipy.interpolate import griddata

TESTNAME = 'turb/driving-decaying'

# Two rows of projections, plus a full-width panel for the energy history
fig = plt.figure(figsize=(12, 9))
gs = fig.add_gridspec(nrows=3, ncols=3, hspace=0.3)
ax = np.array([[fig.add_subplot(gs[i, j]) for j in range(3)] for i in range(2)])
ax_ekin = fig.add_subplot(gs[2, :])

# Load RAMSES output
data = visu_ramses.load_snapshot(2)
x      = data["data"]["x"]
y      = data["data"]["y"]
z      = data["data"]["z"]
dx     = data["data"]["dx"]
rho    = data["data"]["density"]
vx     = data["data"]["velocity_x"]
vy     = data["data"]["velocity_y"]
vz     = data["data"]["velocity_z"]

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
z2 = griddata(points,vx ,(grid_x,grid_y, grid_z),method='nearest')
z3 = griddata(points,vy ,(grid_x,grid_y, grid_z),method='nearest')
z4 = griddata(points,vz ,(grid_x,grid_y, grid_z),method='nearest')

rho_proj2 = np.sum(z1, axis=0) #proj along y-axis
rho_proj1 = np.sum(z1, axis=1) #proj along x-axis
rho_proj3 = np.sum(z1, axis=2) #proj along z-axis
vx_proj = np.sum(z2, axis=1)
vy_proj = np.sum(z3, axis=0)
vz_proj = np.sum(z4, axis=2)

im1 = ax[0,0].imshow(rho_proj1, origin="lower", aspect='equal', extent=[xmin, xmax, ymin, ymax])
im2 = ax[0,1].imshow(rho_proj2, origin="lower", aspect='equal', extent=[xmin, xmax, ymin, ymax])
im3 = ax[0,2].imshow(rho_proj3, origin="lower", aspect='equal', extent=[xmin, xmax, ymin, ymax])
im4 = ax[1,0].imshow(vx_proj, origin="lower", aspect='equal', extent=[xmin, xmax, ymin, ymax], cmap='magma')
im5 = ax[1,1].imshow(vy_proj, origin="lower", aspect='equal', extent=[xmin, xmax, ymin, ymax], cmap='magma')
im6 = ax[1,2].imshow(vz_proj, origin="lower", aspect='equal', extent=[xmin, xmax, ymin, ymax], cmap='magma')

cb = []
cb.append(plt.colorbar(im1, ax=ax[0,0], label='Density'))
cb.append(plt.colorbar(im2, ax=ax[0,1], label='Density'))
cb.append(plt.colorbar(im3, ax=ax[0,2], label='Density'))
cb.append(plt.colorbar(im4, ax=ax[1,0], label='Velocity_x'))
cb.append(plt.colorbar(im5, ax=ax[1,1], label='Velocity_y'))
cb.append(plt.colorbar(im6, ax=ax[1,2], label='Velocity_z'))

for i in [0,1]:
    ax[i,0].set_xlabel('y')
    ax[i,0].set_ylabel('z')
    ax[i,1].set_xlabel('x')
    ax[i,1].set_ylabel('z')
    ax[i,2].set_xlabel('x')
    ax[i,2].set_ylabel('y')

for c in cb:
    c.ax.yaxis.set_label_coords(-1.1, 0.5)

# Energy history. This is only a diagnostic: it is not compared against the
# reference solution, and a failure to read the log must not break the test.
try:
    hist = turb_log.read_turb_history(os.path.join('..', '..', 'test_suite.log'), TESTNAME)
except Exception as e:
    hist = turb_log.read_turb_history('', TESTNAME)
    print("Could not read energy history: %s" % e)

# The rms is identically zero here, so plotting it would only add a flat line at
# zero; it is checked and reported as text instead.
turb_log.plot_history(ax_ekin, hist, show_rms=False,
                      missing='No energy history found in ../../test_suite.log')

rms = hist['rms']
if len(rms) > 0:
    # No forcing should ever be applied in a decaying run
    if np.all(rms == 0.0):
        note = 'turbulent forcing rms = 0 at all %d steps (no driving)' % len(rms)
        colour = 'C2'
    else:
        note = 'WARNING: turbulent forcing rms is not zero (max %.3e)' % np.abs(rms).max()
        colour = 'C3'
    ax_ekin.text(0.99, 0.95, note, ha='right', va='top', fontsize='small',
                 color=colour, transform=ax_ekin.transAxes)

fig.savefig('driving-decaying.pdf',bbox_inches='tight')

# Check results against reference solution
visu_ramses.check_solution(data["data"],'driving-decaying', threshold=1e-30, overwrite=False)
