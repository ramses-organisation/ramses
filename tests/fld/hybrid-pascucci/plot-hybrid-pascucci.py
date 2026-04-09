import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import visu_ramses

fig = plt.figure()
ratio = 0.75
sizex = 12.0
fig.set_size_inches(sizex,ratio*sizex)

ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)

# Load RAMSES output
data = visu_ramses.load_snapshot(2)

au = 1.5e13

#scale_d = data["data"]["unit_d"]
scale_l = data["data"]["unit_l"]
#scale_t = data["data"]["unit_t"]
x    = data["data"]["x"]*scale_l/au
y    = data["data"]["y"]*scale_l/au
z    = data["data"]["z"]*scale_l/au
dx   = data["data"]["dx"]*scale_l/au
T    = np.log10(data["data"]["temperature"])
erad = np.log10(data["data"]["radiative_energy_01"])
phot = np.log10(data["data"]["photon_number_flux_01"])
lev  = data["data"]["level"]

# SLICES =====================================

dx_im = 100.

# Re-centre coordinates
x += -0.5*data["data"]["boxlen"]*scale_l/au #np.amax(x+0.5*dx)
y += -0.5*data["data"]["boxlen"]*scale_l/au #np.amax(y+0.5*dx)
z += -0.5*data["data"]["boxlen"]*scale_l/au #np.amax(z+0.5*dx)

dist = np.sqrt(x**2+y**2+z**2) - np.sqrt(3.0)*0.5*dx

cube = np.where(np.logical_and(np.abs(x) <= 0.5000000001*dx,np.abs(dist) <= dx_im*0.5*np.sqrt(2.0)))
im_x = y[cube]
im_y = z[cube]

xmin = -0.5*dx_im
xmax =  0.5*dx_im
ymin = -0.5*dx_im
ymax =  0.5*dx_im

nx  = 128
dpx = (xmax-xmin)/float(nx)
dpy = (ymax-ymin)/float(nx)
xpx = np.linspace(xmin+0.5*dpx,xmax-0.5*dpx,nx)
ypx = np.linspace(ymin+0.5*dpy,ymax-0.5*dpy,nx)
grid_x, grid_y = np.meshgrid(xpx,ypx)
points = np.transpose([im_x,im_y])
z1 = griddata(points,T[cube]   ,(grid_x,grid_y),method='nearest')
z2 = griddata(points,erad[cube],(grid_x,grid_y),method='nearest')
z3 = griddata(points,phot[cube],(grid_x,grid_y),method='nearest')
z4 = griddata(points,lev[cube] ,(grid_x,grid_y),method='nearest')

nc=21
im1 = ax1.contourf(xpx,ypx,z1,nc,cmap='hot')
im2 = ax2.contourf(xpx,ypx,z2,nc,cmap='jet')
im3 = ax3.contourf(xpx,ypx,z3,nc,cmap='jet')
im4 = ax4.contourf(xpx,ypx,z4,nc,cmap='viridis')

cb1 = plt.colorbar(im1,ax=ax1,label='log(temperature)')
cb2 = plt.colorbar(im2,ax=ax2,label='log(radiative energy)')
cb3 = plt.colorbar(im3,ax=ax3,label='log(photon number flux)')
cb4 = plt.colorbar(im4,ax=ax4,label='level')
cb1.ax.yaxis.set_label_coords(-1.1,0.5)
cb2.ax.yaxis.set_label_coords(-1.1,0.5)
cb3.ax.yaxis.set_label_coords(-1.1,0.5)
cb4.ax.yaxis.set_label_coords(-1.1,0.5)

ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax3.set_xlabel('x')
ax3.set_ylabel('y')
ax4.set_xlabel('x')
ax4.set_ylabel('y')
ax1.set_aspect('equal')
ax2.set_aspect('equal')
ax3.set_aspect('equal')
ax4.set_aspect('equal')
ax1.set_xlim([xmin,xmax])
ax1.set_ylim([ymin,ymax])
ax2.set_xlim([xmin,xmax])
ax2.set_ylim([ymin,ymax])
ax3.set_xlim([xmin,xmax])
ax3.set_ylim([ymin,ymax])
ax4.set_xlim([xmin,xmax])
ax4.set_ylim([ymin,ymax])

fig.subplots_adjust(wspace=0.25)
fig.savefig('hybrid-pascucci.pdf',bbox_inches='tight')

# Check results against reference solution
visu_ramses.check_solution(data["data"],'hybrid-pascucci')
