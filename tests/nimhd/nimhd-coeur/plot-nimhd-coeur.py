import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import visu_ramses

fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(12, 6))
ax = ax.flatten()

# Load RAMSES output
data = visu_ramses.load_snapshot(2)

# Fundamental constants
G = 6.67259e-8 #cm^3 g^-1 s^-2             # gravitational constant
YR = 3.1556926e7 #s                        # 1 year
MSUN = 1.989e33 #g                         # solar mass
MH = 1.6737236e-24 #g                      # hydrogen mass
KB = 1.38064852e-16 #cm^2 g s^-2 K^-1      # Boltzman constant
PC = 3.0857e18 #cm                         # 1 parsec
AU = 1.49597871e13 #cm                     # 1 astronomical unit

au = 1.5e13

data["boxlen"] = 1.9140529831524967e-2

scale_d = 2.31*1.6605390e-24
scale_l = 3.08e18
scale_t = 1.0/np.sqrt(6.67e-8*scale_d)
scale_b = np.sqrt(4.0*np.pi*scale_d*(scale_l/scale_t)**2)
x    = data["data"]["x"]*scale_l/au
y    = data["data"]["y"]*scale_l/au
z    = data["data"]["z"]*scale_l/au
dx   = data["data"]["dx"]*scale_l/au
lev  = data["data"]["level"]
rho  = np.log10(data["data"]["density"]*scale_d)
ux   = data["data"]["velocity_x"]*scale_l/scale_t/1.0e5
uy   = data["data"]["velocity_y"]*scale_l/scale_t/1.0e5
uz   = data["data"]["velocity_z"]*scale_l/scale_t/1.0e5
bx   = 0.5*(data["data"]["B_x_left"]+data["data"]["B_x_right"])*scale_b
by   = 0.5*(data["data"]["B_y_left"]+data["data"]["B_y_right"])*scale_b
bz   = 0.5*(data["data"]["B_z_left"]+data["data"]["B_z_right"])*scale_b
B    = np.log10(np.sqrt(bx**2 + by**2 + bz**2))
cs2 = (data["data"]["pressure"]/data["data"]["density"]) * (scale_l/scale_t)**2
T = np.log10(cs2 * 2.37 * MH /KB)

dmin = -19.5
dmax = -9.0
tmin = 0.8
tmax = 3.0
bmin = -5.0
bmax = 1.0

nx = 129
# Construct some edge specifiers for the histogram2d function call
d_edges = np.linspace(dmin,dmax,nx)
t_edges = np.linspace(tmin,tmax,nx)
b_edges = np.linspace(bmin,bmax,nx)
# Call the numpy histogram2d function
za, yedges1, xedges1 = np.histogram2d(T,rho,bins=(t_edges,d_edges))
zb, yedges1, xedges1 = np.histogram2d(B,rho,bins=(b_edges,d_edges))
with np.errstate(divide="ignore",invalid="ignore"):
    z1 = np.log10(za)
    z2 = np.log10(zb)
# In the contour plots, x and y are the centers of the cells, instead of the edges.
d_mesh = np.zeros([nx-1])
t_mesh = np.zeros([nx-1])
b_mesh = np.zeros([nx-1])
for i in range(nx-1):
    d_mesh[i] = 0.5*(d_edges[i]+d_edges[i+1])
for i in range(nx-1):
    t_mesh[i] = 0.5*(t_edges[i]+t_edges[i+1])
for i in range(nx-1):
    b_mesh[i] = 0.5*(b_edges[i]+b_edges[i+1])
# Plot histograms
nc = 21
cont1 = ax[0].contourf(d_mesh,t_mesh,z1,nc,cmap='Reds')
cont2 = ax[3].contourf(d_mesh,b_mesh,z2,nc,cmap='Blues')
cont3 = ax[0].contour (d_mesh,t_mesh,za,colors='r',levels=[1.0])
cont4 = ax[3].contour (d_mesh,b_mesh,zb,colors='b',levels=[1.0])

ax[0].set_xlabel('log(rho)')
ax[0].set_ylabel('log(T)')
ax[3].set_xlabel('log(rho)')
ax[3].set_ylabel('log(B)')


# SLICES =====================================

dx_im = 400.0

# Re-centre coordinates
x += -0.5*data["boxlen"]*scale_l/au #np.amax(x+0.5*dx)
y += -0.5*data["boxlen"]*scale_l/au #np.amax(y+0.5*dx)
z += -0.5*data["boxlen"]*scale_l/au #np.amax(z+0.5*dx)

dist = np.sqrt(x**2+y**2+z**2) - np.sqrt(3.0)*0.5*dx

cube = np.where(np.logical_and(np.abs(z) <= 0.5000000001*dx,np.abs(dist) <= dx_im*0.5*np.sqrt(2.0)))
im_x = x[cube]
im_y = y[cube]

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
z1 = griddata(points,rho[cube],(grid_x,grid_y),method='nearest')
z2 = griddata(points,T[cube]  ,(grid_x,grid_y),method='nearest')
z3 = griddata(points,ux[cube] ,(grid_x,grid_y),method='nearest')
z4 = griddata(points,uy[cube] ,(grid_x,grid_y),method='nearest')
z5 = np.around(griddata(points,lev[cube],(grid_x,grid_y),method='nearest'))

nc=21
im1 = ax[1].contourf(xpx,ypx,z1,nc,cmap='jet')
im2 = ax[4].contourf(xpx,ypx,z2,nc,cmap='hot')

ctr = ax[1].contour(xpx,ypx,z5,colors='w',levels=range(0,20))
ax[1].clabel(ctr,inline=1,fmt="%i")
vskip = 6
vec = ax[4].quiver(xpx[::vskip],ypx[::vskip],z3[::vskip,::vskip],z4[::vskip,::vskip],color="w")

cb1 = plt.colorbar(im1,ax=ax[1],label='log(Density)')
cb2 = plt.colorbar(im2,ax=ax[4],label='log(T)')
cb1.ax.yaxis.set_label_coords(-1.1,0.5)
cb2.ax.yaxis.set_label_coords(-1.1,0.5)

ax[1].set_xlabel('x')
ax[1].set_ylabel('y')
ax[4].set_xlabel('x')
ax[4].set_ylabel('y')
ax[1].set_aspect('equal')
ax[4].set_aspect('equal')
ax[1].set_xlim([xmin,xmax])
ax[1].set_ylim([ymin,ymax])
ax[4].set_xlim([xmin,xmax])
ax[4].set_ylim([ymin,ymax])

# Side map with B field streamlines
cube = np.where(np.logical_and(np.abs(y) <= 0.5000000001*dx,np.abs(dist) <= dx_im*0.5*np.sqrt(2.0)))
im_x = x[cube]
im_y = z[cube]
points = np.transpose([im_x,im_y])
z1 = griddata(points,rho[cube],(grid_x,grid_y),method='nearest')
z3 = griddata(points,bx[cube] ,(grid_x,grid_y),method='nearest')
z4 = griddata(points,bz[cube] ,(grid_x,grid_y),method='nearest')
im3 = ax[2].contourf(xpx,ypx,z1,nc,cmap='jet')
stm = ax[2].streamplot(xpx,ypx,z3,z4,color="w")
cb3 = plt.colorbar(im3,ax=ax[2],label='log(Density)')
cb3.ax.yaxis.set_label_coords(-1.1,0.5)

ax[2].set_xlabel('x')
ax[2].set_ylabel('y')
ax[2].set_aspect('equal')
ax[2].set_xlim([xmin,xmax])
ax[2].set_ylim([ymin,ymax])


# Full box map
dx_im = 3800.0
cube = np.where(np.abs(y) <= 0.5000000001*dx)
im_x = x[cube]
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
z1 = griddata(points,rho[cube],(grid_x,grid_y),method='nearest')
z3 = griddata(points,ux[cube] ,(grid_x,grid_y),method='nearest')
z4 = griddata(points,uz[cube] ,(grid_x,grid_y),method='nearest')

z3 = z3 / np.sqrt(z3**2 + z4**2)
z4 = z4 / np.sqrt(z3**2 + z4**2)

im4 = ax[5].contourf(xpx,ypx,z1,nc,cmap='jet')
vskip = 6
vec = ax[5].quiver(xpx[::vskip],ypx[::vskip],z3[::vskip,::vskip],z4[::vskip,::vskip],color="w")
cb4 = plt.colorbar(im4,ax=ax[5],label='log(Density)')
cb4.ax.yaxis.set_label_coords(-1.1,0.5)
ax[5].set_xlabel('x')
ax[5].set_ylabel('y')
ax[5].set_aspect('equal')
ax[5].set_xlim([xmin,xmax])
ax[5].set_ylim([ymin,ymax])

fig.savefig('nimhd-coeur.pdf',bbox_inches='tight')

# Check results against reference solution
visu_ramses.check_solution(data["data"],'nimhd-coeur', overwrite=True)
