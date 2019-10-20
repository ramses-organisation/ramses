import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import visu_ramses

# Make figure
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
scale_d = data["data"]["unit_d"]
scale_l = data["data"]["unit_l"]
scale_t = data["data"]["unit_t"]
scale_v = scale_l/scale_t
scale_b = np.sqrt(4.0*np.pi*scale_d*(scale_l/scale_t)**2)
x = data["data"]["x"]
y = data["data"]["y"]
r = np.sqrt(x**2 + y**2)
rho = data["data"]["density"]*scale_d
P = data["data"]["pressure"]*scale_d*scale_v**2
ps1 = data["data"]["scalar_01"]
ps2 = data["data"]["scalar_02"]
ps3 = data["data"]["scalar_03"]
ux  = data["data"]["velocity_x"]
uy  = data["data"]["velocity_y"]

# Bin the data in r to avoid having too many symbols in figure
rmin = 0.0
rmax = 1.0
nr   = 301
r_edges = np.linspace(rmin,rmax,nr)
n_bin, xedges1 = np.histogram(r,bins=(r_edges))
rho_bin, xedges1 = np.histogram(r,bins=(r_edges),weights=rho)
P_bin, xedges1 = np.histogram(r,bins=(r_edges),weights=P)
ps1_bin, xedges1 = np.histogram(r,bins=(r_edges),weights=ps1)
ps2_bin, xedges1 = np.histogram(r,bins=(r_edges),weights=ps2)
ps3_bin, xedges1 = np.histogram(r,bins=(r_edges),weights=ps3)
cube = np.where(n_bin > 0.0)
dens = rho_bin[cube]/n_bin[cube]
pres = P_bin[cube]/n_bin[cube]
ion1 = ps1_bin[cube]/n_bin[cube]
ion2 = ps2_bin[cube]/n_bin[cube]
ion3 = ps3_bin[cube]/n_bin[cube]
rr = np.zeros([nr-1])
for i in range(nr-1):
    rr[i] = 0.5*(r_edges[i]+r_edges[i+1])
r_mesh = rr[cube]

# Radial profiles
ax1.plot(r_mesh,np.log10(dens),'o',mec='b',mfc='None',label='density')
ax5 = ax1.twinx()
ax5.plot(r_mesh,np.log10(pres),'o',mec='r',mfc='None',label='Pressure')

ax2.plot(r_mesh,np.log10(ion1),'o',mec='k',mfc='None',label='Ions 1')
ax2.plot(r_mesh,np.log10(ion2),'o',mec='lime',mfc='None',label='Ions 2')
ax2.plot(r_mesh,np.log10(ion3),'o',mec='magenta',mfc='None',label='Ions 3')
ax1.set_xlabel('Distance (pc)')
ax1.set_ylabel('log(Density)')
ax5.set_ylabel('log(Pressure)')
ax1.legend(loc=(0.68,0.72),fontsize=12)
ax5.legend(loc=(0.65,0.2),fontsize=12)
ax1.set_xlim([0.0,rmax])

ax2.set_xlabel('Distance (pc)')
ax2.set_ylabel('log(Ion)')
ax2.legend(loc=1,fontsize=12)
ax2.set_xlim([0.0,rmax])

# 2D maps
slice_x = x
slice_y = y
slice_d = np.log10(rho)
slice_p = np.log10(P)

nx = 128
xmin = ymin = 0.0
xmax = ymax = rmax
dpx = (xmax-xmin)/float(nx)
dpy = (ymax-ymin)/float(nx)
xpx = np.linspace(xmin+0.5*dpx,xmax-0.5*dpx,nx)
ypx = np.linspace(ymin+0.5*dpy,ymax-0.5*dpy,nx)
grid_x, grid_y = np.meshgrid(xpx,ypx)
points = np.transpose([slice_x,slice_y])
map_d = griddata(points,slice_d,(grid_x,grid_y),method='nearest')
map_p = griddata(points,slice_p,(grid_x,grid_y),method='nearest')
map_ux = griddata(points,ux ,(grid_x,grid_y),method='nearest')
map_uy = griddata(points,uy ,(grid_x,grid_y),method='nearest')

im1 = ax3.contourf(xpx,ypx,map_d,cmap='Blues',levels=np.linspace(np.nanmin(map_d),np.nanmax(map_d),10))
im2 = ax4.contourf(xpx,ypx,map_p,cmap='Reds')
vskip=4
vec = ax4.quiver(xpx[::vskip],ypx[::vskip],
	             map_ux[::vskip,::vskip],
	             map_uy[::vskip,::vskip],color="k",pivot='mid',scale=0.05)
cb1 = plt.colorbar(im1,ax=ax3,label='log(Density)')
cb2 = plt.colorbar(im2,ax=ax4,label='log(Pressure)')
ax3.set_xlabel('Distance x (pc)')
ax3.set_ylabel('Distance z (pc)')
ax3.set_xlim([xmin,xmax])
ax3.set_ylim([ymin,ymax])
ax4.set_xlabel('Distance x (pc)')
ax4.set_ylabel('Distance z (pc)')
ax4.set_xlim([xmin,xmax])
ax4.set_ylim([ymin,ymax])

fig.subplots_adjust(wspace=0.35)
fig.savefig('stromgren2d.pdf',bbox_inches='tight')

# Check results against reference solution
visu_ramses.check_solution(data["data"], 'stromgren2d', tolerance={"all":3.0e-06})
