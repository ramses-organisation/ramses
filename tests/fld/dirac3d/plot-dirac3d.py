import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import visu_ramses

# Compute analytical solution
def analytical_solution(r,t):
    # Compute analytical solution
    chi=1.0e10
    E0=1.0e5
    ana2 = E0/(8.0*(chi*np.pi*t)**1.5)
    E_ana = ana2*np.exp(-r**2/(4.0*t*chi)) + 1.0
    return E_ana

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
data  = visu_ramses.load_snapshot(2)
dxmin = np.amin(data["data"]["dx"])
x = data["data"]["x"]-0.5
y = data["data"]["y"]-0.5
z = data["data"]["z"]-0.5
r = np.sqrt(x**2 + y**2 + z**2)
Er = data["data"]["radiative_energy_01"]
t  = data["data"]["time"]

# Bin the data in r to avoid having too many symbols in figure
rmin = 0.0
rmax = 0.7
nr   = 301
r_edges = np.linspace(rmin,rmax,nr)
E0, xedges1 = np.histogram(r,bins=(r_edges))
E1, xedges1 = np.histogram(r,bins=(r_edges),weights=Er)
cube = np.where(E0 > 0.0)
E2 = E1[cube]/E0[cube]
rr = np.zeros([nr-1])
for i in range(nr-1):
    rr[i] = 0.5*(r_edges[i]+r_edges[i+1])
r_mesh = rr[cube]

# Generate analytical solutions
r_ana  = np.linspace(0.0,1.0,100)
Er_ana = analytical_solution(r_ana,t)
Er_mesh = analytical_solution(r_mesh,t)
error  = abs(E2 - Er_mesh)/Er_mesh

# Er(r) profile
ax1.plot(r_mesh,np.log10(E2),'o',mec='k',mfc='None',label='simulation')
ax1.plot(r_ana,np.log10(Er_ana),color='r',label='analytical')
ax1.plot([0.0,1.0],[-10.0,-10.0],'.',color='grey',label='Error')
ax1.set_xlabel('Distance (cm)')
ax1.set_ylabel('log(Er)')
ax1.legend(loc=1,fontsize=12)
er1 = ax1.twinx()
er1.plot(r_mesh,np.log10(error),'.',color='grey')
er1.set_ylabel('log(Error)')
ax1.set_xlim([0.0,0.7])
ax1.set_ylim([-1.0,8.0])

# Er(x,z) and error maps
cube = np.where(np.abs(y)<=0.51*data["data"]["dx"])
slice_x = x[cube]
slice_z = z[cube]
slice_E = np.log10(Er[cube])

nx = 128
xmin = ymin = -0.5
xmax = ymax =  0.5
dpx = (xmax-xmin)/float(nx)
dpy = (ymax-ymin)/float(nx)
xpx = np.linspace(xmin+0.5*dpx,xmax-0.5*dpx,nx)
ypx = np.linspace(ymin+0.5*dpy,ymax-0.5*dpy,nx)
grid_x, grid_y = np.meshgrid(xpx,ypx)
points = np.transpose([slice_x,slice_z])
map_E = griddata(points,slice_E,(grid_x,grid_y),method='nearest')

im1 = ax2.contourf(xpx,ypx,map_E)
cb1 = plt.colorbar(im1,ax=ax2,label='log(Er)')
ax2.set_xlabel('Distance x (cm)')
ax2.set_ylabel('Distance z (cm)')
ax2.set_xlim([xmin,xmax])
ax2.set_ylim([ymin,ymax])

map_R = np.zeros([nx,nx])
for j in range(nx):
    r = np.sqrt(xpx**2 + ypx[j]**2)
    map_R[:,j] = abs(np.power(10.0,map_E[:,j]) - analytical_solution(r,t))/analytical_solution(r,t)

im2 = ax3.contourf(xpx,ypx,np.log10(map_R),cmap='cubehelix')
cb2 = plt.colorbar(im2,ax=ax3,label='Error')
ax3.set_xlabel('Distance x (cm)')
ax3.set_ylabel('Distance z (cm)')
ax3.set_xlim([xmin,xmax])
ax3.set_ylim([ymin,ymax])


# Er(x,y) map
cube = np.where(np.abs(z)<=0.51*data["data"]["dx"])
slice_x = x[cube]
slice_y = y[cube]
slice_E = np.log10(Er[cube])
points = np.transpose([slice_x,slice_y])
map_E = griddata(points,slice_E,(grid_x,grid_y),method='nearest')

im3 = ax4.contourf(xpx,ypx,map_E)
cb3 = plt.colorbar(im3,ax=ax4,label='log(Er)')
ax4.set_xlabel('Distance x (cm)')
ax4.set_ylabel('Distance y (cm)')
ax4.set_xlim([xmin,xmax])
ax4.set_ylim([ymin,ymax])


fig.subplots_adjust(wspace=0.25)
fig.savefig('dirac3d.pdf',bbox_inches='tight')

# Check results against reference solution
visu_ramses.check_solution(data["data"],'dirac3d',tolerance={"velocity_x":1.5e-12,"velocity_y":1.5e-12,"velocity_z":1.5e-12})
