import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy.special import erf
import visu_ramses

# Compute analytical solution
def analytical_solution(x,dx,t):
    b     = 0.00108967
    B_ana = 0.5*(erf((-x+b/2.0)/np.sqrt(4.0*t))+erf((x+b/2.0)/np.sqrt(4.0*t)))
    return B_ana

# Make figure
fig = plt.figure()
ratio = 0.75
sizex = 12.0
fig.set_size_inches(sizex,ratio*sizex)

ax1 = plt.subplot2grid((2, 2), (0, 0), colspan=2)
ax2 = plt.subplot2grid((2, 2), (1, 0))
ax3 = plt.subplot2grid((2, 2), (1, 1))

# Load RAMSES output
data  = visu_ramses.load_snapshot(2)
dxmin = np.amin(data["dx"])
x = data["x"]-(0.5+0.5*dxmin)
z = data["z"]-(0.5+0.5*dxmin)
rc = np.sqrt(x**2 + z**2)
By = 0.5*(data["B_left_y"] + data["B_right_y"])
t  = data["time"]

# Bin the data in r to avoid having too many symbols in figure
rmin = 0.0
rmax = 0.7
nr   = 301
r_edges = np.linspace(rmin,rmax,nr)
B0, xedges1 = np.histogram(rc,bins=(r_edges))
B1, xedges1 = np.histogram(rc,bins=(r_edges),weights=By)
cube = np.where(B0 > 0.0)
B2 = B1[cube]/B0[cube]
rr = np.zeros([nr-1])
for i in range(nr-1):
    rr[i] = 0.5*(r_edges[i]+r_edges[i+1])
r_mesh = rr[cube]

# Generate analytical solutions
r_ana  = np.linspace(rmin,rmax,100)
By_ana = analytical_solution(r_ana,dxmin,t)
error  = np.where(np.isnan(B2),0.0,abs(B2 - analytical_solution(r_mesh,dxmin,t))/analytical_solution(r_mesh,dxmin,t))

# By(r) profile
#ax1.scatter(r_mesh,np.log10(B2),marker='o',edgecolor='k',color='None',label='simulation')
ax1.plot(r_mesh,np.log10(B2),'o',mec='k',mfc='None',label='simulation')
ax1.plot(r_ana,np.log10(By_ana),color='r',label='analytical')
ax1.plot([0.0,1.0],[10.0,10.0],'.',color='grey',label='Error')
ax1.set_xlabel('Distance (cm)')
ax1.set_ylabel('By')
ax1.legend(loc=6)
er1 = ax1.twinx()
#er1.plot(r_mesh,np.cumsum(error),color='grey',ls='dotted')
er1.plot(r_mesh,error,'.',color='grey')
er1.set_ylabel('Error')
ax1.set_xlim([rmin,rmax])
ax1.set_ylim([-18,0])
er1.set_ylim([-1,9])

# By(x,z) and error maps
cube = np.where(np.abs(data["y"]-0.51)<=0.51*data["dx"])
slice_x = x[cube]
slice_z = z[cube]
slice_B = By[cube]

nx = 128
xmin = ymin = -0.5
xmax = ymax =  0.5
dpx = (xmax-xmin)/float(nx)
dpy = (ymax-ymin)/float(nx)
xpx = np.linspace(xmin+0.5*dpx,xmax-0.5*dpx,nx)
ypx = np.linspace(ymin+0.5*dpy,ymax-0.5*dpy,nx)
grid_x, grid_y = np.meshgrid(xpx,ypx)
points = np.transpose([slice_x,slice_z])
map_B = griddata(points,slice_B,(grid_x,grid_y),method='nearest')

im1 = ax2.contourf(xpx,ypx,map_B)
cb1 = plt.colorbar(im1,ax=ax2,label='By')
cb1.ax.yaxis.set_label_coords(-1.1,0.5)
ax2.set_xlabel('Distance x (cm)')
ax2.set_ylabel('Distance z (cm)')
ax2.set_xlim([-0.25,0.25])
ax2.set_ylim([-0.25,0.25])

map_E = np.zeros([nx,nx])
for j in range(nx):
    r = np.sqrt(xpx**2 + ypx[j]**2)
    map_E[:,j] = (map_B[:,j] - analytical_solution(r,dxmin,t))/analytical_solution(r,dxmin,t)

im2 = ax3.contourf(xpx,ypx,map_E,cmap='RdBu',levels=np.linspace(-0.01,0.01,11),extend='both')
cb2 = plt.colorbar(im2,ax=ax3,label='Error')
cb2.ax.yaxis.set_label_coords(-1.1,0.5)
ax3.set_xlabel('Distance x (cm)')
ax3.set_ylabel('Distance z (cm)')
ax3.set_xlim([-0.25,0.25])
ax3.set_ylim([-0.25,0.25])

fig.subplots_adjust(wspace=0.25)
fig.savefig('nimhd-diffusion-ohm.pdf',bbox_inches='tight')

# Check results against reference solution
visu_ramses.check_solution(data,'nimhd-diffusion-ohm')
