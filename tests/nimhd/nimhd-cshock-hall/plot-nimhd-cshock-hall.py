import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import visu_ramses
from scipy.interpolate import interp1d

# Make figure
fig = plt.figure()
ratio = 0.6
sizex = 12.0
fig.set_size_inches(sizex,ratio*sizex)

ax1 = fig.add_subplot(231)
ax2 = fig.add_subplot(232)
ax3 = fig.add_subplot(233)
ax4 = fig.add_subplot(234)
ax5 = fig.add_subplot(235)
ax6 = fig.add_subplot(236)

# Load RAMSES output
data = visu_ramses.load_snapshot(2)
x    = data["x"]
dens = data["density"]
vx   = data["velocity_x"]
vy   = data["velocity_y"]
vz   = data["velocity_z"]
By   = 0.5*(data["B_left_y"] + data["B_right_y"])
Bz   = 0.5*(data["B_left_z"] + data["B_right_z"])
P    = data["thermal_pressure"]

# Analytical solution
data_ana = np.loadtxt('nimhd-cshock-hall-ana.dat')
x_ana = data_ana[:,0]

# Bin the data in r to avoid having too many symbols in figure
xmin = 0.0
xmax = 1.0
nx   = 101
x_edges = np.linspace(xmin,xmax,nx)
B0, xedges1 = np.histogram(x,bins=(x_edges))
B1, xedges1 = np.histogram(x,bins=(x_edges),weights=dens)
B2, xedges1 = np.histogram(x,bins=(x_edges),weights=vx)
B3, xedges1 = np.histogram(x,bins=(x_edges),weights=vy)
#B4, xedges1 = np.histogram(x,bins=(x_edges),weights=vz)
B4, xedges1 = np.histogram(x,bins=(x_edges),weights=By)
B5, xedges1 = np.histogram(x,bins=(x_edges),weights=Bz)

cube = np.where(B0 > 0.0)

rho  = B1[cube]/B0[cube]
ux   = B2[cube]/B0[cube]
uy   = B3[cube]/B0[cube]
#uz   = B4[cube]/B0[cube]
MagBy= B4[cube]/B0[cube]
MagBz= B5[cube]/B0[cube]

xx = np.zeros([nx-1])
for i in range(nx-1):
    xx[i] = 0.5*(x_edges[i]+x_edges[i+1])
x_mesh = xx[cube]

colors = ['b','r','g','k','magenta']

# Density
ax1.plot(x_mesh,rho,'o',mec=colors[0],mfc='None',label='simulation')
ax1.plot(x_ana,data_ana[:,1],color=colors[0],label='analytical')
ax1.set_xlabel('Distance x (cm)')
ax1.set_ylabel('Density')
ax1.set_xlim([xmin,xmax])

# Vx
ax2.plot(x_mesh,ux,'o',mec=colors[1],mfc='None',label='simulation')
ax2.plot(x_ana,data_ana[:,2],color=colors[1],label='analytical')
ax2.set_xlabel('Distance x (cm)')
ax2.set_ylabel('Velocity x')
ax2.set_xlim([xmin,xmax])

# Vy
ax3.plot(x_mesh,uy,'o',mec=colors[2],mfc='None',label='simulation')
ax3.plot(x_ana,data_ana[:,3],color=colors[2],label='analytical')
ax3.set_xlabel('Distance x (cm)')
ax3.set_ylabel('Velocity y')
ax3.set_xlim([xmin,xmax])

## Vz
#ax4.plot(x_mesh,uz,'o',mec=colors[2],mfc='None',label='simulation')
#ax4.plot(x_ana,data_ana[:,4],color=colors[2],label='analytical')
#ax4.set_xlabel('Distance x (cm)')
#ax4.set_ylabel('Velocity z')
#ax4.set_xlim([xmin,xmax])

# By
ax4.plot(x_mesh,MagBy,'o',mec=colors[3],mfc='None',label='simulation')
ax4.plot(x_ana,data_ana[:,5],color=colors[3],label='analytical')
ax4.set_xlabel('Distance x (cm)')
ax4.set_ylabel('By')
ax4.set_xlim([xmin,xmax])

# Bz
ax5.plot(x_mesh,MagBz,'o',mec=colors[4],mfc='None',label='simulation')
ax5.plot(x_ana,data_ana[:,6],color=colors[4],label='analytical')
ax5.set_xlabel('Distance x (cm)')
ax5.set_ylabel('Bz')
ax5.set_xlim([xmin,xmax])

# Errors
icol = [1,2,3,5,6]
var = [rho,ux,uy,MagBy,MagBz]
names = ['rho','Vx','Vy','By','Bz']
for i in range(5):
    y1 = data_ana[:,icol[i]]
    f1 = interp1d(x_ana,y1)
    err1 = var[i]-f1(x_mesh)
    ax6.plot(x_mesh,err1,color=colors[i],label=names[i])
ax6.legend(loc=1,fontsize=8)
ax6.set_xlabel('Distance x (cm)')
ax6.set_ylabel('Error')
ax6.set_xlim([xmin,xmax])

fig.subplots_adjust(wspace=0.33)
fig.savefig('nimhd-cshock-hall.pdf',bbox_inches='tight')

# Check results against reference solution
visu_ramses.check_solution(data,'nimhd-cshock-hall')
