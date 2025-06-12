import matplotlib as mpl
mpl.use('Agg')
#mpl.use('MacOSX')
#mpl.use('QtAgg')
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import visu_ramses

# Make figure
fig = plt.figure()
ratio = 0.75/2
sizex = 12.0
fig.set_size_inches(sizex,ratio*sizex)

ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

# Load RAMSES output
data = visu_ramses.load_snapshot(2)
scale_d = data["data"]["unit_d"]
scale_l = data["data"]["unit_l"]
scale_t = data["data"]["unit_t"]
scale_v = scale_l/scale_t
x = data["data"]["x"]
y = data["data"]["y"]
rho = data["data"]["density"]*scale_d
P = data["data"]["pressure"]*scale_d*scale_v**2
mu_gas = 1.e0       # mean molecular weight
mH = 1.6737236e-24  #g                      # hydrogen mass
kB = 1.38064852e-16 #cm^2 g s^-2 K^-1       # Boltzman constant
T = P/rho * mu_gas * mH/kB
ion1 = data["data"]["scalar_00"]

# 2D maps
slice_x = x
slice_y = y
slice_T = np.log10(T)
slice_ion1 = np.log10(ion1)

nx = 128
xmin = ymin = 0.0
xmax = ymax = data["data"]["boxlen"]
dpx = (xmax-xmin)/float(nx)
dpy = (ymax-ymin)/float(nx)
xpx = np.linspace(xmin+0.5*dpx,xmax-0.5*dpx,nx)
ypx = np.linspace(ymin+0.5*dpy,ymax-0.5*dpy,nx)
grid_x, grid_y = np.meshgrid(xpx,ypx)
points = np.transpose([slice_x,slice_y])
map_T = griddata(points,slice_T,(grid_x,grid_y),method='nearest')
map_ion1 = griddata(points,slice_ion1,(grid_x,grid_y),method='nearest')

im1 = ax1.contourf(xpx,ypx,map_ion1,cmap='rainbow',levels=np.linspace(-5.,1.,100))
im2 = ax2.contourf(xpx,ypx,map_T,cmap='rainbow',levels=np.linspace(2.,5.,100))
cb1 = plt.colorbar(im1,ax=ax1,label='log(ion1)')
cb2 = plt.colorbar(im2,ax=ax2,label='log(Temperature)')
ax1.set_xlabel('Distance x (pc)')
ax1.set_ylabel('Distance z (pc)')
ax1.set_xlim([xmin,xmax])
ax1.set_ylim([ymin,ymax])
ax2.set_xlabel('Distance x (pc)')
ax2.set_ylabel('Distance z (pc)')
ax2.set_xlim([xmin,xmax])
ax2.set_ylim([ymin,ymax])
ax2.set_title('t='+str(round(data["info"]["time"],2))+' Myr')

fig.subplots_adjust(wspace=0.35)
fig.savefig('crossing-beams.pdf',bbox_inches='tight')

plt.show()

# Check results against reference solution
visu_ramses.check_solution(data["data"], 'crossing-beams', tolerance={"all":3.0e-06},overwrite=False)
