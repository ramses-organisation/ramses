from pylab import *

fig = matplotlib.pyplot.figure()
ratio = 1.3
sizex = 8.0
fig.set_size_inches(sizex,ratio*sizex)

# Read in data
data   = loadtxt('data.dat')
amrlev = data[:,0]
x      = data[:,1] - 1.5
rho    = data[:,2]
By     = data[:,8]

# Read analytical solution
data_ana = loadtxt('imhd-tube-ana.dat')
x_ana    = data_ana[:,0]
rho_ana  = data_ana[:,1]
By_ana   = data_ana[:,6]

# Density
density = subplot(211)
density.plot(x,rho,'o',color='black',markerfacecolor='none')
density.plot(x_ana,rho_ana,color='red')
density.set_xlabel('Distance (cm)')
density.set_ylabel('Density (g/cm3)')
levels = density.twinx()
majorLocatorY = MultipleLocator(1.0)
levels.yaxis.set_major_locator(majorLocatorY)
levels.plot(x,amrlev,color='black',ls='dotted')
levels.set_ylabel('AMR Level')

# Y magnetic field
magy = subplot(212)
magy.plot(x,By,'o',color='black',markerfacecolor='none')
magy.plot(x_ana,By_ana,color='red')
magy.set_xlabel('Distance (cm)')
magy.set_ylabel('By')

savefig('imhd-tube.pdf',bbox_inches='tight')
