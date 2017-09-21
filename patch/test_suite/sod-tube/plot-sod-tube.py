from pylab import *

fig = matplotlib.pyplot.figure()
ratio = 1.5
sizex = 8.0
fig.set_size_inches(sizex,ratio*sizex)

# Read in data
data   = loadtxt('data.dat')
amrlev = data[:,0]
x      = data[:,1]
rho    = data[:,2]
u      = data[:,3]
p      = data[:,4]

# Read analytical solution
data_ana = loadtxt('sod-tube-ana.dat')
x_ana    = data_ana[:,1]
rho_ana  = data_ana[:,3]
u_ana    = data_ana[:,2]
p_ana    = data_ana[:,4]

# Density
density = subplot(311)
density.plot(x,rho,'o',color='black',markerfacecolor='none')
density.plot(x_ana,rho_ana,color='red')
density.set_xlabel('Distance (cm)')
density.set_ylabel('Density (g/cm3)')
levels = density.twinx()
majorLocatorY = MultipleLocator(1.0)
levels.yaxis.set_major_locator(majorLocatorY)
levels.plot(x,amrlev,color='black',ls='dotted')
levels.set_ylabel('AMR Level')

# Velocity
velocity = subplot(312)
velocity.plot(x,u,'o',color='black',markerfacecolor='none')
velocity.plot(x_ana,u_ana,color='red')
velocity.set_xlabel('Distance (cm)')
velocity.set_ylabel('Velocity (cm/s)')

# Pressure
pressure = subplot(313)
pressure.plot(x,p,'o',color='black',markerfacecolor='none')
pressure.plot(x_ana,p_ana,color='red')
pressure.set_xlabel('Distance (cm)')
pressure.set_ylabel('Pressure (g/cm/s2)')

savefig('sod-tube.pdf',bbox_inches='tight')
