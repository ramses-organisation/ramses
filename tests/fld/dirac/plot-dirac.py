import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import visu_ramses

fig = plt.figure()
ratio = 0.8
sizex = 12.0
fig.set_size_inches(sizex,ratio*sizex)

ax1 = plt.subplot(221)
ax2 = plt.subplot(222)
ax3 = plt.subplot(223)
ax4 = plt.subplot(224)
ax5 = ax1.twinx()
ax6 = ax3.twinx()

# Load RAMSES output
data   = visu_ramses.load_snapshot(2)
order  = data["data"]["x"].argsort()
x      = data["data"]["x"][order]
amrlev = data["data"]["level"][order]
er     = data["data"]["radiative_energy_01"][order]
ux     = np.abs(data["data"]["velocity_x"][order])
et     = data["data"]["internal_energy"][order]

# Compute analytical solution
t = data["data"]["time"]
chi=1.0e10
x0=0.5
E0=1.0e5
ana2 = E0/(2.0*(chi*np.pi*t)**.5)
ana = 1.+ana2*np.exp(-(x-x0)**2/(4.0*t*chi))

majorLocatorY = plt.MultipleLocator(1.0)

# Radiative energy
ax1.semilogy(x,er,'o',color='black',markerfacecolor='none')
ax1.semilogy(x,ana,color='red')
ax1.set_xlabel('Distance (cm)')
ax1.set_ylabel('Radiative energy')
ax5.yaxis.set_major_locator(majorLocatorY)
ax5.plot(x,amrlev,color='black',ls='dotted')
ax5.set_ylabel('AMR Level')

# Velocity
ax2.semilogy(x,ux,'o',color='black',markerfacecolor='none')
ax2.set_xlabel('Distance (cm)')
ax2.set_ylabel('Velocity (cm/s)')

# Relative error
ax3.semilogy(x,abs(er-ana)/ana,color='red')
ax3.set_xlabel('Distance (cm)')
ax3.set_ylabel('Relative error')
ax6.yaxis.set_major_locator(majorLocatorY)
ax6.plot(x,amrlev,color='black',ls='dotted')
ax6.set_ylabel('AMR Level')

# Energy ratio
ax4.semilogy(x,er/et,'o',color='black',markerfacecolor='none')
ax4.set_xlabel('Distance (cm)')
ax4.set_ylabel('Radiative energy/Thermal energy')

fig.subplots_adjust(wspace=0.33)
fig.savefig('dirac.pdf',bbox_inches='tight')

# Check results against reference solution
visu_ramses.check_solution(data["data"],'dirac',tolerance={"velocity_x":1.5e-12})
