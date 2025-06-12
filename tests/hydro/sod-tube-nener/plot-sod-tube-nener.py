import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import visu_ramses

fig = plt.figure(figsize=(8, 12))
ax1 = plt.subplot(311)
ax2 = plt.subplot(312)
ax3 = plt.subplot(313)
ax5 = ax1.twinx()

# Load RAMSES output
data = visu_ramses.load_snapshot(2)

order  = data["data"]["x"].argsort()
x      = data["data"]["x"][order]
amrlev = data["data"]["level"][order]
rho    = data["data"]["density"][order]
u      = data["data"]["velocity_x"][order]
p      = data["data"]["pressure"][order]
nener1 = data["data"]["non_thermal_pressure_01"][order]
nener2 = data["data"]["non_thermal_pressure_02"][order]

# Density
ax1.plot(x,rho,'o',color='black',markerfacecolor='none')
ax1.set_xlabel('Distance (cm)')
ax1.set_ylabel('Density (g/cm3)')
ax5.plot(x,amrlev,color='black',ls='dotted')
ax5.set_ylabel('AMR Level')

# Velocity
ax2.plot(x,u,'o',color='black',markerfacecolor='none')
ax2.set_xlabel('Distance (cm)')
ax2.set_ylabel('Velocity (cm/s)')

# Pressure
ax3.plot(x,p,'o',color='green',markerfacecolor='none', label='P')
ax3.plot(x,nener1,'o',color='blue',markerfacecolor='none',label='P_nener1')
ax3.plot(x,nener2,'o',color='yellow',markerfacecolor='none',label='P_nener2')
ax3.plot(x,p+nener1+nener2,'o',color='black',markerfacecolor='none',label='P_tot')
ax3.set_xlabel('Distance (cm)')
ax3.set_ylabel('Pressure (g/cm/s2)')
ax3.legend()

# Read analytical solution
data_ana = np.loadtxt('../sod-tube/sod-tube-ana.dat')
x_ana    = data_ana[:,1]
rho_ana  = data_ana[:,3]
u_ana    = data_ana[:,2]
p_ana    = data_ana[:,4]
ax1.plot(x_ana,rho_ana,color='red')
ax2.plot(x_ana,u_ana,color='red')
ax3.plot(x_ana,p_ana,color='red')

fig.subplots_adjust(wspace=0.3)
fig.savefig('sod-tube-nener.pdf',bbox_inches='tight')

# Check results against reference solution
visu_ramses.check_solution(data["data"],'sod-tube-nener',overwrite=False)
