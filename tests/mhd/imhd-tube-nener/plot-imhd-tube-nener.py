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

# Load RAMSES output
data = visu_ramses.load_snapshot(2)

order  = data["data"]["x"].argsort()
x      = data["data"]["x"][order] - 1.5
amrlev = data["data"]["level"][order]
rho    = data["data"]["density"][order]
u      = data["data"]["velocity_x"][order]
p      = data["data"]["pressure"][order]
nener1 = data["data"]["non_thermal_pressure_01"][order]
nener2 = data["data"]["non_thermal_pressure_02"][order]
By     = 0.5*(data["data"]["B_y_left"][order] + data["data"]["B_y_right"][order])

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

# B field
ax4.plot(x,By,'o',color='black',markerfacecolor='none')
ax4.set_xlabel('Distance (cm)')
ax4.set_ylabel('By (G)')

# Read analytical solution
data_ana = np.loadtxt('../imhd-tube/imhd-tube-ana.dat')
x_ana    = data_ana[:,0]
rho_ana  = data_ana[:,1]
u_ana    = data_ana[:,2]
p_ana    = data_ana[:,8]
By_ana   = data_ana[:,6]
ax1.plot(x_ana,rho_ana,color='red')
ax2.plot(x_ana,u_ana,color='red')
ax3.plot(x_ana,p_ana,color='red')
ax4.plot(x_ana,By_ana,color='red')

fig.subplots_adjust(wspace=0.3)
fig.savefig('imhd-tube-nener.pdf',bbox_inches='tight')

# Check results against reference solution
visu_ramses.check_solution(data["data"],'imhd-tube-nener')
