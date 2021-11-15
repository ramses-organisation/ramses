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
ax5 = ax4.twinx()

# Load RAMSES output
data   = visu_ramses.load_snapshot(2)
time   = data["data"]["time"]
ngr    = data["info"]["ngrp"]
ar     = 7.56591469318689378e-015
xshift = 1013.67

order   = data["data"]["x"].argsort()
x       = data["data"]["x"][order] - xshift
amrlev  = data["data"]["level"][order]
rho     = data["data"]["density"][order]
u       = data["data"]["velocity_x"][order]
p       = data["data"]["pressure"][order]
T       = data["data"]["temperature"][order]
Tr      = (data["data"]["radiative_energy_01"][order]/ar)**0.25

# Read analytical solution
data_ana = np.loadtxt('rshock-mach2-ana.dat')

xmin = -1000.0
xmax =  1000.0
lw=2

# Density
ax1.plot(x,rho,'o',color='black',markerfacecolor='none',label='simulation')
ax1.plot(data_ana[:,0],data_ana[:,1]*5.459E-013,color='r',label='analytical',lw=lw)
ax1.set_xlabel('Distance (cm)')
ax1.set_ylabel('Density (g/cm3)')
ax1.set_xlim([xmin,xmax])
ax1.legend(loc=2,fontsize=12)

# Velocity
ax2.plot(x,u/1.0e5,'o',color='black',markerfacecolor='none')
ax2.set_xlabel('Distance (cm)')
ax2.set_ylabel('Velocity (km/s)')
ax2.set_xlim([xmin,xmax])

# Pressure
ax3.plot(x,p,'o',color='black',markerfacecolor='none')
ax3.set_xlabel('Distance (cm)')
ax3.set_ylabel('Pressure (g/cm/s2)')
ax3.set_xlim([xmin,xmax])

# Temperature
ax4.plot(x,T,'o',markeredgecolor='b',lw=lw,markerfacecolor='none',label='T_simu')
ax4.plot(data_ana[:,2],data_ana[:,3]*100.0,color='orange',lw=lw,label='T_ana')
ax4.plot(x,Tr,'s',markeredgecolor='r',lw=lw,markerfacecolor='none',label='Tr_simu')
ax4.plot(data_ana[:,4],data_ana[:,5]*100.0,color='cyan',lw=lw,label='Tr_ana')
ax4.set_xlabel('Distance (cm)')
ax4.set_ylabel('Temperature (K)')
ax4.set_xlim([xmin,xmax])
ax4.legend(loc=2,fontsize=12)

ax5.plot(x,amrlev,color='grey',ls='dotted',label='AMR Level')
ax5.set_ylabel('AMR Level')
ax5.legend(loc=7,fontsize=12)

#fig.subplots_adjust(wspace=0.3)
fig.savefig('rshock-mach2.pdf',bbox_inches='tight')

# Check results against reference solution
visu_ramses.check_solution(data["data"],'rshock-mach2')
