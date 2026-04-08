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
scale_d = data["data"]["unit_d"]
scale_l = data["data"]["unit_l"]
scale_t = data["data"]["unit_t"]
time    = data["data"]["time"]*scale_t
ngr     = data["info"]["ngrp"]
ar      = 7.56591469318689378e-015
scale_v = scale_l/scale_t

order   = data["data"]["x"].argsort()
x       = data["data"]["x"][order]*scale_l + scale_v*time
amrlev  = data["data"]["level"][order]
rho     = data["data"]["density"][order]*scale_d
u       = data["data"]["velocity_x"][order]*scale_v
p       = data["data"]["pressure"][order]*scale_d*scale_v**2
T       = data["data"]["temperature"][order]
er = np.zeros([data["data"]["ncells"],ngr])
for i in range(ngr):
    er[:,i] = data["data"]["radiative_energy_0"+str(i+1)][order]
Tr      = (er*scale_d*scale_v**2/ar)**0.25
Tr_tot  = (er.sum(axis=1)*scale_d*scale_v**2/ar)**0.25

xmin = 2.0e+10
xmax = 7.0e+10

# Density
ax1.plot(x,rho,'o-',color='black',markerfacecolor='none')
ax1.set_xlabel('Distance (cm)')
ax1.set_ylabel('Density (g/cm3)')
ax1.set_xlim([xmin,xmax])

# Velocity
ax2.plot(x,u/1.0e5,'o-',color='black',markerfacecolor='none')
ax2.set_xlabel('Distance (cm)')
ax2.set_ylabel('Velocity (km/s)')
ax2.set_xlim([xmin,xmax])

# Pressure
ax3.plot(x,p,'o-',color='black',markerfacecolor='none')
ax3.set_xlabel('Distance (cm)')
ax3.set_ylabel('Pressure (g/cm/s2)')
ax3.set_xlim([xmin,xmax])

# Temperature
lw=2
ax4.plot(x,T,color='black',label='T',lw=lw)
ax4.set_xlabel('Distance (cm)')
ax4.set_ylabel('Temperature (K)')
for ig in range(ngr):
    ax4.plot(x,Tr[:,ig],label='Tr('+str(ig+1)+')',lw=lw)
ax4.plot(x,Tr_tot,label='Tr tot',lw=lw)
ax4.set_xlabel('Distance (cm)')
ax4.set_ylabel('Temperature (K)')
ax4.set_xlim([xmin,xmax])
ax4.legend(loc=1,fontsize=12)

ax5.plot(x,amrlev,color='black',ls='dotted',label='AMR Level')
ax5.set_ylabel('AMR Level')
ax5.legend(loc=7,fontsize=12)

fig.subplots_adjust(wspace=0.3)
fig.savefig('radiative-shock.pdf',bbox_inches='tight')

# Check results against reference solution
visu_ramses.check_solution(data["data"],'radiative-shock')
