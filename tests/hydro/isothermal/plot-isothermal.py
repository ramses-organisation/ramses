import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import visu_ramses

# Fundamental constants
G = 6.67259e-8 #cm^3 g^-1 s^-2             # gravitational constant
YR = 3.1556926e7 #s                        # 1 year
MSUN = 1.989e33 #g                         # solar mass
MH = 1.6737236e-24 #g                      # hydrogen mass
KB = 1.38064852e-16 #cm^2 g s^-2 K^-1      # Boltzman constant
PC = 3.0857e18 #cm                         # 1 parsec
AU = 1.49597871e13 #cm                     # 1 astronomical unit

unit_d=1.50492957435e-20
unit_t=3.1556926e13
unit_l=3.0857e18
unit_v=unit_l/unit_t

fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(12, 4))

# Load RAMSES output
data = visu_ramses.load_snapshot(2)
x      = (data["data"]["x"] - 0.5) * unit_l / AU
rho    = data["data"]["density"] * unit_d
p      = data["data"]["pressure"] * unit_d * unit_l**2 / unit_t**2

ax[0].plot(x,rho,'o',color='black',markerfacecolor='none')
ax[0].set_xlabel('Distance (AU)')
ax[0].set_ylabel('Density (g/cm3)')
ax[0].set_xscale('symlog', linthresh=1.*PC/AU/2**16)
ax[0].set_yscale('log')

ax[1].plot(rho,p,'o',color='black',markerfacecolor='none')
ax[1].set_xlabel('Density (g/cm3)')
ax[1].set_ylabel('Pressure (g/cm/s2)')
ax[1].set_xscale('log')
ax[1].set_yscale('log')

cs = np.sqrt(p/rho)
cs2 = p/rho
T = cs2 * 2.37 * MH /KB

ax[2].plot(rho,T,'o',color='black',markerfacecolor='none')
ax[2].plot(rho, np.full(len(rho),10), color='red')
ax[2].set_xlabel('Density (g/cm3)')
ax[2].set_ylabel('Temperature (K)')
#ax[2].set_ylabel('cs (cm/s)')
ax[2].set_xscale('log')
ax[2].set_yscale('log')
#ax[2].set_ylim(1.7e4,2e4)
ax[2].set_ylim(1,max(T)*10)

fig.savefig('isothermal.pdf',bbox_inches='tight')

# Check results against reference solution
tolerance={"density":5.0e-12, "pressure":1.0e-11, "time":2.0e-12, "velocity_x":2.0e-12}
visu_ramses.check_solution(data["data"],'isothermal', tolerance=tolerance, overwrite=False)
