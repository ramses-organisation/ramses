import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import visu_ramses

# Fundamental constants
MH = 1.6605390e-24 #g                      # hydrogen mass
KB = 1.3806490e-16 #cm^2 g s^-2 K^-1      # Boltzman constant
PC = 3.0856776e+18
MYR = 3.15576000e+13

unit_d=2.3247434e-24
unit_t=2.5395079e15
unit_l=3.0857e18

unit_v=unit_l/unit_t

fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10, 5))

# ICs
data = visu_ramses.load_snapshot(2, read_grav=True)
z   = data["data"]["y"] - 500
rho = data["data"]["density"]
phi = data["data"]["phi"]
a_z = data["data"]["a_y"]

# density
z_ana = np.linspace(-500, 500)
rho_ana = 0.66 * np.exp(-z_ana**2 / (2.*150**2))
ax[0].plot(z_ana, rho_ana, color='black', label='IC')
ax[0].scatter(z, rho, s=4, marker='o', label='sim', color='red')
ax[0].set_xlabel('disk height [pc]')
ax[0].set_ylabel('density [c.u]')
ax[0].legend()

# grav acc
a1 = 1.42e-3 * 1e3 * PC / MYR**2 / unit_l * unit_t**2
a2 = 5.49e-4 / MYR**2 * unit_t**2
z0 = 0.18e3 * PC / unit_l
ana_az = -a1 * z_ana / np.sqrt(z_ana**2+z0**2) - a2 *z_ana
ax[1].plot(z_ana, ana_az, color='black', label='analytical')
ax[1].scatter(z, a_z, s=4, marker='o', color='red', label='simulation')
ax[1].set_xlabel('disk height [pc]')
ax[1].set_ylabel('vertical accel. [c.u]')
ax[1].legend()

fig.savefig('ana-disk-potential.pdf',bbox_inches='tight')

# Check results against reference solution
visu_ramses.check_solution(data["data"],'ana-disk-potential', overwrite=False)
