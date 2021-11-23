import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import visu_ramses

# Fundamental constants
MH = 1.6605390e-24 #g                      # hydrogen mass
KB = 1.3806490e-16 #cm^2 g s^-2 K^-1      # Boltzman constant

unit_d=2.3247434e-24 
unit_t=2.5395079e15
unit_l=3.0857e18

unit_v=unit_l/unit_t

fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(15, 10))

num_snaps = 6
for i in range(1, num_snaps+1):
    # Load RAMSES output
    data = visu_ramses.load_snapshot(i, read_grav=True)
    #x      = data["data"]["x"] - 500
    #y      = data["data"]["y"] - 500
    z      = data["data"]["z"] - 500
    vz      = data["data"]["velocity_z"]
    rho    = data["data"]["density"]
    mass   = rho * (1000./2**5)**3
    p      = data["data"]["pressure"]
    T      = (p / rho)*unit_v**2 * 1.4 * MH / KB
    Bx     = 0.5*(data["data"]["B_x_left"]+data["data"]["B_x_right"])
    By     = 0.5*(data["data"]["B_y_left"]+data["data"]["B_y_right"])
    Bz     = 0.5*(data["data"]["B_z_left"]+data["data"]["B_z_right"])
    a_z = data["data"]["a_z"]

    # density
    z_ana = np.linspace(-500, 500)
    rho_ana = 1.5 * np.exp(-z_ana**2 / (2.*150**2))
    ax[0][0].plot(z, rho, 'o', markerfacecolor='none', label='sim '+str(i))
    if i==1:
        ax[0][0].plot(z_ana, rho_ana, color='black', label='ana')
        ax[0][0].set_xlabel('disk height [pc]')
        ax[0][0].set_ylabel('density [c.u]')
    ax[0][0].legend()

    # velocity
    ax[0][1].plot(z, vz, 'o', markerfacecolor='none')
    if i==1:
        ax[0][1].set_xlabel('disk height [pc]')
        ax[0][1].set_ylabel('velocity [c.u.]')

    # grav acc
    ax[0][2].plot(z, a_z, 'o', markerfacecolor='none')


    # temperature
    #P_ana = rho_ana * (KB*8000/(1.4*MH))/unit_v**2
    ax[1][0].plot(z, T, 'o', markerfacecolor='none', label='sim')
    #ax[1][0].plot(z, p, 'o', color='red', markerfacecolor='none', label='sim')
    if i==1:
        ax[1][0].plot([-500, 500], [8000,8000], color='black', label='ana')
        #ax[1][0].plot(z_ana, P_ana, color='black', label='ana')
        ax[1][0].set_xlabel('disk height [pc]')
        ax[1][0].set_ylabel('T [K]')

    # magnetic field
    Bx_ana = 1 * 3.0 * np.exp(-z_ana**2/(2.*150**2))
    ax[1][1].plot(z, Bx, 'X', label='sim Bx')
    #ax[1][1].plot(z, By, '1', label='sim By')
    #ax[1][1].plot(z, Bz, 's', markerfacecolor='none', label='sim Bz')
    if i==1:
        ax[1][1].plot(z_ana, Bx_ana, color='black', label='ana')
        ax[1][1].set_xlabel('disk height [pc]')
        ax[1][1].set_ylabel('magnetic field [c.u.]')
        ax[1][1].legend()

    # phase diagram rho - T
    #hist, xedges, yedges = np.histogram2d(np.log10(rho*unit_d/(1.4*MH)), np.log10(T), weights=mass)
    #imag = ax[2].imshow(hist, interpolation='none', origin='lower',
    #                    extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
    #                    aspect=(xedges[-1]-xedges[0])/(yedges[-1]-yedges[0]))
    #ax[2].set_xlabel('log(density [cm-3])')
    #ax[2].set_ylabel('log(temperature [K])')

fig.savefig('strat_disk.pdf',bbox_inches='tight')

# Check results against reference solution
visu_ramses.check_solution(data["data"],'strat_disk', overwrite=True)
