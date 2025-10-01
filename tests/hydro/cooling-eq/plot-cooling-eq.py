import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import visu_ramses
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import glob

# Some constants
mH = 1.6605390e-24
kB = 1.3807e-16
nx=2 # Number of cells in each dimension

# Routine for showing maps, in one column, of density and temperature from a given snapshot
def show_maps(ax, isnap, icol):
    data = visu_ramses.load_snapshot(isnap)
    scale_d = data["data"]["unit_d"]
    scale_l = data["data"]["unit_l"]
    scale_t = data["data"]["unit_t"]
    scale_v = scale_l/scale_t
    scale_T2 = mH/kB * scale_v**2
    rho    = data["data"]["density"] * scale_d / 1.6605390e-24 * 0.76
    p      = data["data"]["pressure"] * scale_d * scale_v**2
    Tmu    = data["data"]["pressure"]/data["data"]["density"] * scale_T2
    time   = data["data"]["time"] * scale_t / 3.156e13

    nhmin=3e-2 ; nhmax=3e0
    Tmin=1e2 ; Tmax=1e7
    im_nH = ax[0,icol].imshow(np.array(rho).reshape((nx,nx)), origin="lower", aspect='equal', norm=colors.LogNorm(vmin=nhmin, vmax=nhmax), cmap='viridis')
    im_T = ax[1,icol].imshow(np.array(Tmu).reshape((nx,nx)), origin="lower", aspect='equal', norm=colors.LogNorm(vmin=Tmin, vmax=Tmax), cmap='coolwarm')
    for i in range(2):
      ax[i,icol].text(0.02, 0.98, '%d Myr'%(time), transform=ax[i,icol].transAxes, ha='left', va='top', fontsize=9, color='white')
    for a in ax.flat:
      a.set_xticks([])
      a.set_yticks([])
      for spine in a.spines.values():
        spine.set_visible(True)
    plt.colorbar(im_nH, ax=ax[0,icol], label=r'$n_{\rm H} \, {\rm [cm^{-3}]}$', orientation='horizontal', pad=0.05)
    plt.colorbar(im_T, ax=ax[1,icol], label=r'$T/\mu \, {\rm [K]}$', orientation='horizontal', pad=0.05)
# End of show_maps routine

# Find the number of snapshots
sSnaps = glob.glob('./output_*') # Snapshot strings
nSnaps = len(sSnaps)

with PdfPages('cooling-eq.pdf') as pdf:

    # Create a single figure
    fig = plt.figure(figsize=(9, 7))  # Adjust size as needed
    gs = gridspec.GridSpec(3, 4, height_ratios=[1, 1, 1], hspace=0.4)  # 2 rows for maps, 1 row for plot

    # 1) Show maps of density and temperature before and after (density does not change over time) ------------------------------
    ax_maps = np.array([[fig.add_subplot(gs[0, i]) for i in range(4)],
                        [fig.add_subplot(gs[1, i]) for i in range(4)]])
    show_maps(ax_maps, 1,0)
    show_maps(ax_maps, 59,1)
    show_maps(ax_maps, 74,2)
    show_maps(ax_maps, nSnaps,3)

    # 2) Plot the time-evolution of the temperature -----------------------------------------------------------------------------
    ax = fig.add_subplot(gs[2, :])  # Bottom row, spans all columns

    # Collect the time-evolution of Tmu in each cell from all the snapshots
    times = np.zeros(nSnaps)
    Ts = np.zeros((nx**2,nSnaps))
    for isnap in range(nSnaps):
        data = visu_ramses.load_snapshot(isnap+1)
        scale_l = data["data"]["unit_l"]
        scale_t = data["data"]["unit_t"]
        scale_v = scale_l/scale_t
        scale_T2 = mH/kB * scale_v**2
        Ts[:,isnap]    = data["data"]["pressure"]/data["data"]["density"] * scale_T2
        times[isnap] = data["data"]["time"] * scale_t / 3.156e13

    # Plot the data for the four cells in the simulation
    ax.plot(times, Ts[0,:], ls='--', label=r'$n_{\rm H} = 3\times10^{-2} \, {\rm cm^{-3}}$', color='black')
    ax.plot(times, Ts[1,:], ls='-', label=r'$n_{\rm H} = 3 \, {\rm cm^{-3}}$', color='black')
    ax.plot(times, Ts[2,:], ls='--', color='black')
    ax.plot(times, Ts[3,:], ls='-', color='black')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel('time [Myr]')
    ax.set_ylabel(r'$T/\mu \, {\rm [K]}$')
    ax.legend(loc='center left')

    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

visu_ramses.check_solution(data["data"],'cooling-eq',overwrite=False)
