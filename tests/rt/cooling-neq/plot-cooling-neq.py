import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import visu_ramses
#import matplotlib.colors as colors
import glob

# Some constants
mH = 1.6605390e-24
kB = 1.3807e-16
nx=2 # Number of cells in each dimension

# Find the number of snapshots
sSnaps = glob.glob('./output_*') # Snapshot strings
nSnaps = len(sSnaps)

# Collect the time-evolution of Tmu in each cell from all the snapshots
times = np.zeros(nSnaps)
TKs = np.zeros((nx**2,nSnaps))
xHIs = np.zeros((nx**2,nSnaps))
xHIIs = np.zeros((nx**2,nSnaps))
xHeIIs = np.zeros((nx**2,nSnaps))
xHeIIIs = np.zeros((nx**2,nSnaps))

for isnap in range(nSnaps):
    data = visu_ramses.load_snapshot(isnap+1)
    scale_l = data["data"]["unit_l"]
    scale_t = data["data"]["unit_t"]
    scale_v = scale_l/scale_t
    scale_T2 = mH/kB * scale_v**2
    times[isnap]     = data["data"]["time"] * scale_t / 3.156e13
    TKs[:,isnap]     = data["data"]["pressure"]/data["data"]["density"] * scale_T2
    xHIs[:,isnap]    = data["data"]["scalar_00"]
    xHIIs[:,isnap]   = data["data"]["scalar_01"]
    xHeIIs[:,isnap]  = data["data"]["scalar_02"]
    xHeIIIs[:,isnap] = data["data"]["scalar_03"]

TK_label = r'$T/\mu \, {\rm [K]}$'
xion_label = r'$x_{\rm ion}$'
TK_color = 'blue'
xh_color = 'red'
xhe_color = 'orange'
titles = [r'$n_{\rm H} = 3\times10^{-2} \, {\rm cm^{-3}}$',
              r'$n_{\rm H} = 3 \, {\rm cm^{-3}}$',
              r'$n_{\rm H} = 3\times10^{-2} \, {\rm cm^{-3}}$',
              r'$n_{\rm H} = 3 \, {\rm cm^{-3}}$']
Tmin=1e0 ; Tmax=1e8

# Create 4 subplots in one column, sharing the x-axis
fig, axes = plt.subplots(4, 1, sharex=True, figsize=(7, 6))

# Loop through each subplot
for i in range(4):
    ax = axes[i]
    ax_right = ax.twinx()
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylim((Tmin,Tmax))
    ax_right.set_ylim((-0.05,1.05))

    # Left y-axis is temperature
    ax.plot(times, TKs[i,:], color=TK_color)
    ax.set_ylabel(TK_label, color=TK_color)
    ax.tick_params(axis='y', labelcolor=TK_color)
    #ax.grid(True)

    # Right y-axis is for ionization fractions
    line_xhi, = ax_right.plot(
      times, 1.-xHIs[i,:]-xHIIs[i,:],    color=xh_color, ls='-', label=r'$x_{\rm H2}$')
    line_xhii, = ax_right.plot(times, xHIs[i,:],   color=xh_color, ls='--', label=r'$x_{\rm HI}$')
    line_xhei, = ax_right.plot(times, 1.-xHeIIs[i,:]-xHeIIIs[i,:],  color=xhe_color, ls=':', label=r'$x_{\rm HeI}$')
    line_xheii, = ax_right.plot(times, xHeIIs[i,:], color=xhe_color, ls='-.', label=r'$x_{\rm HeII}$')
    ax_right.set_ylabel(xion_label, color=xh_color)
    ax_right.tick_params(axis='y', labelcolor=xh_color)

    ax.set_title(titles[i], fontsize=8)

lines = [line_xhi, line_xhii, line_xhei, line_xheii]
labels = [line.get_label() for line in lines]
axes[3].legend(lines, labels, loc='center left',ncol=2)

axes[-1].set_xlabel('time [Myr]')

plt.tight_layout(pad=0.1)
plt.savefig("cooling-neq.pdf", bbox_inches='tight')

scal_tol =  6.0e-12
tolerance={"scalar_01":scal_tol, "scalar_02":scal_tol, "scalar_03":scal_tol}
visu_ramses.check_solution(data["data"],'cooling-neq',tolerance=tolerance, overwrite=False)
