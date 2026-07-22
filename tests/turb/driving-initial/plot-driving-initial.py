# REMARK: it is normal the reference value for density is negative. This is the sum of the log of the cell densities.
#
# Driving on top of a turbulent initial velocity field: the only test with both
# initial_turb and turb switched on, covering the interaction between the two.
# The box starts as broadband turbulence (initial_turb_spectrum='power_law') and
# is pushed by a static forcing on the single mode k=(1,0,0).
#
# Rows:
#   1. initial velocity (output 1, the initial conditions)
#   2. final velocity
#   3. final density
#   4. driving rms, pinned exactly on target (a static field forces
#      turb_exact_rms) and undisturbed by the initial field, plus the
#      kinetic energy
#
# The velocity rows show the three signed components on a slice through the
# middle of the box, with a diverging colour scale shared between the two rows so
# they can be compared directly. Slices rather than projections, because summing
# a signed quantity along a line of sight largely cancels it out.
#
# The driving is fully compressive (comp_frac=1), so the Helmholtz projection
# leaves only the component along k: vx picks up the mode while vy and vz are
# untouched and simply decay from their initial values. With comp_frac<1 the
# solenoidal part would drive vy and vz as well, still varying along x.
#
# The density row keeps the usual projections along x, y and z. The mode is
# invisible in the projection along x, since it is coherent in that direction.

import matplotlib as mpl
mpl.use('Agg')
import os
import numpy as np
import matplotlib.pyplot as plt
import visu_ramses
import turb_log
from scipy.interpolate import griddata

# Number of dimensions this test is built with (see config.txt)
NDIM = 3
TESTNAME = 'turb/driving-initial'
NPIX = 2**6


def resample(data, values):
    """Interpolate a cell quantity onto a regular NPIX^3 grid."""
    x, y, z, dx = data["x"], data["y"], data["z"], data["dx"]
    bounds = [(np.amin(c-0.5*dx), np.amax(c+0.5*dx)) for c in (x, y, z)]
    axes = []
    for lo, hi in bounds:
        d = (hi-lo)/float(NPIX)
        axes.append(np.linspace(lo+0.5*d, hi-0.5*d, NPIX))
    grid_x, grid_y, grid_z = np.meshgrid(*axes)
    # meshgrid defaults to 'xy' indexing, so axis 0 is y, axis 1 is x, axis 2 is z
    cube = griddata(np.transpose([x, y, z]), values,
                    (grid_x, grid_y, grid_z), method='nearest')
    return cube, [bounds[0][0], bounds[0][1], bounds[1][0], bounds[1][1]]


def velocity_slices(data):
    """Mid-box slice of each velocity component."""
    out = []
    for comp in ("velocity_x", "velocity_y", "velocity_z"):
        cube, extent = resample(data, data[comp])
        out.append(cube[:, :, NPIX//2])
    return out, extent


def projections(data, values):
    """Projection of a quantity along x, y and z."""
    cube, extent = resample(data, values)
    return [np.sum(cube, axis=1), np.sum(cube, axis=0), np.sum(cube, axis=2)], extent


# Load the initial conditions and the final state
initial = visu_ramses.load_snapshot(1)["data"]
data = visu_ramses.load_snapshot(2)
final = data["data"]

v_ini, extent = velocity_slices(initial)
v_fin, _ = velocity_slices(final)
rho_fin, _ = projections(final, final["density"])

# Three rows of maps, plus a full-width panel for the history
fig = plt.figure(figsize=(9.5, 11))
gs = fig.add_gridspec(nrows=4, ncols=3, hspace=0.04, wspace=0.04,
                      height_ratios=[1, 1, 1, 0.85])

# Diverging scale centred on zero, shared by both velocity rows
vmax = max(np.amax(np.abs(p)) for p in v_ini + v_fin)

# The velocity rows show components on a slice, the density row shows
# projections, so each panel is labelled in place rather than on the axes
comps = ('$v_x$', '$v_y$', '$v_z$')
projs = (r'$\int$dx', r'$\int$dy', r'$\int$dz')

rows = ((0, v_ini, 'Initial velocity', 'RdBu_r', -vmax, vmax, comps),
        (1, v_fin, 'Final velocity', 'RdBu_r', -vmax, vmax, comps),
        (2, rho_fin, 'Final density', 'viridis', None, None, projs))

for irow, panels, label, cmap, lo, hi, tags in rows:
    axes = [fig.add_subplot(gs[irow, col]) for col in range(3)]
    for ax, panel, tag in zip(axes, panels, tags):
        im = ax.imshow(panel, origin="lower", aspect='equal', extent=extent,
                       cmap=cmap, vmin=lo, vmax=hi)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.text(0.04, 0.96, tag, transform=ax.transAxes, va='top', ha='left',
                fontsize='small',
                bbox=dict(boxstyle='square,pad=0.15', fc='white', ec='none',
                          alpha=0.7))
    fig.colorbar(im, ax=axes, label=label, fraction=0.025, pad=0.01)

ax_rms = fig.add_subplot(gs[3, :])

# Turbulent rms and kinetic energy history. These are only diagnostics: they are
# not compared against the reference solution, and a failure to read the log must
# not break the test.
try:
    hist = turb_log.read_turb_history(os.path.join('..', '..', 'test_suite.log'), TESTNAME)
    target = turb_log.read_target_rms('driving-initial.nml')
except Exception as e:
    hist, target = turb_log.read_turb_history('', TESTNAME), None
    print("Could not read turbulence history: %s" % e)

turb_log.plot_history(ax_rms, hist, target_rms=target, ndim=NDIM)

fig.savefig('driving-initial.pdf',bbox_inches='tight')

# Check results against reference solution
visu_ramses.check_solution(data["data"],'driving-initial', threshold=1e-30, overwrite=False)
