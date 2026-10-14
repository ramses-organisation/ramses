import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import visu_ramses
from scipy.interpolate import griddata
from mpl_toolkits.axes_grid1 import make_axes_locatable

fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(8, 10))

# Parameters
nu_visco = 0.05
twopi= 6.2831853
amplitude = 0.01

# Load RAMSES output
data = visu_ramses.load_snapshot(2)
x      = data["data"]["x"]
vx     = data["data"]["velocity_x"]
rho    = data["data"]["density"]
y      = data["data"]["y"]
dx     = data["data"]["dx"]
k = twopi / data["data"]["boxlen"]

ymin = np.amin(y-0.5*dx)
ymax = np.amax(y+0.5*dx)

nx  = 2**7
dpy = (ymax-ymin)/float(nx)
ypx = np.linspace(ymin+0.5*dpy,ymax-0.5*dpy,nx)
points = np.transpose([x,y])

# Density
time = data["data"]["time"]

expected_amplitude = amplitude * np.exp(-nu_visco*k**2*time)
solution_vx = lambda y: expected_amplitude*np.sin(k*y)

profile_vx = griddata(y,vx,ypx,method='linear')
profile_dens = griddata(y,rho,ypx,method='linear')


solution = solution_vx(ypx)

ax[0].plot(ypx, solution, color="tab:red", lw=5, alpha=0.6, label="Analytical solution")
ax[0].plot(ypx, profile_vx, color="k", label="Data")
ax[0].set_xlabel("y")
ax[0].set_ylabel(r"$v_x$")
ax[0].legend()

ax[1].plot(ypx, profile_dens, color="k", label="density")
ax[1].set_xlabel("y")
ax[1].set_ylabel(r"$\rho$")

error = np.abs(profile_vx - solution) / expected_amplitude
data["data"]["error_wrt_analytical"] = error


fig.savefig('visco-shear-wave-decay.pdf',bbox_inches='tight')

# Check results against reference solution
visu_ramses.check_solution(data["data"],'visco-shear-wave-decay', tolerance={"all":1e-10}, overwrite=False)
