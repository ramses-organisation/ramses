"""Plot results of sink accretion test"""
import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot as plt
from matplotlib import ticker
import visu_ramses
import numpy as np
from scipy.interpolate import griddata
from scipy.io import FortranFile

out_end = 2

# Fundamental constants
YR = 3.1556926e7 #s                        # 1 year
MH = 1.6737236e-24 #g                      # hydrogen mass
KB = 1.38064852e-16 #cm^2 g s^-2 K^-1      # Boltzman constant
# code units
unit_d=1.66e-24
unit_t=3.004683525921981e15
unit_l=3.08567758128200e+18
unit_v=unit_l/unit_t

data = visu_ramses.load_snapshot(out_end, read_rt=True)
#for key in data["stellars"].keys():
#    data["data"]["stellar_"+key] = data["stellars"][key]

x   = data["data"]["x"]
y   = data["data"]["y"]
z   = data["data"]["z"]
dx  = data["data"]["dx"]
data['data']['temperature'] = data["data"]["pressure"]/ data["data"]["density"]
xmin = np.amin(x-0.5*dx)
xmax = np.amax(x+0.5*dx)
ymin = np.amin(y-0.5*dx)
ymax = np.amax(y+0.5*dx)
zmin = np.amin(z-0.5*dx)
zmax = np.amax(z+0.5*dx)
ext=[xmin, xmax, ymin, ymax]

nx  = 2**7
dpx = (xmax-xmin)/float(nx)
dpy = (ymax-ymin)/float(nx)
dpz = (zmax-zmin)/float(nx)
xpx = np.linspace(xmin+0.5*dpx,xmax-0.5*dpx,nx)
ypx = np.linspace(ymin+0.5*dpy,ymax-0.5*dpy,nx)
zpx = np.linspace(zmin+0.5*dpz,zmax-0.5*dpz,nx)
grid_x, grid_y, grid_z = np.meshgrid(xpx,ypx,zpx)
points = np.transpose([x,y,z])

# Set up each plot
vars=["density", "temperature", "photon_flux_03", "scalar_00", "scalar_01", "scalar_02"]
units=[unit_d/MH,  unit_v**2 * 2.37 * MH /KB, unit_v, 1., 1., 1.]
labels=["log(density [H/cc])", "log(temperature [K])", "log(photon 3 flux [cm-2 s-1])",
             "log(xHII)", "log(xHeII)", "log(xHeIII)"]
axx=[0,0,0,1,1,1]
axy=[0,1,2,0,1,2]
vmin=[None, None, None, -10, -10, -10]
vmax=[None, None, None, 0, 0, 0]
fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(10, 10))

# Do the plots
for i in range(len(vars)):
    ps = data["data"][vars[i]]*units[i]
    gd = griddata(points,ps,(grid_x,grid_y,grid_z),method='nearest')
    im = ax[axx[i],axy[i]].imshow(np.log10(gd[:,:,int(2**7/2.)])
                ,origin="lower", aspect='equal', extent=[xmin, xmax, ymin, ymax],vmin=vmin[i], vmax=vmax[i])

    plt.colorbar(im, ax=ax[axx[i],axy[i]], label=labels[i])
    ax[axx[i],axy[i]].set_axis_off()
    ax[axx[i],axy[i]].scatter(data["sinks"]['x'],data["sinks"]['y'], s = 20, marker='x', color='red')


# Read binary movie frames, show and add to output data for testing ================================
movie_path = './movie1/'
vars = ['temp', 'Fp3', 'xHII']
dvars = ['movie_temp', 'movie_Fp3', 'movie_xHII']
labels = ['log(movie temp [K])', 'log(movie photon 3 flux [cm-2 s-1])', 'log(movie xHII)']
units=[2.37, unit_v, 1]
axx=[2,2,2]
axy=[0,1,2]
vmin=[None, None, -10]
vmax=[None, None, 0]

for i in range(len(vars)):

    fname = '%s_00002.map'%(vars[i])
    ffile = FortranFile(f"{movie_path}{fname}")
    [time, fdw, fdh, fdd] = ffile.read_reals('d')
    [frame_nx, frame_ny] = ffile.read_ints()
    mov_array = np.array(ffile.read_reals('f4'), dtype=np.float64) * units[i]
    data["data"][dvars[i]] = mov_array
    mov_map = mov_array.reshape(frame_nx,frame_ny)
    im = ax[axx[i],axy[i]].imshow(np.log10(mov_map), cmap='viridis', origin="lower", aspect='equal'
            ,vmin=vmin[i], vmax=vmax[i])
    ax[axx[i],axy[i]].axis('off')
    plt.colorbar(im, ax=ax[axx[i],axy[i]], label=labels[i])

fig.savefig('stellar-HII.pdf', bbox_inches="tight")

# Check if test results match reference, both for outputs and movie frames =========================
# Why is this so inaccurate on multiple cores?
red_tol = 3e-6
tolerance={"scalar_00":1e-7, "scalar_01":3.0e-7, "scalar_02":1.0e-7, "velocity_x":1.0e-8, "velocity_y":1.0e-8, "velocity_z":1.0e-8
  ,'movie_Fp3':red_tol, 'movie_xHII':red_tol
  ,'photon_flux_01':red_tol, 'photon_flux_01_x':red_tol, 'photon_flux_01_y':red_tol, 'photon_flux_01_z':red_tol
  ,'photon_flux_02':red_tol, 'photon_flux_02_x':red_tol, 'photon_flux_02_y':red_tol, 'photon_flux_02_z':red_tol
  ,'photon_flux_03':red_tol, 'photon_flux_03_x':red_tol, 'photon_flux_03_y':red_tol, 'photon_flux_03_z':red_tol}

visu_ramses.check_solution(data["data"],'stellar-HII',tolerance=tolerance,overwrite=False)
