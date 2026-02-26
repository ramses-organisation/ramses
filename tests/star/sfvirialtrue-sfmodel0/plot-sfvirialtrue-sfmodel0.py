import matplotlib as mpl
import numpy as np
from scipy.interpolate import griddata
# uncomment for command line use, comment for notebook
#mpl.use("Agg")
import matplotlib.pyplot as plt
import visu_ramses
from matplotlib.colors import LogNorm

# NOTE: this test is designed to run on 2 proc
# on a different number of proc, it will fail unless tolerance are significantly increased
# this is because of the random number generator use in the star formation routine


fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(12, 8))

# Load RAMSES output
data = visu_ramses.load_snapshot(2,read_hydro=True)
xp = data["particle"]["position_x"]
yp = data["particle"]["position_y"]
zp = data["particle"]["position_z"]
mp = data["particle"]["mass"]
x      = data["data"]["x"]
y      = data["data"]["y"]
z      = data["data"]["z"]
dx     = data["data"]["dx"]
rho    = data["data"]["density"]
p      = data["data"]["pressure"]
l      = data["data"]["level"]
metals = data["data"]["metallicity"]

temp=p/rho
# find highest density cell as most likely SN site
itmax=np.argmax(rho)
xtmax=x[itmax]
ytmax=y[itmax]
ztmax=z[itmax]

# override => use center of domain
xtmax=25.
ytmax=25.
ztmax=25.

zoomdx=5.

xmin = xtmax-zoomdx
xmax = xtmax+zoomdx
ymin = ytmax-zoomdx
ymax = ytmax+zoomdx
zmin = ztmax-zoomdx
zmax = ztmax+zoomdx



# make a subcube of half-size zoomdx around max rho
nx  = 2**7
dpx = (xmax-xmin)/float(nx)
dpy = (ymax-ymin)/float(nx)
dpz = (zmax-zmin)/float(nx)
xpx = np.linspace(xmin+0.5*dpx,xmax-0.5*dpx,nx)
ypx = np.linspace(ymin+0.5*dpy,ymax-0.5*dpy,nx)
zpx = np.linspace(zmin+0.5*dpz,zmax-0.5*dpz,nx)
grid_x, grid_y, grid_z = np.meshgrid(xpx,ypx,zpx)
points = np.transpose([x,y,z])
z1 = griddata(points,rho,(grid_x,grid_y, grid_z),method='nearest')
z5 = griddata(points,p,(grid_x,grid_y, grid_z),method='nearest')
zl = griddata(points,l,(grid_x,grid_y, grid_z),method='nearest')
Z = griddata(points,metals,(grid_x,grid_y, grid_z),method='nearest')

temp=z5/z1

rho_proj2=z1[:,int(nx/2.),:]
rho_proj1=z1[int(nx/2.),:,:]
rho_proj3=z1[:,:,int(nx/2.)]

Z_proj2=Z[:,int(nx/2.),:]
Z_proj1=Z[int(nx/2.),:,:]
Z_proj3=Z[:,:,int(nx/2.)]



T_proj2=temp[:,int(nx/2.),:]
T_proj1=temp[int(nx/2.),:,:]
T_proj3=temp[:,:,int(nx/2.)]

im1 = ax[0,0].imshow(np.transpose(np.log10(rho_proj1)), origin="lower",aspect='equal', extent=[ymin, ymax, zmin, zmax],interpolation='none')
im2 = ax[0,1].imshow(np.transpose(np.log10(rho_proj2)), origin="lower", aspect='equal', extent=[xmin, xmax, zmin, zmax],interpolation='none')
im3 = ax[0,2].imshow(np.transpose(np.log10(rho_proj3)), origin="lower", aspect='equal', extent=[xmin, xmax, ymin, ymax],interpolation='none')
im4 = ax[1,0].imshow(np.transpose(np.log10(T_proj1)), origin="lower", aspect='equal', extent=[ymin, ymax, zmin, zmax],interpolation='none')
im5 = ax[1,1].imshow(np.transpose(np.log10(T_proj2)), origin="lower", aspect='equal', extent=[xmin, xmax, zmin, zmax],interpolation='none')
im6 = ax[1,2].imshow(np.transpose(np.log10(T_proj3)), origin="lower", aspect='equal', extent=[xmin, xmax, ymin, ymax],interpolation='none')
im7 = ax[2,0].imshow(np.transpose(np.log10(Z_proj1+1.e-10)), origin="lower", aspect='equal', extent=[ymin, ymax, zmin, zmax],interpolation='none')
im8 = ax[2,1].imshow(np.transpose(np.log10(Z_proj2+1.e-10)), origin="lower", aspect='equal', extent=[xmin, xmax, zmin, zmax],interpolation='none')
im9 = ax[2,2].imshow(np.transpose(np.log10(Z_proj3+1.e-10)), origin="lower", aspect='equal', extent=[xmin, xmax, ymin, ymax],interpolation='none')


cb = []
cb.append(plt.colorbar(im1, ax=ax[0,0], label='Density slice'))
cb.append(plt.colorbar(im2, ax=ax[0,1], label='Density slice'))
cb.append(plt.colorbar(im3, ax=ax[0,2], label='Density slice'))
cb.append(plt.colorbar(im4, ax=ax[1,0], label='Temp slice'))
cb.append(plt.colorbar(im5, ax=ax[1,1], label='Temp slice'))
cb.append(plt.colorbar(im6, ax=ax[1,2], label='Temp slice'))
cb.append(plt.colorbar(im7, ax=ax[2,0], label='Abs metals slice'))
cb.append(plt.colorbar(im8, ax=ax[2,1], label='Abs metals slice'))
cb.append(plt.colorbar(im9, ax=ax[2,2], label='Abs metals slice'))

for i in [0,1,2]:
    ax[i,0].set_xlabel('y')
    ax[i,0].set_ylabel('z')
    ax[i,1].set_xlabel('x')
    ax[i,1].set_ylabel('z')
    ax[i,2].set_xlabel('x')
    ax[i,2].set_ylabel('y')


for c in cb:
    c.ax.yaxis.set_label_coords(-1.1, 0.5)

fig.savefig('sfvirialtrue-sfmodel0.pdf',bbox_inches='tight')
fig.savefig('sfvirialtrue-sfmodel0.png') # => produces a sharper figure on mac

#to_check={}
#to_check = data["particle"]
#to_check["pressure"]=z5.flatten()  # data["data"]["pressure"]

to_check = data["data"]
for key in data["particle"].keys():
    to_check['particle '+key] = data["particle"][key]

# default test nproc is 2
nproc=2

# tweak tolerance to allow for 2p vs 4p/8p/12p deviations:
tolerance={}
if(nproc!=2):
    tolerance={"density":1.0e-4,\
               "pressure":1.e-4,\
               "temperature":1.e-4,\
               "metallicity":1.e-4,\
}

# and this is the "default tolerance"
if(nproc==2):
     tolerance={"all":1.0e-13,\
               }
visu_ramses.check_solution(to_check, 'sfvirialtrue-sfmodel0',tolerance=tolerance)#,overwrite=True)
