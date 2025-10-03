import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
import visu_ramses
import pyvista
from scipy.interpolate import griddata

# Load RAMSES output
data = visu_ramses.load_snapshot(2)
x = data["data"]["x"]
y = data["data"]["y"]
z = data["data"]["z"]
lvl = data["data"]["level"]
dx = data["data"]["dx"]
b2 = 0.25*(data["data"]["B_x_left"]+data["data"]["B_x_right"])**2 + \
0.25*(data["data"]["B_y_left"]+data["data"]["B_y_right"])**2 + \
0.25*(data["data"]["B_z_left"]+data["data"]["B_z_right"])**2

plot_in_3d=False
plot_in_2d=True

if plot_in_3d:
    logb2 = np.log10(0.5*b2)
    lvlmax = max(lvl)
    filt=((logb2>-5.5)&(logb2<-4.)&(lvl==lvlmax))
    #At time t=80, one should rather use:
    #filt=((logb2>19.)&(logb2<21.5)&(lvl==lvlmax))

    plotter = pyvista.Plotter(shape=(1,1), off_screen=True)
    points = np.transpose([x[filt],y[filt],z[filt]])
    cloud = pyvista.PolyData(points)
    cloud['B2/2'] = logb2[filt]
    p = plotter.add_mesh(cloud, cmap="magma", point_size=7,opacity='linear')
    plotter.camera_position = [(3.2, 3.2, 3.2), (0.0, 0.0, 0.0), (0.0, 0.0, 1.0)]
    plotter.show_grid()
    plotter.save_graphic('abc-flow.pdf')
    plotter.close()


if plot_in_2d:
    fig, ax = plt.subplots(tight_layout=True)

    xmin = np.amin(x-0.5*dx)
    xmax = np.amax(x+0.5*dx)
    ymin = np.amin(y-0.5*dx)
    ymax = np.amax(y+0.5*dx)
    zmin = np.amin(z-0.5*dx)
    zmax = np.amax(z+0.5*dx)

    nx  = 64
    dpx = (xmax-xmin)/float(nx)
    dpy = (ymax-ymin)/float(nx)
    dpz = (zmax-zmin)/float(nx)
    xpx = np.linspace(xmin+0.5*dpx,xmax-0.5*dpx,nx)
    ypx = np.linspace(ymin+0.5*dpy,ymax-0.5*dpy,nx)
    zpx = np.linspace(zmin+0.5*dpz,zmax-0.5*dpz,nx)
    grid_x, grid_y, grid_z = np.meshgrid(xpx,ypx,zpx)
    points = np.transpose([x,y,z])
    z1 = griddata(points,0.5*b2,(grid_x,grid_y, grid_z),method='nearest')

    b2_proj = np.sum(z1, axis=0) #proj along y-axis

    im = ax.imshow(np.log10(b2_proj), origin="lower", aspect='equal', extent=[xmin, xmax, ymin, ymax])

    cb = []
    cb.append(plt.colorbar(im, ax=ax, label='log10(B^2/2)'))

    ax.set_xlabel('x')
    ax.set_ylabel('z')

    plt.savefig('abc-flow.pdf')

# Check results against reference solution
visu_ramses.check_solution(data["data"], 'abc-flow') #, tolerance={"all":3.0e-06},overwrite=False)
