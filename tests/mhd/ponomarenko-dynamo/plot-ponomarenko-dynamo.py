import numpy as np
import visu_ramses
import pyvista

# Load RAMSES output
data = visu_ramses.load_snapshot(2)
x = data["data"]["x"]
y = data["data"]["y"]
z = data["data"]["z"]
lvl = data["data"]["level"]
b2 = 0.25*(data["data"]["B_x_left"]+data["data"]["B_x_right"])**2 + \
0.25*(data["data"]["B_y_left"]+data["data"]["B_y_right"])**2 + \
0.25*(data["data"]["B_z_left"]+data["data"]["B_z_right"])**2
logb2 = np.log10(0.5*b2)
lvlmax = max(lvl)

filt=((logb2>-5.5)&(logb2<-3.4)&(lvl==lvlmax))
#At time t=200, one should rather use:
#filt=((logb2>1.7)&(logb2<2.3)&(lvl==lvlmax))
plotter = pyvista.Plotter(shape=(1,1), off_screen=True)
points = np.transpose([x[filt],y[filt],z[filt]])
cloud = pyvista.PolyData(points)
cloud['B2/2'] = logb2[filt]
p = plotter.add_mesh(cloud, cmap="magma", point_size=7, opacity='linear')
plotter.camera.roll += 90
plotter.show_grid()
plotter.save_graphic('ponomarenko-dynamo.pdf')
plotter.close()

# Check results against reference solution
visu_ramses.check_solution(data["data"], 'ponomarenko-dynamo') #, tolerance={"all":3.0e-06},overwrite=False)
