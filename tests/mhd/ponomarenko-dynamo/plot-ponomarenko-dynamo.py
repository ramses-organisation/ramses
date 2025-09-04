import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import visu_ramses
import pyvista

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
b2 = 0.5*b2
lvlmax = max(lvl)

#fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(4.5,5), sharey=True, sharex=True)
#ax.hist(np.log10(b2),bins=np.linspace(-8,-1,30))
#plt.show()

logb2 = np.log10(b2)
filt=((logb2>-5.5)&(logb2<-3.4)&(lvl==lvlmax))
#At time t=200, one should rather use: 
#filt=((logb2>1.7)&(logb2<2.3)&(lvl==lvlmax))
plotter = pyvista.Plotter(shape=(1,1))
points = np.transpose([x[filt],y[filt],z[filt]])
cloud = pyvista.PolyData(points)
cloud['B2/2'] = logb2[filt]
print(cloud)
p = plotter.add_mesh(cloud, cmap="magma", point_size=7,
opacity='linear')
plotter.camera.roll += 90
#plotter.camera.elevation += 45
#plotter.camera.azimuth += 180
plotter.show_grid()
#plotter.show()
plotter.save_graphic('ponomarenko-dynamo.pdf')

# Check results against reference solution
visu_ramses.check_solution(data["data"], 'ponomarenko-dynamo') #, tolerance={"all":3.0e-06},overwrite=False)
