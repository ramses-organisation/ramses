import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import visu_ramses
from scipy.interpolate import griddata

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 8))

# Load RAMSES output
data = visu_ramses.load_snapshot(2)
units_density=2.3247434e-24
units_time = 2.5395079e15
units_length=3.0857e18
mu=1.4
MH = 1.6737236e-24 #g                      # hydrogen mass
rho    = data["data"]["density"]*units_density/(mu*MH)
p      = data["data"]["pressure"]*units_density*(units_length/units_time)**2

ax.scatter(rho, p, s=3, marker='s', color='grey', alpha=0.1, edgecolors='none')
ax.scatter([1*units_density/(mu*MH)], [319266*units_density*(units_length/units_time)**2], marker='*', color='black')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('density [cm-3]')
ax.set_ylabel('pressure [erg/cm3]')

fig.savefig('cooling-frig.pdf',bbox_inches='tight')

# Check results against reference solution
visu_ramses.check_solution(data["data"],'cooling-frig', threshold=1e-30, overwrite=False)
