"""Plot results of sink accretion test"""
import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot as plt
from matplotlib import ticker
import visu_ramses
import numpy as np

out_end = 17

# Fundamental constants
G = 6.67259e-8 #cm^3 g^-1 s^-2             # gravitational constant
YR = 3.1556926e7 #s                        # 1 year
MSUN = 1.989e33 #g                         # solar mass
MH = 1.6737236e-24 #g                      # hydrogen mass
KB = 1.38064852e-16 #cm^2 g s^-2 K^-1      # Boltzman constant
PC = 3.0857e18 #cm                         # 1 parsec
AU = 1.49597871e13 #cm                     # 1 astronomical unit

# code units
unit_d=1.66e-24
unit_t=3.004683525921981e15
unit_l=3.08567758128200e+18
unit_v=unit_l/unit_t
unit_m = unit_d * unit_l**3

m_sink = []
m_stellar = {0:[], 1:[], 2:[]}
t_sink = []
t_stellar = {0:[], 1:[], 2:[]}

for i in range(1,out_end+1):
    data = visu_ramses.load_snapshot(i)
    time = data["info"]["time"] * data["info"]["unit_t"] / YR / 1000.0
    if len(data["stellars"])>1:
        for s in range(data["stellars"]['nstellars']):
            m_stellar[s].append(data['stellars']['mstellar'][s]*unit_m/MSUN)
            t_stellar[s].append(time)
    m_sink.append(data['sinks']['msink']*unit_m/MSUN)
    t_sink.append(time)

data = visu_ramses.load_snapshot(out_end)
for key in data["stellars"].keys():
    data["data"]["stellar_"+key] = data["stellars"][key]
for key in ['nsinks','id','msink','dmfsink','x', 'y', 'z']:
    data["data"]["sink_"+key] = data["sinks"][key]

fig = plt.figure()
axis = fig.add_subplot(111)
axis.plot(t_sink, m_sink, marker='o', label='sink')
axis.plot([0, max(t_sink)], [200,200], color='grey', label='threshold')
axis.plot([0, max(t_sink)], [400,400], color='grey')
axis.plot([0, max(t_sink)], [600,600], color='grey')

for s in range(data["stellars"]['nstellars']):
    axis.plot(t_stellar[s], m_stellar[s], marker='*', label='stellar')
axis.set_xlabel('time [Myr]')
axis.set_ylabel('mass [Msun]')
axis.legend()

fig.savefig('stellar-spawn.pdf', bbox_inches="tight")

visu_ramses.check_solution(data["data"],'stellar-spawn',overwrite=False)
