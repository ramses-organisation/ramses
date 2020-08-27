"""Plot results of sink accretion test"""
import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot as plt
from matplotlib import ticker
import visu_ramses
import numpy as np



def get_accretion_rate(nout=0):

    data = visu_ramses.load_snapshot(nout)
    yr = 86400.0 * 365.25
    time = data["info"]["time"] * data["info"]["unit_t"] / yr / 1000.0

    u_m = data["info"]["unit_d"] * data["info"]["unit_l"]**3
    msun = 1.98847e33
    rate = data["sinks"]["acc_rate"] * u_m * yr / msun

    return time, rate


times = []
rates = []

for i in range(2, 16):
    t, r = get_accretion_rate(i)
    times.append(t)
    rates.append(r)


fig = plt.figure()
axis = fig.add_subplot(111)
axis.semilogy(times, rates, '-+')
axis.set_xlabel('Time [kyr]')
axis.set_ylabel(r'$\dot{M}_{\mathrm{Bondi}}~[M_{\odot}/\mathrm{yr}]$')
fig.savefig('smbh-bondi.pdf', bbox_inches="tight")

data = visu_ramses.load_snapshot(15)
for key in data["sinks"].keys():
    data["data"]["sink_"+key] = data["sinks"][key]
visu_ramses.check_solution(data["data"],'smbh-bondi')
