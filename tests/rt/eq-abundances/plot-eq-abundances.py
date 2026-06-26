import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import visu_ramses

data = np.loadtxt('eq-abundances.txt', skiprows=0)
T2     = data[:,0]
xH2    = data[:,1]
xHI    = data[:,2]
xHII   = data[:,3]
xHeI   = data[:,4]
xHeII  = data[:,5]
xHeIII = data[:,6]

plt.plot(T2, xH2, label=r'$x_{\mathrm{H2}}$', color='blue')
plt.plot(T2, xHI, label=r'$x_{\mathrm{HI}}$', color='cornflowerblue')
plt.plot(T2, xHII, label=r'$x_{\mathrm{HII}}$', color='cyan')
plt.plot(T2, xHeI, label=r'$x_{\mathrm{HeI}}$', color='orange',ls='--')
plt.plot(T2, xHeII, label=r'$x_{\mathrm{HeII}}$', color='red',ls='--')
plt.plot(T2, xHeIII, label=r'$x_{\mathrm{HeIII}}$', color='firebrick',ls='--')
plt.xscale('log')
plt.xlabel(r'${\rm T/\mu \ [k]}$')
plt.ylabel('fraction')
plt.grid(True)
plt.legend()
plt.show()

plt.tight_layout(pad=0.1)
plt.savefig("eq-abundances.pdf", bbox_inches='tight')

data = {"data": {}}

data["data"]['T2']     = T2
data["data"]['xH2']    = xH2
data["data"]['xHI']    = xHI
data["data"]['xHII']   = xHII
data["data"]['xHeI']   = xHeI
data["data"]['xHeII']  = xHeII
data["data"]['xHeIII'] = xHeIII

visu_ramses.check_solution(data["data"],'eq-abundances',overwrite=False)
