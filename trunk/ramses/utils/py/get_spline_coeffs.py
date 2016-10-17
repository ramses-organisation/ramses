#!/usr/bin/env python
"""Finds most massive halo in RAMSES cosmological zoom simulation,
fits its position through cosmic times and returns namelist values
for MOVIE_PARAMS for centres"""

import os
import sys
import glob
import warnings
import numpy as np
from matplotlib import pyplot as plt
from numpy.polynomial.polynomial import polyfit, polyval


warnings.filterwarnings('error')

if len(sys.argv) < 2:
    sys.exit("Not enough arguments given. Call with:\
              get_spline_coeffs.py <path> [min max]")

path = sys.argv[1]
print 'Fitting coefficients for: '+path

# find max output number
min_out = int(min(glob.glob(path+'/output_?????')).split('_')[-1])
max_out = int(max(glob.glob(path+'/output_?????')).split('_')[-1])

skip = []  # can specify which outputs to omit in the fit
num_cell_threshold = 25000  # may need adjusting; excludes too small halos
order = 3  # fit order - from 0 to 3

if len(sys.argv) >= 4:
    custom_min_out = int(sys.argv[2])
    custom_max_out = int(sys.argv[3])
    min_out = max(custom_min_out, min_out)
    max_out = min(custom_max_out, max_out)

# init
aexp = np.zeros(max_out-min_out+1)
x = np.zeros(max_out-min_out+1)
y = np.zeros(max_out-min_out+1)
z = np.zeros(max_out-min_out+1)
m = np.zeros(max_out-min_out+1)

# get positions of most massive halo from PHEW halo catalogues
for j in xrange(max_out-min_out+1):
    if (not os.path.exists(path+'/output_{j:05d}/halo_{j:05d}.txt00001'
                           .format(j=j+min_out)) or j in skip):
        x[j] = np.NaN
        y[j] = np.NaN
        z[j] = np.NaN
        m[j] = np.NaN
        continue
    else:
        max_mass = 0.0
        max_x = np.nan
        max_y = np.nan
        max_z = np.nan

        for halo_file in glob.glob(path+'/output_{j:05d}/halo_{j:05d}.txt*'
                                   .format(j=j+min_out)):
            try:
                hm, hx, hy, hz = np.genfromtxt(halo_file, skip_header=1,
                                               usecols=(1, 2, 3, 4)).T
            except UserWarning:  # catching empty files
                continue

            max_ncell = np.argmax(hm)
            try:
                if ((hm[max_ncell] > max_mass) and
                        (hm[max_ncell] > num_cell_threshold)):
                    max_mass = hm[max_ncell]
                    max_x = hx[max_ncell]
                    max_y = hy[max_ncell]
                    max_z = hz[max_ncell]
            except IndexError:  # catching files with only one value
                if (hm > max_mass) and (hm > num_cell_threshold):
                    max_mass = hm
                    max_x = hx
                    max_y = hy
                    max_z = hz

        x[j] = max_x
        y[j] = max_y
        z[j] = max_z
        m[j] = max_mass


# get aexps
for j in xrange(max_out-min_out+1):
    fp = open(path+'/output_{j:05d}/info_{j:05d}.txt'.format(j=j+min_out))
    for i, line in enumerate(fp):
        if i == 9:
            aexp[j] = float(line.split('=')[-1])
        elif i > 9:
            break
    fp.close()

# remove NaNs
x = x[np.logical_not(np.isnan(x))]
y = y[np.logical_not(np.isnan(y))]
z = z[np.logical_not(np.isnan(z))]
aexp = aexp[-len(x):]

# fit coefficients
coeffx = polyfit(aexp, x, order, full=True, w=aexp)[0]
coeffy = polyfit(aexp, y, order, full=True, w=aexp)[0]
coeffz = polyfit(aexp, z, order, full=True, w=aexp)[0]

# print result
print 'xcentre_frame='+','.join('{:6f}'.format(i) for i in coeffx)
print 'ycentre_frame='+','.join('{:6f}'.format(i) for i in coeffy)
print 'zcentre_frame='+','.join('{:6f}'.format(i) for i in coeffz)

# plotting for a check
if len(sys.argv) >= 3:
    plt.subplot(3, 1, 1)
    plt.plot(aexp, x, 'ro', ms=5)
    plt.plot(aexp, polyval(aexp, coeffx), label='x', lw=3)
    plt.subplot(3, 1, 2)
    plt.plot(aexp, y, 'go', ms=5)
    plt.plot(aexp, polyval(aexp, coeffy), label='x', lw=3)
    plt.subplot(3, 1, 3)
    plt.plot(aexp, z, 'bo', ms=5)
    plt.plot(aexp, polyval(aexp, coeffz), label='x', lw=3)

    plt.savefig('coeffs.pdf')
    plt.show()
