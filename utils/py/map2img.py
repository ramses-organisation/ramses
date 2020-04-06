#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from scipy.io import FortranFile
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("file", help="enter filename dens.map")
parser.add_argument("--log", help="plot log variable",action="store_true")
parser.add_argument("--out", help="output a png image")
args = parser.parse_args()
print("Reading "+args.file)

# path the the file
path_to_output = args.file

# read image data
with FortranFile(path_to_output, 'r') as f:
    t, dx, dy, dz = f.read_reals('f8')
    nx, ny = f.read_ints('i')
    dat = f.read_reals('f4')

print(nx,ny)
# reshape the output
dat = np.array(dat)
dat = dat.reshape(nx, ny)
dat = np.transpose(dat)
# plot the map
my_dpi = 96
fig, ax = plt.subplots(figsize=(512/my_dpi, 512/my_dpi), dpi=my_dpi)

if args.log:
    dat=np.log10(dat)
    
ax.imshow(dat[:, :].T, interpolation='nearest', origin='lower')
ax.set_xlabel("nx")
ax.set_ylabel("ny")
if args.out:
    plt.savefig(args.out)
plt.show()


