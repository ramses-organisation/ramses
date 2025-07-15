import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import visu_ramses
from matplotlib.colors import LogNorm
import glob
import os
import filecmp
import numpy as np


def concatenate_clump_files():
    # Find all files matching the pattern
    pattern = "output_00002/clump_00002.txt00*"
    files = glob.glob(pattern)

    if not files:
        print(f"No files found matching pattern: {pattern}")
        return

    # Sort files to ensure consistent ordering
    files.sort()
    print(f"Found {len(files)} files: {files}")

    output_file = "clump_all.txt"
    header_written = False

    with open(output_file, 'w') as outfile:
        for i, filename in enumerate(files):
            print(f"Processing {filename}...")

            with open(filename, 'r') as infile:
                lines = infile.readlines()

                # Skip empty files
                if not lines:
                    continue

                # Write header only from the first file
                if not header_written and lines:
                    outfile.write(lines[0])  # Write header
                    header_written = True

                # Write data lines (skip header for all files)
                for line in lines[1:]:
                    outfile.write(line)

    print(f"Successfully concatenated {len(files)} files into {output_file}")

    # Show summary
    with open(output_file, 'r') as f:
        total_lines = sum(1 for _ in f)
    print(f"Output file contains {total_lines} lines (including 1 header)")

# Execute the map part of the test

fig = plt.figure(figsize=(12, 3.75))
axes = fig.subplots(nrows=1, ncols=3)

# Load RAMSES output
data = visu_ramses.load_snapshot(2,read_hydro=False)
xp = data["particle"]["position_x"]
yp = data["particle"]["position_y"]
zp = data["particle"]["position_z"]
mp = data["particle"]["mass"]

im = axes[0].hist2d(xp,yp,weights=mp,bins=128,range=[[0, 1], [0, 1]],norm=LogNorm(vmin=8e-6,vmax=8e-4),cmap='bone',edgecolor='face')
im = axes[1].hist2d(xp,zp,weights=mp,bins=128,range=[[0, 1], [0, 1]],norm=LogNorm(vmin=8e-6,vmax=8e-4),cmap='bone',edgecolor='face')
im = axes[2].hist2d(yp,zp,weights=mp,bins=128,range=[[0, 1], [0, 1]],norm=LogNorm(vmin=8e-6,vmax=8e-4),cmap='bone',edgecolor='face')
#plt.colorbar(im[3], ax=axes[2])
for ax in axes:
    ax.axis('equal')
    ax.set_xlim([0,1])
    ax.set_ylim([0,1])
axes[0].set_xlabel('x')
axes[0].set_ylabel('y')
axes[1].set_xlabel('x')
axes[1].set_ylabel('z')
axes[2].set_xlabel('y')
axes[2].set_ylabel('z')

fig.savefig("cosmo.pdf", bbox_inches="tight")

to_check = data["particle"]
to_check['time'] = data["data"]["time"]

#Now  add the clumpfinder output to the solution to be checked
# First concatenate the clumpfinder output per proc into 1 file
# Execute concatenate:
concatenate_clump_files()

# Now add the clump data to the solution to be checked.
# Read the file using numpy.loadtxt, skipping the header row
# Also skip indexes since they may not be consistent
clumps = np.loadtxt("clump_all.txt", skiprows=1, usecols=[2,4,5,6,7,8,9,10,11,12])

clumps_ref = np.loadtxt("clump-ref.txt", skiprows=1, usecols=[2,4,5,6,7,8,9,10,11,12])

print(clumps)

# Sort by ncell
sort = np.argsort(clumps[:,2])
clump_sorted = clumps[sort,:]
sort_ref = np.argsort(clumps_ref[:,2])
clumps_ref = clumps_ref[sort_ref,:]

# Compare all properties
all_same = np.all(clump_sorted - clumps_ref == 0)

if not(all_same):
    print("ERROR: Clumps are differents. Check them manually")
else:
    print("Cosmo tests: Clumps are the same!")
# Extract the columns we need for the overall test
# Based on file format:
# Column indices: 0=lev,  1=ncell, 2=peak_x, 3=peak_y, 4=peak_z, 5=rho-, 6=rho+, 7=rho_av, 8=mass_cl, 9=relevance
ncell_output_array = clumps[:, 1]      # ncell column
mass_cl_output_array = clumps[:, 8]   # mass_cl column

print("ncell values:", ncell_output_array)
print("mass_cl values:", mass_cl_output_array)

# Add the clump ref solution to the to_check dict.
to_check["ncell_clumpfinder"]=ncell_output_array
to_check["mass_cl_clumpfinder"]=mass_cl_output_array

# then run the check_solution method
visu_ramses.check_solution(to_check, 'cosmo',overwrite=False) # True if you want to overwrite the ref solution
