import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import visu_ramses
from matplotlib.colors import LogNorm

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

visu_ramses.check_solution(to_check, 'cosmo')

# AND NOW CHECK THE CLUMPFINDER OUTPUT

# First concatenate the clumpfinder output per proc into 1 file
import glob
import os

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

if __name__ == "__main__":
    concatenate_clump_files()


# Now compare the clump_all with the reference file

import filecmp

def simple_diff():
    file1 = "clump_all.txt"
    file2 = "clump-ref.txt"
    
    # Quick check if files are identical
    if filecmp.cmp(file1, file2):
        print("OK! Clump Files are identical => Passed")
        return
    
    print("ERROR clump Files differ from reference clump-ref.dat")
    
    # Show first differing line
    with open(file1) as f1, open(file2) as f2:
        for i, (line1, line2) in enumerate(zip(f1, f2), 1):
            if line1 != line2:
                print(f"First difference at line {i}:")
                print(f"  {file1}: {line1.rstrip()}")
                print(f"  {file2}: {line2.rstrip()}")
                break
        else:
            # Files have different lengths
            print("Files have different lengths")

if __name__ == "__main__":
    simple_diff()


import pandas as pd
import numpy as np

df = pd.read_csv("clump_all.txt", sep='\s+')
#df_ref = pd.read_csv("clump-ref.txt", sep='\s+')

ncell_output_array = df['ncell'].to_numpy()
mass_cl_output_array = df['mass_cl'].to_numpy()

#ncell_ref_array=df_ref['ncell'].to_numpy()
#print('ncell_output_array',ncell_output_array)

# Add the clump ref solution to the to_check dict.

to_check["ncell_clumpfinder"]=ncell_output_array
to_check["mass_cl_clumpfinder"]=mass_cl_output_array


# then run the check_solution method
visu_ramses.check_solution(to_check, 'cosmo',overwrite=False) # True if you want to overwrite the ref solution



