'''
    Module to write grafic binary files used as initial conditions for
    non-cosmo simulations with RAMSES.

    Used for the decaying turbulence test case in tests/hydro/decaying-turbulence
'''
from scipy.io import FortranFile
import numpy as np
import sys

def write_grafic_header(f, ncells, size, endian):
    ''' Write a simple header for the grafic binary file '''
    i4 = np.dtype(f"{endian}i4").type
    f4 = np.dtype(f"{endian}f4").type
    # variables for the header line
    n1 = n2 = n3 = i4(ncells)
    dx = f4(size/ncells)
    xoff1 = xoff2 = xoff3 = f4(0.0)
    boxlen = f4(size)
    f1 = f2 = f3 = f4(0.0)

    f.write_record(n1, n2, n3, dx, xoff1, xoff2, xoff3, boxlen, f1, f2, f3)


def write_grafic_data(f, data, endian):
    ''' Write data in slices to grafic binary file '''
    ncells = data.shape[0]
    data_temp=data.transpose(2,0,1)

    dtype = np.dtype(f"{endian}f4")

    for i in range(0,ncells):
        # Convert input array into float32 with selected endianness
        a = np.array(data_temp[i,:,:], dtype=dtype)
        # flatten in Fortran order explicitly
        f.write_record(a.ravel(order="F"))


def write_grafic_file(filename, data, size, endian='='):
    ''' Write a grafic binary file with header and data.
        Default endian is the native one. Files will be consistent
        when written and read by the same system. If ICs are generated
        on different systems, use '<' for little-endian or '>' for big-endian,
        to match the endianness of the system running RAMSES. '''
    print('writing grafic file',filename)
    if endian == '=':
        endian_string = sys.byteorder
        if endian_string == 'little':
            endian = '<'
        else:
            endian = '>'

    header_dtype = np.dtype(endian + "i4")
    with FortranFile(filename, "w", header_dtype=header_dtype) as f:
        write_grafic_header(f, data.shape[0], size, endian)
        write_grafic_data(f, data, endian)
