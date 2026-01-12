'''
    Module to write grafic binary files used as initial conditions for 
    non-cosmo simulations with RAMSES.

    Be careful in cases where the order of the array representation matters!
    Little/big endian issues are NOT dealt with!

    Used for the decaying turbulence test case in tests/hydro/decaying-turbulence
'''

import numpy as np
import struct
import sys

def write_grafic_header(filename, ncells, size, endian='<'):
    ''' Write a simple header for the grafic binary file '''
    # variables for the header line
    n1 = n2 = n3 = int(ncells)
    dx = float(size/ncells)
    xoff1 = xoff2 = xoff3 = float(0.0)
    boxlen = float(size)
    f1 = f2 = f3 = float(0.0)

    # pack variables for Fortran record payload
    payload = struct.pack(
        endian + "3i8f",    # 3 int32 + 8 float32
        n1, n2, n3,
        dx, xoff1, xoff2, xoff3,
        boxlen, f1, f2, f3)

    reclen = len(payload)

    with open(filename, "wb") as f:
        f.write(struct.pack(endian+"i", reclen))  # leading record marker
        f.write(payload)
        f.write(struct.pack(endian+"i", reclen))  # trailing record marker


def write_grafic_data(filename, data, endian='<'):
    ''' Write data in slices to grafic binary file '''
    ncells = data.shape[0]
    data_temp=data.transpose(2,0,1)

    # check if the native byte order matches the requested one
    native_byteorder = sys.byteorder
    swap = ( (endian == '<' and native_byteorder == 'big') or
             (endian == '>' and native_byteorder == 'little') )

    with open(filename, "ab") as f:  # append binary
        for i in range(0,ncells):
            # Ensure float32
            a = np.asarray(data_temp[i,:,:], dtype=np.float32, order="F")

            # convert byte order if needed
            if ((a.dtype.byteorder == "=" and swap) or
                (a.dtype.byteorder == ">" and endian == "<") or
                (a.dtype.byteorder == "<" and endian == ">") ):
                a = a.byteswap().newbyteorder()

            raw = a.tobytes(order="F")
            reclen = len(raw)

            f.write(struct.pack(endian + "i", reclen))
            f.write(raw)
            f.write(struct.pack(endian + "i", reclen))


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
    write_grafic_header(filename, data.shape[0], size, endian)
    write_grafic_data(filename, data, endian)
