import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from scipy.io import FortranFile

class Cool:
    def __init__(self,n1,n2):
        self.n1 = n1
        self.n2 = n2
        self.nH = np.zeros([n1])
        self.T2 = np.zeros([n2])
        self.cool = np.zeros([n1,n2])
        self.heat = np.zeros([n1,n2])
        self.spec = np.zeros([n1,n2,6])
        self.xion = np.zeros([n1,n2])

def clean(dat,n1,n2):
    dat = np.array(dat)
    dat = dat.reshape(n2, n1)
    return dat

def clean_spec(dat,n1,n2):
    dat = np.array(dat)
    dat = dat.reshape(6, n2, n1)
    return dat

def rd_cool(filename):
    with FortranFile(filename, 'r') as f:
        n1, n2 = f.read_ints('i')
        c = Cool(n1,n2)
        nH = f.read_reals('f8')
        T2 = f.read_reals('f8')
        cool = f.read_reals('f8')
        heat = f.read_reals('f8')
        cool_com = f.read_reals('f8')
        heat_com = f.read_reals('f8')
        metal = f.read_reals('f8')
        cool_prime = f.read_reals('f8')
        heat_prime = f.read_reals('f8')
        cool_com_prime = f.read_reals('f8')
        heat_com_prime = f.read_reals('f8')
        metal_prime = f.read_reals('f8')
        mu = f.read_reals('f8')
        n_spec = f.read_reals('f8')
        c.nH = nH
        c.T2 = T2
        c.cool = clean(cool,n1,n2)
        c.heat = clean(heat,n1,n2)
        c.spec = clean_spec(n_spec,n1,n2)
        c.xion = c.spec[0]
        for i in range(0,n2):
            c.xion[i,:] = c.spec[0,i,:] - c.nH
        return c

class Map:
    def __init__(self,nx,ny):
        self.nx = nx
        self.ny = ny
        self.data = np.zeros([nx,ny])

def rd_map(filename):
    with FortranFile(filename, 'r') as f:
        t, dx, dy, dz = f.read_reals('f8')
        nx, ny = f.read_ints('i')
        dat = f.read_reals('f4')
    
    dat = np.array(dat)
    dat = dat.reshape(ny, nx)
    m = Map(nx,ny)
    m.data = dat
    m.time = t
    m.nx = nx
    m.ny = ny
    
    return m

class Part:
    def __init__(self,nnp,nndim):
        self.np = nnp
        self.ndim = nndim
        self.xp = np.zeros([nndim,nnp])
        self.vp = np.zeros([nndim,nnp])
        self.mp = np.zeros([nnp])

def rd_part(nout):
    car1=str(nout).zfill(5)
    filename="output_"+car1+"/part_"+car1+".out00001"
    ncpu = 1
    ndim = 1
    with FortranFile(filename, 'r') as f:
        ncpu2 = f.read_ints('i')
        ndim2 = f.read_ints('i')

    ncpu = ncpu2[0]
    ndim = ndim2[0]
    npart = 0
    for icpu in range(0,ncpu):
        car1 = str(nout).zfill(5)
        car2 = str(icpu+1).zfill(5)
        filename="output_"+car1+"/part_"+car1+".out"+car2
        with FortranFile(filename, 'r') as f:
            ncpu2 = f.read_ints('i')
            ndim2 = f.read_ints('i')
            npart2 = f.read_ints('i')
        npart = npart + npart2[0]
    print("Found ",npart," particles")
    print("Reading data...")

    p = Part(npart,ndim)
    p.np = npart
    p.ndim = ndim
    ipart = 0

    for	icpu in	range(0,ncpu):
        car1 = str(nout).zfill(5)
        car2 = str(icpu+1).zfill(5)
        filename = "output_"+car1+"/part_"+car1+".out"+car2

        with FortranFile(filename, 'r') as f:
            ncpu2 = f.read_ints('i')
            ndim2 = f.read_ints('i')
            npart2 = f.read_ints('i')
            npart2 = npart2[0]
            
            dummy1 = f.read_reals('f8')
            dummy2 = f.read_reals('f4')
            dummy3 = f.read_reals('f8')
            dummy4 = f.read_reals('f8')
            dummy5 = f.read_reals('f4')

            for idim in range(0,ndim):
                xp = f.read_reals('f8')
                p.xp[idim,ipart:ipart+npart2] = xp

            for idim in range(0,ndim):
                xp = f.read_reals('f8')
                p.vp[idim,ipart:ipart+npart2] = xp

            xp = f.read_reals('f8')
            p.mp[ipart:ipart+npart2] = xp

        ipart = ipart + npart2
    return p
