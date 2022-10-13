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
    def __init__(self,np,ndim):
        self.np = np
        self.ndim = ndim
	self.xp = np.zeros([np,ndim])
	self.vp = np.zeros([np,ndim])
	self.mp = np.zeros([np])

def rd_part(nout):
    car=str(nout).zfill(5)
    filename="output_"+car+"+/part_"+car+".out00001"
    with FortranFile(filename, 'r') as f:
        ncpu = f.read_reals('i')
        ndim = f.read_reals('i')
        npart = f.read_reals('i')

