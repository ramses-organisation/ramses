import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from scipy.io import FortranFile
from tqdm import tqdm
from astropy.io import ascii

import time

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
    car1 = str(nout).zfill(5)
    filename = "output_"+car1+"/part_"+car1+".out00001"
    with FortranFile(filename, 'r') as f:
        ncpu, = f.read_ints('i')
        ndim, = f.read_ints('i')

    npart = 0
    for icpu in range(0,ncpu):
        car1 = str(nout).zfill(5)
        car2 = str(icpu+1).zfill(5)
        filename="output_"+car1+"/part_"+car1+".out"+car2
        with FortranFile(filename, 'r') as f:
            ncpu2, = f.read_ints('i')
            ndim2, = f.read_ints('i')
            npart2, = f.read_ints('i')
        npart = npart + npart2
        
    txt = "Found "+str(npart)+" particles"
    tqdm.write(txt)
    tqdm.write("Reading particle data...")
    time.sleep(0.5)
    
    p = Part(npart,ndim)
    p.np = npart
    p.ndim = ndim
    ipart = 0

    for	icpu in	tqdm(range(0,ncpu)):
        car1 = str(nout).zfill(5)
        car2 = str(icpu+1).zfill(5)
        filename = "output_"+car1+"/part_"+car1+".out"+car2

        with FortranFile(filename, 'r') as f:
            ncpu2, = f.read_ints('i')
            ndim2, = f.read_ints('i')
            npart2, = f.read_ints('i')
            
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

class Level:
    def __init__(self,iilev,nndim):
        self.level = iilev
        self.ngrid = 0
        self.ndim = nndim
        self.xg = np.empty(shape=(nndim,0))
        self.refined = np.empty(shape=(2**nndim,0),dtype=bool)

def rd_amr(nout):
    car1 = str(nout).zfill(5)
    filename = "output_"+car1+"/amr_"+car1+".out00001"
    with FortranFile(filename, 'r') as f:
        ncpu, = f.read_ints('i')
        ndim, = f.read_ints('i')
        nx,ny,nz = f.read_ints('i')
        nlevelmax, = f.read_ints('i')

    txt = "ncpu="+str(ncpu)+" ndim="+str(ndim)+" nlevelmax="+str(nlevelmax)
    tqdm.write(txt)
    tqdm.write("Reading grid data...")
    time.sleep(0.5)

    amr=[]
    for ilevel in range(0,nlevelmax):
        amr.append(Level(ilevel,ndim))

    for icpu in tqdm(range(0,ncpu)):

        car1 = str(nout).zfill(5)
        car2 = str(icpu+1).zfill(5)
        filename = "output_"+car1+"/amr_"+car1+".out"+car2

        with FortranFile(filename, 'r') as f:
            ncpu2, = f.read_ints('i')
            ndim2, = f.read_ints('i')
            nx2,ny2,nz2 = f.read_ints('i')
            nlevelmax2, = f.read_ints('i')
            ngridmax, = f.read_ints('i')
            nboundary, = f.read_ints('i')
            ngrid_current, = f.read_ints('i')
            boxlen, = f.read_reals('f8')

            noutput,iout,ifout = f.read_ints('i')
            tout = f.read_reals('f8')
            aout = f.read_reals('f8')
            t, = f.read_reals('f8')
            dtold = f.read_reals('f8')
            dtnew = f.read_reals('f8')
            nstep,nstep_coarse = f.read_ints('i')
            einit,mass_tot_0,rho_tot = f.read_reals('f8')
            omega_m,omega_l,omega_k,omega_b,h0,aexp_ini,boxlen_ini = f.read_reals('f8')
            aexp,hexp,aexp_old,epot_tot_int,epot_tot_old = f.read_reals('f8')
            mass_sph, = f.read_reals('f8')

            headl = f.read_ints('i')
            taill = f.read_ints('i')
            numbl = f.read_ints('i')
            numbl = numbl.reshape(nlevelmax,ncpu)

            numbtot = f.read_ints('i')

            xbound=[0,0,0]
            if ( nboundary > 0 ):
                headb = f.read_ints('i')
                tailb = f.read_ints('i')
                numbb = f.read_ints('i')
                numbb = numbb.reshape(nlevelmax,nboundary)
                xbound = [float(nx//2),float(ny//2),float(nz//2)]
                
            headf,tailf,numbf,used_mem,used_mem_tot = f.read_ints('i')

            ordering = f.read_ints("i")

            bound_key = f.read_ints("i8")

            son = f.read_ints("i")
            flag1 = f.read_ints("i")
            cpu_map = f.read_ints("i")

            for ilevel in range(0,nlevelmax):
                for ibound in range(0,nboundary+ncpu):
                    if(ibound<ncpu):
                        ncache=numbl[ilevel,ibound]
                    else:
                        ncache=numbb[ilevel,ibound-ncpu]

                    if (ncache>0):
                        index = f.read_ints("i")
                        nextg = f.read_ints("i")
                        prevg = f.read_ints("i")
                        xg = np.zeros([ndim,ncache])
                        for idim in range(0,ndim):
                            xg[idim,:] = f.read_reals('f8')-xbound[idim]
                        if(ibound == icpu):
                            amr[ilevel].xg=np.append(amr[ilevel].xg,xg,axis=1)
                            amr[ilevel].ngrid=amr[ilevel].ngrid+ncache
                        father = f.read_ints("i")
                        for ind in range(0,2*ndim):
                            nbor = f.read_ints("i")
                        son = np.zeros([2**ndim,ncache])
                        for ind in range(0,2**ndim):
                            son[ind,:] = f.read_ints("i")
                        if(ibound == icpu):
                            ref = np.zeros([2**ndim,ncache],dtype=bool)
                            ref = np.where(son > 0, True, False)
                            amr[ilevel].refined=np.append(amr[ilevel].refined,ref,axis=1)
                        for ind in range(0,2**ndim):
                            cpumap = f.read_ints("i")
                        for ind in range(0,2**ndim):
                            flag1 = f.read_ints("i")                        

    return amr

class Hydro:
    def __init__(self,iilev,nndim,nnvar):
        self.level = iilev
        self.ngrid = 0
        self.ndim = nndim
        self.nvar = nnvar
        self.u = np.empty(shape=(nnvar,2**nndim,0))

def rd_hydro(nout):
    car1 = str(nout).zfill(5)
    filename = "output_"+car1+"/hydro_"+car1+".out00001"
    with FortranFile(filename, 'r') as f:
        ncpu, = f.read_ints('i')
        nvar, = f.read_ints('i')
        ndim, = f.read_ints('i')
        nlevelmax, = f.read_ints('i')
        nboundary, = f.read_ints('i')
        gamma, = f.read_reals('f8')
        
    txt = "ncpu="+str(ncpu)+" ndim="+str(ndim)+" nvar="+str(nvar)+" nlevelmax="+str(nlevelmax)+" gamma="+str(gamma)
    tqdm.write(txt)
    tqdm.write("Reading hydro data...")
    time.sleep(0.5)

    hydro=[]
    for ilevel in range(0,nlevelmax):
        hydro.append(Hydro(ilevel,ndim,nvar))

    for icpu in tqdm(range(0,ncpu)):

        car1 = str(nout).zfill(5)
        car2 = str(icpu+1).zfill(5)
        filename = "output_"+car1+"/hydro_"+car1+".out"+car2

        with FortranFile(filename, 'r') as f:
            ncpu2, = f.read_ints('i')
            nvar2, = f.read_ints('i')
            ndim2, = f.read_ints('i')
            nlevelmax2, = f.read_ints('i')
            nboundary2, = f.read_ints('i')
            gamma2, = f.read_reals('f8')

            for ilevel in range(0,nlevelmax):
                for ibound in range(0,nboundary+ncpu):
                    ilevel2, = f.read_ints('i')
                    ncache, = f.read_ints('i')

                    if (ncache>0):
                        uu = np.zeros([nvar,2**ndim,ncache])
                        for ind in range(0,2**ndim):
                            for ivar in range(0,nvar):
                                uu[ivar,ind,:] = f.read_reals('f8')
                            
                        if(ibound == icpu):
                            hydro[ilevel].u = np.append(hydro[ilevel].u,uu,axis=2)
                            hydro[ilevel].ngrid = hydro[ilevel].ngrid + ncache                        

    return hydro
    
class Cell:
    def __init__(self,nndim,nnvar):
        self.ncell = 0
        self.ndim = nndim
        self.nvar = nnvar
        self.x = np.empty(shape=(nndim,0))
        self.u = np.empty(shape=(nnvar,0))

def rd_cell(nout):
    
    a = rd_amr(nout)
    h = rd_hydro(nout)

    nlevelmax = len(a)
    ndim = a[0].ndim
    nvar = h[0].nvar
    
    offset = np.zeros([ndim,2**ndim])
    offset[0,:]=[-0.5,0.5,-0.5,0.5,-0.5,0.5,-0.5,0.5]
    offset[1,:]=[-0.5,-0.5,0.5,0.5,-0.5,-0.5,0.5,0.5]
    offset[2,:]=[-0.5,-0.5,-0.5,-0.5,0.5,0.5,0.5,0.5]

    ncell = 0
    for ilev in range(0,nlevelmax):
        ncell = ncell + np.count_nonzero(a[ilev].refined == False)

    print("Found",ncell,"leaf cells")
    print("Extracting leaf cells...")

    c = Cell(ndim,nvar)
        
    for ilev in range(0,nlevelmax):
        dx = 1./2**ilev
        for ind in range(0,2**ndim):
            nc = np.count_nonzero(a[ilev].refined[ind] == False)
            if (nc > 0):
                xc = np.zeros([ndim,nc])
                for idim in range(0,ndim):
                    xc[idim,:]= a[ilev].xg[idim,np.where(a[ilev].refined[ind] == False)]+offset[idim,ind]*dx/2
                c.x = np.append(c.x,xc,axis=1)
                uc = np.zeros([nvar,nc])
                for ivar in range(0,nvar):
                    uc[ivar,:]= h[ilev].u[ivar,ind,np.where(a[ilev].refined[ind] == False)]
                c.u = np.append(c.u,uc,axis=1)

    return c

class Info:
    def __init__(self):
        self.nlevelmax = 0
        self.ndim = 0

def rd_info(nout):
    car1 = str(nout).zfill(5)
    filename = "output_"+car1+"/amr_"+car1+".out00001"
    with FortranFile(filename, 'r') as f:
        ncpu, = f.read_ints('i')
        ndim, = f.read_ints('i')
        nx,ny,nz = f.read_ints('i')
        nlevelmax, = f.read_ints('i')

        txt = "ncpu="+str(ncpu)+" ndim="+str(ndim)+" nlevelmax="+str(nlevelmax)
        tqdm.write(txt)
        tqdm.write("Reading info data...")
        time.sleep(0.5)
        
        i = Info()
        
        ngridmax, = f.read_ints('i')
        nboundary, = f.read_ints('i')
        ngrid_current, = f.read_ints('i')
        boxlen, = f.read_reals('f8')
        
        noutput,iout,ifout = f.read_ints('i')
        tout = f.read_reals('f8')
        aout = f.read_reals('f8')
        t, = f.read_reals('f8')
        dtold = f.read_reals('f8')
        dtnew = f.read_reals('f8')
        nstep,nstep_coarse = f.read_ints('i')
        einit,mass_tot_0,rho_tot = f.read_reals('f8')
        omega_m,omega_l,omega_k,omega_b,h0,aexp_ini,boxlen_ini = f.read_reals('f8')
        aexp,hexp,aexp_old,epot_tot_int,epot_tot_old = f.read_reals('f8')
        mass_sph, = f.read_reals('f8')
        
        headl = f.read_ints('i')
        taill = f.read_ints('i')
        numbl = f.read_ints('i')
        numbl = numbl.reshape(nlevelmax,ncpu)
        
        numbtot = f.read_ints('i')
        
        xbound=[0,0,0]
        if ( nboundary > 0 ):
            headb = f.read_ints('i')
            tailb = f.read_ints('i')
            numbb = f.read_ints('i')
            numbb = numbb.reshape(nlevelmax,nboundary)
            xbound = [float(nx//2),float(ny//2),float(nz//2)]
            
        headf,tailf,numbf,used_mem,used_mem_tot = f.read_ints('i')
        
        ordering = f.read_ints("i")
        
        bound_key = f.read_ints("i8")

    filename = "output_"+car1+"/info_"+car1+".txt"

    return i
