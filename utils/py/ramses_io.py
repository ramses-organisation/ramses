import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from scipy.io import FortranFile
from tqdm import tqdm
from astropy.io import ascii

import time

class Cool:
    """
    This is the class for RAMSES cooling table.
    """
    def __init__(self,n1,n2):
        """
        This function initialize the cooling table. 
        Args:
            n1: number of points for the gas density axis
            n2: number of points for the gas temperature axis
        """
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
    """This function reads a RAMSES cooling table file (unformatted Fortran binary) 
    and store it in a cooling object.

    Args:
        filename: the complete path (including the name) of the cooling table file.

    Returns:
        A cooling table (Cool) object.
    """
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
    """This class defines a map object.
    """
    def __init__(self,nx,ny):
        """This function initalize a map object.
        
        Args:
            nx: number of pixels in the x direction
            ny: number of pixels in the y direction
        """
        self.nx = nx
        self.ny = ny
        self.data = np.zeros([nx,ny])

def rd_map(filename):
    """This function reads a RAMSES map file (unformatted Fortran binary) 
    as produced by the RAMSES utilities amr2map or part2map and store it in a map object.

    Args:
        filename: the complete path (including the name) of the map file.

    Returns:
        A map (class Map) object.
    """
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

def rd_part(nout,**kwargs):
    
    car1 = str(nout).zfill(5)
    filename = "output_"+car1+"/part_"+car1+".out00001"
    with FortranFile(filename, 'r') as f:
        ncpu, = f.read_ints('i')
        ndim, = f.read_ints('i')

    center = kwargs.get("center")
    radius = kwargs.get("radius")

    if ( not (center is None)  and not (radius is None) ):
        info = rd_info(nout)
        cpulist = get_cpu_list(info,**kwargs)
        print("Will open only",len(cpulist),"files")
    else:
        cpulist = range(1,ncpu+1)

    npart = 0
    for icpu in cpulist:
        car1 = str(nout).zfill(5)
        car2 = str(icpu).zfill(5)
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

    for	icpu in	tqdm(cpulist):
        car1 = str(nout).zfill(5)
        car2 = str(icpu).zfill(5)
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

    if ( not (center is None)  and not (radius is None) ):
        # Filtering particles
        r = np.sqrt((p.xp[0]-center[0])**2+(p.xp[1]-center[1])**2+(p.xp[2]-center[2])**2)
        p.np = np.count_nonzero(r < radius)
        p.mp = p.mp[r < radius]
        p.xp = p.xp[:,r < radius]
        p.vp = p.vp[:,r < radius]

    return p

class Level:
    def __init__(self,nndim):
        self.level = 0
        self.ngrid = 0
        self.ndim = nndim
        self.xg = np.empty(shape=(nndim,0))
        self.refined = np.empty(shape=(2**nndim,0),dtype=bool)

def rd_amr(nout,**kwargs):
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

    center = kwargs.get("center")
    radius = kwargs.get("radius")

    if ( not (center is None)  and not (radius is None) ):
        info = rd_info(nout)
        cpulist = get_cpu_list(info,**kwargs)
        print("Will open only",len(cpulist),"files")
    else:
        cpulist = range(1,ncpu+1)

    time.sleep(0.5)

    amr=[]
    for ilevel in range(0,nlevelmax):
        amr.append(Level(ndim))
        
    # Reading and computing total AMR grids count
    for icpu in cpulist:

        car1 = str(nout).zfill(5)
        car2 = str(icpu).zfill(5)
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
            amr[0].boxlen = boxlen

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

            for ilevel in range(0,nlevelmax):
                amr[ilevel].ngrid = amr[ilevel].ngrid + numbl[ilevel,icpu-1]

    # Allocating memory
    for ilevel in range(0,nlevelmax):
        amr[ilevel].xg = np.zeros([ndim,amr[ilevel].ngrid],dtype=float)
        amr[ilevel].refined = np.zeros([2**ndim,amr[ilevel].ngrid],dtype=bool)

    iskip = np.zeros(nlevelmax, dtype=int)
    
    # Reading and storing data
    for icpu in tqdm(cpulist):

        car1 = str(nout).zfill(5)
        car2 = str(icpu).zfill(5)
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
                            xg[idim,:] = (f.read_reals('f8')-xbound[idim])*boxlen
                        if(ibound == icpu-1):
                            amr[ilevel].xg[:,iskip[ilevel]:iskip[ilevel]+ncache] = xg
                        father = f.read_ints("i")
                        for ind in range(0,2*ndim):
                            nbor = f.read_ints("i")
                        son = np.zeros([2**ndim,ncache])
                        for ind in range(0,2**ndim):
                            son[ind,:] = f.read_ints("i")
                        if(ibound == icpu-1):
                            ref = np.zeros([2**ndim,ncache],dtype=bool)
                            ref = np.where(son > 0, True, False)
                            amr[ilevel].refined[:,iskip[ilevel]:iskip[ilevel]+ncache] = ref
                        for ind in range(0,2**ndim):
                            cpumap = f.read_ints("i")
                        for ind in range(0,2**ndim):
                            flag1 = f.read_ints("i")
                        if(ibound == icpu-1):
                            iskip[ilevel] = iskip[ilevel] + ncache

    return amr

class Hydro:
    def __init__(self,nndim,nnvar):
        self.level = 0
        self.ngrid = 0
        self.ndim = nndim
        self.nvar = nnvar
        self.u = np.empty(shape=(nnvar,2**nndim,0))

def rd_hydro(nout,**kwargs):
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

    center = kwargs.get("center")
    radius = kwargs.get("radius")

    if ( not (center is None)  and not (radius is None) ):
        info = rd_info(nout)
        cpulist = get_cpu_list(info,**kwargs)
        print("Will open only",len(cpulist),"files")
    else:
        cpulist = range(1,ncpu+1)

    time.sleep(0.5)

    hydro=[]
    for ilevel in range(0,nlevelmax):
        hydro.append(Hydro(ndim,nvar))
        hydro[ilevel].level = ilevel

    for icpu in tqdm(cpulist):

        car1 = str(nout).zfill(5)
        car2 = str(icpu).zfill(5)
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
                            
                        if(ibound == icpu-1):
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
        self.dx = np.empty(shape=(0))

def rd_cell(nout,**kwargs):
    
    a = rd_amr(nout,**kwargs)
    h = rd_hydro(nout,**kwargs)

    center = kwargs.get("center")
    radius = kwargs.get("radius")

    nlevelmax = len(a)
    ndim = a[0].ndim
    nvar = h[0].nvar
    boxlen = a[0].boxlen
    
    offset = np.zeros([ndim,2**ndim])
    if (ndim == 2):
        offset[0,:]=[-0.5,0.5,-0.5,0.5]
        offset[1,:]=[-0.5,-0.5,0.5,0.5]
    if (ndim == 3):
        offset[0,:]=[-0.5,0.5,-0.5,0.5,-0.5,0.5,-0.5,0.5]
        offset[1,:]=[-0.5,-0.5,0.5,0.5,-0.5,-0.5,0.5,0.5]
        offset[2,:]=[-0.5,-0.5,-0.5,-0.5,0.5,0.5,0.5,0.5]

    ncell = 0
    for ilev in range(0,nlevelmax):
        ncell = ncell + np.count_nonzero(a[ilev].refined == False)


    time.sleep(0.5)
    print("Found",ncell,"leaf cells")
    print("Extracting leaf cells...")

    c = Cell(ndim,nvar)
    c.nc = ncell
    
    for ilev in range(0,nlevelmax):
        dx = 0.5*boxlen/2**ilev
        for ind in range(0,2**ndim):
            nc = np.count_nonzero(a[ilev].refined[ind] == False)
            if (nc > 0):
                xc = np.zeros([ndim,nc])
                for idim in range(0,ndim):
                    xc[idim,:]= a[ilev].xg[idim,np.where(a[ilev].refined[ind] == False)]+offset[idim,ind]*dx
                c.x = np.append(c.x,xc,axis=1)
                uc = np.zeros([nvar,nc])
                for ivar in range(0,nvar):
                    uc[ivar,:]= h[ilev].u[ivar,ind,np.where(a[ilev].refined[ind] == False)]
                c.u = np.append(c.u,uc,axis=1)
                dd = np.ones(nc)*dx
                c.dx = np.append(c.dx,dd)

    if ( not (center is None)  and not (radius is None) ):
        # Filtering cells
        r = np.sqrt((c.x[0]-center[0])**2+(c.x[1]-center[1])**2+(c.x[2]-center[2])**2) - dx
        c.nc = np.count_nonzero(r < radius)
        c.u  = c.u[:,r < radius]
        c.x  = c.x[:,r < radius]
        c.dx = c.dx[r < radius]

    return c

def save_cell(c,filename):

    with open(filename,'wb') as f:
        np.save(f,c.nc)
        np.save(f,c.ndim)
        np.save(f,c.nvar)
        np.save(f,c.dx)
        np.save(f,c.x)
        np.save(f,c.u)

def load_cell(filename):

    with open(filename,'rb') as f:
        nc = np.load(f)
        ndim = np.load(f)
        nvar = np.load(f)
        c = Cell(ndim,nvar)
        c.nc = nc
        c.ndim = ndim
        c.nvar = nvar
        c.dx = np.append(c.dx,np.load(f))
        c.x  = np.append(c.x, np.load(f),axis=1)
        c.u  = np.append(c.u, np.load(f),axis=1)

    return c

class Info:
    def __init__(self,nncpu):
        self.bound_key = np.empty(shape=(nncpu+1),dtype=np.double)
        
def rd_info(nout):
    car1 = str(nout).zfill(5)
    filename = "output_"+car1+"/amr_"+car1+".out00001"
    with FortranFile(filename, 'r') as f:
        ncpu, = f.read_ints('i')
        ndim, = f.read_ints('i')
        nx,ny,nz = f.read_ints('i')
        nlevelmax, = f.read_ints('i')

        txt = "ncpu="+str(ncpu)+" ndim="+str(ndim)+" nlevelmax="+str(nlevelmax)
        print(txt)
        print("Reading info data...")
                
        i = Info(ncpu)
        i.nlevelmax = nlevelmax
        i.ndim = ndim
        i.ncpu = ncpu
        
        ngridmax, = f.read_ints('i')
        nboundary, = f.read_ints('i')
        ngrid_current, = f.read_ints('i')
        boxlen, = f.read_reals('f8')

        i.boxlen = boxlen
        
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
        
        i.omega_m = omega_m
        i.omega_l = omega_l
        i.omega_k = omega_k
        i.omega_b = omega_b
        i.h0 = h0
        i.aexp = aexp
        i.t = t
        
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
        
        bound_key = f.read_ints("f8")
        
        i.bound_key[:] = bound_key
        
    filename = "output_"+car1+"/info_"+car1+".txt"
    data = ascii.read(filename, header_start=0, data_start=0, data_end=18, delimiter='=', names=["field","value"])
    name = np.array(data["field"])
    val = np.array(data["value"])
    
    i.ordering, = val[np.where(name=="ordering type")]
    unit_l, = val[np.where(name=="unit_l")]
    unit_d, = val[np.where(name=="unit_d")]
    unit_t, = val[np.where(name=="unit_t")]

    i.unit_l = float(unit_l)
    i.unit_d = float(unit_d)
    i.unit_t = float(unit_t)
    
    return i

def hilbert3d(x,y,z,bit_length):
    
    state_diagram = [ 1, 2, 3, 2, 4, 5, 3, 5,
                      0, 1, 3, 2, 7, 6, 4, 5,
                      2, 6, 0, 7, 8, 8, 0, 7,
                      0, 7, 1, 6, 3, 4, 2, 5,
                      0, 9,10, 9, 1, 1,11,11,
                      0, 3, 7, 4, 1, 2, 6, 5,
                      6, 0, 6,11, 9, 0, 9, 8,
                      2, 3, 1, 0, 5, 4, 6, 7,
                      11,11, 0, 7, 5, 9, 0, 7,
                      4, 3, 5, 2, 7, 0, 6, 1,
                      4, 4, 8, 8, 0, 6,10, 6,
                      6, 5, 1, 2, 7, 4, 0, 3,
                      5, 7, 5, 3, 1, 1,11,11,
                      4, 7, 3, 0, 5, 6, 2, 1,
                      6, 1, 6,10, 9, 4, 9,10,
                      6, 7, 5, 4, 1, 0, 2, 3,
                      10, 3, 1, 1,10, 3, 5, 9,
                      2, 5, 3, 4, 1, 6, 0, 7,
                      4, 4, 8, 8, 2, 7, 2, 3,
                      2, 1, 5, 6, 3, 0, 4, 7,
                      7, 2,11, 2, 7, 5, 8, 5,
                      4, 5, 7, 6, 3, 2, 0, 1,
                      10, 3, 2, 6,10, 3, 4, 4,
                      6, 1, 7, 0, 5, 2, 4, 3]

    state_diagram = np.array(state_diagram)
    state_diagram = state_diagram.reshape((8,2,12),order='F')

    n = len(x)
    order = np.zeros(n,dtype="double")
    x_bit_mask = np.zeros(bit_length  ,dtype="bool")
    y_bit_mask = np.zeros(bit_length  ,dtype="bool")
    z_bit_mask = np.zeros(bit_length  ,dtype="bool")
    i_bit_mask = np.zeros(3*bit_length,dtype=bool)
    
    for ip in  range(0,n):
        
        for i in range(0,bit_length):
            x_bit_mask[i] = x[ip] & (1 << i)
            y_bit_mask[i] = y[ip] & (1 << i)
            z_bit_mask[i] = z[ip] & (1 << i)
            
        for i in range(0,bit_length):
            i_bit_mask[3*i+2] = x_bit_mask[i]
            i_bit_mask[3*i+1] = y_bit_mask[i]
            i_bit_mask[3*i  ] = z_bit_mask[i]
            
        cstate = 0
        for i in range(bit_length-1,-1,-1):
            b2 = 0
            if (i_bit_mask[3*i+2]):
                b2 = 1
            b1 = 0
            if (i_bit_mask[3*i+1]):
                b1 = 1
            b0 = 0
            if (i_bit_mask[3*i  ]):
                b0 = 1
            sdigit = b2*4 + b1*2 + b0
            nstate = state_diagram[sdigit,0,cstate]
            hdigit = state_diagram[sdigit,1,cstate]
            i_bit_mask[3*i+2] = hdigit & (1 << 2)
            i_bit_mask[3*i+1] = hdigit & (1 << 1)
            i_bit_mask[3*i  ] = hdigit & (1 << 0)
            cstate = nstate
            
        order[ip]= 0
        for i in range(0,3*bit_length):
            b0 = 0
            if (i_bit_mask[i]):
                b0 = 1
            order[ip] = order[ip] + float(b0)*2.**i
                
    return order

def hilbert2d(x,y,bit_length):
    
    state_diagram = [ 1, 0, 2, 0, 
                      0, 1, 3, 2, 
                      0, 3, 1, 1, 
                      0, 3, 1, 2, 
                      2, 2, 0, 3, 
                      2, 1, 3, 0, 
                      3, 1, 3, 2, 
                      2, 3, 1, 0 ]
    
    state_diagram = np.array(state_diagram)    
    state_diagram = state_diagram.reshape((4,2,4), order='F')
    
    n = len(x)
    order = np.zeros(n,dtype="double")
    x_bit_mask = np.zeros(bit_length  ,dtype="bool")
    y_bit_mask = np.zeros(bit_length  ,dtype="bool")
    i_bit_mask = np.zeros(2*bit_length,dtype=bool)
    
    for ip in  range(0,n):
        
        for i in range(0,bit_length):
            x_bit_mask[i] = bool(x[ip] & (1 << i))
            y_bit_mask[i] = bool(y[ip] & (1 << i))
            
        for i in range(0,bit_length):
            i_bit_mask[2*i+1] = x_bit_mask[i]
            i_bit_mask[2*i  ] = y_bit_mask[i]
            
        cstate = 0
        for i in range(bit_length-1,-1,-1):
            b1 = 0
            if (i_bit_mask[2*i+1]):
                b1 = 1
            b0 = 0
            if (i_bit_mask[2*i  ]):
                b0 = 1
            sdigit = b1*2 + b0
            nstate = state_diagram[sdigit,0,cstate]
            hdigit = state_diagram[sdigit,1,cstate]
            i_bit_mask[2*i+1] = hdigit & (1 << 1)
            i_bit_mask[2*i  ] = hdigit & (1 << 0)
            cstate = nstate
            
        order[ip]= 0
        for i in range(0,2*bit_length):
            b0 = 0
            if (i_bit_mask[i]):
                b0 = 1
            order[ip] = order[ip] + float(b0)*2.**i
                
    return order

def get_cpu_list(info,**kwargs):

    center = kwargs.get("center")
    radius = kwargs.get("radius")
    center = np.array(center)
    radius = float(radius)
    
    for ilevel in range(0,info.nlevelmax):
        dx = 1/2**ilevel
        if (dx < 2*radius/info.boxlen):
            break

    levelmin = np.max([ilevel,1])
    bit_length = levelmin-1
    nmax = 2**bit_length
    ndim = info.ndim
    ncpu = info.ncpu
    nlevelmax = info.nlevelmax
    dkey = 2**(ndim*(nlevelmax+1-bit_length))
    ibound = [0, 0, 0, 0, 0, 0]
    if(bit_length > 0):
        ibound[0:3] = (center-radius)*nmax/info.boxlen
        ibound[3:6] = (center+radius)*nmax/info.boxlen
        ibound[0:3] = np.array(ibound[0:3]).astype(int)
        ibound[3:6] = np.array(ibound[3:6]).astype(int)
        ndom = 8
        idom = [ibound[0], ibound[3], ibound[0], ibound[3], ibound[0], ibound[3], ibound[0], ibound[3]]
        jdom = [ibound[1], ibound[1], ibound[4], ibound[4], ibound[1], ibound[1], ibound[4], ibound[4]]
        kdom = [ibound[2], ibound[2], ibound[2], ibound[2], ibound[5], ibound[5], ibound[5], ibound[5]]
        order_min = hilbert3d(idom,jdom,kdom,bit_length)
    else:
        ndom = 1
        order_min = np.array([0.])
        
    bounding_min = order_min*dkey
    bounding_max = (order_min+1)*dkey

    cpu_min = np.zeros(ndom, dtype=int)
    cpu_max = np.zeros(ndom, dtype=int)
    for icpu in range(0,ncpu):
        for idom in range(0,ndom):
            if( (info.bound_key[icpu] <= bounding_min[idom]) and (info.bound_key[icpu+1] > bounding_min[idom]) ):
                cpu_min[idom] = icpu+1
            if( (info.bound_key[icpu] < bounding_max[idom]) and (info.bound_key[icpu+1] >= bounding_max[idom]) ):
                cpu_max[idom] = icpu+1


    ncpu_read = 0
    cpu_read = np.zeros(ncpu, dtype=bool)
    cpu_list = []
    for idom in range(0,ndom):
        for icpu in range(cpu_min[idom]-1,cpu_max[idom]):
            if ( not cpu_read[icpu] ):
                cpu_list.append(icpu+1)
                ncpu_read = ncpu_read+1
                cpu_read[icpu] = True

    return cpu_list
    
def visu2d(c,**kwargs):

    ivar = kwargs.get("ivar",0)
    log = kwargs.get("log",None)
    
    px = 1/plt.rcParams['figure.dpi'] 
    fig, ax = plt.subplots(figsize=(1000*px,1000*px))
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
    plt.scatter(c.x[0],c.x[1],s=0.0001)
    ymin, ymax = ax.get_ylim()

    vv = c.u[ivar]
    if( not (log is None)):
        vv = np.log10(c.u[ivar])

    print(np.min(vv),np.max(vv))

    plt.scatter(c.x[0],c.x[1],c=vv,s=(c.dx*800/(ymax-ymin))**2,cmap="viridis",marker="s")


def visu3d(c, **kwargs):

    ivar = kwargs.get("ivar",0)
    log = kwargs.get("log",None)
    axis = kwargs.get("axis","z")
    vmin = kwargs.get("vmin",None)
    vmax = kwargs.get("vmax",None)

    if(ivar>=0) and ivar<c.nvar:
        vv=c.u[ivar]
    elif(ivar == 25):
        vv=c.u[4]/c.u[0]
    elif(ivar == 35):
        vv=c.u[10]/c.u[0]
    else:
        vv=c.dx
    if( not (log is None)):
        vv = np.log10(vv)

    if(axis == "z"):
        i1=0
        i2=1
    if(axis == "y"):
        i1=0
        i2=2
    if(axis == "x"):
        i1=1
        i2=2
                
    print("min=",np.min(vv)," max=",np.max(vv))

    px = 1/plt.rcParams['figure.dpi'] 
    fig, ax = plt.subplots(figsize=(1000*px,1000*px))
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1)
    plt.scatter(c.x[i1],c.x[i2],s=0.0001)
    ymin, ymax = ax.get_ylim()
    ind = np.argsort(np.log10(c.u[0]))
        
    ax.set_aspect("equal")
    plt.scatter(c.x[i1,ind],c.x[i2,ind],c=vv[ind],s=(c.dx[ind]*800/(ymax-ymin))**2,cmap="viridis",marker="s",vmin=vmin,vmax=vmax)
    plt.colorbar(location="bottom",shrink=0.7)
