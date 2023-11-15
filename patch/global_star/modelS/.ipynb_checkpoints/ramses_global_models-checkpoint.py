"""
    ramses_global_models.py

    José Roberto Canivete Cuissa
    06.12.2022, Zurich

This code contains routines to read and analyze Ramses convective global
models.
"""

# Imports and constants
import numpy as np
import osyris
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

kB = 1.380649*10**-16 # erg K^-1
m = 1.672*10**-24 #g
# ------------------


class RamsesOutput:
    """
    An object containing the output of a ramses simulation
    """
    def __init__(self):
        """
        Constructor
        """
        # Variables
        self._variables = {'density':None,
                           'velocity':None, 
                           'B_field':None, 
                           'pressure':None, 
                           'equilibrium_density':None,
                           'equilibrium_pressure':None}
        self._derived_variables = {
                           'entropy':None,
                           'temperature':None,
                           'equilibrium_entropy':None,
                           'equilibrium_temperature':None,
                           'internal_energy':None,
                           'total_energy':None,
                           'equilibrium_internal_energy':None,
                           'equilibrium_total_energy':None,}
        # Info
        self.dx = None
        self.boxlen = None
        self.levelmax = None
        self.ndim = None
        self.meta = None
        self.gamma = None
        self.time = None
        # Other attributes
    # ------------------

    @property
    def d(self):
        return self._variables['density']
    # -----------------------------

    @property
    def p(self):
        return self._variables['pressure']
    # -----------------------------

    @property
    def B(self):
        return self._variables['B_field']
    # -----------------------------

    @property
    def v(self):
        return self._variables['velocity']
    # -----------------------------

    @property
    def deq(self):
        return self._variables['equilibrium_density']
    # -----------------------------

    @property
    def peq(self):
        return self._variables['equilibrium_pressure']
    # -------------------------------

    @property
    def T(self):
        if self._derived_variables['temperature'] is None:
            d = self._variables['density']
            p = self._variables['pressure']
            self._derived_variables['temperature'] = p*m/(d*kB)
        return self._derived_variables['temperature']
        
    # -------------------------------

    @property
    def Teq(self):
        if self._derived_variables['equilibrium_temperature'] is None:
            d = self._variables['equilibrium_density']
            p = self._variables['equilibrium_pressure']
            self._derived_variables['equilibrium_temperature'] = p*m/(d*kB)
        return self._derived_variables['equilibrium_temperature']
    # -------------------------------

    @property
    def s(self):
        if self._derived_variables['entropy'] is None:
            d = self._variables['density']
            p = self._variables['pressure']
            gamma = self.gamma
            self._derived_variables['entropy'] = p/(d**gamma)
        return self._derived_variables['entropy']
    # -------------------------------

    @property
    def seq(self):
        if self._derived_variables['equilibrium_entropy'] is None:
            d = self._variables['equilibrium_density']
            p = self._variables['equilibrium_pressure']
            gamma = self.gamma
            self._derived_variables['equilibrium_entropy'] = p/(d**gamma)
        return self._derived_variables['equilibrium_entropy']
    # -------------------------------
    
    @property
    def e(self):
        if self._derived_variables['internal_energy'] is None:
            p = self._variables['pressure']
            gamma = self.gamma
            self._derived_variables['internal_energy'] = p/(gamma-1.0)
        return self._derived_variables['internal_energy']
    # -------------------------------
    
    @property
    def eeq(self):
        if self._derived_variables['equilibrium_internal_energy'] is None:
            p = self._variables['equilibrium_pressure']
            gamma = self.gamma
            self._derived_variables['equilibrium_internal_energy'] = p/(gamma-1.0)
        return self._derived_variables['equilibrium_internal_energy']
    # -------------------------------
    
    @property
    def etot(self):
        if self._derived_variables['total_energy'] is None:
            e = self.e
            d = self.d
            v_norm = np.sqrt(self.v.x**2 + self.v.y**2)
            B_norm = np.sqrt(self.B.x**2 + self.B.y**2)
            self._derived_variables['total_energy'] = e + 0.5*d*v_norm**2 + 0.5*B_norm**2
        return self._derived_variables['total_energy']
    # --------------------------------------------------------------
    
    @property
    def etoteq(self):
        if self._derived_variables['equilibrium_internal_energy'] is None:
            e = self.eeq
            self._derived_variables['equilibrium_total_energy'] = e 
        return self._derived_variables['equilibrium_internal_energy']
    # -------------------------------

    def read(self, folder, n, vars=[], save=None, norms={}):
        """
        This method reads a binary output and saves the variables
        into the object attributes. One can choose which variables
        to store with vars and if to save the quantities in a h5py
        file.

        Parameters:
        -----------
            folder : string
                The path to the output folder
            n : int
                The number of the output
            vars : list = []
                The list of variables to store. Possible choices:
                ['v', 'd', 'p', 'B', 'deq', 'peq']
            save : string = None
                If save is a string, then it saves the output in 
                a hdf5 file with h5py that can be loaded with 
                this same class.
        """
        # Check variables to fetch
        if vars==[]:
            vars = list(self._variables.keys())
        # Read output
        output, meta = readRamsesOutput(folder, n, vars)
        # Save variables
        for var_name in self._variables.keys():
            self._variables[var_name] = output[var_name]
        self.normalize(norms)
        # Save info
        self.meta = meta
        self.boxlen = self.meta['boxlen']
        self.levelmax = self.meta['levelmax']
        self.ndim = self.meta['ndim']
        self.gamma = self.meta['gamma']
        self.time = self.meta['time']
        self.dx = self.boxlen/(2**self.levelmax)
        
        # Save output 
        # TODO
        if save:
            pass
        
    def circles(self):
        r_center_1 = 4.86990e2
        r_center_2 = 6.60915e2
        heating_dx = 20
        Nx = 2**self.levelmax
        c = Nx//2
        b1 = int(r_center_1/self.dx) 
        b2 = int((r_center_1+heating_dx)/self.dx) 
        b3 = int((r_center_2-heating_dx)/self.dx) 
        b4 = int(r_center_2/self.dx) 
        
        bottom_interface = plt.Circle((c, c), b1, color='k', fill=False)
        heating_interface = plt.Circle((c, c), b2, color='r', fill=False)
        cooling_interface = plt.Circle((c, c), b3, color='b', fill=False)
        upper_interface = plt.Circle((c, c), b4, color='k', fill=False)
        
        return bottom_interface, heating_interface, cooling_interface, upper_interface
    
    def normalize(self, norms):
        if 'density' in norms.keys():
            self._variables['density']*=norms['density']
            self._variables['equilibrium_density']*=norms['density']
        if 'pressure' in norms.keys():
            self._variables['pressure']*=norms['pressure']
            self._variables['equilibrium_pressure']*=norms['pressure']
        if 'velocity' in norms.keys():
            self._variables['velocity'].x*=norms['velocity']
            self._variables['velocity'].y*=norms['velocity'] 
        if 'B_field' in norms.keys():
            self._variables['B_field'].x*=norms['B_field']
            self._variables['B_field'].y*=norms['B_field'] 
        
        
    def mean_profile(self, quantity):
        if quantity in self._variables.keys():
            q = self._variables[quantity]
        else:
            q = self._derived_variables[quantity]
            
        Nx = 2**self.levelmax
        x,y = np.meshgrid(np.arange(Nx), np.arange(Nx))
        x = x-Nx//2
        y = y-Nx//2
        r = np.sqrt(x**2 + y**2)
        r = np.array(r, dtype=int)
        q_mean = []
        for i in np.arange(np.max(r)+1):
            q_mean.append(np.mean(q[r==i]))

        q_mean = np.array(q_mean)
        return q_mean

class Vector:
    """
    A simple vector class
    """
    def __init__(self, l):
        """
        Constructor of the class

        Parameters:
        -----------
            l : list of arrays
                A list containing the x,y,(z) arrays of the vector
        """
        self.ndim = len(l)
        if self.ndim < 2 or self.ndim > 3:
            raise ValueError('The vector must have 2 or 3 components')
        self.x = l[0]
        self.y = l[1]
        if self.ndim == 3:
            self.z = l[2]
        self._components = l
    # ---------------

    @property
    def norm(self):
        n = np.array([i**2 for i in self._components])
        n = np.sqrt(np.sum(n, axis=0))
        return n
    # ------------------
# ----------------------


def readRamsesOutput(folder, n, var_list):
    """
    A routine that uses the Osyris module to read a Ramses binary 
    output and put it in a Ramses class
    """
    # Read data with osyris
    data = osyris.Dataset(n, path=folder).load()
    hydro = data['hydro']
    amr = data['amr']
    grav = data['grav']
    meta = data.meta
    # Prepare grid 
    grid_position = amr['position']/amr['dx']
    grid_position = grid_position - np.min(grid_position)
    # Prepare dictionary
    grid_data = {}
    # Check variable names
    variables_hydro = ['density', 'pressure', 'equilibrium_density', 'equilibrium_pressure', 'velocity', 'B_field']
    if not all(var in variables_hydro for var in var_list):
        raise ValueError('Variable in list non found')
    # Loop over variables and store the array
    for var_name in var_list:
        var_amr = hydro[var_name]
        # Check if it is a vector or a scalar and read it to a grid
        if isinstance(var_amr, osyris.core.vector.Vector):
            grid_var = grid_vector(var_amr, grid_position, amr['level'])
        elif isinstance(var_amr, osyris.core.array.Array):
            grid_var = grid_scalar(var_amr, grid_position, amr['level'])
        grid_data[var_name] = grid_var
    del data
    return grid_data, meta
# ------------------------


def fill_2D(scalar, x, y, level):
    """ 
    Fill a 2D array with AMR array and information 
    
    Parameters:
    -----------
        scalar : osyris.array
            The 1D array of the scalar field 
        x : osyris.array
            The 1D array of x coordinates 
        y : osyris.array
            The 1D array of y coordinates
        level : osyris.array
            The 1D array of amr levels
        
    Return:
    -------
        grid_var : array
            The 2D map of the scalar
    """
    lmax = np.max(level.values)
    # Prepare array
    grid_var = np.zeros((2**lmax, 2**lmax))
    # loop over positions
    for i in np.arange(len(scalar.values)):
        xi = int(x[i])
        yi = int(y[i])
        grid_var[xi,yi] = scalar.values[i]
        # Fill cubes if not at the max level
        if level.values[i] != lmax:
            for l1 in np.arange(2**(lmax-level.values[i])):
                for l2 in np.arange(2**(lmax-level.values[i])):
                    offset =  2**(level.values[i]-1) - 2**(lmax-1)
                    xi = int((x[i]+offset)*2**(lmax-level.values[i]) + l1)
                    yi = int((y[i]+offset)*2**(lmax-level.values[i]) + l2)
                    grid_var[xi,yi] = scalar.values[i]

    return grid_var
# -----------------

def fill_3D(scalar, x, y, z, level):
    """ 
    Fill a 3D array with AMR array and information 
    
    Parameters:
    -----------
        scalar : osyris.array
            The 1D array of the scalar field 
        x : osyris.array
            The 1D array of x coordinates 
        y : osyris.array
            The 1D array of y coordinates
        z : osyris.array
            The 1D array of z coordinates
        level : osyris.array
            The 1D array of amr levels
        
    Return:
    -------
        grid_var : array
            The 3D map of the scalar
    """
    lmax = np.max(level.values)
    # Prepare empty array
    grid_var = np.zeros((2**lmax, 2**lmax, 2**lmax))
    # loop over positions
    for i in np.arange(len(scalar.values)):
        xi = int(x[i])
        yi = int(y[i])
        zi = int(z[i])
        grid_var[xi,yi,zi] = scalar.values[i]
        # Fill cubes if not at the max level
        if level.values[i] != lmax:
            for l1 in np.arange(2**(lmax-level.values[i])):
                for l2 in np.arange(2**(lmax-level.values[i])):
                    for l3 in np.arange(2**(lmax-level.values[i])):
                        offset =  2**(level.values[i]-1) - 2**(lmax-1)
                        xi = int((x[i]+offset)*2**(lmax-level.values[i]) + l1)
                        yi = int((y[i]+offset)*2**(lmax-level.values[i]) + l2)
                        zi = int((z[i]+offset)*2**(lmax-level.values[i]) + l3)
                        grid_var[xi,yi,zi] = scalar.values[i]

    return grid_var
# -----------------

def grid_scalar(scalar, grid, level):
    """ 
    Puts a 1D array into a 2D/3D array with size 2**l_max 
    according to positions given by grid and levels
    
    Parameters:
    -----------
        scalar : osyris.array
            The 1D array of the scalar field 
        grid : osyris.vector
            The array of amr coordinates 
        level : osyris.array
            The 1D array of amr levels
        
    Return:
    -------
        grid_var : array
            The 2D/3D map of the scalar
    """
    # Gather variables
    ndim = grid.nvec
    x = np.round(grid.x.values)
    y = np.round(grid.y.values)
    # Fill grids
    if ndim == 2:
        grid_var = fill_2D(scalar, x, y, level)
    elif ndim == 3:
        z = np.round(grid.z.values)
        grid_var = fill_3D(scalar, x, y, z, level)
    return grid_var
# -----------------

def grid_vector(vector, grid, level):
    """ 
    Puts a vector into a 2D/3D array with size 2**l_max
    according to positions given by grid and levels
    
    Parameters:
    -----------
        vector : osyris.vector
            The osyris vector of the quantity 
        grid : osyris.vector
            The array of amr coordinates 
        level : osyris.array
            The 1D array of amr levels
        
    Return:
    -------
        grid_vector : Vector
            The vector with 2D/3D maps
    """
    
    ndim = grid.nvec
    # Extract components of the vector
    components = [vector.x, vector.y]
    if ndim == 3:
        components.append(vector.z)
    # Each component is a scalar, so it can be filled with the previous function
    grid_vector = []
    for component in components:
        grid_component = grid_scalar(component, grid, level)
        grid_vector.append(grid_component)
    grid_vector = Vector(grid_vector)
    return grid_vector
# --------------------

# def simple_2D_plot(q, figsize=(5,5), vmax=None, vmin=None, title=None, cmap=None):
#     """ Simple plotting routine for the quantity q """
    
#     fig, ax = plt.subplots(figsize=figsize)
#     if vmax is None:
#         vmax = np.max(q)
#     if vmin is None:
#         vmin = np.min(q)
#     if cmap is None:
#         cmap = 'nipy_spectral'
#     cs = ax.imshow(q.T, 
#                    origin='lower',
#                    vmax=vmax,
#                    vmin=vmin,
#                    cmap=cmap)
#     divider = make_axes_locatable(ax)
#     cax = divider.new_vertical(size = '3%', pad = 0.03)
#     fig.add_axes(cax)
#     fig.colorbar(cs, cax = cax, orientation = 'horizontal', label=title)
#     cax.xaxis.set_ticks_position('top')
#     cax.xaxis.set_label_position('top')
    
#     return fig, ax
# # ----------------


# def simple_radial_plot(q, figsize=(5,5), xlim=None, ylim=None, yscale=None):
#     """ Simple plotting routine for the quantity q """
    
#     fig, ax = plt.subplots(figsize=figsize)
#     ax.plot(q, color='g')
#     if xlim:
#         ax.set_xlim(xlim)
#     if ylim:
#         ax.set_ylim(ylim)
#     if yscale:
#         ax.set_yscale(yscale)
    
#     r_center_1 = 4.86990e2
#     r_center_2 = 6.60915e2
#     heating_dx = 20
#     dx = 1.39140e3/2**7
    
#     ax.axvline(r_center_1/dx, color='k')
#     ax.axvline((r_center_1+heating_dx)/dx, color='r')
#     ax.axvline((r_center_2-heating_dx)/dx, color='b')
#     ax.axvline(r_center_2/dx, color='k')
    
#     return fig, ax
# # ----------------

# def compute_energy(output):
#     d = output.d
#     p = output.p
#     v2 = output.v.norm**2
#     B2 = output.B.norm**2
#     e = p/(output.gamma - 1.0) + 0.5*d*v2 + 0.5*B2
#     return e

# def readRamsesOutput(folder, n, vars):
#     """
#     A routine that uses the Osyris module to read a Ramses binary 
#     output and put it in a Ramses class
#     """
#     # Read data with osyris
#     mydata = osyris.Dataset(n, path=folder).load()
#     # Get important variables     
#     ndim = mydata.info.get('ndim')
#     lmax = mydata.info.get('levelmax')
#     dx = mydata.dx.values
#     level = mydata.level.values
#     # Get coordinates
#     x = mydata.x.values
#     y = mydata.y.values
#     if ndim==3:
#         z = mydata.z.values

#     # Compute grid coordinates
#     x = x/dx
#     x = x - np.min(x)
#     x = np.round(x)

#     y = y/dx
#     y = y - np.min(y)
#     y = np.round(y)
    
#     if ndim==3:
#         z = z/dx
#         z = z - np.min(z)
#         z = np.round(z)

#     mhd = False
#     if mydata.info.get('nvar_hydro')>=11:
#         mhd = True

#     if ndim==2:
#         rho = np.zeros((2**lmax,2**lmax))
#         p = np.zeros((2**lmax, 2**lmax))
#         u = np.zeros((2**lmax, 2**lmax))
#         v = np.zeros((2**lmax, 2**lmax))
#         w = np.zeros((2**lmax, 2**lmax)) 
#         deq = np.zeros((2**lmax, 2**lmax))
#         peq = np.zeros((2**lmax, 2**lmax))
#         Bx = np.zeros((2**lmax, 2**lmax))
#         By = np.zeros((2**lmax, 2**lmax))
#         Bz = np.zeros((2**lmax, 2**lmax))

#         N = len(dx)
#         for i in np.arange(N):
#             xx = int(x[i])
#             yy = int(y[i])
#             rho[xx,yy] = mydata.density.values[i]
#             p[xx,yy] = mydata.pressure.values[i]
#             u[xx,yy] = mydata.velocity_x.values[i]
#             v[xx,yy] = mydata.velocity_y.values[i]
#             deq[xx,yy] = mydata.equilibrium_density.values[i]
#             peq[xx,yy] = mydata.equilibrium_pressure.values[i]
#             if mhd:
#                 Bx[xx,yy] = 0.5*(mydata.B_x_left.values[i] + mydata.B_x_right.values[i]) 
#                 By[xx,yy] = 0.5*(mydata.B_y_left.values[i] + mydata.B_y_right.values[i])
#                 Bz[xx,yy] = 0.5*(mydata.B_z_left.values[i] + mydata.B_z_right.values[i])
#             if level[i] != lmax:
#                 for l1 in np.arange(2**(lmax-level[i])):
#                     for l2 in np.arange(2**(lmax-level[i])):
#                         offset =  2**(level[i]-1) - 2**(lmax-1)
#                         xx = int((x[i]+offset)*2**(lmax-level[i]) + l1)
#                         yy = int((y[i]+offset)*2**(lmax-level[i]) + l2)
#                         rho[xx,yy] = mydata.density.values[i]
#                         p[xx,yy] = mydata.pressure.values[i]
#                         u[xx,yy] = mydata.velocity_x.values[i]
#                         v[xx,yy] = mydata.velocity_y.values[i]
#                         deq[xx,yy] = mydata.equilibrium_density.values[i]
#                         peq[xx,yy] = mydata.equilibrium_pressure.values[i]
#                         if mhd:
#                             Bx[xx,yy] = 0.5*(mydata.B_x_left.values[i] + mydata.B_x_right.values[i]) 
#                             By[xx,yy] = 0.5*(mydata.B_y_left.values[i] + mydata.B_y_right.values[i])
#                             Bz[xx,yy] = 0.5*(mydata.B_z_left.values[i] + mydata.B_z_right.values[i])

#     elif ndim==3:
#         rho = np.zeros((2**lmax, 2**lmax, 2**lmax))
#         p = np.zeros((2**lmax, 2**lmax, 2**lmax))
#         u = np.zeros((2**lmax, 2**lmax, 2**lmax))
#         v = np.zeros((2**lmax, 2**lmax, 2**lmax))
#         w = np.zeros((2**lmax, 2**lmax, 2**lmax))
#         deq = np.zeros((2**lmax, 2**lmax, 2**lmax))
#         peq = np.zeros((2**lmax, 2**lmax, 2**lmax))
#         Bx = np.zeros((2**lmax, 2**lmax, 2**lmax))      
#         By = np.zeros((2**lmax, 2**lmax, 2**lmax))
#         Bz = np.zeros((2**lmax, 2**lmax, 2**lmax))

#         N = len(dx)
#         for i in np.arange(N):

#             xx = int(x[i])
#             yy = int(y[i])
#             zz = int(z[i])
#             rho[xx,yy,zz] = mydata.density.values[i]
#             p[xx,yy,zz] = mydata.pressure.values[i]
#             u[xx,yy,zz] = mydata.velocity_x.values[i]
#             v[xx,yy,zz] = mydata.velocity_y.values[i]
#             w[xx,yy,zz] = mydata.velocity_z.values[i]
#             deq[xx,yy,zz] = mydata.equilibrium_density.values[i]
#             peq[xx,yy,zz] = mydata.equilibrium_pressure.values[i]
#             if mhd:
#                 Bx[xx,yy,zz] = 0.5*(mydata.B_x_left.values[i] + mydata.B_x_right.values[i]) 
#                 By[xx,yy,zz] = 0.5*(mydata.B_y_left.values[i] + mydata.B_y_right.values[i])
#                 Bz[xx,yy,zz] = 0.5*(mydata.B_z_left.values[i] + mydata.B_z_right.values[i])
#             if level[i] != lmax:
#                 for l1 in np.arange(2**(lmax-level[i])):
#                     for l2 in np.arange(2**(lmax-level[i])):
#                         for l3 in np.arange(2**(lmax-level[i])):
#                             offset =  2**(level[i]-1) - 2**(lmax-1)
#                             xx = int((x[i]+offset)*2**(lmax-level[i]) + l1)
#                             yy = int((y[i]+offset)*2**(lmax-level[i]) + l2)
#                             zz = int((z[i]+offset)*2**(lmax-level[i]) + l3)
#                             rho[xx,yy,zz] = mydata.density.values[i]
#                             p[xx,yy,zz] = mydata.pressure.values[i]
#                             u[xx,yy,zz] = mydata.velocity_x.values[i]
#                             v[xx,yy,zz] = mydata.velocity_y.values[i]
#                             w[xx,yy,zz] = mydata.velocity_z.values[i]
#                             deq[xx,yy,zz] = mydata.equilibrium_density.values[i]
#                             peq[xx,yy,zz] = mydata.equilibrium_pressure.values[i]
#                             if mhd:
#                                 Bx[xx,yy,zz] = 0.5*(mydata.B_x_left.values[i] + mydata.B_x_right.values[i]) 
#                                 By[xx,yy,zz] = 0.5*(mydata.B_y_left.values[i] + mydata.B_y_right.values[i])
#                                 Bz[xx,yy,zz] = 0.5*(mydata.B_z_left.values[i] + mydata.B_z_right.values[i])
    
#     vel = Vector([u,v,w])
#     B = Vector([Bx,By,Bz])
#     data = {'d':rho, 
#             'v':vel, 
#             'p':p, 
#             'B':B,
#             'deq':deq,
#             'peq':peq,
#             'info':mydata.info}
    
#     del mydata

#     return data
# -------------