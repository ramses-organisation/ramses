'''
    This script generated initial conditions for turbulent box simulations.

    A random velocity field with specified power spectrum is generated for each velocity component.
    This field is then rescaled to the desired Mach number.
    The density and pressure fields are set to uniform values as specified in the namelist.
'''

from matplotlib import pyplot as plt
import numpy as np
import f90nml
import os
# For the current script to find the grafic module in utils,
# it needs to be run from the ramses root folder
from utils.py.write_grafic import write_grafic_file


def make_3dfield(ncells: int, *, powerlaw: int, randomseed: int):
    # Author: Corentin Cadiou
    """Generate a 3D random field with a power-law power spectrum.

    Parameters
    ----------
    ncells : int
        The number of cells along each dimension of the 3D field.
    powerlaw : int
        The exponent of the power-law power spectrum. For example, a value of
        2 would correspond to a power spectrum that scales as k^-2.
    randomseed : int
        The seed for the random number generator, ensuring reproducibility of
        the generated field.

    Returns
    -------
    field : numpy.ndarray
        A 3D array of shape (ncells, ncells, ncells) containing the generated
        random field values.
    """
    shape = (ncells, ncells, ncells)
    shape_fourier = (ncells, ncells, ncells // 2 + 1)

    # Generate the Fourier space whitenoise
    generator = np.random.default_rng(randomseed)
    args = {"scale": 1 / np.sqrt(2), "size": shape_fourier}
    whitenoise = generator.normal(**args) + 1j * generator.normal(**args)

    # Create the k-space grid
    kx = np.fft.fftfreq(ncells) * ncells
    ky = np.fft.fftfreq(ncells) * ncells
    kz = np.fft.rfftfreq(ncells) * ncells
    kx, ky, kz = np.meshgrid(kx, ky, kz, indexing="ij")
    k = np.sqrt(kx**2 + ky**2 + kz**2)

    # Compute power spectrum
    Pk = np.zeros_like(k)
    Pk[k > 0] = k[k > 0] ** -powerlaw

    # Apply the power spectrum to the whitenoise
    field_fourier = whitenoise * np.sqrt(Pk)

    return np.fft.irfftn(field_fourier, s=shape)

''' ---- Parameters ---- '''

# Mach number
mach3d = 8

# powerlaw index
# empirically Gaussian fluctuations power law (modelled by the fBm) is about 3.4 - 3.6 (J.-F. Robitaille)
powerlaw=3.5

# random number seeds
seed1 = 824329
seed2 = 129862
seed3 = 786276

''' ---- Get further information from the namelist ----'''

mypath = os.path.dirname(__file__)
nml = f90nml.read(mypath + '/decaying-turbulence.nml')

# refinement level (assume uniform)
level = nml['amr_params']['levelmin']
ncells = 2**level
num_cells = ncells**3

# box length
size = nml['amr_params']['boxlen']

# temperature and molecular weight
T = nml['cooling_params']['T_eos']
mu = nml['cooling_params']['mu_gas']

# code units for velocity
scale_v = nml['units_params']['units_length']/nml['units_params']['units_time']

# calculate sound speed in code units from temperature
kb_cgs = 1.38064852e-16 # cm2 g s-2 K-1   # Boltzman constant
mH_cgs = 1.6737236e-24 # g                # hydrogen mass
sound_speed = (kb_cgs*T/(mu*mH_cgs))**0.5 / scale_v # code units

''' ---- Generate ICs ---- '''

# generate fields for the three components of the velocity
vx = make_3dfield(ncells, powerlaw=powerlaw, randomseed=seed1)
vy = make_3dfield(ncells, powerlaw=powerlaw, randomseed=seed2)
vz = make_3dfield(ncells, powerlaw=powerlaw, randomseed=seed3)

# correct for possible bulk motion
vx = vx - vx.mean()
vy = vy - vy.mean()
vz = vz - vz.mean()

# calculate current velocity dispersion
sigma_x = np.sqrt(np.sum(vx**2)/num_cells)
sigma_y = np.sqrt(np.sum(vy**2)/num_cells)
sigma_z = np.sqrt(np.sum(vz**2)/num_cells)

# rescale velocity to Mach requested Mach number
sigma_new = (mach3d/np.sqrt(3))*sound_speed # code units
vx = vx * sigma_new/sigma_x
vy = vy * sigma_new/sigma_y
vz = vz * sigma_new/sigma_z

check = False
if check:
    # check new velocity dispersion
    sigma_x_new = np.sqrt(np.sum(vx**2)/num_cells)
    sigma_y_new = np.sqrt(np.sum(vy**2)/num_cells)
    sigma_z_new = np.sqrt(np.sum(vz**2)/num_cells)
    print(sigma_x_new,sigma_y_new,sigma_z_new)

    # make plots to verify the fields look OK
    for seed, v in zip([seed1, seed2, seed3], [vx,vy,vz]):
        fig, ax = plt.subplots(nrows=1, ncols=3, figsize=[10, 3])
        ax[0].imshow(v.mean(0), origin='lower')
        ax[1].imshow(v.mean(1), origin='lower')
        ax[2].imshow(v.mean(2), origin='lower')
        plt.tight_layout()
        plt.savefig('seed{}_pl{}.png'.format(seed, powerlaw))
        plt.close()


''' ---- Write IC files in a folder in the current test directory ---- '''

ic_dir = mypath + '/ic/'
if not os.path.exists(ic_dir):
    os.makedirs(ic_dir)

# write velocity grafic files
for v, file in zip([vx,vy,vz],['ic_u','ic_v','ic_w']):
    write_grafic_file(ic_dir+file, v, size)

# write uniform density and pressure to make ramses happy
# (it first tries to read ic_d to get the header and exits when this file is not found)
d = np.full((ncells,ncells,ncells),nml['init_params']['d_region'],dtype='f4')
write_grafic_file(ic_dir+'ic_d', d, size)
p = np.full((ncells,ncells,ncells),nml['init_params']['p_region'],dtype='f4')
write_grafic_file(ic_dir+'ic_p', p, size)
