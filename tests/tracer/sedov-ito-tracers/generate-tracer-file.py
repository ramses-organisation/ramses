import numpy as np
tracer_mass = 1e-6
ndim = 2
d = 1
boxlen = 0.5
total_mass = d * boxlen**ndim
n_tracers = int(total_mass / tracer_mass)
n_tracers_per_dim = int(np.ceil(n_tracers**(1/ndim)))

half_spacing = 1 / (2 * n_tracers_per_dim)

xyz_dict = {
    1: np.mgrid[half_spacing:1-half_spacing:n_tracers_per_dim*1j, 0:0:1j, 0:0:1j],
    2: np.mgrid[half_spacing:1-half_spacing:n_tracers_per_dim*1j,
                 half_spacing:1-half_spacing:n_tracers_per_dim*1j, 0:0:1j],
    3: np.mgrid[half_spacing:1-half_spacing:n_tracers_per_dim*1j,
                    half_spacing:1-half_spacing:n_tracers_per_dim*1j,
                    half_spacing:1-half_spacing:n_tracers_per_dim*1j],
}

xyz = xyz_dict[ndim]

list_tracers = xyz.reshape(3,n_tracers).T
print(f"Created a file with {list_tracers.shape[0]} tracers")
np.savetxt('tracers.txt', list_tracers, fmt='%.6f')
