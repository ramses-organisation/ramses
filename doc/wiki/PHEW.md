

# Clumpfinder (PHEW) #

The block named `&CLUMPFIND_PARAMS` contains the parameters related to the built-in RAMSES clump finder (PHEW). The parameters are described only briefly, for more backround information read the [PHEW paper](http://www.comp-astrophys-cosmol.com/content/pdf/s40668-015-0009-7.pdf).

| Variable name, syntax, default value | Fortran type  | Description               |
|:---------------------------- |:------------- |:------------------------- |
| `ivar_clump=1`               | `integer`     | Control which density field is used for clump finding (1: gas density, 0: particle density)
| `density_threshold=-1.d0`              | `float`     | Density threshold for the clump finder (code units)
| `rho_clfind=-1.d0`                     | `float`     | Density threshold for the clump finder (g/cc). Not recommended - use `density_threshold` instead.
| `n_clfind=-1.d0`                       | `float`     | Density threshold for the clump finder (H/cc). Not recommended -use `density_threshold` instead.
| `relevance_threshold=2.d0`             | `float`     | Relevance (peak-to-saddle ratio) threshold for a clump to be considered real (instead of noise). 
| `saddle_threshold=-1.d0`               | `float`     | Saddle density threshold for sub-structure merging (code units). A negative value turns sub-structure merging off (PHEW will only merge noise).
| `mass_threshold=0.d0`                  | `float`     | When set to a value > 0, the properties of those clumps/haloes with `mass > mass_threshold * particle_mass` are written to disk. `particle_mass` is the smallest strictly positive particle mass in the simulation. Setting this parameter does NOT affect the merging of noise/clumps.
| `age_cut_clfind=0.d0`                  | `float`     | When set to a value > 0, only the stars with an age lower than `age_cut_clfind` are used by the clump finder. Stellar and DM particles are used otherwise.