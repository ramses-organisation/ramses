

Particle Unbinding
===========================================================

 Quickstart Guide
-------------------------------------

To activate particle unbinding, you need to set two runtime parameters in the `&RUN_PARAMS` namelist (in addition to whatever you have in there):

```
&RUN_PARAMS
clumpfind=.true.
unbind=.true.
/
```

The code will create `output_XXXXX/unbinding_XXXXX.outYYYYY` files following the same logic as the other particle
output files `output_XXXXX/part_XXXXX.outYYYYY` files.
They will contain the associated clump ID of each particle.
In only works when particles, i.e. dark matter, is present in your simulation.






Important notes
------------------------

 * You can't do particle unbinding without doing clump/halo finding.

 * Some parts of the code (e.g. binning particles in mass profiles of halos) rely on consistent floating-point  operations. 
 [The (intel) fortran compiler however doesn't necessarily use value-safe optimisations][1], which may lead to **errors resulting in warnings**,
 but the code doesn't crash.  The error should be small  (~1e-16), and you may choose to ignore it. 
 Otherwise, you might want to compile the code with the  `-fp-model precise` flag for intel, or the 
 appropriate flag for the compiler you'd like to use.

 * The previously available `unbinding_formatted_output` namelist parameter has been removed. Make sure you remove it from your namelists if you used it!

 * For anything else regarding particle unbinding, feel free to contact Mladen Ivkovic (mladen.ivkovic [at] hotmail DOT com)












 Documentation
-------------------------------------

### What it does

The purpose of particle unbinding is to identify unbound particles in clumps as identified by the clumpfinder and pass 
them on to the parent clumps, until the halo-namegiver clumps are reached (where there are no more parent
structures to pass the particles on to.)

It will write unformatted output in `output_XXXXX/unbinding_XXXXX.outYYYYY` files the same way it is done for
any other backup files in `ramses`, containing the assigned clump IDs of every particle after unbinding. The
clump IDs correspond to the clump IDs as used in the `halo_XXXXX.txtYYYYY` and `clump_XXXXX.txtYYYYY` files. 
If a particle has clump ID 0, it wasn't found to be in any clump.



### How it works

First all particles that are in a clump are gathered and assigned the corresponding clump ID. 
Simultaneously, linked lists of these particles are created for each clump that is not a halo-namegiver, i.e. 
that is not a clump whose ID will be the halo ID. Such clumps are never merged into another clump, but may have
arbitrarily many clumps merged into them. Then clump properties such as centre of mass, the mass profile and bulk 
velocity are acquired using the linked lists. Starting with the lowest clump level, particles are checked if they
are bound to the clump they are assigned to. By default, this unbinding is done iteratively: The clump properties
are recomputed using only the remaining particles, and then all particles checked again. Furthermore, by default
particles are not allowed to leave the boundaries of the clump they're assigned to in order to be considered
bound; For this, the potential at the closest saddle from the clump's centre of mass to a neighbouring clump is
subtracted from the particle's energy. (This default behaviour can however be changed with the namelist 
parameters `saddle_pot` and `iter_properties`). Also note that the bulk velocity and centre of mass of a clump are
recovered using only the bound particles per iteration, the mass profile however always uses all included particles,
including the substructure particles.
The iteration per clump level stops when the bulk density of each clump of that level has converged, i.e.
`v_clump_old/v_clump_new < conv_limit` (or when a maximal number of iterations is reached). Particles that are
found to be not bound are passed to the parent structure for examination, provided such a structure exists, and the 
iterations repeated for the next clump level, provided there are clumps of a higher level.
There are two possibilities implemented to define the clump's centre. By default, the centre will be the position of the
associated density peak. However using the `-DUNBINDINGCOM` preprocessing flag, you can make the centre be the 
centre of mass, which will also be determined iteratively like the bulk velocity.

More details can be found [here][2].



### Namelist Parameters for unbinding

Can be set in the `UNBINDING_PARAMS` block




|   Name                        |   default                 |   type   |   function                                        |
|-------------------------------|---------------------------|----------|---------------------------------------------------|
| `particlebased_clump_output=` | `.false.`                 | logical  | write resulting clump properties based on particles after unbinding, not default cell-based properties  |
| `nmassbins=`                  | `50`                      | integer  | Number of bins for the mass binning of the cumulative mass profile. Any integer > 1.         |
| `logbins=`                    | `.true.`                  | logical  | use logarithmic binning distances for cumulative mass profiles (and gravitational potential of clumps). If false, the code  will use linear binning distances. |
| `saddle_pot=`                 | `.true.`                  | logical  | Take neighbouring structures into account; Cut  potential off at closest saddle. |
| `iter_properties= `           | `.true.`                  | logical  | whether to unbind multiple times with updated clump properties determined by earlier unbindings |
| `conv_limit =`                | `0.01`                    | real     | convergence limit. If `v_clump_old/v_clump_new < conv_limit`, stop iterating for this clump. (only used when `iter_properties=.true.`)                         |
| `repeat_max =`                | `100`                     | integer  | maximal number of loops per level for iterative unbinding (in case a clump doesn't converge) (shouldn't happen) (only used when `iter_properties=.true.`)         |











[1]: https://www.nccs.nasa.gov/images/FloatingPoint_consistency.pdf
[2]: https://drive.google.com/file/d/0B7IyoMUxCr-3V3NFSVFjY1lMbk0
[3]: https://arxiv.org/pdf/1207.6105.pdf 
[4]: https://github.com/ramses-organisation/ramses/blob/stable/doc/wiki/PHEW.md
[5]: https://drive.google.com/open?id=1q0RSMeTIF7gQ7s2DXzYZUSG6cQ1LstkF
