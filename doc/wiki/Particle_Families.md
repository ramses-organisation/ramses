

New particle types [WORK IN PROGRESS]
=====================================

We are implementing particle types as two extra arrays (grouped in a derived type `part_t`, stored in the variable `typep`). The first array contains the actual family of the particle (`typep(:)%family`), while the second is the tag of a specific subtype such as pop II and pop III for stars (`typep(:)%tag`).

As of now, the code support 127 different families (+ their 127 attached tracers) and 256 different tags (each is a 1 byte array).
We are trying to foresee most of plausible uses, but me might have overlooked some. Feel free to suggest any other!

For tracers, we have in mind the Monte Carlo tracers of Corentin Cadiou, with "attached" tracers being attached to particles and have a `family` equal to `-family` of the particle, and "detached" tracers (advected with the gas) with type 0.

|   Particle type    |  `family`  |
| ------------------ | ---------- |
|         DM         |     1      |
|        Stars       |     2      |
| Clouds (for sinks) |     3      |
|       Debris       |     4      |
|         TBD        |  5 to 127  |
|  Attached tracers  | -127 to -1 |
|  Detached tracers  |     0      |

Note that **this is not a backward compatible change**: the particle output will be changed, and more importantly, the selection of particles e.g. for SN feedback, BH accretion, etc.

**If particle has no type defined, then by default undef=127 is assigned.**

## Note for users ##
A new output file `part_file_descriptor.txt` has been introduced. It contains information about the fields in the `part_XXXXX.outYYYYY` files.
The `header_XXXXX.txt` has also been updated so that it counts the number of particles and output them. Please note that the fields provided in this file *should* be the same as the ones in `part_file_descriptor.txt`, but you should trust the `part_file_descriptor.txt` file for this.

## Note for developers ##
In `pm_commons.f90`, there is now a group of functions designed to check the particle type, with names like `is_star` or `is_tracer`. This is the preferred method to check for particle type and should be used instead of matching the id, mass or formation time of the particles. 

For implementing a new subtype of particle(e.g. popIII and popII stars), we recommend that users create a new tag `TAG_<tagname>` in pm_commons and a function that matches the tag (similarly to e.g. `is_star`).
Here is an example:
```fortran
module pm_commons
    [...]
    ! Customize here for particle tags within particle types (e.g. different kind of stars).
    [...]
    integer(1) :: TAG_POPIII=1, TAG_POPII=2
    [...]
contains
    [...]

    logical pure function is_star_popIII(typep)
        type(part_t), intent(in) :: typep
        is_star_popIII = (typep%family == FAM_STAR) .and. (typep%family == TAG_POPIII)
    end function

    logical pure function is_star_popII(typep)
        type(part_t), intent(in) :: typep
        is_star_popII = (typep%family == FAM_STAR) .and. (typep%family == TAG_POPII)
    end function    
end module pm_commons
```

For implementing any new kind  particle (e.g. pop III stars), we recommend that users create similar functions, such as `is_popIII_star`.

Authors: Corentin Cadiou, Maxime Trebitsch, Rebekka Bieri, Hoseung Choi, Owain Snaith