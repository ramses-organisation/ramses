---
orphan: true
---

# 1.2 Compiling the code

You need to go first to the `bin/`directory:
```
$ cd trunk/ramses/bin
$ ls -F
Makefile	Makefile.rt
```

We will use the first `Makefile` to compile the code.  The first thing
to do is to edit the `Makefile` and modify the two variables `F90` and
`FFLAGS`.   Several  examples   corresponding  to   different  Fortran
compilers are given. The default values are:

```
F90 = gfortran -O3 -frecord-marker=4 -fbacktrace -ffree-line-length-none
FFLAGS = -x f95-cpp-input -DWITHOUTMPI $(DEFINES)
```

The first variable is obviously the command used to invoke the Fortran
compiler.  In this  case, this is the [GNU  Fortran compiler][1].  The
second variable  contains Fortran  compilation flags  and preprocessor
directives.  The first directive, `-DWITHOUTMPI`, switches off all MPI
routines.  On  the other hand,  if you  don't use this  directive, the
code must  be linked to  the MPI library.  We will discuss  this point
later.

[1]: http://gcc.gnu.org/fortran

Other preprocessor directives are defined in variable `DEFINES` in the 
`Makefile`:

```
# Compilation time parameters
NVECTOR = 64
NDIM = 3
NPRE = 8
NVAR = 8
SOLVER = hydro
PATCH =
EXEC = ramses
DEFINES = -DNVECTOR=$(NVECTOR) -DNDIM=$(NDIM) -DNPRE=$(NPRE) -DNVAR=$(NVAR) -DSOLVER$(SOLVER)
```

These   additional    directives   are   called    _Compilation   Time
Parameters_. They should be defined  in the `Makefile` and the code must 
be recompiled entirely using:

```
$ make clean
$ make
```

We list now the definitions of these parameters.

`NVECTOR=64`
: This parameter is used   to    set   the   vector    size   for 
computation-intensive   operations.   It    must   be   determined
experimentally on each new hardware.

`NPRE=4`
: This parameter sets the precision of the floating point operations. 
`NPRE=4`stands for single precision arithmetics, while `NPRE=8` is 
for double precision.

`NENER=0`
: This parameter sets the number of energy variables used in the hydro or mhd solver.

`NDIM=3`
: This parameter sets the dimensionality of the problem.
The value `NDIM=1`is for 1D, plan-parallel flows. `NDIM=2` and 
`NDIM=3` are resp. for 2D and 3D flows.

`SOLVER=hydro`
: This parameter selects the type of hyperbolic solver used. 
Possible values are: `hydro` for the adiabatic Euler equations,
`mhd`, the Constrained Transport scheme for the ideal MHD equations.
and `rhd` for relativistic hydro.

`NVAR=8`
: This parameters defines the number of variables in the hyperbolic solver. 
For `SOLVER=hydro`, `NVAR>=NDIM+2`. 
For `SOLVER=mhd`, `NVAR>=8` and for `SOLVER=rhd`, one has `NVAR>=5`. 

Our goal is now to compile the code for a simple one-dimensional problem.
You need to modify the `Makefile` so that:

```
NDIM=1
SOLVER=hydro
NVAR=3
```

Then type:

```
$ make
```

If everything goes well, all source files will be compiled and 
linked into an executable called `ramses1d`.

## 1.2.1 Additional compilation preprocessor flags

`-DTSC`
: This parameter sets the triangular shape cloud approximation; only works with `NDIM=3`

`-DOUTPUT_PARTICLE_POTENTIAL`
: This parameter forces the code to output particle potentials at snapshots

`-DQUADHILBERT`
: This parameter sets longer Hilbert curve necessary if `levelmax>19`

`-DLONGINT`
: This parameter switches to long ints (necessary when one has lots of particles)

`-DNOSYSTEM`
: This parameter handles operating system commands

## [Next step: executing the test case !](Start3)