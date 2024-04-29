---
orphan: true
---

# 1.1 Obtaining the package

The package can be downloaded from the GitHub repository using git:

```
$ git clone https://github.com/ramses-organisation/ramses
```

This will create a new repository called `ramses/`
In this directory, you will see:

```
$ ls -F
README		bin/		mhd/		patch/		rt/
amr/		doc/		namelist/	pm/		utils/
aton/		hydro/		pario/		poisson/
```

Each directory  contains a set of  files with a given  common purpose.
For example, `amr/`  contains all Fortran 90 routines  dealing with the
AMR grid  management and  MPI communications, while  `hydro/` obviously
contains  all Fortran  90  routines dealing  with hydrodynamics.   The
first directory you are interested in is the `bin/` directory, in which
the code will be compiled.

## [Next step: compiling the code !](Start2)