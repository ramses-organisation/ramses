---
orphan: true
---

# Executing the test case

To test the compilation, you need to execute a simple test case. Go up
one level and type the following command:

```
$ cd trunk/ramses
$ bin/ramses1d namelist/tube1d.nml
```

The first part of the command is the executable we have just compiled.
The second part, the only command line argument, is an input file 
containing the _run time parameters_. Several examples of such parameter 
files are stored in the `namelist/` directory. The namelist file we 
have just used `tube1d.nml` is the Sod test, a simple shock tube simulation 
in 1D. For comparison, we now show the last 14 lines of the standard output: 

```
 Mesh structure
 Level  1 has          1 grids (       1,       1,       1,)
 Level  2 has          2 grids (       2,       2,       2,)
 Level  3 has          4 grids (       4,       4,       4,)
 Level  4 has          8 grids (       8,       8,       8,)
 Level  5 has         16 grids (      16,      16,      16,)
 Level  6 has         27 grids (      27,      27,      27,)
 Level  7 has         37 grids (      37,      37,      37,)
 Level  8 has         17 grids (      17,      17,      17,)
 Level  9 has         16 grids (      16,      16,      16,)
 Level 10 has         13 grids (      13,      13,      13,)
 Main step=    43 mcons=-1.97E-16 econs= 1.61E-16 epot= 0.00E+00 ekin= 1.38E+00
 Fine step=   688 t= 2.45047E-01 dt= 3.561E-04 a= 1.000E+00 mem= 7.6%
 Run completed
```

If  your  execution looks  similar,  it  means your  installation  was
successfull.   Users are  encouraged to  redirect the  standard output
into  a _log  file_. This  log  file contains  all simulation  control
variables, as well as output variables, but for 1D simulations only.

```
$ bin/ramses1d namelist/tube1d.nml > tube.log
```

## [Next step: Reading the log file !](Start4)




