
# Testing 

## 1. Running the automatic test suite

To run the automatic tests, navigate to the [tests](https://github.com/ramses-organisation/ramses/tree/stable/tests) directory, and run the `run_test_suite.sh` script:
```
>$ cd tests
>$ ./run_test_suite.sh
```
The tests will begin and the output should look like:
```
############################################
#   Running RAMSES automatic test suite    #
############################################
Will perform the following tests:
 [ 1] hydro/implosion
 [ 2] hydro/sod-tube
 [ 3] mhd/imhd-tube
 [ 4] mhd/orszag-tang
 [ 5] rt/stromgren2d
 [ 6] sink/smbh-bondi
--------------------------------------------
Test 1/6: hydro/implosion
Cleanup
Compiling source
```
and so on.
Once the tests have completed, a report is generated in a `.pdf` file named `test_results.pdf`, alongside a log file `test_suite.log`.

### Options

- Run the suite in parallel (on 4 cpus):
```
./run_test_suite.sh -p 4
```

- Do not delete results data:
```
./run_test_suite.sh -d
```

- Run in verbose mode:
```
./run_test_suite.sh -v
```

- Select individual tests (for tests 3 to 5, and 10):
```
./run_test_suite.sh -t 3-5,10
```

- Run all tests in `mhd` directory:
```
./run_test_suite.sh -t mhd
```

## 2. Creating a new test

The following steps describe how to add a new test to the test suite. In this example, the test will be named `sedov-3d`.

The first step is to create a new directory `sedov-3d` in one of the `hydro`, `mhd`, `rt`, or `sinks` directories. No need to modify the `run_test_suite.sh` script, the new test will automatically be picked up and added to the list. We will choose to place it inside the `hydro` directory. Please use hyphens (`-`) in your test names instead of underscores (`_`) as `latex` does not like underscores.

```
>$ cd hydro
>$ mkdir sedov-3d
```

**Note**: use one directory per test. If you want to run a 2D and a 3D sedov test, create separate `sedov-2d` and `sedov-3d` directories.

In that directory, you will need:

- A `config.txt` file: usually just contains the Makefile flags, e.g. `FLAGS: NDIM=3 PATCH= SOLVER=hydro`

- A namelist: `sedov-3d.nml` (the name needs to be the same as the test directory)

- A file for plotting and checking the solution against a reference: `plot-sedov-3d.py`. It is advised to copy a file from the other directories to see how to write this. **Note that this file needs to contain at least one call to `visu_ramses.check_solution(data["data"], 'sedov-3d')`**.

- A reference solution: `sedov-3d-ref.txt`. To create it, run your test and once the final output (number 2 in this case) has been created, do the following:
```
import visu_ramses
data = visu_ramses.load_snapshot(2)
visu_ramses.check_solution(data["data"], 'sedov-3d', overwrite=True)
```

- A `Readme.md` containing a short description of the test

Optional files:

- `condinit.f90`: you can have your own initial setup if it's not entirely definable in a namelist. **REMEMBER** to set the correct `PATCH` in the `config.txt` file! (e.g. `PATCH=../tests/hydro/sedov-3d`)

- `before-test.sh`: if this file is present in the test directory, it will be run before the test begins (useful for e.g. creating symbolic links to libraries...)

- `after-test.sh`: if this file is present in the test directory, it will be run after the test begins (useful for e.g. cleaning up symbolic links to libraries...)

### Tuning tolerances for solution verification

By default, relative differences between the sums of all the variables inside all leaf cells in the domain and the reference solution cannot exceed `3.0e-13`.
Sometimes, some variables are more volatile than others when running simulations on different numbers of CPUs, and this limit is too low, leading to false failed tests.
The `check_solution` method in the `visu/visu_ramses.py` module can be tuned to work for your test using the following options:

- `tolerance`: a dictionary listing the allowed relative difference between the sum over all leaf cells and reference value. The default is `{"all":3.0e-13}`. To make the check on `density` less restrictive, use for instance `tolerance={"density":1.0e-10}`.

- `threshold`: relative value below which a vector component is set to zero. Default is `2.0e-14`.

- `norm_min`: minimum value for the norm of a vector, to protect against null vectors. Default is `1.0e-30`.

- `min_variance`: if the data differs by less than this value from the average value, it is set to the average. Default is `1.0e-14`.


## 3. Creating a new group of tests

If your test does not fall under the categories already present in the `tests` directory (`hydro`, `mhd`, `rt`, `sinks`), you can create a new directory and put your tests in there. You will then have to edit the `run_test_suite.sh` file to ensure your new tests will be picked up.

Say you want to create 3 new tests, `sedov-1d`, `sedov-2d`, and `sedov-3d` inside a new `sedov` directory, you have to find the line describing the list of directories to be scanned at the top of the `run_test_suite.sh` file:

```
# List of directories to scan
testlist="hydro,mhd,rt,sink";
```
and add your new directory separated from the previous one by a comma, i.e.
```
# List of directories to scan
testlist="hydro,mhd,rt,sink,sedov";
```