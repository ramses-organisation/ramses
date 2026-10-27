# Testing

## 1. Running the automatic test suite

### Requirements

The test suite compiles and runs RAMSES, then analyses the outputs with Python. You will need:

- a Fortran compiler and, to run in parallel, MPI (`openmpi` for instance). Some tests (`turb`) also need FFTW.
- Python 3 with the packages listed in [`.github/workflows/requirements.txt`](https://github.com/ramses-organisation/ramses/tree/stable/.github/workflows/requirements.txt) (`numpy`, `matplotlib`, `scipy`, `f90nml`, `pyvista`):
```
>$ pip install -r .github/workflows/requirements.txt
```
- `pdflatex` or, if it is not installed, `tectonic`, to compile the PDF report.

### A first try

The full suite contains about forty tests and takes a long time to run. For a first try, run only a few fast 1D tests: a hydro advection test, a hydro shock tube and an MHD shock tube. Navigate to the `tests` directory, and run the `run_test_suite.sh` script with the `-t` option to select the tests:
```
>$ cd tests
>$ ./run_test_suite.sh -t 1,12,16
```
The tests will begin and the output should look like:
```
############################################
#   Running RAMSES automatic test suite    #
############################################
Will perform the following tests:
 [ 1] hydro/advect1d
 [12] hydro/sod-tube-nener
 [16] mhd/imhd-tube-nener
--------------------------------------------
Test 1/3: hydro/advect1d
Cleanup
Compiling source
Running test:
Restart: step 1 ...
Restart: step 2 ...
Restart: step 3 ...
 0h0m1s
Plotting and analysing results
Test hydro/advect1d passed            [ OK ]
--------------------------------------------
Test 2/3: hydro/sod-tube-nener
...
All tests were completed successfully
```
Each test runs in a few seconds; most of the time is spent compiling RAMSES, which is done again for each test. The whole run takes a few minutes.

If the numbers do not match the test names above (because tests were added since this page was written), use the list printed at the start of a run to find them again.

### Running the tests

To run the whole test suite, run `run_test_suite.sh` without options:
```
>$ cd tests
>$ ./run_test_suite.sh
```
The tests are numbered in the order of the directories listed at the top of `run_test_suite.sh` (`hydro`, `mhd`, `poisson`, `rt`, `sink`, `star`, `turb`, `tracer`), then alphabetically within each directory. Use the list printed at the start of a run to find the number of a test.

For each test, the script compiles RAMSES with the flags of that test, runs it, and calls the test's plotting script, which compares the result to a reference solution. The test passes if the plotting script prints `PASSED`.

Once the tests have completed, a report is generated in a `.pdf` file named `test_results.pdf`, alongside a log file `test_suite.log`. The script exits with status 1 if any test failed, so it can be used in scripts and CI.

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

- Run all tests in `mhd` directory (several directories can be given, separated by commas):
```
./run_test_suite.sh -t mhd
```

- Run test suite with coverage (see Section 3, Code coverage):
```
./run_test_suite.sh -s
```

- Control which tests are run through the restart mechanism (see Restart tests below):
```
./run_test_suite.sh            # restart only the tests that opt in
./run_test_suite.sh -r         # restart all tests
./run_test_suite.sh -r all     # same as -r
./run_test_suite.sh -r none    # restart no test
```

- Delete the results of all tests, clean the `bin` directory, and exit without running any test:
```
./run_test_suite.sh -c
```

Options can be combined, e.g. `./run_test_suite.sh -p 2 -v -t hydro,mhd`.

### Restart tests

A restart test checks that a run split in two, with a restart in the middle, gives exactly the same final output as an uninterrupted run. The run is stopped half-way, RAMSES is restarted from the last output, and the final output is compared to the same reference solution as usual.

The splitting is done by `tests/run_with_restart.py`, which edits the namelist before each of the two runs and restores it afterwards. By default, a test is run this way only if its `config.txt` contains the line
```
RESTART: true
```
Most tests opt in. Use `-r` to force the restart on all tests, or `-r none` to turn it off for all tests.

### Continuous integration

On GitHub, the workflow `.github/workflows/tests.yaml` runs the test suite on every pull request, one job per test directory, with
```
bash run_test_suite.sh -p 2 -v -t <directory>
```
The `test_results.pdf` and `test_suite.log` of each job can be downloaded from the artifacts of the workflow run.

## 2. Creating a new test

The following steps describe how to add a new test to the test suite. In this example, the test will be named `sedov-3d`.

The first step is to create a new directory `sedov-3d` in one of the existing category directories. No need to modify the `run_test_suite.sh` script, the new test will automatically be picked up and added to the list. We will choose to place it inside the `hydro` directory. Please use hyphens (`-`) in your test names instead of underscores (`_`) as `latex` does not like underscores.

```
>$ cd hydro
>$ mkdir sedov-3d
```

**Note**: use one directory per test. If you want to run a 2D and a 3D sedov test, create separate `sedov-2d` and `sedov-3d` directories.

In that directory, you will need:

- A `config.txt` file, containing the Makefile flags and, if the test should be run through the restart mechanism, the `RESTART` line:
```
FLAGS: NDIM=3 PATCH= SOLVER=hydro
RESTART: true
```

- A namelist: `sedov-3d.nml` (the name needs to be the same as the test directory)

:::{warning}
For the restart mechanism to work, the end of the run and its outputs must be set in `&OUTPUT_PARAMS` with `tout`/`aout` (a single value or a list), `tend`/`aend`, and/or `delta_tout`/`delta_aout`. For instance
```
&OUTPUT_PARAMS
tout=0.620
/
```
or, for cosmological runs,
```
&OUTPUT_PARAMS
aout=1.355E-01
/
```
`foutput` is only supported with `foutput=1`, and a run stopped with `nstepmax` cannot be restarted. If the namelist schedules outputs before the end of the run, the restart happens at one of them. Otherwise, an extra output is written half-way and the outputs are renumbered afterwards, so the final output keeps the number it would have without the restart.
:::

- A file for plotting and checking the solution against a reference: `plot-sedov-3d.py`. It is advised to copy a file from the other directories to see how to write this. The script should save its figure as `sedov-3d.pdf`, which is included in the report. **Note that this file needs to contain at least one call to `visu_ramses.check_solution(data["data"], 'sedov-3d')`**, which prints `PASSED` or `FAILED`.

- A reference solution: `sedov-3d-ref.dat`. To create it, run your test (with `-d`, so that the outputs are kept) and once the final output (number 2 in this case) has been created, do the following in the test directory:
```
import visu_ramses
data = visu_ramses.load_snapshot(2)
visu_ramses.check_solution(data["data"], 'sedov-3d', overwrite=True)
```

- A `Readme.md` containing a short description of the test

Optional files:

- `condinit.f90`: you can have your own initial setup if it's not entirely definable in a namelist. **REMEMBER** to set the correct `PATCH` in the `config.txt` file! (e.g. `PATCH=../tests/hydro/sedov-3d`). Another option is to add a custom condinit routine option in the existing condinit.f90 file.

- `before-test.sh`: if this file is present in the test directory, it will be run before the test begins (useful for e.g. creating symbolic links to libraries, generating initial conditions, or running extra checks such as the parameter study of `hydro/sod-tube`, see below)

- `after-test.sh`: if this file is present in the test directory, it will be run when the results of the test are deleted at the end of the suite (useful for e.g. cleaning up symbolic links to libraries). It is not run when the suite is run with `-d`.

### Tuning tolerances for solution verification

By default, relative differences between the sums of all the variables inside all leaf cells in the domain and the reference solution cannot exceed `3.0e-13`.
Sometimes, some variables are more volatile than others when running simulations on different numbers of CPUs, and this limit is too low, leading to false failed tests.
The `check_solution` method in the `visu/visu_ramses.py` module can be tuned to work for your test using the following options:

- `tolerance`: a dictionary listing the allowed relative difference between the sum over all leaf cells and reference value. The default is `{"all":3.0e-13}`. To make the check on `density` less restrictive, use for instance `tolerance={"density":1.0e-10}`.

- `threshold`: relative value below which a vector component is set to zero. Default is `2.0e-14`.

- `norm_min`: minimum value for the norm of a vector, to protect against null vectors. Default is `1.0e-30`.

- `min_variance`: if the data differs by less than this value from the average value, it is set to the average. Default is `1.0e-14`.

### Parameter studies

To test many combinations of namelist parameters without adding one test per combination, a test can call `tests/run_analytical_test.py` from its `before-test.sh`. The script runs RAMSES for every combination of the parameters listed in `ana_parameters.py` and compares each result with a reference stored in `<test>-parameter-study-ref.csv`, using the `check_solution.py` of the test. See `hydro/sod-tube` for an example:
```
python ../../run_analytical_test.py "../../../bin/test_exe_1d sod-tube.nml" -t sod-tube
```

### Creating a new group of tests

If your test does not fall under the categories already present in the `tests` directory, you can create a new directory and put your tests in there. You will then have to edit the `run_test_suite.sh` file to ensure your new tests will be picked up.

Say you want to create 3 new tests, `sedov-1d`, `sedov-2d`, and `sedov-3d` inside a new `sedov` directory, you have to find the line describing the list of directories to be scanned at the top of the `run_test_suite.sh` file:

```
# List of directories to scan
testlist="hydro,mhd,poisson,rt,sink,star,turb,tracer";
```
and add your new directory separated from the previous one by a comma, i.e.
```
# List of directories to scan
testlist="hydro,mhd,poisson,rt,sink,star,turb,tracer,sedov";
```

To have the new group run by the continuous integration, also add it to the `test-family` list in `.github/workflows/tests.yaml`.


## 3. Code coverage

Running the suite with `-s` compiles RAMSES with `GCOV=1` and measures which lines of the source code each test executes:
```
./run_test_suite.sh -s
```
At the end of the run, the results are collected in a directory named after the branch and commit that were tested, `tests/coverage_<branch>_<commit date>_<short hash>/`, which contains:

- `coverage_report.txt`: the coverage of each source file and of the whole code.
- `coverage_notbuilt.txt`: the lines that no build of the run compiled, because they sit behind a preprocessor directive (`#ifdef ...`) that none of the tests enabled. These lines count as uncovered: covering them needs another build of RAMSES rather than another namelist.
- `gcov_files/`: one annotated `.gcov` file per source file, summed over all tests. Each line is marked with the number of times it ran, `#####` (compiled but never run), `-----` (never compiled by any build) or `-` (not executable).
- `coverage_tests.json`: which tests reached which line.
- `coverage_metadata.txt`, `build_records/`: how the run was made, and the flags with which each test was built and run.
- `gcov_per_test/`: the raw `.gcov` files of each test, one subdirectory per run.
- `test_results_<time>.pdf` and `test_suite_<time>.log`: the report and log of each run.

Coverage can be collected over several runs of the same commit, for instance to run some tests in parallel and others in serial, or with different builds. If the coverage directory of the commit already exists, the new run adds to it instead of replacing it.

The reports can be regenerated from the kept per-test files without re-running the tests:
```
python3 multi_gcov_aggregator.py <dir>/gcov_per_test/*/ <outdir> --build-records <dir>/build_records
```
