import os
import argparse
import shutil
import glob
import numpy as np
import csv


try:
    import f90nml
except ImportError:
    print("Please install f90nml to use this script.")
    exit(1)

from ana_parameters import params
from check_solution import check_solution

from itertools import product

def build_param_grid(d):
    for value in product(*d.values()):
        yield dict(zip(d.keys(), value))

class NamelistRecursive:
    def __init__(self, namelist):
        self.data = namelist

    def get_nml_value(self, nml_key):
        res = self.data
        for key in nml_key.split("/"):
            if key in res:
                res = res[key]
            elif key == nml_key.split("/")[-1]:
                res = None
            else:
                raise KeyError(key)
        return res

    def set_nml_value(self, nml_key, value):
        res = self.data
        nml_keys = nml_key.split("/")
        for key in nml_keys[:-1]:
            if key in res:
                res = res[key]
            else:
                raise KeyError(key)
        res[nml_keys[-1]] = value

    def __getitem__(self, key):
        return self.get_nml_value(key)

    def __setitem__(self, key, value):
        return self.set_nml_value(key, value)

    def __repr__(self):
        return self.data.__repr__()

    def __str__(self):
        return self.data.__str__()

def load_namelist(path):
    return NamelistRecursive(f90nml.read(path))


"""
Usage:
python run_analytical_test.py  <command> -t <test_name>

Where:
<command> is the command to launch the test
-t <test_name> is the name of the test (without .nml extension)

The script uses the f90nml library to read and write Fortran namelist files.
"""


def backup_namelist(test_name):
    nml_path = f"{test_name}.nml"

    # backup the original namelist file
    if os.path.exists(nml_path):
        shutil.copyfile(nml_path, nml_path + "_backup")
    else:
        print(f"Warning: {nml_path} does not exist")
        return

def restore_namelist(test_name):

    nml_path = f"{test_name}.nml"
    os.rename(nml_path + "_backup", nml_path)


def set_params(test_name, params):
    """
    Step 2: Modify the namelist file to
    add a second output time and set nrestart to 2.
    """

    nml_path = f"{test_name}.nml"
    nml = load_namelist(nml_path)


    for key in params:
        nml[key] = params[key]


    f90nml.write(nml=nml.data, nml_path=nml_path, force=True)

if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Run analytical test')
    parser.add_argument("command", help="The command to launch the test", type=str)
    parser.add_argument("-t", "--test_name", help="Test name", type=str)

    args = parser.parse_args()
    test_name = args.test_name

    param_grid = build_param_grid(params)
    results = {}
    errors = {}
    for param_set in param_grid:

        if "refine_params/interpol_var" in param_set and \
            "refine_params/interpol_type" in param_set and \
            param_set["refine_params/interpol_var"] == 4 and \
            param_set["refine_params/interpol_type"] != 2:
            # Skip this combination
            continue

        # Create a unique test name based on the parameters
        test_name_combi = "_".join(f"{k.split('/')[1]}_{v}" for k, v in param_set.items())
        print(f"Running test: {test_name_combi}")

        # Backup the original namelist file
        backup_namelist(test_name)

        try:
            # Set parameters in the namelist file
            set_params(test_name, param_set)

            # Run the analytical test (assuming a script or command exists)
            os.system(args.command)

            errors[test_name_combi] = check_solution()

        finally:
            # Restore the original namelist file
            restore_namelist(test_name)

    # Store results as plain text
    out_file = f"{test_name}-parameter-study.csv"
    with open(out_file, "w", newline="", encoding="utf‑8") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["test_name_combi", "error"])      # header row
        for combo, err in errors.items():
            writer.writerow([combo, "%.16e"%err])
