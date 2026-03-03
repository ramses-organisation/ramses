try:
    import f90nml
except ImportError:
    print("Please install f90nml to use this script.")
    exit(1)

import os
import argparse
import shutil
import glob

"""
This script modifies the namelist file of a test to prepare for a restart test.
It performs the following steps:
1. Modify the namelist file to divide the output time by 2 and create a backup.
2. Modify the namelist file to add a second output time and set nrestart to 2 (or stop at the middle if more that one output time).
3. Clean up by recovering the original namelist file and, if needed, copy output_00003 to output_00002.

Usage:
python run_with_restart.py -s <step> -t <test_name>

Where:
-s <step> is the step to run (1, 2, or 3)
-t <test_name> is the name of the test (without .nml extension)

The script uses the f90nml library to read and write Fortran namelist files.
"""

def apply_output_factor(nml, factor, add_extra_output_time=False):
    """
    Apply a factor to the output time in the namelist.
    This function modifies the namelist in place.

    Parameters
    ---------
    nml: namelist
        Namelist object to modify in place
    factor: float
        Factor to apply to the existing output
    add_extra_output_time: bool
        Add the modified output time after the existing one instead of
        modifying it.
    """
    if "tout" in nml["output_params"]:
        tout = nml["output_params"]["tout"]
        if add_extra_output_time:
            nml["output_params"]["tout"] = [tout, tout * factor]
        else:
            nml["output_params"]["tout"] = tout * factor
    elif "aout" in nml["output_params"]:
        aout = nml["output_params"]["aout"]
        if add_extra_output_time:
            nml["output_params"]["aout"] = [aout, aout * factor]
        else:
            nml["output_params"]["aout"] = aout * factor
    else:
        print("ERROR: noutput found but tout or aout not found in output_params")

    return nml


def step_1(test_name):
    """
    Step 1: Modify the namelist file to
    divide the output time by 2 and create a backup.
    """

    nml_path = f"{test_name}.nml"

    # backup the original namelist file
    if os.path.exists(nml_path):
        shutil.copyfile(nml_path, nml_path + "_backup")
    else:
        print(f"Warning: {nml_path} does not exist")
        return

    nml = f90nml.read(nml_path)

    if "noutput" in nml["output_params"]:
        assert(nml["output_params"]["noutput"] == 1)
        nml = apply_output_factor(nml, 0.5)
    elif "tout" in nml["output_params"]:
        nout = len(nml["output_params"]["tout"])
        nml["output_params"]["tout"] = nml["output_params"]["tout"][0:int(nout/2)]
    else:
        assert("tend" in nml["output_params"])
        assert("foutput" in nml["output_params"])
        tend = nml["output_params"]["tend"]
        nml["output_params"]["tend"] = tend / 2

    #nml["output_params"]["write_conservative"] = True

    f90nml.write(nml=nml, nml_path=nml_path, force=True)

def step_2(test_name):
    """
    Step 2: Modify the namelist file to
    add a second output time and set nrestart to 2.
    """

    nml_path = f"{test_name}.nml"
    nml = f90nml.read(nml_path)
    nml_orig = f90nml.read(nml_path + "_backup")

    # Find the last output time
    all_outputs = sorted(glob.glob("output_*"))
    last_output = int(all_outputs[-1].split("_")[-1])
    print(f"Restarting from output {last_output}")
    nml["run_params"]["nrestart"] = last_output

    if "noutput" in nml["output_params"]:
        nml["run_params"]["nrestart"] = 2
        nml["output_params"]["noutput"] = 2
        nml = apply_output_factor(nml, 2, add_extra_output_time=True)
    elif "tout" in nml["output_params"]:
        # recover original output list from namelist
        nml["output_params"]["tout"] = nml_orig["output_params"]["tout"]
    else:
        tend = nml["output_params"]["tend"]
        nml["output_params"]["tend"] = tend * 2

    #nml["output_params"]["read_conservative"] = True
    ##if "write_conservative" in nml_orig["output_params"]:
    #nml["output_params"]["write_conservative"] = False #nml_orig["output_params"]["write_conservative"]

    f90nml.write(nml=nml, nml_path=nml_path, force=True)

def step_3(test_name):
    """
    Step 3: Cleaning: recover the original namelist file
    and, if needed, copy output_00003 to output_00002.
    """

    nml_path = f"{test_name}.nml"
    os.rename(nml_path + "_backup", nml_path)

    nml = f90nml.read(nml_path)

    if "noutput" in nml["output_params"] and os.path.exists("output_00003"):
        try:
            # Remove output 2
            shutil.rmtree(f"output_00002")
            # Move output 3 into output 2
            os.rename(f"output_00003", f"output_00002")

            for file in os.listdir("output_00002"):
                # remove extension
                name, ext = os.path.splitext(file)
                # rename file
                if name.endswith("00003"):
                    new_name = name[:-5] + "00002" + ext
                    os.rename(os.path.join("output_00002", file), os.path.join("output_00002", new_name))
        except FileNotFoundError:
            print(f"Warning: output_00003 or output_00002 does not exist")
            return -1

if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Prepare the test for restarts.')
    parser.add_argument("-s", "--step", help="Step to run", type=int, default=0)
    parser.add_argument("-t", "--test_name", help="Test name", type=str)

    args = parser.parse_args()

    steps = {
        1: step_1,
        2: step_2,
        3: step_3
    }

    # Run the specified step
    if args.step in steps:
        steps[args.step](args.test_name)
    else:
        print(f"Step {args.step} not recognized. Available steps are: {list(steps.keys())}")
