import f90nml
import os
import argparse
import shutil

""" 
This script modifies the namelist file of a test to prepare for a restart test.
It performs the following steps:
1. Modify the namelist file to divide the output time by 2 and create a backup.
2. Modify the namelist file to add a second output time and set nrestart to 2.
3. Clean up by recovering the original namelist file and copying output_00003 to output_00002.

Usage:
python run_with_restart.py -s <step> -t <test_name>

Where:
-s <step> is the step to run (1, 2, or 3)
-t <test_name> is the name of the test (without .nml extension)

The script uses the f90nml library to read and write Fortran namelist files.
"""


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
    
    if "tout" in nml["output_params"]:
        tout = nml["output_params"]["tout"]
        nml["output_params"]["tout"] = tout / 2
    elif "aout" in nml["output_params"]:
        aout = nml["output_params"]["aout"]
        nml["output_params"]["aout"] = aout / 2
    else:
        print("ERROR: tout or aout not found in output_params")
        return
    f90nml.write(nml=nml, nml_path=nml_path, force=True)
    
def step_2(test_name):
    """ 
    Step 2: Modify the namelist file to 
    add a second output time and set nrestart to 2.
    """ 

    nml_path = f"{test_name}.nml"
    nml = f90nml.read(nml_path)
    
    nml["output_params"]["noutput"] = 2
    nml["run_params"]["nrestart"] = 2

    if "tout" in nml["output_params"]:
        tout = nml["output_params"]["tout"]
        nml["output_params"]["tout"] = [tout, 2*tout]
    elif "aout" in nml["output_params"]:
        aout = nml["output_params"]["aout"]
        nml["output_params"]["aout"] = [aout, 2*aout]
    else:
        print("ERROR: tout or aout not found in output_params")
        return

    f90nml.write(nml=nml, nml_path=nml_path, force=True)
    
def step_3(test_name):
    """
    Step 3: Cleaning: recover the original namelist file
    and copy output_00003 to output_00002.
    """

    nml_path = f"{test_name}.nml"
    os.rename(nml_path + "_backup", nml_path)
    
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
        return
        
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
    
