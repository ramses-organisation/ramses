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
Prepare a RAMSES test for a restart test and clean up afterwards.

The idea of a restart test is to reproduce, in two chunks joined by a restart,
exactly the same final output(s) that a single uninterrupted run would produce.
The test suite invokes this script in three steps around two runs of the code:

    python run_with_restart.py -s 1 -t <test_name>   # prepare namelist for shortened run
    <run the code>                                   # run 1 -> stops mid-way
    python run_with_restart.py -s 2 -t <test_name>   # prepare namelist for restart
    <run the code>                                   # run 2 -> reached original end time
    python run_with_restart.py -s 3 -t <test_name>   # restore namelist + renumber outputs

Where:
  -s <step> is the step to run (1, 2 or 3)
  -t <test_name> is the test name

Supported &OUTPUT_PARAMS options
--------------------------------
The output schedule can be driven by any combination of:
  tout / aout            explicit output times / expansion factors (scalar or list)
  tend / aend            end time / end expansion factor (appended as final tout/aout)
  delta_tout/delta_aout  periodic outputs in time / expansion factor
  noutput                number of predefined outputs (recomputed by RAMSES)
  write_conservative     conservative outputs -> read_conservative is set on restart

NOT supported: foutput (unless 1)

These mirror read_output_params() in amr/read_params.f90 and the output logic in
amr/amr_step.f90 / amr/update_time.f90.

Strategy
--------
Whenever possible, we choose the restart point from the options in the default
namelist. When there is no such candidate, an additional output is created and
output numbers are renumbered to match the original in the cleanup phase.
"""


def _as_list(value):
    """Return a namelist entry as a flat list of finite values."""
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        return [x for x in value if x is not None]
    return [value]


def analyze_schedule(output_params):
    """
    Reconstruct, from the &OUTPUT_PARAMS block, the information needed to split
    the run. Mirrors the termination logic of read_output_params/update_time.

    Returns a dict with:
      var    : 'time' or 'expansion' -- the controlling variable
      end    : value of that variable at which the run terminates
      end_type   : how the end is imposed ('tend'/'aend'/'tout'/'aout')
      candidates : sorted intermediate reference outputs strictly in (0, end)
      has_noutput    : whether noutput is written explicitly in the namelist
    """
    tout = _as_list(output_params.get("tout"))
    aout = _as_list(output_params.get("aout"))
    tend = output_params.get("tend", 0) or 0
    aend = output_params.get("aend", 0) or 0
    delta_tout = output_params.get("delta_tout")
    delta_aout = output_params.get("delta_aout")
    foutput = output_params.get("foutput")
    if (foutput not in [None,1]):
        print("ERROR: foutput!=1 not supported for restarts")
        return

    # Store information about which namelist param is controling the end of the run
    # Remark that aend/tend take precedence over tout/aout
    if aend > 0:
        var, end, end_type = "expansion", aend, "aend"
    elif tend > 0:
        var, end, end_type = "time", tend, "tend"
    elif tout:
        var, end, end_type = "time", max(tout), "tout"
    elif aout:
        var, end, end_type = "expansion", max(aout), "aout"
    else:
        raise ValueError(
            "Could not determine the simulation end from &OUTPUT_PARAMS "
            "(need one of tout, aout, tend or aend)."
        )
    # TODO: support for nstepmax

    # Determine possible restart points from expected outputs
    candidates = set()

    # list explicit intermediate output points
    if end_type in ['tout','aout']:
        output_points = _as_list(output_params.get(end_type))
        output_points = sorted(output_points)
        output_points = [i for i in output_points if i>0] #remove 0 if present
        for t in output_points[:-1]:
            candidates.add(t)

    # predicted output times when using delta
    if (var == "time"):
        delta = delta_tout
    else:
        delta = delta_aout
    if delta is not None:
        k = 1
        while k * delta < end:
            candidates.add(k * delta)
            k += 1

    return dict(
        end=end,
        end_type=end_type,
        candidates=sorted(candidates),
        has_noutput=("noutput" in output_params),
        foutput=foutput)


def decide_split(schedule):
    """Pick where run 1 should stop and if that creates an extra output."""
    candidates = schedule["candidates"]
    offset = 0
    if candidates:
        # that the scheduled output in the middle of the list of candidates
        # -> continuous numbering.
        split = candidates[(len(candidates) - 1) // 2]
    else:
        # No intermediate reference output: stop at the mid-point
        # if outputting at each coarse step, then the new end time is a pre-existing output
        split = schedule["end"] / 2.0
        if (schedule["foutput"]!=1): offset = 1
    return split, offset


def impose_end(output_params, schedule, split):
    """Modify the &OUTPUT_PARAMS block in place so the run terminates at `split`."""
    end_type = schedule["end_type"]

    if end_type in ("tend", "aend"):
        output_params[end_type] = split
    elif end_type in ("tout", "aout"):
        # insert new end time in list if not present (this can be the case when using delta)
        updated_params = _as_list(output_params[end_type])
        if split not in updated_params:
            updated_params.append(split)
        updated_params = sorted(updated_params)
        # cut the output list at the new end time
        i_end = updated_params.index(split)
        updated_params = updated_params[0:i_end+1]
        output_params[end_type] = updated_params
        # update noutput if set
        if schedule["has_noutput"]:
            output_params["noutput"] = len(updated_params)


def _list_output_numbers():
    """Sorted list of existing output_XXXXX directory numbers."""
    nums = []
    for path in glob.glob("output_*"):
        tail = path.split("_")[-1]
        if os.path.isdir(path) and tail.isdigit():
            nums.append(int(tail))
    return sorted(nums)


def step_1(test_name):
    """Shorten the run so it stops mid-way, backing up the original namelist."""
    nml_path = f"{test_name}.nml"

    if not os.path.exists(nml_path):
        print(f"Warning: {nml_path} does not exist")
        return

    shutil.copyfile(nml_path, nml_path + "_backup")

    nml = f90nml.read(nml_path)
    schedule = analyze_schedule(nml["output_params"])
    split, offset = decide_split(schedule)
    impose_end(nml["output_params"], schedule, split)

    f90nml.write(nml, nml_path, force=True)


def step_2(test_name):
    """Restore the full output schedule and restart from the last output written."""
    nml_path = f"{test_name}.nml"

    # Start from the original schedule so the second run reproduces it exactly.
    nml = f90nml.read(nml_path + "_backup")
    output_params = nml["output_params"]

    outputs = _list_output_numbers()
    if not outputs:
        print("ERROR: no output_XXXXX directory found to restart from")
        return
    last_output = outputs[-1]

    nml["run_params"]["nrestart"] = last_output

    # Conservative outputs must be read back as conservative on restart.
    if bool(output_params.get("write_conservative", False)):
        output_params["read_conservative"] = True

    f90nml.write(nml, nml_path, force=True)


def step_3(test_name):
    """Restore the original namelist and, if needed, renumber the outputs."""
    nml_path = f"{test_name}.nml"
    backup = nml_path + "_backup"

    # Determine (from the untouched original) whether run 1 created an extra
    # output that must be removed to match a normal run's numbering.
    offset = 0
    if os.path.exists(backup):
        sched = analyze_schedule(f90nml.read(backup)["output_params"])
        _, offset = decide_split(sched)
        os.rename(backup, nml_path)

    if offset == 0:
        # nothing to do, an intermediate snapshot was outputted by the original namelist
        return

    # offset == 1: the original run only produced the start and end snapshot,
    # and the restart produced an extra mid-way output. We remove output 2
    # and rename output_00003 -> output_00002
    if os.path.exists("output_00003"):
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
    parser = argparse.ArgumentParser(description="Prepare a test for restarts.")
    parser.add_argument("-s", "--step", help="Step to run (1, 2 or 3)", type=int, default=0)
    parser.add_argument("-t", "--test_name", help="Test name", type=str)
    args = parser.parse_args()

    steps = {1: step_1, 2: step_2, 3: step_3}

    # Run the specified step
    if args.step in steps:
        steps[args.step](args.test_name)
    else:
        print(f"Step {args.step} not recognized. Available steps are: {list(steps.keys())}")
