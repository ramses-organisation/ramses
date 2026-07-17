try:
    import f90nml
except ImportError:
    print("Please install f90nml to use this script.")
    exit(1)

import os
import re
import argparse
import shutil
import glob

"""
Prepare a RAMSES test for a restart test and clean up afterwards.

The idea of a restart test is to reproduce, in two chunks joined by a restart,
exactly the same final output(s) that a single uninterrupted run would produce.
The test suite invokes this script in three steps around two runs of the code:

    python run_with_restart.py -s 1 -t <test_name>   # shorten the run
    <run the code>                                   # run 1 -> stops mid-way
    python run_with_restart.py -s 2 -t <test_name>   # restore + set nrestart
    <run the code>                                   # run 2 -> completes
    python run_with_restart.py -s 3 -t <test_name>   # restore + renumber

Where:
  -s <step> is the step to run (1, 2 or 3)
  -t <test_name> is the test name

Supported &OUTPUT_PARAMS options
--------------------------------
The output schedule can be driven by any combination of:
  tout / aout            explicit output times / expansion factors (scalar or list)
  tend / aend            end time / end expansion factor (appended as final tout/aout)
  delta_tout/delta_aout  periodic outputs in time / expansion factor
  foutput                outputs every foutput coarse steps
  noutput                number of predefined outputs (recomputed by RAMSES)
  write_conservative     conservative outputs -> read_conservative is set on restart

These mirror read_output_params() in amr/read_params.f90 and the output logic in
amr/amr_step.f90 / amr/update_time.f90.

Strategy
--------
Whenever possible, we choose the restart point from the options in the default
namelist. When there is no such candidate, an additional output is created and
output numbers are renumbered to match the original in the cleanup phase.
"""

# foutput default in amr/amr_parameters.f90 (effectively "off").
FOUTPUT_OFF = 1000000


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
      end_type   : how the end is imposed ('tend'/'aend'/'tout_list'/'aout_list')
      candidates : sorted intermediate reference outputs strictly in (0, end)
      foutput_active : True if an output is written every coarse step-ish
      has_noutput    : whether noutput is written explicitly in the namelist
      write_conservative : whether outputs store conservative variables
    """
    tout = _as_list(output_params.get("tout"))
    aout = _as_list(output_params.get("aout"))
    tend = output_params.get("tend", 0) or 0
    aend = output_params.get("aend", 0) or 0
    delta_tout = output_params.get("delta_tout")
    delta_aout = output_params.get("delta_aout")
    foutput = output_params.get("foutput")
    foutput_active = foutput is not None and 0 < foutput < FOUTPUT_OFF

    # Controlling variable and end value. aend/tend are appended by RAMSES as the
    # final aout/tout and thus set the termination, taking precedence over the
    # explicit lists.
    if aend > 0:
        var, end, end_type = "expansion", aend, "aend"
    elif tend > 0:
        var, end, end_type = "time", tend, "tend"
    elif tout:
        var, end, end_type = "time", max(tout), "tout_list"
    elif aout:
        var, end, end_type = "expansion", max(aout), "aout_list"
    else:
        raise ValueError(
            "Could not determine the simulation end from &OUTPUT_PARAMS "
            "(need one of tout, aout, tend or aend)."
        )
    # TODO: support for nstepmax

    # Intermediate reference outputs strictly inside the run.
    candidates = set()
    if var == "time":
        explicit, delta = tout, delta_tout
    else:
        explicit, delta = aout, delta_aout
    for x in explicit:
        if 0 < x < end:
            candidates.add(x)
    if delta is not None and 0 < delta < FOUTPUT_OFF:
        k = 1
        while k * delta < end and k < FOUTPUT_OFF:
            candidates.add(k * delta)
            k += 1

    return dict(
        var=var,
        end=end,
        end_type=end_type,
        candidates=sorted(candidates),
        foutput_active=foutput_active,
        has_noutput=("noutput" in output_params),
        write_conservative=bool(output_params.get("write_conservative", False)),
    )


def decide_split(schedule):
    """Pick where run 1 should stop and how many extra outputs that creates."""
    candidates = schedule["candidates"]
    if candidates:
        # Stop on the middle scheduled output -> continuous numbering.
        split = candidates[(len(candidates) - 1) // 2]
        offset = 0
    else:
        # No intermediate reference output: stop at the mid-point.
        split = schedule["end"] / 2.0
        # With foutput active every step is an output, so the mid-point dump is
        # a regular output and numbering stays continuous.
        offset = 0 if schedule["foutput_active"] else 1
    return split, offset


def impose_end(output_params, schedule, split):
    """Modify the &OUTPUT_PARAMS block in place so the run terminates at `split`."""
    end_type = schedule["end_type"]
    single = len(schedule["candidates"]) == 0

    if end_type in ("tend", "aend"):
        output_params[end_type] = split
    elif end_type in ("tout_list", "aout_list"):
        key = "tout" if end_type == "tout_list" else "aout"
        if single:
            output_params[key] = split
            if schedule["has_noutput"]:
                output_params["noutput"] = 1
        else:
            new = sorted(x for x in _as_list(output_params.get(key)) if x <= split)
            if split not in new:
                new.append(split)
                new.sort()
            output_params[key] = new
            if schedule["has_noutput"]:
                output_params["noutput"] = len(new)


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
    """Restore the full schedule and restart from the last output written."""
    nml_path = f"{test_name}.nml"

    # Start from the original schedule so the second run reproduces it exactly.
    nml = f90nml.read(nml_path + "_backup")
    output_params = nml["output_params"]

    outputs = _list_output_numbers()
    if not outputs:
        print("ERROR: no output_XXXXX directory found to restart from")
        return
    last_output = outputs[-1]
    #print(f"[restart] step 2: restarting from output_{last_output:05d}")

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
        # nothing to do
        return

    # offset == 1: run 1 produced output_00001 (start) and output_00002 (extra
    # mid-point snapshot). Remove the extra and shift everything above down by 1.
    extra = 2
    outputs = _list_output_numbers()
    if not outputs:
        return
    last = outputs[-1]
    if not os.path.isdir(f"output_{extra:05d}") or last <= extra:
        return

    shutil.rmtree(f"output_{extra:05d}")
    for k in range(extra + 1, last + 1):
        src = f"output_{k:05d}"
        dst = f"output_{k - 1:05d}"
        os.rename(src, dst)
        for fname in os.listdir(dst):
            # Replace the (first) output-number token; leave the .outNNNNN cpu
            # index untouched.
            new = re.sub(f"{k:05d}", f"{k - 1:05d}", fname, count=1)
            if new != fname:
                os.rename(os.path.join(dst, fname), os.path.join(dst, new))
    #print(f"[restart] step 3: removed extra snapshot, renumbered {extra + 1}..{last} down by 1")


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
