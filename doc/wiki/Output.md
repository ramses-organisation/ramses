

# Output Parameters

This namelist block, called `&OUTPUT_PARAMS`, is used to set up the frequency and properties of data output to disk.

| Variable name, syntax, default value | Fortran type | Description |
|:---------------------------- |:------------- |:------------------------- |
| `noutput=1`  | `integer` | Number of specified output times.  At least one output time should be given, corresponding to the end of the simulation. |
| `foutput=1000000` | `integer` | Frequency of additional outputs in units of coarse time steps. `foutput=1` means one output at each time step. Specified outputs (see above) will not be superseded by this parameter. |
| `tout=0.0,0.0,0.0,` | `real array` | Value of specified output times. |
| `aout=1.1,1.1,1.1,` | `real array` | Value of specified output expansion factor (for cosmology runs only). `aout=1.0` means "present epoch" or "zero redshift". |
| `delta_tout=0` | `real` | Frequency of outputs in user time units. |
| `delta_aout=0` | `real` | Frequency of outputs in expansion factor (for cosmology runs only). |
| `tend=10` | `real` | End time of the simulation, in user time units. |
| `aend=0` | `real` | End time of the simulation, in expansion factor (for cosmology runs only). |
| `walltime_hrs=-1` | `real` | Wallclock time given in ramses job submission, used for dumping an output at the end. Default value of -1 means this is not used.|
| `minutes_dump=1` | `real` | Dump an output this many minutes before walltime_hrs |



## Restart from previous output
A simulation which has been terminated during run time can be restarted from the last (or any) snapshot output, by setting
```
nrestart=64
```
in the namelist file to the output number. If you don't want to change the namelist file, simply append the restart output number to the command execution, e.g.
```
./ramses3d parameters.nml 64
```

## Saving progress before job limit termination
SLURM jobs on computer clusters often have a time limit, after which the process will be terminated. If you don't want to lose the computation progress since your last regular output, you can instruct the SLURM scheduler to send a "warning" signal to the process <n> seconds before killing it with **--signal=10@<n>**. An example sbatch script with <n> set to 120 seconds would look like this:

```
#!/bin/bash
#SBATCH -J simulation
#SBATCH -p normal
#SBATCH -n 128
#SBATCH --time=24:00:00
#SBATCH --output=logfile-%j.txt
#SBATCH --error=error-%j.txt
#SBATCH --signal=10@120

aprun -B ./ramses3d parameters.nml 64
```

RAMSES will catch this signal and dump the current simulation state to a new output, which can be used to restart the simulation from.

:::{warning}
The signalling does not work on all machines. Sometimes the signal 10 is accompanied by a kill signal and the job is dead before it can perform an output. In this case, there are a couple of useful parameters in the output_params namelist: `walltime_hrs` can be used to specify the walltime given to a job in hours, and `minutes_dump` can then be used to tell RAMSES to dump an output a few minutes before the walltime runs out.
:::

## Dump immediate output
The above mechanism can be used to force an output be written to the disk immediately at any time during the simulation by sending signal 10 to the process:
```
scancel --signal=10 <jobid>
```
or, if you run without SLURM:
```
kill --signal=10 <processid>
```
