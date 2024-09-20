

# Restart from previous output #
A simulation which has been terminated during run time can be restarted from the last (or any) snapshot output, by setting 
```
nrestart=64
```
in the namelist file to the output number. If you don't want to change the namelist file, simply append the restart output number to the command execution, e.g.
```
./ramses3d parameters.nml 64
```

# Saving progress before job limit termination #
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

The signalling does not work on all machines. Sometimes the signal 10 is accompanied by a kill signal and the job is dead before it can perform an output. In this case, there are a couple of useful parameters in the output_params namelist: walltime_hrs can be used to specify the walltime given to a job in hours, and minutes_dump can then be used to tell RAMSES to dump an output a few minutes before the walltime runs out.

# Dump immediate output #
The above mechanism can be used to force an output be written to the disk immediately at any time during the simulation by sending signal 10 to the process:
```
scancel --signal=10 <jobid>
```
or, if you run without SLURM:
```
kill --signal=10 <processid>
```