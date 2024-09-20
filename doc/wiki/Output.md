

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