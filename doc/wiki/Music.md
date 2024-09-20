

# Generate Zoom Initial Conditions with MUSIC

After having obtained the final output of the unigrid run, it can be analyzed for possible halo candidates for example with the `HOP` halo finder. Once the halo of interest is found and its position at z=0 determined, the high resolution region is in principle identified. This holds only for z=0, or likewise for the final output. The particles belonging to this region have of course been spread over a larger volume at earlier stages before they collapsed into the halo. Hence those particles’ positions in the first output (at z=zmax) need to be determined as well. Given this information an ellipsoid can be found containing all relevant particles, but with its volume being as small as possible. This connection from the final output to the initial conditions is implemented in the `get_mucic_refmask.f90` code. For input center of mass coordinates and radius of the considered region, this routine produces an output file containing the relevant information which serves as source for MUSIC, as the latter is run again to create initial conditions for the zoomed run.
 
It follows a description of the individual steps of this procedure:

(1) Create initial conditions with MUSIC for a unigrid simulation of given box size, by choosing levelmax=levelmin in the MUSIC configuration file.

(2) Run the unigrid RAMSES simulation with these initial conditions as described above.

(3) Apply the `HOP` halofinder to the final output, to obtain the list of halos contained in the simulated box.
```
$ ~/ramses/utils/f90/hop_ramses/hop -in output_00100/part_00100.out -p 1. -o hop_outputname
$ ~/ramses/utils/f90/hop_ramses/regroup -root hop_outputname -douter 80. -dsaddle 200. -dpeak 240. -f77 -o regroup_outputname
$ ~/ramses/utils/f90/hop_ramses/poshalo -inp output_00100 -pre regroup_outputname 
```
(4) Select a halo from the list
```
  #   npart       mass  cont.frac         xc         yc         zc         uc         vc         wc
  .     .           .       .              .          .          .          .          .          .
  .     .           .       .              .          .          .          .          .          .
  .     .           .       .              .          .          .          .          .          .
 88   67456  5.026E-04  0.000E+00  7.063E-01  4.010E-01  6.728E-01 -1.625E-02 -2.879E-04  1.586E-02
  .     .           .       .              .          .          .          .          .          .
  .     .           .       .              .          .          .          .          .          .
  .     .           .       .              .          .          .          .          .          .

```
and extract the bounding ellipsoid for the initial condition by running `get_mucic_refmask.f90`.

```
$ ~/MUSIC/get_music_refmask -ini output_00001 -inf output_00100 -out chosen_halo.part -xc 0.7063 -yc 0.4010 -zc 0.6728 -rad 0.008  
```

(5) Create new initial conditions for the zoom run by choosing the desired resolution of the zoomed region with levelmax>levelmin in the MUSIC configuration file, and by specifying the input file containing the information about the ellipsoid region, through adding the following line to the configuration file:

```
[setup]
region 			= ellipsoid
region_point_file 	= chosen_halo.part

```

(6) Once initial conditions for the zoom run are in place, re-run RAMSES with the specifications for the higher resolution. For this step the folders containing the initial conditions need to be listed in the namelist.   

```
&INIT_PARAMS
filetype='grafic'
initfile(1)='ics_88_zoom_2_r_4/level_007'
initfile(2)='ics_88_zoom_2_r_4/level_008'
initfile(3)='ics_88_zoom_2_r_4/level_009'
/
```

Further modifications to the ramses namelist are

```
&AMR_PARAMS
nexpand=5,1,1
/

&REFINE_PARAMS
mass_cut_refine=1.49012e-08
/
```
All three modifications are already set in the ramses namelist created by MUSIC. `nexpand` hereby depends on the padding value which is chosen in the MUSIC configuration file. The value for `mass_cut_refine` depends on the choice of levelmax in the initial conditions (since the resolution in the refined volume is higher, the effective mass of the particles therein is smaller).

After the completion of the zoomed ramses run, the chosen halo can be analyzed. However one should check first if it really consist only of high resolution particles. For this purpose `HOP` can be used. If it is applied to the final output of the zoom run, just as described above it will again list all halos sorted by their particle numbers. The difference is that this time the selected halo will have a higher particle number since it consist of the refined particles, which have smaller mass, therefore the halo will move up in the list.    
To check if there are no low resolution particles contained in it, after the zoom run, one can use the option 'contamination' of `HOP`.
Therefore the mass_cut value needs to be specified in executing `poshalo.f90`. To distinguish between the two values of particle masses the cut needs to be set slightly higher than the value of mass_cut_refine in the ramses namelist.  

```
$ ~/ramses/utils/f90/hop_ramses/poshalo -inp output_00100 -pre regroup_zoom_outputname -cut 1.5e-08
```
```
  #   npart       mass  cont.frac         xc         yc         zc         uc         vc         wc
  1   67757  5.026E-04  5.157E-02  7.063E-01  4.010E-01  6.728E-01 -1.625E-02 -2.879E-04  1.586E-02
  .     .           .       .              .          .          .          .          .          .
  .     .           .       .              .          .          .          .          .          .
  .     .           .       .              .          .          .          .          .          .

```
If the contamination of the selected halo is higher than desired, a possibility is to extract the bounding ellipsoid from the zoom run output itself, and to re-do the zoom run with initial conditions using this `get_music_refmask.f90` output. In principle this method can be applied more than one time. 

However, one needs to know the following: when MUSIC creates the initial conditions for a zoom run, it will shift the coordinate system so that the ellipsoid’s center of mass will be the center of the box (0.5, 0.5, 0.5) in the new initial conditions. Hence extra care needs to be taken, if a zoom run is undertaken on halo coordinates extracted from a previous zoom run. It is necessary to specify the shift of the coordinate system in the new MUSIC configuration file, through the levelmin value and displacement vector of the previous zoom run. 

```
[setup]
region 			= ellipsoid
region_point_file 	= chosen_halo.part
region_point_shift 	= -29,13,-18
region_point_levelmin 	= 7
```

This information is contained within the MUSIC logfile of the first zoom run.    
```
12:35:50 | info    |    Domain shifted by      (  -29,   13,  -18)
```