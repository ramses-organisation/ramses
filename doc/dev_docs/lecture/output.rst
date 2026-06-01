Outputting
------------

In this chapter we cover

- How outputs are written, structured, and extended.
- The link between I/O routines and the AMR (adaptive mesh refinement) hierarchy and MPI domain decomposition.

.. contents::


1. The different output files in a RAMSES snapshot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. admonition:: **Exercise**

   Look in a ramses output. What files are there?
   (If you have no outputs laying around,
   you can run one of the test suite cases with ``-d``, for example
   ``./run_test_suite.sh -t 1 -p 2 -d`` and look in the individual test
   folders, e.g. *tests/hydro/advect*)

   .. admonition:: **Solution**
      :class: dropdown

      Several types of files can be found in a RAMSES snapshot *output_xxxxx/*
      (non-exhaustive):

      * General simulation info:

         * info_xxxxxx.txt: units, boxlen, output time, cosmological constants,…
         * info_rt_xxxxxx.txt: info on photon groups and chemistry
         * header_xxxxxx.txt: number of particles per family

      * Data outputted by each MPI process (yyyyy):

         * amr_xxxxx.outyyyyyy: the grid linked list variables, grid center position, octree variables, loadbalancing and refinement map
         * hydro_xxxxx.outyyyyyy: all hydro variables
         * rt_xxxxx.outyyyyyy:
         * grav_xxxxx.outyyyyyy: gravitational potential and acceleration \* part_xxxxx.outyyyyyy: particle fields

      * File descriptors: these files list which variables are outputted in the corresponding data files and in which order they are:

         * hydro_file_descriptor.txt
         * rt_file_descriptor.txt
         * part_file_descriptor.txt

      * Data files outputted by cpu 0 only:

         * sink_xxxxxx.csv
         * stellar_xxxxxx.csv

      * Information about the execution:

         * namelist.txt: copy of the input namelist file
         * timer_xxxxxx.txt: execution time per module

      * Information about the compilation:

         * compilation.txt: date, commit hash, etc
         * makefile.txt: copy of the Makefile used
         * patches.txt: copy of the content of any patch files

2. The main output routine: ``dump_all``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

RAMSES outputs simulation data through a centralized routine,
``dump_all`` in *amr/output_amr.f90*, which orchestrates the writing of
all output files. It is executed by each MPI process, and so each MPI
process produces its own set of files in the snapshot.

Each module typically has its own ``output_xxxx.f90`` file / subroutine
to deal with outputting, which is then called by ``dump_all``.

3. Outputting the AMR structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The outputting of the AMR structure is handled by the routine
``backup_amr`` which can be found in ``amr/output_amr.f90``. It writes:

* the grid linked list variables,
* the grid center position, and the octree variables
* the loadbalancing map
* refinement map.

Normally, there is no need to alter something here.

4. Outputting the variables on the grid (hydro, rt, gravity)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Outputting the hydro variables, which are stored in ``uold`` (see the
chapter on Hydrodynamics `Section
2.2 <https://codimd.math.cnrs.fr/_iN5oHYnTJyNxOQL4biZuw#22-Hydro-variables-in-RAMSES>`__),
is done by the routine ``backup_hydro`` found in the file
*hydro/output_hydro.f90*. The structure of the routine is as follows:

.. admonition:: Click to show the code
   :class: dropdown

   .. code:: fortran

      ! Open the file and the file descriptor
      ...
      ! write some information to the header
      write(unit_out) ncpu
      ...
      ! Loop over refinement levels
      do ilevel = 1, nlevelmax
         ! Loop over physical boundaries and MPI domains
         do ibound = 1, nboundary+ncpu
            ! Get the head of grid linked list and number of elements in it
            if (ibound <= ncpu) then
               ncache = numbl(ibound, ilevel)
               istart = headl(ibound, ilevel)
            else
               ncache = numbb(ibound-ncpu, ilevel)
               istart = headb(ibound-ncpu, ilevel)
            end if
            write(unit_out) ilevel
            write(unit_out) ncache
            if (ncache > 0) then
               ! allocate work arrays to gather grid indices and variable values
               allocate(ind_grid(1:ncache), xdp(1:ncache))
               ! Gather the grid index ind_grid by following the linked list
               igrid = istart
               do i = 1, ncache
                  ind_grid(i) = igrid
                  igrid = next(igrid)
               end do
               ! Loop over cells in the grids
               do ind = 1, twotondim
                  iskip = ncoarse+(ind-1)*ngridmax
                  ! Write the data of all fields in order
                  ! starting with the Euler variables
                  do ivar = 1, neul-1
                     if (ivar == 1) then
                        ! Write density
                        do i = 1, ncache
                           xdp(i) = uold(ind_grid(i)+iskip, 1)
                        end do
                        field_name = 'density'
                     else
                        ! Write velocity field
                        do i = 1, ncache
                           xdp(i) = uold(ind_grid(i)+iskip, ivar)/max(uold(ind_grid(i)+iskip, 1), smallr)
                        end do
                        field_name = 'velocity_' // dim_keys(ivar - 1)
                     end if
                     call generic_dump(field_name, info_var_count, xdp, unit_out, dump_info_flag, unit_info)
                  end do
                  ...
            end do
            deallocate(ind_grid, xdp)
         end if
      end do
      end do
      close(unit_out)

An equivalent exists for outputting the content of ``rt_old``, done by
the subroutine ``rt_backup_hydro`` in *rt/rt_output_hydro.f90*.
Similarly, the outputting of the gravity variables ``phi``, ``f`` and
optionally ``rho`` is handled by the routine ``backup_poisson`` in
*poisson/output_poisson.f90*.

Generalized, the part of the inner loop that does the outputting of a
variable looks like this:

.. code:: fortran

   ! Gather the values of the variable ivar into the work array xdp
   do i = 1, ncache
      xdp(i) = uold(ind_grid(i)+iskip, ivar)
   end do
   ! Specify the name of the variable
   field_name = 'my_variable_name'
   ! Pass data to dump utils for writing
   call generic_dump(field_name, info_var_count, xdp, unit_out,&
                   & dump_info_flag, unit_info)

The routine ``generic_dump`` (defined in the file *io/dump_utils.f90*)
is performing the ``write`` instruction for any type. When adding new
variables on the grid, the corresponding backup routine needs to be
updated to include a new block following this structure. The name you
choose for your new variable will appear in the file descriptor, so make
sure to choose something informative.

5. Outputting particle fields
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The routine for outputting particles, ``backup_part``, can be found in
*pm/output_part.f90*. The code is simpler than for grid variables, since
we do not follow any linked list, but simply loop over all allocated
particle memory and check whether a particle exists at each memory
location. The structure of the code is as follows:

.. admonition:: Click to show the code
   :class: dropdown

   .. code:: fortran

      ! write some information to the header
      write(unit_out) ncpu
      write(unit_out) ndim
      write(unit_out) npart
      ...
      ! allocate float work array
      allocate(xdp(1:npart))
      ! Write float properties in order, starting with the position
      do idim = 1, ndim
         ! gather data in work array
         ipart = 0
         do i = 1, npartmax
            if (levelp(i) > 0) then ! check if particle exists
               ipart = ipart+1 ! count the number of particles
               xdp(ipart) = xp(i, idim)
            end if
         end do
         ! write data
         call generic_dump("position_"//dim_keys(idim), ivar, xdp, unit_out, dump_info, unit_info)
      end do
      ...
      ! Write integer properties
      ...

Remark that the counter ``ivar`` is updated automatically in
``generic_dump``.

Generalized, new particle fields can be added to the output by adding a
new block of the form:

.. code:: fortran

     ! allocate work array of the correcto type if needed
     allocate(work_array(1:npart))

     ! Gather data into work array by looping over all allocated space
     ipart = 0
     do i = 1, npartmax
        ! check if particle exists
        if (levelp(i) > 0) then
           ipart = ipart+1
           ! retrieve value
           work_array(ipart) = my_particle_array(i)
        end if
     end do
     ! send the data to be outputted
     call generic_dump("my_field_name", ivar, work_array, unit_out, dump_info, unit_info)

     ! deallocate work array if needed
     deallocate(work_array)

Be careful to choose a work array of the type corresponding to the type
of your new particle array. Give the variable an informative name.

6. Other files
~~~~~~~~~~~~~~~

Some other output files exist that do not follow the traditional rules.
For example:

* clumpfinder: ``write_clump_field`` in *pm/output_clump.f90* (writes clump field)
* sinks: ``output_sink_csv`` in *pm/output_sink.f90*
* stellars: ``output_stellar_csv`` in *pm/output_stellar.f90*

7. Restarts
~~~~~~~~~~~~

When restarting a simulation, ramses will read the information that was
previously outputted. It is thus important that the parts of the code
that read the input for restarting, match those of the outputting.

The restart parts are not encapsulated in their own routines, but can be
found in the various init routines when searching for
``if(nrestart>0)``.

Exercises
~~~~~~~~~~~~

.. admonition:: **Exercise**

   Add a new variable to uold that will keep track of the
   temperature.

.. admonition:: **Exercise**

   Add a file descriptor for the gravity output.
