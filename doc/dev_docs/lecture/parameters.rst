Interlude: Parameters, initialisation, outputs
==============================================

-  Go back to the code and review it from scratch starting with commons,
   namelist blocs, initialisations, etc… , and outputs.
-  **EXERCISE:** add a parameter to the namelist, read it, convert it to
   some units, and print it out. (Try adding parameters to different
   blocs).
-  **EXERCISE:** set the units of a passive scalar in the namelist and
   use that to change the units of some outputs. (e.g. output
   metallicity in solar units, with :math:`Z_\odot` provided in the
   namelist).

Input/Output and User Interaction
=================================

In this chapter we cover - How users control simulations via namelist
parameters. - How initial conditions are prepared and read. - How
outputs are written, structured, and extended. - The link between I/O
routines and the AMR (adaptive mesh refinement) hierarchy and MPI domain
decomposition.

[TOC]

1. Namelist parameters
----------------------

.. container:: success

   **Topics**: processing and adding user parameters

1.1 What is a namelist
~~~~~~~~~~~~~~~~~~~~~~

In RAMSES, namelist parameters form the primary user interface for
configuring a simulation. They allow users to specify runtime settings
without recompiling the code — such as solver options, output frequency,
cosmological constants, or model flags. At runtime, the user provides an
input file (typically *simulation_name.nml*) containing these
parameters. An example namelist file and overview of all runtime
parameters can be found in the `User
Documentation <https://ramses-organisation.readthedocs.io/en/latest/wiki/Runtime_Parameters.html>`__.
For example, the namelist block setting the AMR-related parameters can
look like this:

::

   &AMR_PARAMS
   levelmin=3
   levelmax=10
   ngridmax=2000
   nexpand=1
   boxlen=1.0
   /

Namelists are Fortran constructs that group configuration variables.
They are defined using the keyword ``namelist``, followed by a chosen
name and a list of the associated variables. For example the
``amr_params`` namelist in RAMSES is defined as follows:

.. code:: fortran

      namelist/amr_params/levelmin,levelmax,ngridmax,ngridtot &
      & ,npartmax,nparttot,nexpand,boxlen,nlevel_collapse

In RAMSES, the variables associated to namelist parameters are usually
declared in the one of the parameter or common modules (found in the
*module_parameter.f90* or *module_commons.f90* files). They are globally
accessible by importing the module, for example:

.. code:: fortran

   use amr_parameters

1.2 Reading and processing namelist parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Namelists are read once during program initialization. The keywords in
the namelist file provided by the user are mapped to the variables names
in the code. Their value is altered from the default value specified in
the code to the value provided by the user. In RAMSES, this is handled
by the ``read_params`` subroutine (found in *amr/read_params.f90*),
which is called at the very beginning of the simulation. The code to
process a namelist will go through the following steps:

.. code:: fortran=

      ! Definition of the namelist
      namelist/something_params/param1,param2,param3

      ! Go to the beginning of the namelist file
      rewind(namelist_unit)

      ! Read namelist
      read(namelist_unit,NML=something_params,IOSTAT=nml_err)

      ! Checks in whether the namelist was present in the file
      if(nml_err<0)then
         ! EOF reached before namelist was found
      elseif(nml_err>0)then
         ! Problem with formatting in the file
      endif

      ! Do some check on the read values
      if(param1<0) nml_ok=.false.
      ...

Older parts of the code may use more unorthodox syntax, but effectively
do the same thing. For example:

.. code:: fortran=

     ! Read namelist file
     rewind(1)
     read(1,NML=init_params,END=121)
     goto 122
   121 write(*,*)' You need to set up namelist &INIT_PARAMS in parameter file'
     call clean_stop
   122 rewind(1)

The identifier number for the namelist file is 1.

Remark that each MPI process will read the namelist file, and so each
process has access to all parameters.

1.3 Where are namelists defined in RAMSES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each module or subsystem (e.g. AMR, hydrodynamics, gravity) is
controlled by its own namelist definition and associated variables. The
naming convention for namelists is ``<module>_params``. The namelist
definitions can be found in one of the following places: -
*amr/read_params.f90*, which contains: - ``run_params`` - ``amr_params``
- ``output_params`` - ``movie_params`` - ``lightcone_params`` -
``tracer_params`` - ``poisson_params`` - *hydro/read_hydro_params.f90*,
which contains: - ``init_params`` - ``hydro_params`` - ``refine_params``
- ``boundary_params`` - ``feedback_params`` - ``cooling_params`` -
``sf_params`` - ``units_params`` - ``grackle_params`` -
``physics_params`` (legacy) - the corresponding code module: -
``clumpfind_params`` in *pm/clump_finder.f90* - ``mergertree_params`` in
*pm/merger_tree.f90* - ``stellar_params`` in
*pm/read_sink_feedback_params.f90* - ``sink_params`` in
*pm/sink_particle.f90* - ``unbinding_params`` in *pm/unbinding.f90* -
``rt_params`` and ``rt_groups`` in *rt/rt_init.f90* - ``turb_params`` in
*turb/read_turb_params.f90*

Remark that recently dedicated subroutines have been created in
*amr/read_params.f90* to handle the namelist defined in this file.

1.4 Adding a new namelist parameter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When introducing a new algorithm or making small modification to the
code, you may need to add a new parameter to an existing namelist. The
procedure is straightforward: - Step 1: Determine to which namelist the
new variable belongs and add it to the namelist declaration. - Step 2:
Identify the module to which to add the new variable and declare it. A
default value should be set, which will be used in case the parameter is
not specified in the namelist file. This also ensures backward
compatibility with existing namelist files. A comment should be added to
describe what the variable represents. - Step 3: Since the variable is
added to an existing namelist, it will be automatically read when the
corresponding namelist is parsed. Optionally, you can add checks to
verify whether the user provided a sensible value. - Step 4: Document
the new parameter. Add an entry to the markdown file listing the
parameters of the altered namelist, which can be found in the folder
*doc/wiki*.

.. container:: info

   **Exercise:** Add a new fictional parameter called
   ``tuto_heating_model``. Imagine that setting this parameter enables
   an additional heating source for which three models exist in the
   literature. By default, we want this source to be turned off. Test
   your code by printing the value of ``tuto_heating_model``. :::spoiler
   **Solution**

   **Step 1** Since the new parameter deals adds a heating source, it
   should be added to the ``cooling_params`` namelist, which deals with
   cooling/heating and chemistry. This namelist is defined in
   *hydro/read_hydro_params.f90*. We add it at the end, since the order
   in which the parameters are listed is irrelevant:

   .. code:: fortran

        ! Cooling / basic chemistry parameters
        namelist/cooling_params/cooling,metal,isothermal,haardt_madau,J21 &
             & ,barotropic_eos,barotropic_eos_form,polytrope_rho,polytrope_index,T_eos,mu_gas &
             & ,a_spec,self_shielding,z_ave,z_reion,ind_rsink,T2max,neq_chem &
             & ,cooling_ism,tuto_heating_model

   **Step 2** Searching for the other variables which are associated
   with this namelist, we find that they are defined in
   *amr/amr_parameters.f90*, rather than in
   *hydro/hydro_parameters.f90*. Legacy codes are often plagued by this
   sort of inconsistencies. We declare the new variable as an integer,
   since it will take the values 0,1,2 or 3 and we set its default value
   to 0, which turns the feature off:

   .. code:: fortran

        logical ::cooling_ism = .false.      ! Use cooling module from Audit & Hennebelle 2005 (non-RT)
                                              ! instead of ramses classical cooling
        integer ::tuto_heating_model=0       ! Model for the tutorial heating (1=Cleopatra+1996, 2=Thutmose+1998, 3=Seti+2004)
                                             ! 0 turns off tutorial heating

   **Step 3** We add some a check to verify that the user chose an
   existing model:

   .. code:: fortran

        !--------------------------------------------------
        ! Check tutorial heating model
        !--------------------------------------------------
        if(tuto_heating_model<0.or.tuto_heating_model>3)then
           if(myid==1)write(*,*)'Error: unknown tuto_heating_model. Choose 1,2,3 or 0.'
           nml_ok=.false.
        endif
        if(myid==1)write(*,*)'TUTORIAL: tuto_heating_model is set to',tuto_heating_model

   Several other parameters in the cooling namelist have checks. We can
   place our block of code after the checks on ``barotropic_eos`` for
   example.

   **Step 4** We find the description of the ``cooling_params`` namelist
   in the file *doc/wiki/Physics.md*. We add a line to the Cooling
   parameters table:

   ::

      | `tuto_heating_model`     | `integer` | `0`  | Model for the tutorial heating. 1=Cleopatra+1996, 2=Thutmose+1998, 3=Seti+2004. 0 turns off tutorial heating.

   Don’t forget to test your code. This can be done, for example, by
   adding the new parameter to one of the namelists in the test suite.

1.5 Defining a new namelist block
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When adding major new features, it is appropriate to define a new
namelist to group the parameters that control the new part of the code.
When developing bigger modules, we encourage to follow the example of
the turbulence module. In short: - Define module-specific variables in a
*module/module_parameters.f90* file - Write a routine
``read_module_params`` in which you define the namelist, read it and
process its parameters. Place it in a file *read_module_params.f90* file
included in the module’s directory. - Add the call to
``read_module_params`` to ``read_params``. - Update the documentation of
your module.

Remember that when creating new directories and files, the Makefile
needs to be updated.

.. container:: info

   **Exercise:** In *amr/read_params.f90*, add a subroutine that will
   read and process a new fictional namelist called ``tuto_params``,
   which contains two parameters: ``tuto_efficiency`` and
   ``tuto_timescale``. Both of these parameters have to be positive and
   cannot be zero. For simplicity, you can define the variables inside
   the new subroutine. :::spoiler **Solution** At the bottom of
   *amr/read_params.f90*, we add:

   .. code:: fortran=

      subroutine read_tuto_params(namelist_unit,nml_ok)
         use amr_parameters, only:dp
         use amr_commons, only:myid
         implicit none
         integer,intent(in)::namelist_unit
         logical,intent(inout)::nml_ok
         integer::nml_err
         real(dp) :: tuto_efficiency=1
         real(dp) :: tuto_timescale=1

         namelist/tuto_params/tuto_efficiency,tuto_timescale

         ! Go to the beginning of the file
         rewind(namelist_unit)

         ! Read namelist
         read(namelist_unit,NML=tuto_params,IOSTAT=nml_err)

         if(nml_err>0)then
            if(myid==1)write(*,*)'Error reading namelist &TUTO_PARAMS. Check formatting.'
            nml_ok=.false.
         endif

         if(tuto_efficiency<=0)then
            if(myid==1)write(*,*)'Error in the namelist: tuto_efficiency must be larger than 0'
            nml_ok=.false.
         endif

         if(tuto_timescale<=0)then
            if(myid==1)write(*,*)'Error in the namelist: tuto_timescale must be larger than 0'
            nml_ok=.false.
         endif

         if(myid==1)write(*,*)'TUTO: tuto_efficiency=',tuto_efficiency,', tuto_timescale=',tuto_timescale

      end subroutine read_tuto_params

   In the main routine of *amr/read_params.f90*, we add the call:

   .. code:: fortan

        call read_poisson_params(1,nml_ok)
        call read_tuto_params(1,nml_ok)

   Remark that because we defined the variables directly in the
   subroutine, we avoided having to update the Makefile to add new
   module files.

1.6 Bonus: namelist and Python
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you want to read/write namelist with python, check out the ``f90nml``
package.

2 Initial conditions
--------------------

There are several ways to provide the starting state of a RAMSES
simulation. For basic geometries, the user can use the parameters in the
``init_params`` namelist block. For more advanced setups that can be
described analytically, there are the ``condinit`` subroutines. Finally,
RAMSES also supports input from files with specified formats.

2.1 Analytical initial conditions on the grid with ``condinit``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

RAMSES users may already have experience with implementing initial
conditions using the ``condinit`` routine. This routine is called by
``init_flow_fine`` during the setup phase of the simulation. At the time
of writing, there are two versions available in the public version: one
for hydro and one for mhd. Remark that recently, the system has been
reworked to support various default setups (rather than replying on the
patch system).

As input, the ``condinit`` routine receives the cell center positions
(in code units) of the ``nn`` cells in the vector sweep, as well as the
cell size of the current level. From this information, the primitive
variables are then calculated. Several prescriptions are available by
default and can be selected by setting ``condinit_kind`` in the
namelist. One can easily add their own analytical initial conditions as
a new subroutine following the existing examples. Finally, the primitive
variables are converted to the conservative ones. The conservative
variables are then returned to ``init_flow_fine`` through the array
``u``, where they are written to ``uold``.

.. container:: info

   **Exercise:** Implement a new ``condinit_type`` that adds a
   sinusoidal perturbation on a uniform density background in 1D:
   :math:`\rho(x) = \rho_0 [1 + A \cos(\frac{2\pi x}{\lambda})]` The
   pressure is set to the same value as the density. :::spoiler
   **Solution**

   .. code:: fortran=

      ...
        case('jeans_instability_cos')
           call jeans_instability_cos_condinit(x, q, dx, nn)
      ...
      !================================================================
      subroutine jeans_instability_cos_condinit(x,q,dx,nn)
        use amr_parameters
        use hydro_parameters
        use constants, only:pi

        implicit none
        integer ::nn                            ! Number of cells
        real(dp)::dx                            ! Cell size
        real(dp),dimension(1:nvector,1:nvar)::q ! Primitive variables
        real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
        !================================================================
        ! sinusoidal perturbation of amplitude A
        !================================================================
        integer::i
        real(dp),parameter::A=1d-4, lambda=0.5

        ! Call built-in initial condition generator to init the fields
        call region_condinit(x,q,dx,nn)

        do i=1,nn
          ! density
          q(i,1) = 1+A*cos(2*pi*x(i,1)/lambda)
          ! pressure
          q(i,3)=q(i,1)
        end do

      end subroutine jeans_instability_cos_condinit

2.2 Input file formats
~~~~~~~~~~~~~~~~~~~~~~

Another way to provide the initial conditions is through files, by
setting the namelist parameters ``initfile`` and ``filetype``. For the
variables on the grid, supported formats are ascii and grafic (see
*init_flow_fine.f90*), while for particles ascii, grafic, and gadget are
available (see *init_part.f90*, and `Section
1.3 <https://codimd.math.cnrs.fr/z55fgvBcTjiGxrWt7VBUVw?both#13-Initialising-particles>`__
in the chapter on particles).

If you want to add your own input file, good luck.

3 Outputting
------------

3.1 The different output files in a RAMSES snapshot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. container:: info

   **Exercise:** Run the test suite with ``./run_`` Look in a ramses
   output. What files are there? (If you have no outputs laying around,
   you can run one of the test suite cases with ``-d``, for example
   ``./run_test_suite.sh -t 1 -p 2 -d`` and look in the individual test
   folders, e.g. *tests/hydro/advect*) :::spoiler **Solution** Several
   types of files can be found in a RAMSES snapshot *output_xxxxx/*
   (non-exhaustive): \* General simulation info: \* info_xxxxxx.txt:
   units, boxlen, output time, cosmological constants,… \*
   info_rt_xxxxxx.txt: info on photon groups and chemistry \*
   header_xxxxxx.txt: number of particles per family \* Data outputted
   by each MPI process (yyyyy): \* amr_xxxxx.outyyyyyy: the grid linked
   list variables, grid center position, octree variables, loadbalancing
   and refinement map \* hydro_xxxxx.outyyyyyy: all hydro variables \*
   rt_xxxxx.outyyyyyy: \* grav_xxxxx.outyyyyyy: gravitational potential
   and acceleration \* part_xxxxx.outyyyyyy: particle fields \* File
   descriptors: these files list which variables are outputted in the
   corresponding data files and in which order they are: \*
   hydro_file_descriptor.txt \* rt_file_descriptor.txt \*
   part_file_descriptor.txt \* Data files outputted by cpu 0 only: \*
   sink_xxxxxx.csv \* stellar_xxxxxx.csv \* Information about the
   execution: \* namelist.txt: copy of the input namelist file \*
   timer_xxxxxx.txt: execution time per module \* Information about the
   compilation: \* *compilation.txt*: date, commit hash, etc \*
   *makefile.txt*: copy of the Makefile used \* *patches.txt*: copy of
   the content of any patch files

3.2 The main output routine: ``dump_all``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

RAMSES outputs simulation data through a centralized routine,
``dump_all`` in *amr/output_amr.f90*, which orchestrates the writing of
all output files. It is executed by each MPI process, and so each MPI
process produces its own set of files in the snapshot.

Each module typically has its own ``output_xxxx.f90`` file / subroutine
to deal with outputting, which is then called by ``dump_all``.

3.3 Outputting the AMR structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The outputting of the AMR structure is handled by the routine
``backup_amr`` which can be found in ``amr/output_amr.f90``. It writes:
\* the grid linked list variables, \* the grid center position, and the
octree variables \* the loadbalancing map \* refinement map.

Normally, there is no need to alter something here.

3.4 Outputting the variables on the grid (hydro, rt, gravity)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Outputting the hydro variables, which are stored in ``uold`` (see the
chapter on Hydrodynamics `Section
2.2 <https://codimd.math.cnrs.fr/_iN5oHYnTJyNxOQL4biZuw#22-Hydro-variables-in-RAMSES>`__),
is done by the routine ``backup_hydro`` found in the file
*hydro/output_hydro.f90*. The structure of the routine is as follows:

:::spoiler Click to show the code

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

:::

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

3.5 Outputting particle fields
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The routine for outputting particles, ``backup_part``, can be found in
*pm/output_part.f90*. The code is simpler than for grid variables, since
we do not follow any linked list, but simply loop over all allocated
particle memory and check whether a particle exists at each memory
location. The structure of the code is as follows:

   :::spoiler Click to show the code

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

3.6 Other files
~~~~~~~~~~~~~~~

Some other output files exist that do not follow the traditional rules.
For example: \* clumpfinder: ``write_clump_field`` in
*pm/output_clump.f90* (writes clump field) \* sinks: ``output_sink_csv``
in *pm/output_sink.f90* \* stellars: ``output_stellar_csv`` in
*pm/output_stellar.f90*

3.7 Restarts
~~~~~~~~~~~~

When restarting a simulation, ramses will read the information that was
previously outputted. It is thus important that the parts of the code
that read the input for restarting, match those of the outputting.

The restart parts are not encapsulated in their own routines, but can be
found in the various init routines when searching for
``if(nrestart>0)``.

Exercise
--------

.. container:: info

   **Exercise:** Add a new variable to uold that will keep track of the
   temperature.

.. container:: info

   **Exercise**: Add a file descriptor for the gravity output.
