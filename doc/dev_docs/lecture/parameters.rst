Namelist parameters
----------------------

In this chapter, we cover how users control simulations via namelist parameters.

.. contents::


1. What is a namelist
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

2. Reading and processing namelist parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Namelists are read once during program initialization. The keywords in
the namelist file provided by the user are mapped to the variables names
in the code. Their value is altered from the default value specified in
the code to the value provided by the user. In RAMSES, this is handled
by the ``read_params`` subroutine (found in *amr/read_params.f90*),
which is called at the very beginning of the simulation. The code to
process a namelist will go through the following steps:

.. code:: fortran

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

.. code:: fortran

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

3. Where are namelists defined in RAMSES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each module or subsystem (e.g. AMR, hydrodynamics, gravity) is
controlled by its own namelist definition and associated variables. The
naming convention for namelists is ``<module>_params``. The namelist
definitions can be found in one of the following places:

* *amr/read_params.f90*, which contains:

   - ``run_params``
   - ``amr_params``
   - ``output_params``
   - ``movie_params``
   - ``lightcone_params``
   - ``tracer_params``
   - ``poisson_params``

* *hydro/read_hydro_params.f90*, which contains:

   - ``init_params``
   - ``hydro_params``
   - ``refine_params``
   - ``boundary_params``
   - ``feedback_params``
   - ``cooling_params``
   - ``sf_params``
   - ``units_params``
   - ``grackle_params``
   - ``physics_params`` (legacy)

* the corresponding code module:

   - ``clumpfind_params`` in *pm/clump_finder.f90*
   - ``mergertree_params`` in *pm/merger_tree.f90*
   - ``stellar_params`` in *pm/read_sink_feedback_params.f90*
   - ``sink_params`` in *pm/sink_particle.f90*
   - ``unbinding_params`` in *pm/unbinding.f90*
   - ``rt_params`` and ``rt_groups`` in *rt/rt_init.f90*
   - ``turb_params`` in *turb/read_turb_params.f90*

Remark that recently dedicated subroutines have been created in
*amr/read_params.f90* to handle the namelist defined in this file.

4. Adding a new namelist parameter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When introducing a new algorithm or making small modification to the
code, you may need to add a new parameter to an existing namelist. The
procedure is straightforward:

- Step 1: Determine to which namelist the new variable belongs and add it to the namelist declaration.
- Step 2: Identify the module to which to add the new variable and declare it. A default value should be set, which will be used in case the parameter is not specified in the namelist file. This also ensures backward compatibility with existing namelist files. A comment should be added to describe what the variable represents.
- Step 3: Since the variable is added to an existing namelist, it will be automatically read when the corresponding namelist is parsed. Optionally, you can add checks to verify whether the user provided a sensible value.
- Step 4: Document the new parameter. Add an entry to the markdown file listing the parameters of the altered namelist, which can be found in the folder
*doc/wiki*.

.. admonition:: **Exercise**

   Add a new fictional parameter called
   ``tuto_heating_model``. Imagine that setting this parameter enables
   an additional heating source for which three models exist in the
   literature. By default, we want this source to be turned off. Test
   your code by printing the value of ``tuto_heating_model``.

   .. admonition:: **Solution**
      :class: dropdown

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

5. Defining a new namelist block
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When adding major new features, it is appropriate to define a new
namelist to group the parameters that control the new part of the code.
When developing bigger modules, we encourage to follow the example of
the turbulence module. In short:

- Define module-specific variables in a *module/module_parameters.f90* file
- Write a routine ``read_module_params`` in which you define the namelist, read it and process its parameters. Place it in a file *read_module_params.f90* file included in the module’s directory.
- Add the call to ``read_module_params`` to ``read_params``.
- Update the documentation of your module.

Remember that when creating new directories and files, the Makefile
needs to be updated.

.. admonition:: **Exercise**

   In *amr/read_params.f90*, add a subroutine that will
   read and process a new fictional namelist called ``tuto_params``,
   which contains two parameters: ``tuto_efficiency`` and
   ``tuto_timescale``. Both of these parameters have to be positive and
   cannot be zero. For simplicity, you can define the variables inside
   the new subroutine.


   .. admonition:: **Solution**
      :class: dropdown

      At the bottom of *amr/read_params.f90*, we add:

      .. code:: fortran

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

6. Bonus: namelist and Python
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you want to read/write namelist with python, check out the ``f90nml``
package.
