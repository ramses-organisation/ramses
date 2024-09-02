subroutine read_params
    use amr_commons
    use mond_parameters
    use pm_parameters
    use poisson_parameters
    use hydro_parameters
    implicit none
#ifndef WITHOUTMPI
    include 'mpif.h'
#endif
    ! --------------------------------------------------
    ! Local variables
    ! --------------------------------------------------
    integer :: i, narg, iargc, ierr, levelmax
    character(LEN=80) :: infile
    character(LEN=80) :: cmdarg
    integer(kind=8) :: ngridtot=0
    integer(kind=8) :: nparttot=0
    real(kind=8) :: delta_tout=0, tend=0
    real(kind=8) :: delta_aout=0, aend=0
    logical :: nml_ok

    ! --------------------------------------------------
    ! Namelist definitions
    ! --------------------------------------------------

    ! ~~~~~~~~~ begin ~~~~~~~~~
    ! Added the "mond" parameter
    namelist / run_params / clumpfind, cosmo, pic, sink, lightcone, poisson, hydro, rt, verbose, debug &
        & , nrestart, ncontrol, nstepmax, nsubcycle, nremap, ordering &
        & , bisec_tol, static, overload, cost_weighting, aton &
        & , mond
    ! ~~~~~~~~~~ end ~~~~~~~~~~
    namelist / output_params / noutput, foutput, aout, tout &
        & , tend, delta_tout, aend, delta_aout, gadget_output
    namelist / amr_params / levelmin, levelmax, ngridmax, ngridtot &
        & , npartmax, nparttot, nexpand, boxlen
    ! ~~~~~~~~~ begin ~~~~~~~~~
    ! Added the "a0" parameter
    namelist / poisson_params / epsilon, gravity_type, gravity_params &
        & , cg_levelmin, cic_levelmax, a0, a0_ms2
    ! ~~~~~~~~~~ end ~~~~~~~~~~
    namelist / lightcone_params / thetay_cone, thetaz_cone, zmax_cone
    namelist / movie_params / levelmax_frame, nw_frame, nh_frame, ivar_frame &
        & , xcentre_frame, ycentre_frame, zcentre_frame &
        & , deltax_frame, deltay_frame, deltaz_frame, movie &
        & , imovout, imov, tendmov, aendmov, proj_axis

    ! MPI initialization
#ifndef WITHOUTMPI
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, ncpu, ierr)
    myid = myid + 1 ! Careful with this...
#endif
#ifdef WITHOUTMPI
    ncpu = 1
    myid = 1
#endif
    ! --------------------------------------------------
    ! Advertise RAMSES
    ! --------------------------------------------------

    ! ~~~~~~~~~ begin ~~~~~~~~~
    ! Advertise also the QUMOND patch to ensure the user
    ! is aware of the fact that this executable has been
    ! patched
    !
    if (myid == 1) then
        write(*, *) ''
        write(*, *) '                    ~  The Phantom of  ~                        '
        write(*, *) ''
        write(*, *) '_/_/_/         .-.     _/    _/    _/_/_/   _/_/_/_/    _/_/_/  '
        write(*, *) '_/    _/     _/ ..\    _/_/_/_/   _/    _/  _/         _/    _/ '
        write(*, *) '_/    _/    ( \  v/__  _/ _/ _/   _/        _/         _/       '
        write(*, *) '_/_/_/       \     \   _/    _/     _/_/    _/_/_/       _/_/   '
        write(*, *) '_/    _/     /     |   _/    _/         _/  _/               _/ '
        write(*, *) '_/    _/  __/       \  _/    _/   _/    _/  _/         _/    _/ '
        write(*, *) '_/    _/ (   _._.-._/  _/    _/    _/_/_/   _/_/_/_/    _/_/_/  '
        write(*, *) '          `-`                                        '
        write(*, *) '                     RAMSES  Version 3               '
        write(*, *) '       written by Romain Teyssier (CEA/DSM/IRFU/SAP) '
        write(*, *) '                     (c) CEA 1999-2007               '
        write(*, *) '                                                     '
        write(*, *) '                  with  MONDifications by            '
        write(*, *) '                 F. Lueghausen  (Uni Bonn)           '
        write(*, *) '                                                     '
        !
        ! ~~~~~~~~~~ end ~~~~~~~~~~
        write(*, '(" Working with nproc = ",I4," for ndim = ",I1)') ncpu, ndim
        ! Check nvar is not too small


        ! ~~~~~~~~~ begin ~~~~~~~~~
        ! Ensure NDIM=3, otherwise throw an error msg and stop
#ifndef NDIM == 3
        if (mond) then
            write(*, '(" ERROR: The QUMOND Poisson solver requires NDIM=3")')
            call clean_stop
        end if
#endif
        ! ~~~~~~~~~~ end ~~~~~~~~~~
#ifdef SOLVERhydro
        write(*, '(" Using the hydro solver with nvar = ",I2)') nvar
        if (nvar < ndim + 2) then
            write(*, *) 'You should have: nvar>=ndim+2'
            write(*, '(" Please recompile with -DNVAR=",I2)') ndim + 2
            call clean_stop
        end if
#endif
#ifdef SOLVERmhd
        write(*, '(" Using the mhd solver with nvar = ",I2)') nvar
        if (nvar < 8) then
            write(*, *) 'You should have: nvar>=8'
            write(*, '(" Please recompile with -DNVAR=8")')
            call clean_stop
        end if
#endif
        write(*, *) ' '

        ! Read namelist filename from command line argument
        narg = command_argument_count()
        IF (narg < 1) THEN
            write(*, *) 'You should type: ramses3d input.nml [nrestart]'
            write(*, *) 'File input.nml should contain a parameter namelist'
            write(*, *) 'nrestart is optional'
            call clean_stop
        END IF
        CALL getarg(1, infile)
    end if
#ifndef WITHOUTMPI
    call MPI_BCAST(infile, 80, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
#endif

    ! -------------------------------------------------
    ! Read the namelist
    ! -------------------------------------------------
    INQUIRE(file = infile, exist = nml_ok)
    if (.not. nml_ok) then
        if (myid == 1) then
            write(*, *) 'File '// TRIM(infile) //' does not exist'
        end if
        call clean_stop
    end if

    open(1, file=infile)
    rewind(1)
    read(1, NML = run_params)
    rewind(1)
    read(1, NML = output_params)
    rewind(1)
    read(1, NML = amr_params)
    rewind(1)
    read(1, NML = lightcone_params, END = 83)
83  continue
    rewind(1)
    read(1, NML = movie_params, END = 82)
82  continue
    rewind(1)
    read(1, NML = poisson_params, END = 81)
81  continue

    ! -------------------------------------------------
    ! Read optional nrestart command-line argument
    ! -------------------------------------------------
    if (myid == 1 .and. narg == 2) then
        CALL getarg(2, cmdarg)
        read(cmdarg,*) nrestart
    end if

#ifndef WITHOUTMPI
    call MPI_BCAST(nrestart, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif

    ! -------------------------------------------------
    ! Compute time step for outputs
    ! -------------------------------------------------
    if (tend > 0) then
        if (delta_tout == 0) delta_tout = tend
        noutput = MIN(int(tend / delta_tout), MAXOUT)
        do i = 1, noutput
            tout(i) = dble(i) * delta_tout
        end do
    else if (aend > 0) then
        if (delta_aout == 0) delta_aout = aend
        noutput = MIN(int(aend / delta_aout), MAXOUT)
        do i = 1, noutput
            aout(i) = dble(i) * delta_aout
        end do
    end if
    noutput = MIN(noutput, MAXOUT)
    if (imovout > 0) then
        allocate(tmovout(1:imovout))
        allocate(amovout(1:imovout))
        tmovout = 1d100
        amovout = 1d100
        if (tendmov > 0) then
            do i = 1, imovout
                tmovout(i) = tendmov * dble(i) / dble(imovout)
            end do
        end if
        if (aendmov > 0) then
            do i = 1, imovout
                amovout(i) = aendmov * dble(i) / dble(imovout)
            end do
        end if
        if (tendmov == 0 .and. aendmov == 0) movie = .false.
    end if
    ! --------------------------------------------------
    ! Check for errors in the namelist so far
    ! --------------------------------------------------
    levelmin = MAX(levelmin, 1)
    nlevelmax = levelmax
    nml_ok = .true.
    if (levelmin < 1) then
        if (myid == 1) write(*, *) 'Error in the namelist:'
        if (myid == 1) write(*, *) 'levelmin should not be lower than 1 !!!'
        nml_ok = .false.
    end if
    if (nlevelmax < levelmin) then
        if (myid == 1) write(*, *) 'Error in the namelist:'
        if (myid == 1) write(*, *) 'levelmax should not be lower than levelmin'
        nml_ok = .false.
    end if
    if (ngridmax == 0) then
        if (ngridtot == 0) then
            if (myid == 1) write(*, *) 'Error in the namelist:'
            if (myid == 1) write(*, *) 'Allocate some space for refinements !!!'
            nml_ok = .false.
        else
            ngridmax = ngridtot / int(ncpu, kind = 8)
        end if
    end if
    if (npartmax == 0) then
        npartmax = nparttot / int(ncpu, kind = 8)
    end if
    if (myid > 1) verbose = .false.
    if (sink .and. (.not. pic)) then
        pic = .true.
    end if
    if (clumpfind .and. (.not. pic)) then
        pic = .true.
    end if
    ! if(pic.and.(.not.poisson))then
    ! poisson=.true.
    ! endif

    call read_hydro_params(nml_ok)
#ifdef RT
    call rt_read_hydro_params(nml_ok)
#endif
    if (sink) call read_sink_params
    if (clumpfind .or. sink) call read_clumpfind_params


    close(1)

    ! -----------------
    ! Max size checks
    ! -----------------
    if (nlevelmax > MAXLEVEL) then
        write(*, *) 'Error: nlevelmax>MAXLEVEL'
        call clean_stop
    end if
    if (nregion > MAXREGION) then
        write(*, *) 'Error: nregion>MAXREGION'
        call clean_stop
    end if

    ! -----------------------------------
    ! Rearrange level dependent arrays
    ! -----------------------------------
    do i = nlevelmax, levelmin, - 1
        nexpand   (i) = nexpand   (i - levelmin + 1)
        nsubcycle (i) = nsubcycle (i - levelmin + 1)
        r_refine  (i) = r_refine  (i - levelmin + 1)
        a_refine  (i) = a_refine  (i - levelmin + 1)
        b_refine  (i) = b_refine  (i - levelmin + 1)
        x_refine  (i) = x_refine  (i - levelmin + 1)
        y_refine  (i) = y_refine  (i - levelmin + 1)
        z_refine  (i) = z_refine  (i - levelmin + 1)
        m_refine  (i) = m_refine  (i - levelmin + 1)
        exp_refine(i) = exp_refine(i - levelmin + 1)
        initfile  (i) = initfile  (i - levelmin + 1)
    end do
    do i = 1, levelmin - 1
        nexpand   (i) = 1
        nsubcycle (i) = 1
        r_refine  (i) = - 1.0
        a_refine  (i) = 1.0
        b_refine  (i) = 1.0
        x_refine  (i) = 0.0
        y_refine  (i) = 0.0
        z_refine  (i) = 0.0
        m_refine  (i) = - 1.0
        exp_refine(i) = 2.0
        initfile  (i) = ' '
    end do

    if (.not. nml_ok) then
        if (myid == 1) write(*, *) 'Too many errors in the namelist'
        if (myid == 1) write(*, *) 'Aborting...'
        call clean_stop
    end if

#ifndef WITHOUTMPI
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif

end subroutine read_params
