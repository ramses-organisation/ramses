module dice_commons
  use amr_commons
  use hydro_commons

  ! particle data
  character(len=512)::ic_file, ic_format
  ! misc
  real(dp)::IG_rho         = 1.0D-5
  real(dp)::IG_T2          = 1.0D7
  real(dp)::IG_metal       = 0.01
  real(dp)::ic_scale_pos   = 1.0
  real(dp)::ic_scale_vel   = 1.0
  real(dp)::ic_scale_mass  = 1.0
  real(dp)::ic_scale_u     = 1.0
  real(dp)::ic_scale_age   = 1.0
  real(dp)::ic_scale_metal = 1.0
  real(dp)::ic_t_restart   = 0.0D0
  integer::ic_mask_ivar    = 0
  real(dp)::ic_mask_min    = 1d40
  real(dp)::ic_mask_max    = -1d40
  integer::ic_mask_ptype   = -1
  integer::ic_ifout        = 1
  integer::ic_nfile        = 1
  integer,dimension(1:6)::ic_skip_type        = -1
  integer,dimension(1:6)::cosmo_add_gas_index = -1
  real(dp),dimension(1:3)::ic_mag_const = (/ 0.0, 0.0, 0.0 /)
  real(dp),dimension(1:3)::ic_center    = (/ 0.0, 0.0, 0.0 /)
  character(len=4)::ic_head_name  = 'HEAD'
  character(len=4)::ic_pos_name   = 'POS '
  character(len=4)::ic_vel_name   = 'VEL '
  character(len=4)::ic_id_name    = 'ID  '
  character(len=4)::ic_mass_name  = 'MASS'
  character(len=4)::ic_u_name     = 'U   '
  character(len=4)::ic_metal_name = 'Z   '
  character(len=4)::ic_age_name   = 'AGE '
  ! Gadget units in cgs
  real(dp)::gadget_scale_l = 3.085677581282D21
  real(dp)::gadget_scale_v = 1.0D5
  real(dp)::gadget_scale_m = 1.9891D43
  real(dp)::gadget_scale_t = 1.0D6*365*24*3600
  real(dp),allocatable,dimension(:)::up
  real(dp),allocatable,dimension(:)::maskp
  logical::dice_init       = .false.
  logical::amr_struct      = .false.
  ! magnetic
  integer,parameter::MAXGAL= 32
  real(dp),dimension(1:MAXGAL)::ic_mag_center_x = 0.0
  real(dp),dimension(1:MAXGAL)::ic_mag_center_y = 0.0
  real(dp),dimension(1:MAXGAL)::ic_mag_center_z = 0.0
  real(dp),dimension(1:MAXGAL)::ic_mag_axis_x   = 0.0
  real(dp),dimension(1:MAXGAL)::ic_mag_axis_y   = 0.0
  real(dp),dimension(1:MAXGAL)::ic_mag_axis_z   = 1.0
  real(dp),dimension(1:MAXGAL)::ic_mag_scale_R  = 1.0
  real(dp),dimension(1:MAXGAL)::ic_mag_scale_H  = 1.0
  real(dp),dimension(1:MAXGAL)::ic_mag_scale_B  = 0.0

end module dice_commons

subroutine read_params
  use amr_commons
  use hydro_parameters, only:nvar,nhydro
  use mpi_mod
  use buildinfo
  use iso_fortran_env, ONLY: output_unit !standard output

  implicit none
  !--------------------------------------------------
  ! Local variables
  !--------------------------------------------------
  integer::i,narg
  character(LEN=80)::infile, info_file
  character(LEN=80)::cmdarg
  character(LEN=5)::nchar
  logical::nml_ok, info_ok
  integer,parameter::tag=1134
#ifndef WITHOUTMPI
  integer::dummy_io,ierr,info2
#endif

#ifdef LIGHT_MPI_COMM
   ! RAMSES legacy communicator (from amr_commons.f90)
   type communicator_legacy
      integer                            ::ngrid_legacy
      integer                            ::npart_legacy
      integer     ,dimension(:)  ,pointer::igrid_legacy
      integer     ,dimension(:,:),pointer::f_legacy
      real(kind=8),dimension(:,:),pointer::u_legacy
      integer(i8b),dimension(:,:),pointer::fp_legacy
      real(kind=8),dimension(:,:),pointer::up_legacy
   end type communicator_legacy
   real(kind=8)::mem_used_legacy_buff, mem_used_new_buff,mem_used_legacy_buff_mg, mem_used_new_buff_mg
   type(communicator_legacy),allocatable,dimension(:,:)::emission_reception_legacy  ! 2D (ncpu,nlevelmax) data emission/reception/active_mg/emission_mg "heavy" buffer
#endif

  !-----------------------------------------------------
  ! Advertise RAMSES and print some general information
  !-----------------------------------------------------

  if(myid==1)then
  write(*,*)'_/_/_/       _/_/     _/    _/    _/_/_/   _/_/_/_/    _/_/_/  '
  write(*,*)'_/    _/    _/  _/    _/_/_/_/   _/    _/  _/         _/    _/ '
  write(*,*)'_/    _/   _/    _/   _/ _/ _/   _/        _/         _/       '
  write(*,*)'_/_/_/     _/_/_/_/   _/    _/     _/_/    _/_/_/       _/_/   '
  write(*,*)'_/    _/   _/    _/   _/    _/         _/  _/               _/ '
  write(*,*)'_/    _/   _/    _/   _/    _/   _/    _/  _/         _/    _/ '
  write(*,*)'_/    _/   _/    _/   _/    _/    _/_/_/   _/_/_/_/    _/_/_/  '
  write(*,*)'                        Version 3.0                            '
  write(*,*)'       written by Romain Teyssier (Princeton University)       '
  write(*,*)'           (c) CEA 1999-2007, UZH 2008-2021, PU 2022           '
  write(*,*)' '
  write(*,'(" Working with nproc = ",I5," for ndim = ",I1)')ncpu,ndim
  ! Check nvar is not too small
#ifdef SOLVERhydro
  write(*,'(" Using solver = hydro with nvar = ",I2)')nvar
  if(nvar<nhydro)then
     write(*,*)'You should have: nvar>=ndim+2'
     write(*,'(" Please recompile with -DNVAR=",I2)')nhydro
     call clean_stop
  endif
#endif
#ifdef SOLVERmhd
  write(*,'(" Using solver = mhd with nvar = ",I2)')nvar
  if(nvar<nhydro)then
     write(*,*)'You should have: nvar>=8'
     write(*,'(" Please recompile with -DNVAR=8")')
     call clean_stop
  endif
#endif

#ifdef TSC
  if (ndim/=3)then
     write(*,*)'TSC not supported for ndim neq 3'
     call clean_stop
  end if
#endif

  !Write I/O group size information
  if(IOGROUPSIZE>0.or.IOGROUPSIZECONE>0.or.IOGROUPSIZEREP>0)write(*,*)' '
  if(IOGROUPSIZE>0) write(*,*)'IOGROUPSIZE=',IOGROUPSIZE
  if(IOGROUPSIZECONE>0) write(*,*)'IOGROUPSIZECONE=',IOGROUPSIZECONE
  if(IOGROUPSIZEREP>0) write(*,*)'IOGROUPSIZEREP=',IOGROUPSIZEREP
  if(IOGROUPSIZE>0.or.IOGROUPSIZECONE>0.or.IOGROUPSIZEREP>0)write(*,*)' '

  ! Write information about git version
  write(*,*)' '
  call write_gitinfo(output_unit)
  write(*,*)' '

  ! Read namelist filename from command line argument
  narg = command_argument_count()
  IF(narg .LT. 1)THEN
     write(*,*)'You should type: ramses3d input.nml [nrestart]'
     write(*,*)'File input.nml should contain a parameter namelist'
     write(*,*)'nrestart is optional'
     call clean_stop
  END IF
  CALL getarg(1,infile)
  endif
#ifndef WITHOUTMPI
  call MPI_BCAST(infile,80,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
#endif

  !-------------------------------------------------
  ! Read the namelist
  !-------------------------------------------------

  ! Wait for the token
#ifndef WITHOUTMPI
     if(IOGROUPSIZE>0) then
        if (mod(myid-1,IOGROUPSIZE)/=0) then
           call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
                & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
        end if
     endif
#endif

  ! Check if namelist file exists
  namelist_file=TRIM(infile)
  INQUIRE(file=infile,exist=nml_ok)
  if(.not. nml_ok)then
     if(myid==1)write(*,*)'File '//TRIM(infile)//' does not exist'
     call clean_stop
  end if

  ! Open namelist
  open(1,file=infile)

  ! Read parameter blocks
  call read_run_params(1,nml_ok)  !should be read first
  call read_amr_params(1,nml_ok)
  call read_output_params(1,nml_ok)
  call read_movie_params(1,nml_ok)
  call read_lightcone_params(1,nml_ok)
  call read_tracer_params(1,nml_ok)
  call read_poisson_params(1,nml_ok)

  call read_hydro_params(nml_ok)
#ifdef RT
  call read_rt_params(nml_ok)
#endif
#if NDIM==3
  if (sink)call read_sink_params
  if (clumpfind .or. sink)call read_clumpfind_params
  if (stellar)call read_stellar_params
  if (unbind)call read_unbinding_params
  if (make_mergertree)call read_mergertree_params
#if USE_TURB==1
  call read_turb_params(nml_ok)
#endif
#endif

  ! DEV INFO: add here your call for new namelist blocks
  call read_dice_params(1,nml_ok)

  ! Close namelist
  close(1)

  if(.not. nml_ok)then
     if(myid==1)write(*,*)'Too many errors in the namelist'
     if(myid==1)write(*,*)'Aborting...'
     call clean_stop
  end if

  !-------------------------------------------------
  ! LIGHT MPI memory usage message
  !-------------------------------------------------

#ifndef WITHOUTMPI
#ifdef LIGHT_MPI_COMM
  if(myid==1 .and. ncpu .gt. 100) then
    write(*,*) "--------------------------------------------------------------------------------------------------------------"
    write(*,*) "> Using Light MPI Communicator data structures to reduce memory footprint advocated by P. Wautelet (IDRIS) in"
    write(*,*) "  http://www.idris.fr/docs/docu/support-avance/ramses.html"
    write(*,*) ""

    allocate(emission_reception_legacy(1:100, 1:nlevelmax))
    mem_used_legacy_buff = dble(sizeof(emission_reception_legacy)*2)*ncpu/100.0
    deallocate(emission_reception_legacy)
    write(*,*) "  * Old MPI communication structures (emission+reception) would have allocated : ", mem_used_legacy_buff/1.0e6," MB"
    write(*,*) "      - reception(1:ncpu,1:nlevelmax) : ", mem_used_legacy_buff/2.0e6," MB"
    write(*,*) "      - emission(1:ncpu,1:nlevelmax)  : ", mem_used_legacy_buff/2.0e6," MB"
    if (poisson) then
        allocate(emission_reception_legacy(1:100, 1:nlevelmax-1))
        mem_used_legacy_buff_mg = dble(sizeof(emission_reception_legacy)*2)*ncpu/100.0
        deallocate(emission_reception_legacy)
        write(*,*) "  * Old Poisson-related MPI communication structures (active_mg+emission_mg) would have allocated : ", mem_used_legacy_buff_mg/1.0e6," MB"
        write(*,*) "      - active_mg(1:ncpu,1:nlevelmax-1) : ", mem_used_legacy_buff_mg/2.0e6," MB"
        write(*,*) "      - emission_mg(1:ncpu,1:nlevelmax-1) : ", mem_used_legacy_buff_mg/2.0e6," MB"
        mem_used_legacy_buff = mem_used_legacy_buff + mem_used_legacy_buff_mg
    endif

    allocate(reception(1:100, 1:nlevelmax))
    allocate(emission(1:nlevelmax))
    allocate(emission_part(1:nlevelmax))
    mem_used_new_buff = dble(sizeof(emission)) + dble(sizeof(reception))*ncpu/100.0 + dble(sizeof(emission_part))
    deallocate(reception)
    deallocate(emission)
    deallocate(emission_part)
    write(*,*) "  * New MPI communication structures (emission+reception) use : ", mem_used_new_buff/1.0e6," MB"
    write(*,*) "       - emission(1:nlevelmax)         : ", dble(sizeof(emission))/1.0e6," MB"
    write(*,*) "       - emission_part(1:nlevelmax)    : ", dble(sizeof(emission_part))/1.0e6," MB"
    write(*,*) "       - reception(1:ncpu,1:nlevelmax) : ", dble(sizeof(reception))*ncpu/1.0e8," MB"
    if (poisson) then
        allocate(reception(1:100, 1:nlevelmax-1)) ! active_mg
        allocate(emission(1:nlevelmax-1)) ! emission_mg
        mem_used_new_buff_mg = dble(sizeof(emission)) + dble(sizeof(reception))*ncpu/100.0
        deallocate(reception)
        deallocate(emission)
        write(*,*) "  * New Poisson-related MPI communication structures (emission_mg+active_mg) use : ", mem_used_new_buff_mg/1.0e6," MB"
        write(*,*) "       - emission_mg(1:nlevelmax-1)         : ", dble(sizeof(emission))/1.0e6," MB"
        write(*,*) "       - active_mg(1:ncpu,1:nlevelmax-1) : ", dble(sizeof(reception))*ncpu/1.0e8," MB"
        mem_used_new_buff = mem_used_new_buff + mem_used_new_buff_mg
    endif
    write(*,*) "    => Overall memory economy : ", (mem_used_legacy_buff-mem_used_new_buff)/1.0e6,"MB"
    write(*,*) "--------------------------------------------------------------------------------------------------------------"
  endif
#endif
#endif

  !-------------------------------------------------
  ! Read optional nrestart command-line argument
  !-------------------------------------------------

  ! Will overwrite whatever was set in the namelist
  if (myid==1 .and. narg == 2) then
     CALL getarg(2,cmdarg)
     read(cmdarg,*) nrestart
  endif

  ! Check if info file of restart output exists,
  ! otherwise look for earlier outputs
  if (myid==1 .and. nrestart .gt. 0) then
     call title(nrestart,nchar)
     info_file='output_'//TRIM(nchar)//'/info_'//TRIM(nchar)//'.txt'
     inquire(file=info_file, exist=info_ok)
     do while(.not. info_ok .and. nrestart .gt. 1)
        nrestart = nrestart - 1
        call title(nrestart,nchar)
        info_file='output_'//TRIM(nchar)//'/info_'//TRIM(nchar)//'.txt'
        inquire(file=info_file, exist=info_ok)
     enddo
  endif

#ifndef WITHOUTMPI
  call MPI_BCAST(info_ok,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
#endif

  if (nrestart .gt. 0 .and. .not. info_ok) then
     if (myid==1) then
         write(*,*) "Error: Could not find restart file"
     endif
     call clean_stop
  endif

#ifndef WITHOUTMPI
  call MPI_BCAST(nrestart,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif


  ! Send the token
#ifndef WITHOUTMPI
  if(IOGROUPSIZE>0) then
     if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
        dummy_io=1
        call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
             & MPI_COMM_WORLD,info2)
     end if
  endif
#endif

  !-----------------------------------
  ! Rearrange level dependent arrays
  !-----------------------------------
  do i=nlevelmax,levelmin,-1
     nexpand   (i)=nexpand   (i-levelmin+1)
     nsubcycle (i)=nsubcycle (i-levelmin+1)
     r_refine  (i)=r_refine  (i-levelmin+1)
     a_refine  (i)=a_refine  (i-levelmin+1)
     b_refine  (i)=b_refine  (i-levelmin+1)
     x_refine  (i)=x_refine  (i-levelmin+1)
     y_refine  (i)=y_refine  (i-levelmin+1)
     z_refine  (i)=z_refine  (i-levelmin+1)
     m_refine  (i)=m_refine  (i-levelmin+1)
     exp_refine(i)=exp_refine(i-levelmin+1)
     initfile  (i)=initfile  (i-levelmin+1)
  end do
  do i=1,levelmin-1
     nexpand   (i)= 1
     nsubcycle (i)= 1
     r_refine  (i)=-1
     a_refine  (i)= 1
     b_refine  (i)= 1
     x_refine  (i)= 0
     y_refine  (i)= 0
     z_refine  (i)= 0
     m_refine  (i)=-1
     exp_refine(i)= 2
     initfile  (i)= ' '
  end do

#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

end subroutine read_params
!###############################################################
!###############################################################
!###############################################################
subroutine read_run_params(namelist_unit,nml_ok)
   use amr_parameters
   use amr_commons
   use hydro_parameters, only:nhydro,nvar
   implicit none
   integer,intent(in)::namelist_unit
   logical,intent(inout)::nml_ok
   integer::nml_err
#if NVAR>NHYDRO
   integer::i
#endif

   ! Namelist definition for general run parameters
   namelist/run_params/clumpfind,cosmo,pic,sink,tracer,lightcone,poisson,hydro,rt,verbose,debug &
   & ,nrestart,ncontrol,nstepmax,nsubcycle,nremap,ordering &
   & ,bisec_tol,static,overload,cost_weighting,aton,nrestart_quad,restart_remap &
   & ,static_dm,static_gas,static_stars,convert_birth_times,use_proper_time,remap_pscalar &
   & ,unbind,make_mergertree,stellar

  ! Default passive scalar map
#if NVAR>NHYDRO
   allocate(remap_pscalar(1:nvar-nhydro))
   do i=1,nvar-nhydro
      remap_pscalar(i) = i+nhydro
   enddo
#endif

   ! Go to the beginning of the file
   rewind(namelist_unit)

   ! Read namelist
   read(namelist_unit,NML=run_params,IOSTAT=nml_err)

   if(nml_err<0)then
      ! EOF reached before namelist was found
      if(myid==1)write(*,*)'You need to set up namelist &RUN_PARAMS in parameter file.'
      nml_ok=.false.
   elseif(nml_err>0)then
      if(myid==1)write(*,*)'Error reading namelist &RUN_PARAMS. Check formatting.'
      nml_ok=.false.
   endif

   ! Verify input
   if(myid>1)verbose=.false.

   if(stellar.and.(.not.sink))then
      if(myid==1)write(*,*)'Error in the namelist: sink=.true. is needed if stellar=.true. !'
      nml_ok=.false.
   endif
   if(sink.and.(.not.pic))then
      pic=.true.
   endif
   !if(clumpfind.and.(.not.pic))then
   !   pic=.true.
   !endif
   !if(pic.and.(.not.poisson))then
   !   poisson=.true.
   !endif

   if(.not.cosmo)then
      use_proper_time=.false.
      convert_birth_times=.false.
   endif

end subroutine read_run_params
!###############################################################
!###############################################################
!###############################################################
subroutine read_amr_params(namelist_unit,nml_ok)
   use amr_parameters
   use amr_commons, only:myid,ncpu
   use pm_parameters, only:npartmax
   implicit none
   integer,intent(in)::namelist_unit
   logical,intent(inout)::nml_ok
   integer::nml_err
   integer::levelmax
   integer(kind=8)::ngridtot=0 !total number of grid among all MPI processes
   integer(kind=8)::nparttot=0 !total number of particles
#if NDIM==1
   integer, parameter :: max_level_wout_quadhilbert = 61
#elif NDIM==2
   integer, parameter :: max_level_wout_quadhilbert = 29
#elif NDIM==3
   integer, parameter :: max_level_wout_quadhilbert = 19
#endif

   !TODO replace nlevelmax by levelmax everywhere

   ! AMR grid parameters
   namelist/amr_params/levelmin,levelmax,ngridmax,ngridtot &
   & ,npartmax,nparttot,nexpand,boxlen,nlevel_collapse

   ! Go to the beginning of the file
   rewind(namelist_unit)

   ! Read namelist
   read(namelist_unit,NML=amr_params,IOSTAT=nml_err)

   if(nml_err<0)then
      ! EOF reached before namelist was found
      if(myid==1)write(*,*)'You need to set up namelist &AMR_PARAMS in parameter file.'
      nml_ok=.false.
   elseif(nml_err>0)then
      if(myid==1)write(*,*)'Error reading namelist &INIT_PARAMS. Check formatting.'
      nml_ok=.false.
   endif

   ! Verify input
   levelmin=MAX(levelmin,1)
   nlevelmax=levelmax
   if(nlevelmax<levelmin)then
      if(myid==1)write(*,*)'Error in the namelist:'
      if(myid==1)write(*,*)'levelmax should not be lower than levelmin'
      nml_ok=.false.
   end if
   if(ngridmax==0)then
      if(ngridtot==0)then
         if(myid==1)write(*,*)'Error in the namelist: ngridmax and ngridtot are 0!'
         if(myid==1)write(*,*)'Allocate some space for refinements !!!'
         nml_ok=.false.
      else
         ngridmax=int(ngridtot/int(ncpu,kind=8),kind=4)
      endif
   end if

   if(pic .and. npartmax==0) then
      npartmax=int(nparttot/int(ncpu,kind=8),kind=4)
   else if (.not. pic) then
      npartmax=0
   end if

   if(nlevelmax>MAXLEVEL)then
      write(*,*) 'Error: nlevelmax>MAXLEVEL'
      call clean_stop
   endif

#ifndef QUADHILBERT
   if(nlevelmax>=max_level_wout_quadhilbert) then
      if (myid == 1) then
         write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
         write(*,"(a,i2,a)")"WARNING: running with nlevelmax>=", max_level_wout_quadhilbert, " will likely fail."
         write(*,*)"It is recommended to compiling with -DQUADHILBERT"
         write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
      end if
   end if
#endif

end subroutine read_amr_params
!###############################################################
!###############################################################
!###############################################################
subroutine read_output_params(namelist_unit,nml_ok)
   use amr_parameters
   use amr_commons, only:myid
   implicit none
   integer,intent(in)::namelist_unit
   logical,intent(inout)::nml_ok
   integer::nml_err
   real(kind=8)::tend=0  ! end time for the simulation
   real(kind=8)::aend=0  ! end expansion factor

   ! Parameters to specify when to write an output
   namelist/output_params/noutput,foutput,aout,tout &
   & ,tend,delta_tout,aend,delta_aout,gadget_output,walltime_hrs,minutes_dump,output_to_log

   ! Go to the beginning of the file
   rewind(namelist_unit)

   ! Read namelist
   read(namelist_unit,NML=output_params,IOSTAT=nml_err)

   if(nml_err>0)then
      if(myid==1)write(*,*)'Error reading namelist &OUTPUT_PARAMS. Check formatting.'
      nml_ok=.false.
   endif

   !-------------------------------------------------
   ! Compute time step for outputs
   !-------------------------------------------------
   ! check how many predetermined output times are listed (either give tout or aout)
   if(noutput==0.and..not.all(tout==HUGE(1.0D0)))then
      do while(tout(noutput+1)<HUGE(1.0D0))
         noutput = noutput+1
      enddo
   endif
   if(noutput==0.and..not.all(aout==HUGE(1.0D0)))then
      do while(aout(noutput+1)<HUGE(1.0D0))
         noutput = noutput+1
      enddo
   endif
   ! add final time and expansion factor at the back of the predetermined output list
   if(tend>0)then
      noutput=noutput+1
      tout(noutput)=tend
   endif
   if(aend>0)then
      noutput=noutput+1
      aout(noutput)=aend
   endif
   ! set periodic output params
   tout_next=delta_tout
   aout_next=delta_aout

end subroutine read_output_params
!###############################################################
!###############################################################
!###############################################################
subroutine read_movie_params(namelist_unit,nml_ok)
   use amr_parameters
   use amr_commons, only:myid
   implicit none
   integer,intent(in)::namelist_unit
   logical,intent(inout)::nml_ok
   integer::nml_err,i

   ! Namelist parameters for making movies on the fly
   namelist/movie_params/levelmax_frame,nw_frame,nh_frame,ivar_frame &
   & ,xcentre_frame,ycentre_frame,zcentre_frame &
   & ,deltax_frame,deltay_frame,deltaz_frame,movie,zoom_only_frame &
   & ,imovout,imov,tstartmov,astartmov,tendmov,aendmov,proj_axis,movie_vars_txt &
   & ,theta_camera,phi_camera,dtheta_camera,dphi_camera,focal_camera,dist_camera,ddist_camera &
   & ,perspective_camera,smooth_frame,shader_frame,tstart_theta_camera,tstart_phi_camera &
   & ,tend_theta_camera,tend_phi_camera,method_frame,varmin_frame,varmax_frame

   ! Go to the beginning of the file
   rewind(namelist_unit)

   ! Read namelist
   read(namelist_unit,NML=movie_params,IOSTAT=nml_err)

   if(nml_err>0)then
      if(myid==1)write(*,*)'Error reading namelist &MOVIE_PARAMS. Check formatting.'
      nml_ok=.false.
   endif

   ! Setup up movie times
   if(imovout>0) then
      allocate(tmovout(0:imovout))
      allocate(amovout(0:imovout))
      tmovout=1d100
      amovout=1d100
      if(tendmov>0)then
         do i=0,imovout
            tmovout(i)=(tendmov-tstartmov)*dble(i)/dble(imovout)+tstartmov
         enddo
      endif
      if(aendmov>0)then
         do i=0,imovout
            amovout(i)=(aendmov-astartmov)*dble(i)/dble(imovout)+astartmov
         enddo
      endif
      if(tendmov==0.and.aendmov==0)movie=.false.
   endif

   if (movie)call set_movie_vars

end subroutine read_movie_params
!###############################################################
!###############################################################
!###############################################################
subroutine read_lightcone_params(namelist_unit,nml_ok)
   use amr_parameters
   use amr_commons, only:myid
   implicit none
   integer,intent(in)::namelist_unit
   logical,intent(inout)::nml_ok
   integer::nml_err

   namelist/lightcone_params/thetay_cone,thetaz_cone,zmax_cone

   ! Go to the beginning of the file
   rewind(namelist_unit)

   ! Read namelist
   read(namelist_unit,NML=lightcone_params,IOSTAT=nml_err)

   if(nml_err>0)then
      if(myid==1)write(*,*)'Error reading namelist &LIGHTCONE_PARAMS. Check formatting.'
      nml_ok=.false.
   endif

end subroutine read_lightcone_params
!###############################################################
!###############################################################
!###############################################################
subroutine read_tracer_params(namelist_unit,nml_ok)
   use amr_parameters, only:tracer,mc_tracer,pic,nlevelmax
   use amr_commons, only:myid
   use pm_parameters
   implicit none
   integer,intent(in)::namelist_unit
   logical,intent(inout)::nml_ok
   integer::nml_err

   namelist/tracer_params/MC_tracer,tracer_feed,tracer_feed_fmt &
   & ,tracer_mass,tracer_first_balance_part_per_cell &
   & ,tracer_first_balance_levelmin

   ! Go to the beginning of the file
   rewind(namelist_unit)

   ! Read namelist
   read(namelist_unit,NML=tracer_params,IOSTAT=nml_err)

   if(nml_err>0)then
      if(myid==1)write(*,*)'Error reading namelist &TRACER_PARAMS. Check formatting.'
      nml_ok=.false.
   endif

   ! Verify input
   if (tracer_first_balance_levelmin <= 0) tracer_first_balance_levelmin = nlevelmax + 1

   if(MC_tracer .and. (.not. tracer))then
      if(myid==1)write(*,*)'Error: you have activated the MC tracer but not the tracers in RUN_PARAMS.'
      nml_ok=.false.
   end if

   if(MC_tracer .and. (.not. pic)) then
      if(myid==1)write(*,*)'Error: you have activated the MC tracer but pic is false.'
      nml_ok=.false.
   end if

end subroutine read_tracer_params
!###############################################################
!###############################################################
!###############################################################
subroutine read_poisson_params(namelist_unit,nml_ok)
   use amr_commons, only:myid
   use poisson_parameters
   implicit none
   integer,intent(in)::namelist_unit
   logical,intent(inout)::nml_ok
   integer::nml_err

   namelist/poisson_params/epsilon,gravity_type,gravity_params &
   & ,cg_levelmin,cic_levelmax

   ! Go to the beginning of the file
   rewind(namelist_unit)

   ! Read namelist
   read(namelist_unit,NML=poisson_params,IOSTAT=nml_err)

   if(nml_err>0)then
      if(myid==1)write(*,*)'Error reading namelist &POISSON_PARAMS. Check formatting.'
      nml_ok=.false.
   endif

end subroutine read_poisson_params
!###############################################################
!###############################################################
!###############################################################
subroutine read_dice_params(namelist_unit,nml_ok)
   use amr_commons, only:myid
   use dice_commons
   implicit none
   integer,intent(in)::namelist_unit
   logical,intent(inout)::nml_ok
   integer::nml_err

   namelist/dice_params/ ic_file,ic_nfile,ic_format,IG_rho,IG_T2,IG_metal &
       & ,ic_head_name,ic_pos_name,ic_vel_name,ic_id_name,ic_mass_name &
       & ,ic_u_name,ic_metal_name,ic_age_name &
       & ,gadget_scale_l, gadget_scale_v, gadget_scale_m ,gadget_scale_t &
       & ,ic_scale_pos,ic_scale_vel,ic_scale_mass,ic_scale_u,ic_scale_age &
       & ,ic_scale_metal,ic_center,ic_ifout,amr_struct,ic_t_restart,ic_mag_const &
       & ,ic_mag_center_x,ic_mag_center_y,ic_mag_center_z &
       & ,ic_mag_axis_x,ic_mag_axis_y,ic_mag_axis_z &
       & ,ic_mag_scale_R,ic_mag_scale_H,ic_mag_scale_B,cosmo_add_gas_index,ic_skip_type &
       & ,ic_mask_ivar,ic_mask_min,ic_mask_max,ic_mask_ptype

   ! Go to the beginning of the file
   rewind(namelist_unit)

   ! Read namelist
   read(namelist_unit,NML=dice_params,IOSTAT=nml_err)

   if(nml_err>0)then
      if(myid==1)write(*,*)'Error reading namelist &DICE_PARAMS. Check formatting.'
      nml_ok=.false.
   endif

end subroutine read_dice_params
