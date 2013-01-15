subroutine read_params
  use amr_commons
  use pm_parameters
  use poisson_parameters
  use openmp_support
#ifndef WITHOUTMPI
  use mpi
#endif
!$ use omp_lib

  implicit none
  !--------------------------------------------------
  ! Local variables
  !--------------------------------------------------
  integer::i,narg,iargc,ierr,levelmax
  character(LEN=80)::infile
  integer(kind=8)::ngridtot=0
  integer(kind=8)::nparttot=0
  real(kind=8)::delta_tout=0,tend=0
  real(kind=8)::delta_aout=0,aend=0
  logical::nml_ok

#ifndef WITHOUTMPI
  integer::mpi_provided
#endif

  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------
  namelist/run_params/cosmo,pic,sink,lightcone,poisson,hydro,verbose,debug &
       & ,nrestart,ncontrol,nstepmax,nsubcycle,nremap,ordering &
       & ,bisec_tol,static,geom,overload,cost_weighting
  namelist/output_params/noutput,foutput,fbackup,aout,tout,output_mode &
       & ,tend,delta_tout,aend,delta_aout
  namelist/amr_params/levelmin,levelmax,ngridmax,ngridtot &
       & ,npartmax,nparttot,nsinkmax,nexpand,boxlen
  namelist/poisson_params/epsilon,gravity_type,gravity_params &
       & ,cg_levelmin,cic_levelmax
  namelist/lightcone_params/thetay_cone,thetaz_cone,zmax_cone
!!$  namelist/movie_params/f_frame,nx_frame,ny_frame,ivar_frame &
!!$       & ,xmin_frame,xmax_frame,ymin_frame,ymax_frame,zmin_frame,zmax_frame

  ! MPI initialization
#ifndef WITHOUTMPI

#ifdef _OPENMP
  call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE,mpi_provided,ierr)
#else
  call MPI_INIT(ierr)
#endif

  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu,ierr)
  myid=myid+1 ! Carefull with this...

#ifdef _OPENMP
  if (mpi_provided==MPI_THREAD_SINGLE.and.myid==1) then
    print *, 'MPI implementation provides support for MPI calls outside openmp regions'
    print *, 'This is the minimum requirement for this version of ramses'
  endif
  if (mpi_provided==MPI_THREAD_FUNNELED.and.myid==1) then
    print *, 'MPI implementation provides support for MPI calls in master regions'
    print *, 'This is more than what ramses requires in this version.'
  endif
  if (mpi_provided==MPI_THREAD_SERIALIZED.and.myid==1) then
    print *, 'MPI implementation provides support for serial MPI calls by any thread'
    print *, 'This is more than what ramses requires in this version.'
  endif
  if (mpi_provided==MPI_THREAD_MULTIPLE .and. myid==1) then
    print *, 'MPI implementation provides support for simultaneous MPI calls by different threads'
    print *, 'This is more than what ramses requires in this version.'
 endif
 if (myid==1) print*,' ' 
#endif

#else  ! WITHOUTMPI
  ncpu=1
  myid=1
#endif

  ! OpenMP initialization
  call init_openmp

  !--------------------------------------------------
  ! Advertise RAMSES
  !--------------------------------------------------
  if(myid==1)then

  write(*,*)'_/_/_/       _/_/     _/    _/    _/_/_/   _/_/_/_/    _/_/_/  '
  write(*,*)'_/    _/    _/  _/    _/_/_/_/   _/    _/  _/         _/    _/ '
  write(*,*)'_/    _/   _/    _/   _/ _/ _/   _/        _/         _/       '
  write(*,*)'_/_/_/     _/_/_/_/   _/    _/     _/_/    _/_/_/       _/_/   '
  write(*,*)'_/    _/   _/    _/   _/    _/         _/  _/               _/ '
  write(*,*)'_/    _/   _/    _/   _/    _/   _/    _/  _/         _/    _/ '
  write(*,*)'_/    _/   _/    _/   _/    _/    _/_/_/   _/_/_/_/    _/_/_/  '
  write(*,*)'                        Version 3.0                            '
  write(*,*)'       written by Romain Teyssier (CEA/DSM/IRFU/SAP)           '
  write(*,*)'                     (c) CEA 1999-2007                         '
  write(*,*)' '
  write(*,'(" Working in ",I1," dimensions")') ndim
  if (ncpu > 1)   write(*,'("  # of MPI tasks           :",I8)') ncpu
  if (omp_nthreads > 1) then
                  write(*,'("  # of OpenMP threads      :",I8)') omp_nthreads
    if (ncpu > 1) write(*,'("  # of cores used in total :",I8)') ncore
  endif
  write(*,*)' '

  ! Read namelist filename from command line argument
  narg = iargc()
  IF(narg .LT. 1)THEN
     write(*,*)'You should type: hydro3d input.nml'
     write(*,*)'File input.nml should contain a parameter namelist'
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
  INQUIRE(file=infile,exist=nml_ok)
  if(.not. nml_ok)then
     if(myid==1)then
        write(*,*)'File '//TRIM(infile)//' does not exist'
     endif
     call clean_stop
  end if

  open(1,file=infile)
  rewind(1)
  read(1,NML=run_params)
  rewind(1)
  read(1,NML=output_params)
  rewind(1)
  read(1,NML=amr_params)
  rewind(1)
  read(1,NML=lightcone_params,END=83)
83 continue
  rewind(1)
!!$  read(1,NML=movie_params,END=82)
!!$82 continue
!!$  rewind(1)
  read(1,NML=poisson_params,END=81)
81 continue
  !-------------------------------------------------
  ! Compute time step for outputs
  !-------------------------------------------------
  if(tend>0)then
     if(delta_tout==0)delta_tout=tend
     noutput=int(tend/delta_tout)
     do i=1,noutput
        tout(i)=dble(i)*delta_tout
     end do
  else if(aend>0)then
     if(delta_aout==0)delta_aout=aend
     noutput=int(aend/delta_aout)
     do i=1,noutput
        aout(i)=dble(i)*delta_aout
     end do
  endif
  !--------------------------------------------------
  ! Check for errors in the namelist so far
  !--------------------------------------------------
  levelmin=MAX(levelmin,1)
  nlevelmax=levelmax
  nml_ok=.true.
  if(levelmin<1)then
     if(myid==1)write(*,*)'Error in the namelist:'
     if(myid==1)write(*,*)'levelmin should not be lower than 1 !!!'
     nml_ok=.false.
  end if
  if(nlevelmax<levelmin)then
     if(myid==1)write(*,*)'Error in the namelist:'
     if(myid==1)write(*,*)'levelmax should not be lower than levelmin'
     nml_ok=.false.
  end if
  if(ngridmax==0)then
     if(ngridtot==0)then
        if(myid==1)write(*,*)'Error in the namelist:'
        if(myid==1)write(*,*)'Allocate some space for refinements !!!'
        nml_ok=.false.
     else
        ngridmax=ngridtot/int(ncpu,kind=8)
     endif
  end if
  if(npartmax==0)then
     npartmax=nparttot/int(ncpu,kind=8)
  endif
  if(myid>1)verbose=.false.
  if(sink.and.(.not.pic))then
     pic=.true.
  endif
  if(pic.and.(.not.poisson))then
     poisson=.true.
  endif

  call read_hydro_params(nml_ok)

  close(1)

  !-----------------
  ! Max size checks
  !-----------------
  if(nlevelmax>MAXLEVEL)then
     write(*,*) 'Error: nlevelmax>MAXLEVEL'
     call clean_stop
  end if
  if(nregion>MAXREGION)then
     write(*,*) 'Error: nregion>MAXREGION'
     call clean_stop
  end if
  
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
     r_refine  (i)=-1.0
     a_refine  (i)= 1.0
     b_refine  (i)= 1.0
     x_refine  (i)= 0.0
     y_refine  (i)= 0.0
     z_refine  (i)= 0.0
     m_refine  (i)=-1.0
     exp_refine(i)= 2.0
     initfile  (i)= ' '
  end do
     
  if(.not. nml_ok)then
     if(myid==1)write(*,*)'Too many errors in the namelist'
     if(myid==1)write(*,*)'Aborting...'
     call clean_stop
  end if

#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

end subroutine read_params

