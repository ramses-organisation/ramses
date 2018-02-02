program part2birth
  use utils
  !--------------------------------------------------------------------------
  ! This program converts birth times of star particles to Gyr and saves
  ! as `birth` files for an output of RAMSES simulation 
  ! Version F90 par R. Teyssier le 01/04/01.
  !--------------------------------------------------------------------------
  implicit none
  integer::ncpu,ndim,npart,i,j,k,icpu,ipos,nstar
  integer::ncpu2,npart2,ndim2,levelmin,levelmax,iii
  integer::ncpu_read,n_frw
  real(KIND=8)::mtot,time,time_tot,time_simu
  real(KIND=8)::age
  integer::npart_actual
  real(KIND=8)::boxlen
  real(KIND=8)::aexp,t,omega_m,omega_l,omega_b,omega_k,h0,unit_l,unit_t,unit_d,unit_m
  real(KIND=8),dimension(:),allocatable::aexp_frw,hexp_frw,tau_frw,t_frw
  real(KIND=8),dimension(:,:),allocatable::x
  real(KIND=8),dimension(:)  ,allocatable::m,birth,birth_date
  integer,dimension(:)  ,allocatable::id
  character(LEN=5)::nchar,ncharcpu
  character(LEN=80)::ordering
  character(LEN=128)::nomfich,repository
  logical::ok
  integer::impi
  real(kind=8),dimension(:),allocatable::bound_key
  logical,dimension(:),allocatable::cpu_read
  integer,dimension(:),allocatable::cpu_list
  logical::cosmo=.true.
  integer(kind=1),dimension(:),allocatable::family,tag

  call read_params

  !-----------------------------------------------
  ! Lecture du fichier particules au format RAMSES
  !-----------------------------------------------
  ipos=INDEX(repository,'output_')
  nchar=repository(ipos+7:ipos+13)
  nomfich=TRIM(repository)//'/part_'//TRIM(nchar)//'.out00001'
  inquire(file=nomfich, exist=ok) ! verify input file
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     stop
  endif

  nomfich=TRIM(repository)//'/info_'//TRIM(nchar)//'.txt'
  inquire(file=nomfich, exist=ok) ! verify input file
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     stop
  endif
  open(unit=10,file=nomfich,form='formatted',status='old')
  read(10,'("ncpu        =",I11)')ncpu
  read(10,'("ndim        =",I11)')ndim
  read(10,'("levelmin    =",I11)')levelmin
  read(10,'("levelmax    =",I11)')levelmax
  read(10,*)
  read(10,*)
  read(10,*)

  read(10,'("boxlen      =",E23.15)')boxlen
  read(10,'("time        =",E23.15)')t
  read(10,'("aexp        =",E23.15)')aexp
  read(10,'("H0          =",E23.15)')h0
  read(10,'("omega_m     =",E23.15)')omega_m
  read(10,'("omega_l     =",E23.15)')omega_l
  read(10,'("omega_k     =",E23.15)')omega_k
  read(10,'("omega_b     =",E23.15)')omega_b
  read(10,'("unit_l      =",E23.15)')unit_l
  read(10,'("unit_d      =",E23.15)')unit_d
  read(10,'("unit_t      =",E23.15)')unit_t
  unit_m=unit_d*unit_l**3
  read(10,*)

  if(aexp.eq.1.and.h0.eq.1)cosmo=.false.

  read(10,'("ordering type=",A80)'),ordering
  write(*,'(" ordering type=",A20)'),TRIM(ordering)
  read(10,*)
  allocate(cpu_list(1:ncpu))
  if(TRIM(ordering).eq.'hilbert')then
     allocate(bound_key(0:ncpu))
     allocate(cpu_read(1:ncpu))
     cpu_read=.false.
     do impi=1,ncpu
        read(10,'(I8,1X,E23.15,1X,E23.15)')i,bound_key(impi-1),bound_key(impi)
     end do
  endif
  close(10)

  !-----------------------
  ! Cosmological model
  !-----------------------
  if(cosmo)then
     ! Allocate look-up tables
     n_frw=1000
     allocate(aexp_frw(0:n_frw),hexp_frw(0:n_frw))
     allocate(tau_frw(0:n_frw),t_frw(0:n_frw))

     ! Compute Friedman model look up table
     write(*,*)'Computing Friedman model'
     call friedman(dble(omega_m),dble(omega_l),dble(omega_k), &
          & 1.d-6,1.d-3,aexp_frw,hexp_frw,tau_frw,t_frw,n_frw,time_tot)

     ! Find neighboring expansion factors
     i=1
     do while(aexp_frw(i)>aexp.and.i<n_frw)
        i=i+1
     end do
     ! Interploate time
     time_simu=t_frw(i)*(aexp-aexp_frw(i-1))/(aexp_frw(i)-aexp_frw(i-1))+ &
          & t_frw(i-1)*(aexp-aexp_frw(i))/(aexp_frw(i-1)-aexp_frw(i))
     write(*,*)'Time simu=',(time_tot+time_simu)/(h0*1d5/3.08d24)/(365.*24.*3600.*1d9)
     write(*,*)'Hubble time=',(time_tot)/(h0*1d5/3.08d24)/(365.*24.*3600.*1d9)
  else
     time_simu=t
     write(*,*)'Time simu=',time_simu*unit_t/(365.*24.*3600.*1d9)
  endif

  !-----------------------
  ! Get cpu list
  !-----------------------
  ncpu_read=ncpu
  do j=1,ncpu
     cpu_list(j)=j
  end do

  npart=0
  do k=1,ncpu_read
     icpu=cpu_list(k)
     call title(icpu,ncharcpu)
     nomfich=TRIM(repository)//'/part_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
     open(unit=1,file=nomfich,status='old',form='unformatted')
     read(1)ncpu2
     read(1)ndim2
     read(1)npart2
     read(1)
     read(1)nstar
     close(1)
     npart=npart+npart2
  end do
  write(*,*)'Found ',npart,' particles.'
  if(nstar>0)then
     write(*,*)'Keeping star particles.'
  else
     write(*,*)'No star particles.'
     stop
  endif

  !-----------------------------------------------
  ! Compute SFH using histograming of birth dates
  !-----------------------------------------------
  npart_actual=0
  mtot=0.0d0
  do k=1,ncpu_read
     icpu=cpu_list(k)
     call title(icpu,ncharcpu)
     nomfich=TRIM(repository)//'/part_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
     open(unit=1,file=nomfich,status='old',form='unformatted')
     read(1)ncpu2
     read(1)ndim2
     read(1)npart2
     read(1)
     read(1)
     read(1)
     read(1)
     read(1)
     allocate(m(1:npart2))
     allocate(birth(1:npart2))
     allocate(birth_date(1:npart2))
     allocate(id(1:npart2))
     allocate(x(1:npart2,1:ndim2))
     allocate(family(1:npart2))
     allocate(tag(1:npart2))
     ! Read position
     do i=1,ndim
        read(1)m
        x(1:npart2,i)=m/boxlen
     end do
     ! Skip velocity
     do i=1,ndim
        read(1)m
     end do
     ! Read mass
     read(1)m
     read(1)id
     read(1) ! Skip level
     read(1)family
     read(1)tag
     read(1)birth
     close(1)

     do i=1,npart2
        if(family(i)==2)then ! birth would allow the cloud particles to pass
           if(cosmo)then
              iii=1
              do while(tau_frw(iii)>birth(i).and.iii<n_frw)
                 iii=iii+1
              end do
              ! Interploate time
              time=t_frw(iii)*(birth(i)-tau_frw(iii-1))/(tau_frw(iii)-tau_frw(iii-1))+ &
                   & t_frw(iii-1)*(birth(i)-tau_frw(iii))/(tau_frw(iii-1)-tau_frw(iii))
              age=(time_simu-time)/(h0*1d5/3.08d24)/(365.*24.*3600.*1d9)
              birth_date(i)=(time_tot+time)/(h0*1d5/3.08d24)/(365.*24.*3600.*1d9)
           else
              age=(time_simu-birth(i))*unit_t/(365.*24.*3600.*1d9)
              birth_date(i)=birth(i)*unit_t/(365.*24.*3600.*1d9)
           endif
        else
           birth_date(i)=0d0
        endif
     end do
     nomfich=TRIM(repository)//'/birth_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
     open(unit=1,file=nomfich,status='new',form='unformatted')
     write(1)birth_date
     close(1)
     deallocate(x,m)
     deallocate(birth,id,birth_date)
     deallocate(family,tag)
  end do

contains

  subroutine read_params

      implicit none

      integer       :: i,n
      
      character(len=4)   :: opt
      character(len=128) :: arg

      n = command_argument_count()
      if (n < 2) then
         print *, 'usage: part2birth  -inp  input_dir'
         print *, 'ex: part2birth -inp output_00001'
         stop
      end if

      do i = 1,n,2
         call get_command_argument(i,opt)
         if (i == n) then
            print '("option ",a2," has no argument")', opt
            stop 2
         end if
         call get_command_argument(i+1,arg)
         select case (opt)
         case ('-inp')
            repository = trim(arg)
         case default
            print '("unknown option ",a2," ignored")', opt
         end select
      end do

      return

    end subroutine read_params

  end program part2birth
