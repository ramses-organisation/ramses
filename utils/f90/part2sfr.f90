program part2sfr
  use utils
  !--------------------------------------------------------------------------
  ! This program calculates SFR rate of star particls in a given number of
  ! bins from an output of RAMSES simulation 
  ! Version F90 par R. Teyssier le 01/04/01.
  !--------------------------------------------------------------------------
  implicit none
  integer::ncpu,ndim,npart,ngrid,n,i,j,k,icpu,ipos,nstar
  integer::ncpu2,npart2,ndim2,levelmin,levelmax,ilevel,iii
  integer::nx=0,ny=0,ix,iy,ixp1,iyp1,idim,jdim,ncpu_read,n_frw
  real(KIND=8)::mtot,ddx,ddy,dex,dey,time,time_tot,time_simu,weight
  real(KIND=8)::xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1,mmax=1d10
  real(KIND=8)::xcenter,ycenter,zcenter
  real(KIND=8)::age,birth_date
  real(KIND=8)::jxin=0,jyin=0,jzin=0,jx,jy,jz
  real(KIND=8)::kxin,kyin,kzin,kx,ky,kz
  real(KIND=8)::lxin,lyin,lzin,lx,ly,lz
  integer::imin,imax,jmin,jmax,kmin,kmax,lmin,npart_actual
  real(KIND=8)::xxmin,xxmax,yymin,yymax,dx,dy,deltax,boxlen
  real(KIND=8)::aexp,t,omega_m,omega_l,omega_b,omega_k,h0,unit_l,unit_t,unit_d,unit_m
  real(KIND=4),dimension(:,:),allocatable::toto
  real(KIND=4),dimension(:),allocatable::density
  real(KIND=8),dimension(:),allocatable::aexp_frw,hexp_frw,tau_frw,t_frw
  real(KIND=8),dimension(:),allocatable::sfr
  real(KIND=8),dimension(:,:),allocatable::x
  real(KIND=8),dimension(:)  ,allocatable::m,birth
  integer,dimension(:)  ,allocatable::id
  character(LEN=5)::nchar,ncharcpu
  character(LEN=80)::ordering,format_grille
  character(LEN=128)::nomfich,repository,outfich,filedens
  logical::ok,ok_part,periodic=.false.,star=.false.,ageweight=.false.,do_density=.false.
  integer::impi,ndom,bit_length,maxdom
  integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max
  real(KIND=8),dimension(1:8)::bounding_min,bounding_max
  real(KIND=8)::dkey,order_min,dmax
  real(kind=8)::xx,yy,zz
  real(kind=8),dimension(:),allocatable::bound_key
  logical,dimension(:),allocatable::cpu_read
  integer,dimension(:),allocatable::cpu_list
  logical::cosmo=.true.

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
  ! SFR parameters
  !-----------------------
  if(nx==0)nx=1000
  write(*,*)'Working array =',nx
  allocate(sfr(0:nx))
  sfr=0.0d0

  !-----------------------
  ! Get cpu list
  !-----------------------
  if(TRIM(ordering).eq.'hilbert')then

     dmax=max(xmax-xmin,ymax-ymin,zmax-zmin)
     do ilevel=1,levelmax
        deltax=0.5d0**ilevel
        if(deltax.lt.dmax)exit
     end do
     lmin=ilevel
     bit_length=lmin-1
     maxdom=2**bit_length
     imin=0; imax=0; jmin=0; jmax=0; kmin=0; kmax=0
     if(bit_length>0)then
        imin=int(xmin*dble(maxdom))
        imax=imin+1
        jmin=int(ymin*dble(maxdom))
        jmax=jmin+1
        kmin=int(zmin*dble(maxdom))
        kmax=kmin+1
     endif

     dkey=(dble(2**(levelmax+1)/dble(maxdom)))**ndim
     ndom=1
     if(bit_length>0)ndom=8
     idom(1)=imin; idom(2)=imax
     idom(3)=imin; idom(4)=imax
     idom(5)=imin; idom(6)=imax
     idom(7)=imin; idom(8)=imax
     jdom(1)=jmin; jdom(2)=jmin
     jdom(3)=jmax; jdom(4)=jmax
     jdom(5)=jmin; jdom(6)=jmin
     jdom(7)=jmax; jdom(8)=jmax
     kdom(1)=kmin; kdom(2)=kmin
     kdom(3)=kmin; kdom(4)=kmin
     kdom(5)=kmax; kdom(6)=kmax
     kdom(7)=kmax; kdom(8)=kmax

     do i=1,ndom
        if(bit_length>0)then
           call hilbert3d(idom(i),jdom(i),kdom(i),order_min,bit_length,1)
        else
           order_min=0.0d0
        endif
        bounding_min(i)=(order_min)*dkey
        bounding_max(i)=(order_min+1.0D0)*dkey
     end do
     cpu_min=0; cpu_max=0
     do impi=1,ncpu
        do i=1,ndom
           if (   bound_key(impi-1).le.bounding_min(i).and.&
                & bound_key(impi  ).gt.bounding_min(i))then
              cpu_min(i)=impi
           endif
           if (   bound_key(impi-1).lt.bounding_max(i).and.&
                & bound_key(impi  ).ge.bounding_max(i))then
              cpu_max(i)=impi
           endif
        end do
     end do
     ncpu_read=0
     do i=1,ndom
        do j=cpu_min(i),cpu_max(i)
           if(.not. cpu_read(j))then
              ncpu_read=ncpu_read+1
              cpu_list(ncpu_read)=j
              cpu_read(j)=.true.
           endif
        enddo
     enddo
  else
     ncpu_read=ncpu
     do j=1,ncpu
        cpu_list(j)=j
     end do
  end  if

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
     allocate(id(1:npart2))
     allocate(x(1:npart2,1:ndim2))
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
     read(1) ! Skip family
     read(1) ! Skip tag
     read(1)birth
     close(1)

     do i=1,npart2
        ok_part=(x(i,1)>=xmin.and.x(i,1)<=xmax.and. &
             &   x(i,2)>=ymin.and.x(i,2)<=ymax.and. &
             &   x(i,3)>=zmin.and.x(i,3)<=zmax)
        ok_part=ok_part.and.(birth(i).ne.0.0d0)

        if(ok_part)then

           if(cosmo)then
              iii=1
              do while(tau_frw(iii)>birth(i).and.iii<n_frw)
                 iii=iii+1
              end do
              ! Interploate time
              time=t_frw(iii)*(birth(i)-tau_frw(iii-1))/(tau_frw(iii)-tau_frw(iii-1))+ &
                   & t_frw(iii-1)*(birth(i)-tau_frw(iii))/(tau_frw(iii-1)-tau_frw(iii))
              age=(time_simu-time)/(h0*1d5/3.08d24)/(365.*24.*3600.*1d9)
              birth_date=(time_tot+time)/(h0*1d5/3.08d24)/(365.*24.*3600.*1d9)
           else
              age=(time_simu-birth(i))*unit_t/(365.*24.*3600.*1d9)
              birth_date=birth(i)*unit_t/(365.*24.*3600.*1d9)
           endif

           npart_actual=npart_actual+1
           mtot=mtot+m(i)
           if(birth_date>0.and.birth_date<15.)then
              ix=int(birth_date/15.*dble(nx))
              sfr(ix)=sfr(ix)+m(i)*unit_m/2d33
           endif

        end if

     end do
     deallocate(x,m)
     deallocate(birth,id)
  end do
  write(*,*)'Total mass=',mtot*unit_m/2d33
  write(*,*)'npart tot=',npart_actual

  ! Output file
  nomfich=TRIM(outfich)
  write(*,*)'Ecriture des donnees du fichier '//TRIM(nomfich)
  open(unit=10,file=nomfich,form='formatted')
  do i=0,nx
     xx=(dble(i)+0.5)/dble(nx)*15.
     write(10,*)xx,sfr(i)/(15./dble(nx))/1d9
  end do
  close(10)

contains

  subroutine read_params

      implicit none

      integer       :: i,n
      
      character(len=4)   :: opt
      character(len=128) :: arg
      LOGICAL       :: bad, ok

      n = command_argument_count()
      if (n < 4) then
         print *, 'usage: part2sfr  -inp  input_dir'
         print *, '                 -out  output_file'
         print *, '                 [-xmi xmin] '
         print *, '                 [-xma xmax] '
         print *, '                 [-ymi ymin] '
         print *, '                 [-yma ymax] '
         print *, '                 [-zmi zmin] '
         print *, '                 [-zma zmax] '
         print *, '                 [-nx  nx  ] '
         print *, 'ex: part2sfr -inp output_00001 -out sfr.dat'// &
              &   ' -xmi 0.1 -xma 0.7'
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
         case ('-out')
            outfich = trim(arg)
         case ('-xmi')
            read (arg,*) xmin
         case ('-xma')
            read (arg,*) xmax
         case ('-ymi')
            read (arg,*) ymin
         case ('-yma')
            read (arg,*) ymax
         case ('-zmi')
            read (arg,*) zmin
         case ('-zma')
            read (arg,*) zmax
         case ('-nx')
            read (arg,*) nx
         case default
            print '("unknown option ",a2," ignored")', opt
         end select
      end do

      return

    end subroutine read_params

  end program part2sfr

!=======================================================================
subroutine title(n,nchar)
!=======================================================================
  implicit none
  integer::n
  character*5::nchar

  character*1::nchar1
  character*2::nchar2
  character*3::nchar3
  character*4::nchar4
  character*5::nchar5

  if(n.ge.10000)then
     write(nchar5,'(i5)') n
     nchar = nchar5
  elseif(n.ge.1000)then
     write(nchar4,'(i4)') n
     nchar = '0'//nchar4
  elseif(n.ge.100)then
     write(nchar3,'(i3)') n
     nchar = '00'//nchar3
  elseif(n.ge.10)then
     write(nchar2,'(i2)') n
     nchar = '000'//nchar2
  else
     write(nchar1,'(i1)') n
     nchar = '0000'//nchar1
  endif

end subroutine title
