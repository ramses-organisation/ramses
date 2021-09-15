program part2cube
  !--------------------------------------------------------------------------
  ! This software computes projected density cubes.
  ! Version F90 by R. Teyssier le 01/04/01.
  ! Update for new tracers by C. Cadiou, M. Trebitsch, H. Choi, O. Snaith and R. Bieri

  ! Modified from part2cube.f90 by Eric Moseley. Deposits particle momentum onto a cube.
  !--------------------------------------------------------------------------

  use utils
  implicit none
  integer::ncpu,ndim,npart,i,j,k,icpu,ipos,n_frw,nstar,type
  integer::ncpu2,npart2,ndim2,levelmin,levelmax,ilevel,iii
  integer::nx=0,ny=0,nz=0,ix,iy,iz,ixp1,iyp1,izp1,idim,jdim,kdim,ncpu_read
  real(KIND=8)::ddx,ddy,ddz,dex,dey,dez,t,time,time_tot,time_simu,weight
  real(KIND=8)::ptot=0.0 !momentum total
  real(KIND=8)::xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1
  real(KIND=8)::aexp,omega_m,omega_l,omega_b,omega_k,h0,unit_l,unit_t,unit_d
  integer::imin,imax,jmin,jmax,kmin,kmax,lmin
  real(KIND=8)::xxmin,xxmax,yymin,yymax,zzmin,zzmax,dx,dy,dz,deltax
  real(KIND=4),dimension(:,:,:),allocatable::toto
  real(KIND=8),dimension(:),allocatable::aexp_frw,hexp_frw,tau_frw,t_frw
  real(KIND=8),dimension(:,:,:),allocatable::cube
  real(KIND=8),dimension(:,:),allocatable::x,v
  real(KIND=8),dimension(:)  ,allocatable::m,age
  character(LEN=5)::nchar,ncharcpu
  character(LEN=80)::ordering
  character(LEN=80)::GMGM
  character(LEN=128)::nomfich,repository,outfich
  logical::ok,ok_part,periodic=.false.,star=.false.,ageweight=.false.
  integer::impi,ndom,bit_length,maxdom
  integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max
  real(KIND=8),dimension(1:8)::bounding_min,bounding_max
  real(KIND=8)::dkey,order_min,dmax
  real(kind=8),dimension(:),allocatable::bound_key
  logical,dimension(:),allocatable::cpu_read
  integer,dimension(:),allocatable::cpu_list
  integer(kind=1),dimension(:),allocatable::family,tag

  call read_params

  !-----------------------------------------------
  ! Reading particle files in RAMSES format
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
  ! read(10,'("ncpu        =",I11)')ncpu
  read(10,'(A13,I11)')GMGM,ncpu
  ! read(10,'("ndim        =",I11)')ndim
  read(10,'(A13,I11)')GMGM,ndim
  ! read(10,'("levelmin    =",I11)')levelmin
  read(10,'(A13,I11)')GMGM,levelmin
  ! read(10,'("levelmax    =",I11)')levelmax
  read(10,'(A13,I11)')GMGM,levelmax
  read(10,*)
  read(10,*)
  read(10,*)

  read(10,*)
  ! read(10,'("time        =",E23.15)')t
  read(10,'(A13,E23.15)')GMGM,t

  !  read(10,'("aexp        =",E23.15)')aexp
  read(10,'(A13,E23.15)')GMGM,aexp

  !  read(10,'("H0          =",E23.15)')h0
  read(10,'(A13,E23.15)')GMGM,h0

  !read(10,'("omega_m     =",E23.15)')omega_m
  read(10,'(A13,E23.15)')GMGM,omega_m

  !read(10,'("omega_l     =",E23.15)')omega_l
  read(10,'(A13,E23.15)')GMGM,omega_l

  !read(10,'("omega_k     =",E23.15)')omega_k
  read(10,'(A13,E23.15)')GMGM,omega_k

  !read(10,'("omega_b     =",E23.15)')omega_b
  read(10,'(A13,E23.15)')GMGM,omega_b

  !read(10,'("unit_l      =",E23.15)')unit_l
  read(10,'(A13,E23.15)')GMGM,unit_l

  !read(10,'("unit_d      =",E23.15)')unit_d
  read(10,'(A13,E23.15)')GMGM,unit_d

  !read(10,'("unit_t      =",E23.15)')unit_t
  read(10,'(A13,E23.15)')GMGM,unit_t

  read(10,*)

  ! read(10,'("ordering type=",A80)'),ordering
  read(10,'(A14,A80)')GMGM,ordering
  write(*,'(" ordering type=",A20)')TRIM(ordering)
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
  ! Allocate look-up tables
  n_frw=1000
  allocate(aexp_frw(0:n_frw),hexp_frw(0:n_frw))
  allocate(tau_frw(0:n_frw),t_frw(0:n_frw))

  ! Compute Friedman model look up table
  write(*,*)'Computing Friedman model'
  call friedman(dble(omega_m),dble(omega_l),dble(omega_k), &
       & 1d-6,1d-3,aexp_frw,hexp_frw,tau_frw,t_frw,n_frw,time_tot)

  ! Find neighboring expansion factors
  i=1
  do while(aexp_frw(i)>aexp.and.i<n_frw)
     i=i+1
  end do
  ! Interploate time
  time_simu=t_frw(i)*(aexp-aexp_frw(i-1))/(aexp_frw(i)-aexp_frw(i-1))+ &
       & t_frw(i-1)*(aexp-aexp_frw(i))/(aexp_frw(i-1)-aexp_frw(i))
  write(*,*)'Age simu=',(time_tot+time_simu)/(h0*1d5/3.08d24)/(365.*24.*3600.*1d9)

  !-----------------------
  ! Map parameters
  !-----------------------
  if(nx==0)then
     nx=2**levelmin
  endif
  if(ny==0)then
     ny=nx
  end if
  if(nz==0)then
     nz=nx
  end if
  write(*,*)'time=',t
  write(*,*)'Working cube =',nx,ny,nz
  allocate(cube(0:nx,0:ny,0:nz))
  cube=0.0d0
  idim=1
  jdim=2
  kdim=3
  xxmin=xmin ; xxmax=xmax
  yymin=ymin ; yymax=ymax
  zzmin=zmin ; zzmax=zmax
  dx=(xxmax-xxmin)/dble(nx)
  dy=(yymax-yymin)/dble(ny)
  dz=(zzmax-zzmin)/dble(nz)

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
  write(*,*)'Found ',npart,' particles'
  if(nstar>0)then
     if(star)then
        write(*,*)'Keeping star particles.'
     else
        write(*,*)'Discard star particles.'
     endif
  endif

  !-----------------------------------------------
  ! Compute projected mass using CIC smoothing
  !----------------------------------------------


  ptot=0.0d0

  do k=1,ncpu_read
     icpu=cpu_list(k)
     call title(icpu,ncharcpu)
     nomfich=TRIM(repository)//'/part_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
     open(unit=1,file=nomfich,status='old',form='unformatted')
     write(*,*)'Processing file '//TRIM(nomfich)
     read(1)ncpu2
     read(1)ndim2
     read(1)npart2
     read(1)
     read(1)
     read(1)
     read(1)
     read(1)
     allocate(m(1:npart2))
     if(nstar>0)allocate(age(1:npart2))
     allocate(x(1:npart2,1:ndim2))
     allocate(v(1:npart2,1:ndim2))
     allocate(family(1:npart2))
     allocate(tag(1:npart2))
     ! Read position
     do i=1,ndim
        read(1)m
        x(1:npart2,i)=m
     end do
     ! Skip velocity
     do i=1,ndim
        read(1)m
        v(1:npart2,i)=m
     end do
     ! Read mass
     read(1)m
     if(nstar>0)then
        read(1) ! Skip identity
        read(1) ! Skip level
        read(1)family
        read(1)tag
        read(1)age
     endif

     close(1)
     if(periodic)then
        do i=1,npart2
           ok_part=(x(i,1)>=xmin.and.x(i,1)<=xmax.and. &
                &   x(i,2)>=ymin.and.x(i,2)<=ymax.and. &
                &   x(i,3)>=zmin.and.x(i,3)<=zmax)

           if(nstar>0)then
              if(star)then
                 ! Keep only stars
                 ok_part=ok_part.and.(family(i)==2)
                 if(ageweight)then
                    iii=1
                    do while(tau_frw(iii)>age(i).and.iii<n_frw)
                       iii=iii+1
                    end do
                    ! Interpolate time
                    time=t_frw(iii)*(age(i)-tau_frw(iii-1))/(tau_frw(iii)-tau_frw(iii-1))+ &
                         & t_frw(iii-1)*(age(i)-tau_frw(iii))/(tau_frw(iii-1)-tau_frw(iii))
                    time=(time_simu-time)/(h0*1d5/3.08d24)/(365.*24.*3600.*1d9)
                    if(time>0.01)then
                       weight=(time/0.01)**(-0.7)
                    endif
                 endif
              else
                 ! Keep only DM
                 ok_part=ok_part.and.(family(i).eq.1)
              endif
           endif

           if(ok_part)then
              ddx=(x(i,idim)-xxmin)/dx
              ddy=(x(i,jdim)-yymin)/dy
              ddz=(x(i,kdim)-zzmin)/dz
              ix=ddx
              iy=ddy
              iz=ddz
              ddx=ddx-ix
              ddy=ddy-iy
              ddz=ddz-iz
              dex=1.0-ddx
              dey=1.0-ddy
              dez=1.0-ddz
              if(ix<0)ix=ix+nx
              if(ix>=nx)ix=ix-nx
              if(iy<0)iy=iy+ny
              if(iy>=ny)iy=iy-ny
              if(iz<0)iz=iz+nz
              if(iz>=nz)iz=iz-nz
              ixp1=ix+1
              iyp1=iy+1
              izp1=iz+1
              if(ixp1<0)ixp1=ixp1+nx
              if(ixp1>=nx)ixp1=ixp1-nx
              if(iyp1<0)iyp1=iyp1+ny
              if(iyp1>=ny)iyp1=iyp1-ny
              if(izp1<0)izp1=izp1+nz
              if(izp1>=nz)izp1=izp1-nz

              cube(ix  ,iy  ,iz  )=cube(ix  ,iy  ,iz  )+m(i)*dex*dey*dez*v(i,type)
              cube(ix  ,iyp1,iz  )=cube(ix  ,iyp1,iz  )+m(i)*dex*ddy*dez*v(i,type)
              cube(ixp1,iy  ,iz  )=cube(ixp1,iy  ,iz  )+m(i)*ddx*dey*dez*v(i,type)
              cube(ixp1,iyp1,iz  )=cube(ixp1,iyp1,iz  )+m(i)*ddx*ddy*dez*v(i,type)
              cube(ix  ,iy  ,izp1)=cube(ix  ,iy  ,izp1)+m(i)*dex*dey*ddz*v(i,type)
              cube(ix  ,iyp1,izp1)=cube(ix  ,iyp1,izp1)+m(i)*dex*ddy*ddz*v(i,type)
              cube(ixp1,iy  ,izp1)=cube(ixp1,iy  ,izp1)+m(i)*ddx*dey*ddz*v(i,type)
              cube(ixp1,iyp1,izp1)=cube(ixp1,iyp1,izp1)+m(i)*ddx*ddy*ddz*v(i,type)
              ptot=ptot+m(i)*v(i,type) ! A momentum total

           end if
        end do
     else
        do i=1,npart2
           weight=1.0
           ok_part=(x(i,1)>=xmin.and.x(i,1)<=xmax.and. &
                &   x(i,2)>=ymin.and.x(i,2)<=ymax.and. &
                &   x(i,3)>=zmin.and.x(i,3)<=zmax)

           if(nstar>0)then
              if(star)then
                 ok_part=ok_part.and.(family(i).eq.2)
                 if(ageweight)then
                    iii=1
                    do while(tau_frw(iii)>age(i).and.iii<n_frw)
                       iii=iii+1
                    end do
                    ! Interpolate time
                    time=t_frw(iii)*(age(i)-tau_frw(iii-1))/(tau_frw(iii)-tau_frw(iii-1))+ &
                         & t_frw(iii-1)*(age(i)-tau_frw(iii))/(tau_frw(iii-1)-tau_frw(iii))
                    time=(time_simu-time)/(h0*1d5/3.08d24)/(365.*24.*3600.*1d9)
                    if(time>0.01)then
                       weight=(time/0.01)**(-0.7)
                    endif
                 endif
              else
                 ok_part=ok_part.and.(family(i).eq.1)
              endif
           endif

           if(ok_part)then
              ddx=(x(i,idim)-xxmin)/dx
              ddy=(x(i,jdim)-yymin)/dy
              ddz=(x(i,kdim)-zzmin)/dz
              ix=ddx
              iy=ddy
              iz=ddz
              ddx=ddx-ix
              ddy=ddy-iy
              ddz=ddz-iz
              dex=1.0-ddx
              dey=1.0-ddy
              dez=1.0-ddz
              ixp1=ix+1
              iyp1=iy+1
              izp1=iz+1
              if(ix>=0.and.ix<nx.and.iy>=0.and.iy<ny.and.iz>=0.and.iz<nz)then

              cube(ix  ,iy  ,iz  )=cube(ix  ,iy  ,iz  )+m(i)*dex*dey*dez*weight*v(i,type)
              cube(ix  ,iyp1,iz  )=cube(ix  ,iyp1,iz  )+m(i)*dex*ddy*dez*weight*v(i,type)
              cube(ixp1,iy  ,iz  )=cube(ixp1,iy  ,iz  )+m(i)*ddx*dey*dez*weight*v(i,type)
              cube(ixp1,iyp1,iz  )=cube(ixp1,iyp1,iz  )+m(i)*ddx*ddy*dez*weight*v(i,type)
              cube(ix  ,iy  ,izp1)=cube(ix  ,iy  ,izp1)+m(i)*dex*dey*ddz*weight*v(i,type)
              cube(ix  ,iyp1,izp1)=cube(ix  ,iyp1,izp1)+m(i)*dex*ddy*ddz*weight*v(i,type)
              cube(ixp1,iy  ,izp1)=cube(ixp1,iy  ,izp1)+m(i)*ddx*dey*ddz*weight*v(i,type)
              cube(ixp1,iyp1,izp1)=cube(ixp1,iyp1,izp1)+m(i)*ddx*ddy*ddz*weight*v(i,type)
              ptot=ptot+m(i)*v(i,type)

              endif
           end if
        end do
     endif
     deallocate(x,v,m)
     if(nstar>0)deallocate(age)
     deallocate(family,tag)
  end do
  write(*,*)'Total momentum=',ptot

  ! Output file
  nomfich=TRIM(outfich)
  write(*,*)'Ecriture des donnees du fichier '//TRIM(nomfich)
  open(unit=10,file=nomfich,form='unformatted')
  if(periodic)then
     write(10)nx,ny,nz
     allocate(toto(nx,ny,nz))
     toto=cube(0:nx-1,0:ny-1,0:nz-1)
     write(10)toto
  else
     write(10)nx+1,ny+1,nz+1,ndim
     allocate(toto(nx+1,ny+1,nz+1))
     toto=cube(0:nx,0:ny,0:nz)
     write(10)toto
  endif
  close(10)

contains

  subroutine read_params

    implicit none

    integer       :: i,n

    character(len=4)   :: opt
    character(len=128) :: arg

    n = command_argument_count()
    if (n < 4) then
       print *, 'usage: part2cube  -inp  input_dir'
       print *, '                  -out  output_file'
       print *, '                 [-xmi xmin] '
       print *, '                 [-xma xmax] '
       print *, '                 [-ymi ymin] '
       print *, '                 [-yma ymax] '
       print *, '                 [-zmi zmin] '
       print *, '                 [-zma zmax] '
       print *, '                 [-nx  nx  ] '
       print *, '                 [-ny  ny  ] '
       print *, '                 [-nz  nz  ] '
       print *, '                 [-per flag] '
       print *, '                 [-str flag] '
       print *, '                 [-typ type] '
       print *, 'ex: part2cube -inp output_00001 -out cube.dat'// &
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
       case ('-age')
          read (arg,*) ageweight
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
       case ('-ny')
          read (arg,*) ny
       case ('-nz')
          read (arg,*) nz
       case ('-per')
          read (arg,*) periodic
       case ('-str')
          read (arg,*) star
       case ('-typ')
          read (arg,*) type ! specifies which momentum component you'd like.
       case default
          print '("unknown option ",a2," ignored")', opt
       end select
    end do

    return

  end subroutine read_params

end program part2cube
