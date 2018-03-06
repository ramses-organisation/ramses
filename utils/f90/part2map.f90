program part2map
  !--------------------------------------------------------------------------
  ! Ce programme calcule la carte de densite surfacique projetee
  ! des particules de matiere noire d'une simulation RAMSES.
  ! Version F90 par R. Teyssier le 01/04/01.
  !--------------------------------------------------------------------------
  use utils
  implicit none
  integer::ncpu,ndim,npart,ngrid,n,i,j,k,icpu,ipos,nstar
  integer::ncpu2,npart2,ndim2,levelmin,levelmax,ilevel,iii
  integer::nx=0,ny=0,ix,iy,ixp1,iyp1,idim,jdim,ncpu_read,n_frw
  real(KIND=8)::mtot,ddx,ddy,dex,dey,time,time_tot,time_simu,weight
  real(KIND=8)::xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1,mmax=1d10
  real(KIND=8)::xcenter,ycenter,zcenter
  real(KIND=8)::x_coord,y_coord,z_coord
  real(KIND=8)::jxin=0,jyin=0,jzin=0,jx,jy,jz
  real(KIND=8)::kxin,kyin,kzin,kx,ky,kz
  real(KIND=8)::lxin,lyin,lzin,lx,ly,lz
  integer::imin,imax,jmin,jmax,kmin,kmax,lmin,npart_actual
  real(KIND=8)::xxmin,xxmax,yymin,yymax,zzmin,zzmax,dx,dy,deltax,boxlen
  real(KIND=8)::aexp,t,omega_m,omega_l,omega_b,omega_k,h0,unit_l,unit_t,unit_d
  real(KIND=4),dimension(:,:),allocatable::toto
  real(KIND=4),dimension(:),allocatable::density
  real(KIND=8),dimension(:),allocatable::aexp_frw,hexp_frw,tau_frw,t_frw
  real(KIND=8),dimension(:,:),allocatable::map
  real(KIND=8),dimension(:,:),allocatable::x
  real(KIND=8),dimension(:)  ,allocatable::m,age
  integer,dimension(:)  ,allocatable::id
  character(LEN=1)::proj='z'
  character(LEN=5)::nchar,ncharcpu
  character(LEN=80)::GMGM
  character(LEN=80)::ordering,format_grille
  character(LEN=128)::nomfich,repository,outfich,filedens,filetype='bin'
  logical::ok,ok_part,periodic=.false.,star=.false.,ageweight=.false.,do_density=.false.
  integer::impi,ndom,bit_length,maxdom
  integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max
  real(KIND=8),dimension(1:8)::bounding,bounding_min,bounding_max
  real(KIND=8)::dkey,order_min,dmax
  real(kind=8)::xx,yy,zz
  real(kind=8),dimension(:),allocatable::bound_key
  logical,dimension(:),allocatable::cpu_read
  integer,dimension(:),allocatable::cpu_list
  logical::cosmo=.true.
  logical::rotation=.false.
  logical::sideon=.false.
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
  read(10,'(A13,I11)')GMGM,ncpu
  read(10,'(A13,I11)')GMGM,ndim
  read(10,'(A13,I11)')GMGM,levelmin
  read(10,'(A13,I11)')GMGM,levelmax
  read(10,*)
  read(10,*)
  read(10,*)

  read(10,'(A13,E23.15)')GMGM,boxlen
  read(10,'(A13,E23.15)')GMGM,t
  read(10,'(A13,E23.15)')GMGM,aexp
  read(10,'(A13,E23.15)')GMGM,h0
  read(10,'(A13,E23.15)')GMGM,omega_m
  read(10,'(A13,E23.15)')GMGM,omega_l
  read(10,'(A13,E23.15)')GMGM,omega_k
  read(10,'(A13,E23.15)')GMGM,omega_b
  read(10,'(A13,E23.15)')GMGM,unit_l
  read(10,'(A13,E23.15)')GMGM,unit_d
  read(10,'(A13,E23.15)')GMGM,unit_t

  read(10,*)

  if(aexp.eq.1.and.h0.eq.1)cosmo=.false.

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

  ! Read a density file from hop
  if(do_density)then
     open(unit=1,file=filedens,form='unformatted',status='old',access='direct',recl=1)
     read(1,rec=1)npart
     write(*,*)npart
     allocate(density(npart))
     do j=1,npart
        read(1,rec=j+1)density(j)
!        if(mod(j,10000).eq.0)write(*,*)j,density(j)
     end do
     close(1)
     write(*,*)'Min Max density'
     write(*,*)minval(density),maxval(density)
  endif

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
     write(*,*)'Age simu=',(time_tot+time_simu)/(h0*1d5/3.08d24)/(365.*24.*3600.*1d9)
  else
     time_simu=t
     write(*,*)'Age simu=',time_simu*unit_t/(365.*24.*3600.*1d9)
  endif
  !-----------------------
  ! Rotation parameters
  !-----------------------
  if(abs(jxin)+abs(jyin)+abs(jzin)>0)then
     rotation=.true.
     write(*,*)'Performing rotation'
     write(*,*)jxin,jyin,jzin
  endif

  kxin=0.
  kyin=-jzin
  kzin=jyin

  lxin=jyin*kzin-jzin*kyin
  lyin=jzin*kxin-jxin*kzin
  lzin=jxin*kyin-jyin*kxin

  jx=jxin/sqrt(jxin**2+jyin**2+jzin**2)
  jy=jyin/sqrt(jxin**2+jyin**2+jzin**2)
  jz=jzin/sqrt(jxin**2+jyin**2+jzin**2)

  kx=kxin/sqrt(kxin**2+kyin**2+kzin**2)
  ky=kyin/sqrt(kxin**2+kyin**2+kzin**2)
  kz=kzin/sqrt(kxin**2+kyin**2+kzin**2)

  lx=lxin/sqrt(lxin**2+lyin**2+lzin**2)
  ly=lyin/sqrt(lxin**2+lyin**2+lzin**2)
  lz=lzin/sqrt(lxin**2+lyin**2+lzin**2)

  xcenter=0.5*(xmin+xmax)
  ycenter=0.5*(ymin+ymax)
  zcenter=0.5*(zmin+zmax)
  !-----------------------
  ! Map parameters
  !-----------------------
  if(nx==0)then
     nx=2**levelmin
  endif
  if(ny==0)then
     ny=nx
  end if
  write(*,*)'time=',t
  write(*,*)'Working map =',nx,ny
  allocate(map(0:nx,0:ny))
  map=0.0d0

  if(rotation)then
     proj='z'
     if(sideon)then
        proj='y'
     endif
  endif

  if (proj=='x')then
     idim=2
     jdim=3
     xxmin=ymin ; xxmax=ymax
     yymin=zmin ; yymax=zmax
     zzmin=xmin ; zzmax=xmax
  else if (proj=='y') then
     idim=1
     jdim=3
     xxmin=xmin ; xxmax=xmax
     yymin=zmin ; yymax=zmax
     zzmin=ymin ; zzmax=ymax
  else
     idim=1
     jdim=2
     xxmin=xmin ; xxmax=xmax
     yymin=ymin ; yymax=ymax
     zzmin=zmin ; zzmax=zmax
  end if
  dx=(xxmax-xxmin)/dble(nx)
  dy=(yymax-yymin)/dble(ny)

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
           call hilbert3d(idom(i),jdom(i),kdom(i),bounding(1),bit_length,1)
           order_min=bounding(1)
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

  if(do_density)then
     ncpu_read=ncpu
     do j=1,ncpu
        cpu_list(j)=j
     end do
  endif

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
     if(star)then
        write(*,*)'Keeping star particles.'
     else
        write(*,*)'Discard star particles.'
     endif
  endif

  !-----------------------------------------------
  ! Compute projected mass using CIC smoothing
  !----------------------------------------------
  npart_actual=0
  mtot=0.0d0
  do k=1,ncpu_read
     icpu=cpu_list(k)
     call title(icpu,ncharcpu)
     nomfich=TRIM(repository)//'/part_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
     open(unit=1,file=nomfich,status='old',form='unformatted')
!     write(*,*)'Processing file '//TRIM(nomfich)
     read(1)ncpu2
     read(1)ndim2
     read(1)npart2
     write(*,*)k,npart2
     read(1)
     read(1)
     read(1)
     read(1)
     read(1)
     allocate(m(1:npart2))
     if(nstar>0)then
        allocate(age(1:npart2))
        allocate(id(1:npart2))
     endif
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
     if(nstar>0)then
        read(1)id
        read(1) ! Skip level
        read(1)family
        read(1)tag
        read(1)age
     endif
     close(1)

     if(do_density)then
        do i=1,npart2
           ok_part=(x(i,1)>=xmin.and.x(i,1)<=xmax.and. &
                &   x(i,2)>=ymin.and.x(i,2)<=ymax.and. &
                &   x(i,3)>=zmin.and.x(i,3)<=zmax.and. &
                &   m(i)<mmax )

           if(nstar>0)then
              ! Star
              if(star)then
                 if(family(i)==2) then
                    npart_actual=npart_actual+1
                 else
                    ok_part=.false.
                 endif
              else
                 ! Dark matter
                 if(family(i)==1)then
                    npart_actual=npart_actual+1
                 else
                    ok_part=.false.
                 endif
              endif
           else
              npart_actual=npart_actual+1
           endif

           if(ok_part)then
              ddx=(x(i,idim)-xxmin)/dx
              ddy=(x(i,jdim)-yymin)/dy
              ix=ddx
              iy=ddy
              map(ix,iy)=MAX(map(ix,iy),dble(density(npart_actual)))
           endif
        end do
     else
        do i=1,npart2
           weight=1d0
           ok_part=(x(i,1)>=xmin.and.x(i,1)<=xmax.and. &
                &   x(i,2)>=ymin.and.x(i,2)<=ymax.and. &
                &   x(i,3)>=zmin.and.x(i,3)<=zmax.and. &
                &   m(i)<mmax)

           if(nstar>0)then
              if(star)then
                 ok_part=ok_part.and.family(i)==2
                 if(ageweight)then
                    if(cosmo)then
                       iii=1
                       do while(tau_frw(iii)>age(i).and.iii<n_frw)
                          iii=iii+1
                       end do
                       ! Interploate time
                       time=t_frw(iii)*(age(i)-tau_frw(iii-1))/(tau_frw(iii)-tau_frw(iii-1))+ &
                            & t_frw(iii-1)*(age(i)-tau_frw(iii))/(tau_frw(iii-1)-tau_frw(iii))
                       time=(time_simu-time)/(h0*1d5/3.08d24)/(365.*24.*3600.*1d9)
                    else
                       time=(time_simu-age(i))*unit_t/(365.*24.*3600.*1d9)
                    endif
                    if(time>0.01)then
                       weight=(time/0.01)**(-0.7)
                    endif
                 endif
              else
                 ok_part=ok_part.and.family(i)==1
              endif
           endif

           if(ok_part)then
              npart_actual=npart_actual+1
              if(rotation)then
                 xx=x(i,1)-xcenter
                 yy=x(i,2)-ycenter
                 zz=x(i,3)-zcenter
                 z_coord=xx*jx+yy*jy+zz*jz+zcenter
                 y_coord=xx*kx+yy*ky+zz*kz+ycenter
                 x_coord=xx*lx+yy*ly+zz*lz+xcenter
                 ddx=(x_coord-xmin)/dx
                 ddy=(y_coord-ymin)/dy
                 if(sideon)then
                    ddx=(x_coord-xmin)/dx
                    ddy=(z_coord-zmin)/dy
                 endif
              else
                 ddx=(x(i,idim)-xxmin)/dx
                 ddy=(x(i,jdim)-yymin)/dy
              endif
              ix=ddx
              iy=ddy
              ddx=ddx-ix
              ddy=ddy-iy
              dex=1.0-ddx
              dey=1.0-ddy
              if(periodic)then
                 if(ix<0)ix=ix+nx
                 if(ix>=nx)ix=ix-nx
                 if(iy<0)iy=iy+ny
                 if(iy>=ny)iy=iy-ny
              endif
              ixp1=ix+1
              iyp1=iy+1
              if(periodic)then
                 if(ixp1<0)ixp1=ixp1+nx
                 if(ixp1>=nx)ixp1=ixp1-nx
                 if(iyp1<0)iyp1=iyp1+ny
                 if(iyp1>=ny)iyp1=iyp1-ny
              endif
              if(ix>=0.and.ix<nx.and.iy>=0.and.iy<ny.and.ddx>0.and.ddy>0)then
                 map(ix  ,iy  )=map(ix  ,iy  )+m(i)*dex*dey*weight
                 map(ix  ,iyp1)=map(ix  ,iyp1)+m(i)*dex*ddy*weight
                 map(ixp1,iy  )=map(ixp1,iy  )+m(i)*ddx*dey*weight
                 map(ixp1,iyp1)=map(ixp1,iyp1)+m(i)*ddx*ddy*weight
                 mtot=mtot+m(i)
              endif
           end if
        end do
     end if
     deallocate(x,m)
     if(nstar>0)deallocate(age,id)
     deallocate(family,tag)
  end do
  write(*,*)'Total mass=',mtot
  if(.not. star)write(*,*)'npart tot=',npart_actual
  write(*,*)MINVAL(map),MAXVAL(map)

  ! Output file
  nomfich=TRIM(outfich)
  write(*,*)'Writing data to '//TRIM(nomfich)
  if(do_density.or.periodic)then
     allocate(toto(0:nx-1,0:ny-1))
     toto=map(0:nx-1,0:ny-1)
  else
     allocate(toto(0:nx,0:ny))
     toto=map(0:nx,0:ny)
  endif
  ! Binary format (to be read by ramses utilities)
  if(TRIM(filetype).eq.'bin')then
     open(unit=10,file=nomfich,form='unformatted')
     if(do_density.or.periodic)then
        write(10)t, xxmax-xxmin, yymax-yymin, zzmax-zzmin
        write(10)nx,ny
        write(10)toto
        write(10)xxmin,xxmax
        write(10)yymin,yymax
     else
        write(10)t, xxmax-xxmin, yymax-yymin, zzmax-zzmin
        write(10)nx+1,ny+1
        write(10)toto
        write(10)xxmin,xxmax
        write(10)yymin,yymax
     endif
     close(10)
  endif
  ! Ascii format (to be read by gnuplot)
  if(TRIM(filetype).eq.'ascii')then
     open(unit=10,file=nomfich,form='formatted')
     if(do_density.or.periodic)then
        do j=0,ny-1
           do i=0,nx-1
              xx=xxmin+(dble(i)+0.5)/dble(nx)*(xxmax-xxmin)
              yy=yymin+(dble(j)+0.5)/dble(ny)*(yymax-yymin)
              write(10,*)xx,yy,toto(i,j)
           end do
           write(10,*) " "
        end do
     else
        do j=0,ny
           do i=0,nx
              xx=xxmin+dble(i)/dble(nx)*(xxmax-xxmin)
              yy=yymin+dble(j)/dble(ny)*(yymax-yymin)
              write(10,*)xx,yy,toto(i,j)
           end do
           write(10,*) " "
        end do
     endif
     close(10)
  endif

contains

  subroutine read_params

      implicit none

      integer       :: i,n
      
      character(len=4)   :: opt
      character(len=128) :: arg
      LOGICAL       :: bad, ok

    n = command_argument_count()
      if (n < 4) then
         print *, 'usage: part2map  -inp  input_dir'
         print *, '                 -out  output_file'
         print *, '                 [-dir axis] '
         print *, '                 [-xmi xmin] '
         print *, '                 [-xma xmax] '
         print *, '                 [-ymi ymin] '
         print *, '                 [-yma ymax] '
         print *, '                 [-zmi zmin] '
         print *, '                 [-zma zmax] '
         print *, '                 [-jx  jxin] '
         print *, '                 [-jy  jyin] '
         print *, '                 [-jz  jzin] '
         print *, '                 [-cut mmax] '
         print *, '                 [-nx  nx  ] '
         print *, '                 [-ny  ny  ] '
         print *, '                 [-per flag] '
         print *, '                 [-sid flag] '
         print *, '                 [-str flag] '
         print *, '                 [-fil filetype: bin/ascii] '
         print *, '                 [-den filedens] '
         print *, 'ex: part2map -inp output_00001 -out map.dat'// &
              &   ' -dir z -xmi 0.1 -xma 0.7'
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
         case ('-dir')
            proj = trim(arg)
         case('-jx')
            read (arg,*) jxin
         case('-jy')
            read (arg,*) jyin
         case('-jz')
            read (arg,*) jzin
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
         case ('-cut')
            read (arg,*) mmax
         case ('-zma')
            read (arg,*) zmax
         case ('-nx')
            read (arg,*) nx
         case ('-ny')
            read (arg,*) ny
         case ('-per')
            read (arg,*) periodic
         case ('-sid')
            read (arg,*) sideon
         case ('-str')
            read (arg,*) star
         case ('-age')
            read (arg,*) ageweight
         case ('-fil')
            filetype = trim(arg)
         case ('-den')
            filedens = trim(arg)
            do_density=.true.
         case default
            print '("unknown option ",a2," ignored")', opt
         end select
      end do

      return

    end subroutine read_params

  end program part2map
