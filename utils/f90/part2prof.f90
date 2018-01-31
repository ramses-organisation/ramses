program part2prof
  !--------------------------------------------------------------------------
  ! This software computes profiles.
  ! Version F90 by R. Teyssier le 01/04/01.
  ! Update for new tracers by C. Cadiou, M. Trebitsch, H. Choi, O. Snaith and R. Bieri
  !--------------------------------------------------------------------------
  use utils
  implicit none
  integer::ncpu,ndim,npart,ngrid,n,i,j,k,icpu,ipos,nstar
  integer::ncpu2,npart2,ndim2,levelmin,levelmax,ilevel,iii
  integer::nx=0,ny=0,ix,iy,ixp1,iyp1,idim,jdim,ncpu_read,n_frw
  integer::nprof=28,irad,ivar,nrad=100
  integer(kind=8)::nread
  real(KIND=8)::mtot,ddx,ddy,dex,dey,time,time_tot,time_simu,weight
  real(KIND=8)::xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1,rmax=0.5
  real(KIND=8)::xcen=0.5,ycen=0.5,zcen=0.5
  real(KIND=8)::ucen=0.0,vcen=0.0,wcen=0.0
  real(KIND=8)::xx,yy,zz,uu,vv,ww
  integer::imin,imax,jmin,jmax,kmin,kmax,lmin,npart_actual
  real(KIND=8)::xxmin,xxmax,yymin,yymax,dx,dy,deltax,boxlen
  real(KIND=8)::aexp,t,omega_m,omega_l,omega_b,omega_k,h0
  real(KIND=8)::unit_l,unit_t,unit_d,unit_m,unit_v
  real(KIND=8)::rad2,vol,rprev
  real(KIND=8)::mcumstar,ucumstar,vcumstar,wcumstar,lxcumstar,lycumstar,lzcumstar
  real(KIND=8)::mcumcdm,ucumcdm,vcumcdm,wcumcdm,lxcumcdm,lycumcdm,lzcumcdm

  real(KIND=4),dimension(:,:),allocatable::toto
  real(KIND=8),dimension(:),allocatable::aexp_frw,hexp_frw,tau_frw,t_frw
  real(KIND=8),dimension(:,:),allocatable::x,v,prof
  real(KIND=8),dimension(:),allocatable::m,age,r
  integer,dimension(:),allocatable::id
  character(LEN=1)::proj='z'
  character(LEN=5)::nchar,ncharcpu
  character(LEN=80)::ordering,format_grille
  character(LEN=128)::nomfich,repository,outfich,filedens,filetype='bin'
  logical::ok,ok_part,periodic=.false.
  integer::impi,ndom,bit_length,maxdom
  integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max
  real(KIND=8),dimension(1:8)::bounding,bounding_min,bounding_max
  real(KIND=8)::dkey,order_min,dmax
  real(kind=8),dimension(:),allocatable::bound_key
  logical,dimension(:),allocatable::cpu_read
  integer,dimension(:),allocatable::cpu_list
  logical::cosmo=.true.
  integer::idcdm=1,iucdm=2,ivcdm=3,iwcdm=4,ilxcdm=5,ilycdm=6,ilzcdm=7
  integer::imcumcdm=8,iucumcdm=9,ivcumcdm=10,iwcumcdm=11,ilxcumcdm=12,ilycumcdm=13,ilzcumcdm=14
  integer::idstar=15,iustar=16,ivstar=17,iwstar=18,ilxstar=19,ilystar=20,ilzstar=21
  integer::imcumstar=22,iucumstar=23,ivcumstar=24,iwcumstar=25,ilxcumstar=26,ilycumstar=27,ilzcumstar=28
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
  if(h0.eq.1.0)cosmo=.false.
  read(10,'("omega_m     =",E23.15)')omega_m
  read(10,'("omega_l     =",E23.15)')omega_l
  read(10,'("omega_k     =",E23.15)')omega_k
  read(10,'("omega_b     =",E23.15)')omega_b
  read(10,'("unit_l      =",E23.15)')unit_l
  read(10,'("unit_d      =",E23.15)')unit_d
  read(10,'("unit_t      =",E23.15)')unit_t
  unit_m=unit_d*unit_l**3
  unit_v=unit_l/unit_t
  read(10,*)

  read(10,'("ordering type=",A80)')ordering
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

  if(cosmo)then
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
  endif

  !-----------------------
  ! Profile parameters
  !-----------------------
  xmin=MAX(xcen-rmax,0.0d0)
  xmax=MIN(xcen+rmax,1.0d0)
  ymin=MAX(ycen-rmax,0.0d0)
  ymax=MIN(ycen+rmax,1.0d0)
  zmin=MAX(zcen-rmax,0.0d0)
  zmax=MIN(zcen+rmax,1.0d0)

  write(*,*)'time=',t
  write(*,*)'Working array =',nrad
  allocate(r(1:nrad))
  do i=1,nrad
     r(i)=dble(i)*rmax/dble(nrad)
  end do
  allocate(prof(1:nrad,1:nprof))
  prof=0.0d0

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
     read(1)
     read(1)
     read(1)
     read(1)
     read(1)
     allocate(m(1:npart2))
     allocate(age(1:npart2))
     allocate(id(1:npart2))
     allocate(family(1:npart2))
     allocate(tag(1:npart2))
     age=0d0
     id=1
     allocate(x(1:npart2,1:ndim2))
     allocate(v(1:npart2,1:ndim2))
     ! Read position
     do i=1,ndim
        read(1)m
        x(1:npart2,i)=m/boxlen
     end do
     ! Skip velocity
     do i=1,ndim
        read(1)m
        v(1:npart2,i)=m
     end do
     ! Read mass
     read(1)m
     if(nstar>0)then
        read(1)id ! Read identity
        read(1)family
        read(1)tag
        read(1)   ! Skip level
        read(1)age
     endif
     close(1)

     do i=1,npart2
        rad2=(x(i,1)-xcen)**2+(x(i,2)-ycen)**2+(x(i,3)-zcen)**2
        ok_part=(rad2<rmax**2)
        if(ok_part)then
           irad=int(dble(nrad)*sqrt(rad2)/rmax)+1
           xx=x(i,1)-xcen
           yy=x(i,2)-ycen
           zz=x(i,3)-zcen
           uu=v(i,1)-ucen/(unit_v/1d5)
           vv=v(i,2)-vcen/(unit_v/1d5)
           ww=v(i,3)-wcen/(unit_v/1d5)
           if(family(i).eq.2)then
              if(cosmo)then
                 iii=1
                 do while(tau_frw(iii)>age(i).and.iii<n_frw)
                    iii=iii+1
                 end do
                 ! Interpolate time
                 time=t_frw(iii)*(age(i)-tau_frw(iii-1))/(tau_frw(iii)-tau_frw(iii-1))+ &
                      & t_frw(iii-1)*(age(i)-tau_frw(iii))/(tau_frw(iii-1)-tau_frw(iii))
                 time=(time_simu-time)/(h0*1d5/3.08d24)/(365.*24.*3600.*1d9)
              else
                 time=(time_simu-age(i))*unit_t
              end if
              prof(irad,idstar)=prof(irad,idstar)+m(i)
              prof(irad,iustar)=prof(irad,iustar)+m(i)*v(i,1)
              prof(irad,ivstar)=prof(irad,ivstar)+m(i)*v(i,2)
              prof(irad,iwstar)=prof(irad,iwstar)+m(i)*v(i,3)
              prof(irad,ilxstar)=prof(irad,ilxstar)+m(i)*(yy*ww-zz*vv)
              prof(irad,ilystar)=prof(irad,ilystar)-m(i)*(xx*ww-zz*uu)
              prof(irad,ilzstar)=prof(irad,ilzstar)+m(i)*(xx*vv-yy*uu)
           else if(id(i)>0) then
              prof(irad,idcdm)=prof(irad,idcdm)+m(i)
              prof(irad,iucdm)=prof(irad,iucdm)+m(i)*v(i,1)
              prof(irad,ivcdm)=prof(irad,ivcdm)+m(i)*v(i,2)
              prof(irad,iwcdm)=prof(irad,iwcdm)+m(i)*v(i,3)
              prof(irad,ilxcdm)=prof(irad,ilxcdm)+m(i)*(yy*ww-zz*vv)
              prof(irad,ilycdm)=prof(irad,ilycdm)-m(i)*(xx*ww-zz*uu)
              prof(irad,ilzcdm)=prof(irad,ilzcdm)+m(i)*(xx*vv-yy*uu)
           end if
        end if
     end do
     deallocate(x,m,v,age,id)
  end do

  ! Compute cumulated profiles
  mcumstar=0d0; mcumcdm=0d0
  ucumstar=0d0; vcumstar=0d0; wcumstar=0d0
  ucumcdm=0d0; vcumcdm=0d0; wcumcdm=0d0
  lxcumstar=0d0; lycumstar=0d0; lzcumstar=0d0
  lxcumcdm=0d0; lycumcdm=0d0; lzcumcdm=0d0
  do irad=1,nrad
     mcumstar=mcumstar+prof(irad,idstar)
     ucumstar=ucumstar+prof(irad,iustar)
     vcumstar=vcumstar+prof(irad,ivstar)
     wcumstar=wcumstar+prof(irad,iwstar)
     lxcumstar=lxcumstar+prof(irad,ilxstar)
     lycumstar=lycumstar+prof(irad,ilystar)
     lzcumstar=lzcumstar+prof(irad,ilzstar)

     mcumcdm=mcumcdm+prof(irad,idcdm)
     ucumcdm=ucumcdm+prof(irad,iucdm)
     vcumcdm=vcumcdm+prof(irad,ivcdm)
     wcumcdm=wcumcdm+prof(irad,iwcdm)
     lxcumcdm=lxcumcdm+prof(irad,ilxcdm)
     lycumcdm=lycumcdm+prof(irad,ilycdm)
     lzcumcdm=lzcumcdm+prof(irad,ilzcdm)

     prof(irad,imcumstar)=mcumstar
     prof(irad,iucumstar)=ucumstar
     prof(irad,ivcumstar)=vcumstar
     prof(irad,iwcumstar)=wcumstar
     prof(irad,ilxcumstar)=lxcumstar
     prof(irad,ilycumstar)=lycumstar
     prof(irad,ilzcumstar)=lzcumstar

     prof(irad,imcumcdm)=mcumcdm
     prof(irad,iucumcdm)=ucumcdm
     prof(irad,ivcumcdm)=vcumcdm
     prof(irad,iwcumcdm)=wcumcdm
     prof(irad,ilxcumcdm)=lxcumcdm
     prof(irad,ilycumcdm)=lycumcdm
     prof(irad,ilzcumcdm)=lzcumcdm

  end do

  ! Convert profiles into proper astro units
  rprev=0d0
  do irad=1,nrad
     r(irad)=r(irad)*boxlen
     vol=4./3.*3.1415926*(r(irad)**3-rprev**3)
     ! Stars
     if(prof(irad,idstar)>0.0)then
        prof(irad,iustar)=prof(irad,iustar)/prof(irad,idstar)*unit_v/1d5
        prof(irad,ivstar)=prof(irad,ivstar)/prof(irad,idstar)*unit_v/1d5
        prof(irad,iwstar)=prof(irad,iwstar)/prof(irad,idstar)*unit_v/1d5
        prof(irad,ilxstar)=prof(irad,ilxstar)/prof(irad,idstar)*unit_v/1d5/r(irad)
        prof(irad,ilystar)=prof(irad,ilystar)/prof(irad,idstar)*unit_v/1d5/r(irad)
        prof(irad,ilzstar)=prof(irad,ilzstar)/prof(irad,idstar)*unit_v/1d5/r(irad)
     endif
     prof(irad,idstar)=prof(irad,idstar)/vol*unit_d/1.66d-24
     if(prof(irad,imcumstar)>0.0)then
        prof(irad,iucumstar)=prof(irad,iucumstar)/prof(irad,imcumstar)*unit_v/1d5
        prof(irad,ivcumstar)=prof(irad,ivcumstar)/prof(irad,imcumstar)*unit_v/1d5
        prof(irad,iwcumstar)=prof(irad,iwcumstar)/prof(irad,imcumstar)*unit_v/1d5
        prof(irad,ilxcumstar)=prof(irad,ilxcumstar)/prof(irad,imcumstar)*unit_v/1d5/r(irad)
        prof(irad,ilycumstar)=prof(irad,ilycumstar)/prof(irad,imcumstar)*unit_v/1d5/r(irad)
        prof(irad,ilzcumstar)=prof(irad,ilzcumstar)/prof(irad,imcumstar)*unit_v/1d5/r(irad)
     endif
     prof(irad,imcumstar)=sqrt(6.67e-8*prof(irad,imcumstar)*unit_m/r(irad)/unit_l)/1d5
     ! Dark matter
     if(prof(irad,idcdm)>0.0)then
        prof(irad,iucdm)=prof(irad,iucdm)/prof(irad,idcdm)*unit_v/1d5
        prof(irad,ivcdm)=prof(irad,ivcdm)/prof(irad,idcdm)*unit_v/1d5
        prof(irad,iwcdm)=prof(irad,iwcdm)/prof(irad,idcdm)*unit_v/1d5
        prof(irad,ilxcdm)=prof(irad,ilxcdm)/prof(irad,idcdm)*unit_v/1d5/r(irad)
        prof(irad,ilycdm)=prof(irad,ilycdm)/prof(irad,idcdm)*unit_v/1d5/r(irad)
        prof(irad,ilzcdm)=prof(irad,ilzcdm)/prof(irad,idcdm)*unit_v/1d5/r(irad)
     endif
     prof(irad,idcdm)=prof(irad,idcdm)/vol*unit_d/1.66d-24
     if(prof(irad,imcumcdm)>0.0)then
        prof(irad,iucumcdm)=prof(irad,iucumcdm)/prof(irad,imcumcdm)*unit_v/1d5
        prof(irad,ivcumcdm)=prof(irad,ivcumcdm)/prof(irad,imcumcdm)*unit_v/1d5
        prof(irad,iwcumcdm)=prof(irad,iwcumcdm)/prof(irad,imcumcdm)*unit_v/1d5
        prof(irad,ilxcumcdm)=prof(irad,ilxcumcdm)/prof(irad,imcumcdm)*unit_v/1d5/r(irad)
        prof(irad,ilycumcdm)=prof(irad,ilycumcdm)/prof(irad,imcumcdm)*unit_v/1d5/r(irad)
        prof(irad,ilzcumcdm)=prof(irad,ilzcumcdm)/prof(irad,imcumcdm)*unit_v/1d5/r(irad)
     endif
     prof(irad,imcumcdm)=sqrt(6.67e-8*prof(irad,imcumcdm)*unit_m/r(irad)/unit_l)/1d5
     rprev=r(irad)
  end do

  ! Output file
  nomfich=TRIM(outfich)
  write(*,*)'Ecriture des donnees du fichier '//TRIM(nomfich)
  open(unit=10,file=TRIM(nomfich)//".dark",form='formatted')
  write(10,'(A130)')" r(kpc)      n_d(H/cc)   vc_d(H/cc)  mc_d(Msol)  cu_d(km/s)  cv_d(km/s)  cw_d(km/s)  cl_d(km/s)  lx_d        ly_d        lz_d      "
  do irad=1,nrad
     write(10,999)r(irad)*unit_l/3.08d21,prof(irad,idcdm),prof(irad,imcumcdm),prof(irad,imcumcdm)**2*r(irad)*unit_l/6.67e-8*1d10/2d33 &
          & ,prof(irad,iucumcdm),prof(irad,ivcumcdm),prof(irad,iwcumcdm),sqrt(prof(irad,ilxcumcdm)**2+prof(irad,ilycumcdm)**2+prof(irad,ilzcumcdm)**2) &
          & ,prof(irad,ilxcumcdm)/sqrt(prof(irad,ilxcumcdm)**2+prof(irad,ilycumcdm)**2+prof(irad,ilzcumcdm)**2+1d-30) &
          & ,prof(irad,ilycumcdm)/sqrt(prof(irad,ilxcumcdm)**2+prof(irad,ilycumcdm)**2+prof(irad,ilzcumcdm)**2+1d-30) &
          & ,prof(irad,ilzcumcdm)/sqrt(prof(irad,ilxcumcdm)**2+prof(irad,ilycumcdm)**2+prof(irad,ilzcumcdm)**2+1d-30)
  end do
  close(10)
  open(unit=10,file=TRIM(nomfich)//".star",form='formatted')
  write(10,'(A130)')" r(kpc)      n_*(H/cc)   vc_*(H/cc)  mc_*(Msol)  cu_*(km/s)  cv_*(km/s)  cw_*(km/s)  cl_*(km/s)  lx_*        ly_*        lz_*      "
  do irad=1,nrad
     write(10,999)r(irad)*unit_l/3.08d21,prof(irad,idstar),prof(irad,imcumstar),prof(irad,imcumstar)**2*r(irad)*unit_l/6.67e-8*1d10/2d33 &
          & ,prof(irad,iucumstar),prof(irad,ivcumstar),prof(irad,iwcumstar),sqrt(prof(irad,ilxcumstar)**2+prof(irad,ilycumstar)**2+prof(irad,ilzcumstar)**2) &
          & ,prof(irad,ilxcumstar)/sqrt(prof(irad,ilxcumstar)**2+prof(irad,ilycumstar)**2+prof(irad,ilzcumstar)**2+1e-30) &
          & ,prof(irad,ilycumstar)/sqrt(prof(irad,ilxcumstar)**2+prof(irad,ilycumstar)**2+prof(irad,ilzcumstar)**2+1e-30) &
          & ,prof(irad,ilzcumstar)/sqrt(prof(irad,ilxcumstar)**2+prof(irad,ilycumstar)**2+prof(irad,ilzcumstar)**2+1e-30)
  end do
  close(10)
999 format(30(1PE10.3,2X))

contains

  subroutine read_params

    implicit none

    integer       :: i,n
    
    character(len=4)   :: opt
    character(len=128) :: arg
    LOGICAL       :: bad, ok

    n = command_argument_count()
    if (n < 4) then
       print *, 'usage: part2prof -inp  input_dir'
       print *, '                 -out  output_file'
       print *, '                 [-xce xcen] '
       print *, '                 [-yce ycen] '
       print *, '                 [-zce zcen] '
       print *, '                 [-uce ucen] '
       print *, '                 [-vce vcen] '
       print *, '                 [-wce wcen] '
       print *, '                 [-rma rmax] '
       print *, '                 [-nra nrad] '
       print *, '                 [-per flag] '
       print *, 'ex: part2prof -inp output_00001 -out map.dat'// &
            &   ' -xce 0.1 -yce 0.2 -zce 0.2 -rma 0.1 -nra 100'
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
       case ('-xce')
          read (arg,*) xcen
       case ('-yce')
          read (arg,*) ycen
       case ('-zce')
          read (arg,*) zcen
       case ('-uce')
          read (arg,*) ucen
       case ('-vce')
          read (arg,*) vcen
       case ('-wce')
          read (arg,*) wcen
       case ('-nra')
          read (arg,*) nrad
       case ('-rma')
          read (arg,*) rmax
       case ('-per')
          read (arg,*) periodic
       case default
          print '("unknown option ",a2," ignored")', opt
       end select
    end do

    return

  end subroutine read_params

end program part2prof
