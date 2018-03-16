program amr2map
  !--------------------------------------------------------------------------
  ! Ce programme calcule la carte de densite projetee pour les
  ! variables hydro d'une simulation RAMSES.
  ! Version F90 par R. Teyssier le 01/04/01.
  !--------------------------------------------------------------------------
  use utils
  implicit none
  integer::ndim,n,i,j,k,twotondim,ncoarse,type=0,action=0
  integer::ivar,nvar,ncpu,ncpuh,lmax=0,levelmin
  integer::nx,ny,nz,ilevel,iidim,idim,jdim,kdim,icell
  integer::nlevelmax,ilevel1,ngrid1
  integer::nlevelmaxs,nlevel,iout
  integer::ind,ipos,ngrida,ngridh,ilevela,ilevelh
  integer::ngridmax,nstep_coarse,icpu,ncpu_read
  integer::nhx,nhy,ihx,ihy,ivar1,ivar2
  real::gamma,smallr,smallc,gammah
  real::boxlen,boxlen2
  real(kind=8)::t,aexp
  real::hexp,t2,aexp2,hexp2
  real::omega_m,omega_l,omega_k,omega_b
  real::omega_m2,omega_l2,omega_k2,omega_b2
  real::scale_l,scale_d,scale_t
  real(kind=8)::metmax=0d0

  integer::nx_sample=0,ny_sample=0
  integer::ngrid,imin,imax,jmin,jmax,kmin,kmax
  integer::ncpu2,npart2,ndim2,nlevelmax2,nstep_coarse2
  integer::nx2,ny2,nz2,ngridmax2,nvarh,ndimh,nlevelmaxh
  integer::nx_full,ny_full,lmin,nboundary,ngrid_current
  integer::ix,iy,iz,ndom,impi,bit_length,maxdom
  integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max
  real(KIND=8),dimension(1:8)::bounding,bounding_min,bounding_max
  real(KIND=8)::dkey,order_min,dmax,ddx,dxline,ddy,dex,dey,weight
  real(KIND=8)::xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1
  real(KIND=8)::xxmin,xxmax,yymin,yymax,zzmin,zzmax,dx,dy,xx,yy
  real(KIND=8),dimension(:,:),allocatable::x,xg
  real(KIND=8),dimension(:,:,:),allocatable::var
  real(KIND=4),dimension(:,:),allocatable::toto
  real(KIND=8),dimension(:)  ,allocatable::rho,map
  logical,dimension(:)  ,allocatable::ref
  integer,dimension(:)  ,allocatable::isp
  integer,dimension(:,:),allocatable::son,ngridfile,ngridlevel,ngridbound
  real(KIND=8),dimension(1:8,1:3)::xc
  real(KIND=8),dimension(1:3)::xbound=(/0d0,0d0,0d0/)
  character(LEN=5)::nchar,ncharcpu
  character(LEN=80)::ordering
  character(LEN=80)::GMGM
  character(LEN=128)::nomfich,repository,outfich,filetype='bin'
  logical::ok,ok_part,ok_cell,do_max
  real(kind=8),dimension(:),allocatable::bound_key
  logical,dimension(:),allocatable::cpu_read
  integer,dimension(:),allocatable::cpu_list
  character(LEN=1)::proj='z'

  type level
     integer::ilevel
     integer::ngrid
     real(KIND=8),dimension(:,:),pointer::map
     real(KIND=8),dimension(:,:),pointer::rho
     integer::imin
     integer::imax
     integer::jmin
     integer::jmax
  end type level

  type(level),dimension(1:100)::grid

  call read_params

  !-----------------------------------------------
  ! Lecture du fichier hydro au format RAMSES
  !-----------------------------------------------
  ipos=INDEX(repository,'output_')
  nchar=repository(ipos+7:ipos+13)
  nomfich=TRIM(repository)//'/hydro_'//TRIM(nchar)//'.out00001'
  inquire(file=nomfich, exist=ok) ! verify input file
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     stop
  endif
  nomfich=TRIM(repository)//'/amr_'//TRIM(nchar)//'.out00001'
  inquire(file=nomfich, exist=ok) ! verify input file
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     stop
  endif

  nomfich=TRIM(repository)//'/amr_'//TRIM(nchar)//'.out00001'
  open(unit=10,file=nomfich,status='old',form='unformatted')
  read(10)ncpu
  read(10)ndim
  read(10)nx,ny,nz
  read(10)nlevelmax
  read(10)ngridmax
  read(10)nboundary
  read(10)ngrid_current
  read(10)boxlen
  close(10)
  twotondim=2**ndim
  xbound=(/dble(nx/2),dble(ny/2),dble(nz/2)/)

  allocate(ngridfile(1:ncpu+nboundary,1:nlevelmax))
  allocate(ngridlevel(1:ncpu,1:nlevelmax))
  if(nboundary>0)allocate(ngridbound(1:nboundary,1:nlevelmax))

  nomfich=TRIM(repository)//'/info_'//TRIM(nchar)//'.txt'
  inquire(file=nomfich, exist=ok) ! verify input file
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     stop
  endif
  open(unit=10,file=nomfich,form='formatted',status='old')
  read(10,*)
  read(10,*)
  read(10,'(A13,I11)')GMGM,levelmin
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)

  read(10,*)
  read(10,'(A13,E23.15)')GMGM,t
  read(10,'(A13,E23.15)')GMGM,aexp
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,'(A13,E23.15)')GMGM,scale_l
  read(10,'(A13,E23.15)')GMGM,scale_d
  read(10,'(A13,E23.15)')GMGM,scale_t
  read(10,*)

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
  ! Map parameters
  !-----------------------
  if(lmax==0)then
     lmax=nlevelmax
  endif
  write(*,*)'time=',t
  write(*,*)'Working resolution =',2**lmax
  zzmax=1.0
  zzmin=0.0
  if(ndim>2)then
  if (proj=='x')then
     idim=2
     jdim=3
     kdim=1
     xxmin=ymin ; xxmax=ymax
     yymin=zmin ; yymax=zmax
     zzmin=xmin ; zzmax=xmax
  else if (proj=='y') then
     idim=1
     jdim=3
     kdim=2
     xxmin=xmin ; xxmax=xmax
     yymin=zmin ; yymax=zmax
     zzmin=ymin ; zzmax=ymax
  else
     idim=1
     jdim=2
     kdim=3
     xxmin=xmin ; xxmax=xmax
     yymin=ymin ; yymax=ymax
     zzmin=zmin ; zzmax=zmax
  end if
  else
     idim=1
     jdim=2
     xxmin=xmin ; xxmax=xmax
     yymin=ymin ; yymax=ymax
  end if

  if(TRIM(ordering).eq.'hilbert')then

     dmax=max(xmax-xmin,ymax-ymin,zmax-zmin)
     do ilevel=1,lmax
        dx=0.5d0**ilevel
        if(dx.lt.dmax)exit
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

     dkey=(dble(2**(nlevelmax+1)/dble(maxdom)))**ndim
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

  !-----------------------------
  ! Compute hierarchy
  !-----------------------------
  do ilevel=1,lmax
     nx_full=2**ilevel
     ny_full=2**ilevel
     imin=int(xxmin*dble(nx_full))+1
     imax=int(xxmax*dble(nx_full))+1
     jmin=int(yymin*dble(ny_full))+1
     jmax=int(yymax*dble(ny_full))+1
     allocate(grid(ilevel)%map(imin:imax,jmin:jmax))
     allocate(grid(ilevel)%rho(imin:imax,jmin:jmax))
     grid(ilevel)%map(:,:)=0.0
     grid(ilevel)%rho(:,:)=0.0
     grid(ilevel)%imin=imin
     grid(ilevel)%imax=imax
     grid(ilevel)%jmin=jmin
     grid(ilevel)%jmax=jmax
  end do

  !-----------------------------------------------
  ! Compute projected variables
  !----------------------------------------------

  ! Loop over processor files
  do k=1,ncpu_read
     icpu=cpu_list(k)
     call title(icpu,ncharcpu)

     ! Open AMR file and skip header
     nomfich=TRIM(repository)//'/amr_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
     open(unit=10,file=nomfich,status='old',form='unformatted')
     write(*,*)'Processing file '//TRIM(nomfich)
     do i=1,21
        read(10)
     end do
     ! Read grid numbers
     read(10)ngridlevel
     ngridfile(1:ncpu,1:nlevelmax)=ngridlevel
     read(10)
     if(nboundary>0)then
        do i=1,2
           read(10)
        end do
        read(10)ngridbound
        ngridfile(ncpu+1:ncpu+nboundary,1:nlevelmax)=ngridbound
     endif
     read(10)
! ROM: comment the single follwing line for old stuff
     read(10)
     if(TRIM(ordering).eq.'bisection')then
        do i=1,5
           read(10)
        end do
     else
        read(10)
     endif
     read(10)
     read(10)
     read(10)

     ! Open HYDRO file and skip header
     nomfich=TRIM(repository)//'/hydro_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
     open(unit=11,file=nomfich,status='old',form='unformatted')
     read(11)
     read(11)nvarh
     read(11)
     read(11)
     read(11)
     read(11)

     ! Loop over levels
     do ilevel=1,lmax

        ! Geometry
        dx=0.5**ilevel
        dxline=1
        if(ndim==3)dxline=dx
        nx_full=2**ilevel
        ny_full=2**ilevel
        do ind=1,twotondim
           iz=(ind-1)/4
           iy=(ind-1-4*iz)/2
           ix=(ind-1-2*iy-4*iz)
           xc(ind,1)=(dble(ix)-0.5D0)*dx
           xc(ind,2)=(dble(iy)-0.5D0)*dx
           xc(ind,3)=(dble(iz)-0.5D0)*dx
        end do

        ! Allocate work arrays
        ngrida=ngridfile(icpu,ilevel)
        grid(ilevel)%ngrid=ngrida
        if(ngrida>0)then
           allocate(xg(1:ngrida,1:ndim))
           allocate(son(1:ngrida,1:twotondim))
           allocate(var(1:ngrida,1:twotondim,1:nvarh))
           allocate(x  (1:ngrida,1:ndim))
           allocate(rho(1:ngrida))
           allocate(map(1:ngrida))
           allocate(ref(1:ngrida))
        endif

        ! Loop over domains
        do j=1,nboundary+ncpu

           ! Read AMR data
           if(ngridfile(j,ilevel)>0)then
              read(10) ! Skip grid index
              read(10) ! Skip next index
              read(10) ! Skip prev index
              ! Read grid center
              do iidim=1,ndim
                 if(j.eq.icpu)then
                    read(10)xg(:,iidim)
                 else
                    read(10)
                 endif
              end do
              read(10) ! Skip father index
              do ind=1,2*ndim
                 read(10) ! Skip nbor index
              end do
              ! Read son index
              do ind=1,twotondim
                 if(j.eq.icpu)then
                    read(10)son(:,ind)
                 else
                    read(10)
                 end if
              end do
              ! Skip cpu map
              do ind=1,twotondim
                 read(10)
              end do
              ! Skip refinement map
              do ind=1,twotondim
                 read(10)
              end do
           endif

           ! Read HYDRO data
           read(11)
           read(11)
           if(ngridfile(j,ilevel)>0)then
              ! Read hydro variables
              do ind=1,twotondim
                 do ivar=1,nvarh
                    if(j.eq.icpu)then
                       read(11)var(:,ind,ivar)
                    else
                       read(11)
                    end if
                 end do
              end do
           end if
        end do

        ! Compute map
        if(ngrida>0)then

           ! Loop over cells
           do ind=1,twotondim

              ! Compute cell center
              do i=1,ngrida
                 x(i,1)=(xg(i,1)+xc(ind,1)-xbound(1))
                 x(i,2)=(xg(i,2)+xc(ind,2)-xbound(2))
                 if(ndim>2)x(i,3)=(xg(i,3)+xc(ind,3)-xbound(3))
              end do
              ! Check if cell is refined
              do i=1,ngrida
                 ref(i)=son(i,ind)>0.and.ilevel<lmax
              end do
              ! Extract variable
              rho = var(:,ind,1)
              select case (type)
              case (-1)
                 map = icpu
                 metmax=max(metmax,dble(icpu))
              case (0)
                 map = ilevel
                 metmax=max(metmax,dble(ilevel))
              case (12) !! This is for H2 using HI and HII (ramses_rt patch mol)
                 if(action==0)then
                    map = (1.0-var(:,ind,8)-var(:,ind,9))*var(:,ind,1)
                 else
                    map = 1.0-var(:,ind,8)-var(:,ind,9)
                 endif
                 metmax=max(metmax,maxval(1.0-var(:,ind,8)-var(:,ind,9)))
              case (5) !! This is for temperature
                 if(action==0)then
                    map = var(:,ind,5)
                 else
                    map = var(:,ind,5)/var(:,ind,1)
                 endif
                 metmax=max(metmax,maxval(var(:,ind,type)))
              case default ! Hydro variable
                 if(action==0)then
                    map = var(:,ind,type)*var(:,ind,1)
                 else
                    map = var(:,ind,type)
                 endif
                 metmax=max(metmax,maxval(var(:,ind,type)))
              end select
              ! Store data map
              do i=1,ngrida
                 ok_cell= .not.ref(i)
                 if(ok_cell)then
                    ix=int(x(i,idim)*dble(nx_full))+1
                    iy=int(x(i,jdim)*dble(ny_full))+1
                    if(ndim==3)then
                       weight=(min(x(i,kdim)+dx/2.,zzmax)-max(x(i,kdim)-dx/2.,zzmin))/dx
                       weight=min(1.0d0,max(weight,0.0d0))
                    else
                       weight=1.0
                    endif
                    if(    ix>=grid(ilevel)%imin.and.&
                         & iy>=grid(ilevel)%jmin.and.&
                         & ix<=grid(ilevel)%imax.and.&
                         & iy<=grid(ilevel)%jmax)then
                       if(action==2)then
                          if(weight>0.0d0)then
                             grid(ilevel)%map(ix,iy)=max(grid(ilevel)%map(ix,iy),map(i))
                             grid(ilevel)%rho(ix,iy)=max(grid(ilevel)%rho(ix,iy),rho(i))
                          endif
                       else
                          grid(ilevel)%map(ix,iy)=grid(ilevel)%map(ix,iy)+map(i)*dxline*weight/(zzmax-zzmin)
                          grid(ilevel)%rho(ix,iy)=grid(ilevel)%rho(ix,iy)+rho(i)*dxline*weight/(zzmax-zzmin)
                       endif
                    endif
                 end if
              end do

           end do
           ! End loop over cell

           deallocate(xg,son,var,ref,rho,map,x)
        end if

     end do
     ! End loop over levels

     close(10)
     close(11)

  end do
  ! End loop over cpu

  if(type==0.OR.type==1.OR.type>4)then
     write(*,*)'max val=',metmax
  endif

  nx_full=2**lmax
  ny_full=2**lmax
  imin=int(xxmin*dble(nx_full))+1
  imax=int(xxmax*dble(nx_full))
  jmin=int(yymin*dble(ny_full))+1
  jmax=int(yymax*dble(ny_full))

  do ix=imin,imax
     xmin=((ix-0.5)/2**lmax)
     do iy=jmin,jmax
        ymin=((iy-0.5)/2**lmax)
        do ilevel=1,lmax-1
           ndom=2**ilevel
           i=int(xmin*ndom)+1
           j=int(ymin*ndom)+1
           if(action==2) then
              grid(lmax)%map(ix,iy)=max(grid(lmax)%map(ix,iy),grid(ilevel)%map(i,j))
              grid(lmax)%rho(ix,iy)=max(grid(lmax)%rho(ix,iy),grid(ilevel)%rho(i,j))
           else
              grid(lmax)%map(ix,iy)=grid(lmax)%map(ix,iy) + grid(ilevel)%map(i,j)
              grid(lmax)%rho(ix,iy)=grid(lmax)%rho(ix,iy) + grid(ilevel)%rho(i,j)
           endif
        end do
     end do
  end do

  write(*,*)'Norm=',sum(grid(lmax)%rho(imin:imax,jmin:jmax))/(imax-imin+1)/(jmax-jmin+1)

  ! Output file
  nomfich=TRIM(outfich)
  write(*,*)'Ecriture des donnees du fichier '//TRIM(nomfich)

  if (filetype=='bin')then
     open(unit=20,file=nomfich,form='unformatted')
     if(nx_sample==0)then
        write(20)t, xxmax-xxmin, yymax-yymin, zzmax-zzmin
        write(20)imax-imin+1,jmax-jmin+1
        allocate(toto(imax-imin+1,jmax-jmin+1))
        if(action==0)then
           toto=grid(lmax)%map(imin:imax,jmin:jmax)/grid(lmax)%rho(imin:imax,jmin:jmax)
        else 
           toto=grid(lmax)%map(imin:imax,jmin:jmax)
        endif
        write(20)toto
     else
        if(ny_sample==0)ny_sample=nx_sample
        write(20)t, xxmax-xxmin, yymax-yymin, zzmax-zzmin
        write(20)nx_sample+1,ny_sample+1
        allocate(toto(0:nx_sample,0:ny_sample))
        do i=0,nx_sample
           ix=int(dble(i)/dble(nx_sample)*dble(imax-imin+1))+imin
           ix=min(ix,imax)
           do j=0,ny_sample
              iy=int(dble(j)/dble(ny_sample)*dble(jmax-jmin+1))+jmin
              iy=min(iy,jmax)
              if(action==0)then
                 toto(i,j)=grid(lmax)%map(ix,iy)/grid(lmax)%rho(ix,iy)
              else
                 toto(i,j)=grid(lmax)%map(ix,iy)
              endif
           end do
        end do
        write(20)toto
     endif
     close(20)
  endif
  if (filetype=='ascii')then
     open(unit=20,file=nomfich,form='formatted')
     if(nx_sample==0)then
        do j=jmin,jmax
           do i=imin,imax
              xx=xxmin+(dble(i-imin)+0.5)/dble(imax-imin+1)*(xxmax-xxmin)
              yy=yymin+(dble(j-jmin)+0.5)/dble(jmax-jmin+1)*(yymax-yymin)
              if(action==0)then
                 write(20,*)xx,yy,grid(lmax)%map(i,j)/grid(lmax)%rho(i,j)
              else
                 write(20,*)xx,yy,grid(lmax)%map(i,j)
              endif
           end do
           write(20,*) " "
        end do
     else
        if(ny_sample==0)ny_sample=nx_sample
        allocate(toto(0:nx_sample,0:ny_sample))
        do i=0,nx_sample
           ix=int(dble(i)/dble(nx_sample)*dble(imax-imin+1))+imin
           ix=min(ix,imax)
           do j=0,ny_sample
              iy=int(dble(j)/dble(ny_sample)*dble(jmax-jmin+1))+jmin
              iy=min(iy,jmax)
              if(action==0)then
                 toto(i,j)=grid(lmax)%map(ix,iy)/grid(lmax)%rho(ix,iy)
              else
                 toto(i,j)=grid(lmax)%map(ix,iy)
           endif
           end do
        end do
        do j=0,ny_sample
           do i=0,nx_sample
              xx=xxmin+dble(i)/dble(nx_sample)*(xxmax-xxmin)
              yy=yymin+dble(j)/dble(ny_sample)*(yymax-yymin)
              write(20,*)xx,yy,toto(i,j)
           end do
           write(20,*) " "
        end do
     endif
     close(20)
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
       print *, 'usage: amr2map   -inp  input_dir'
       print *, '                 -out  output_file'
       print *, '                 [-dir axis] '
       print *, '                 [-xmi xmin] '
       print *, '                 [-xma xmax] '
       print *, '                 [-ymi ymin] '
       print *, '                 [-yma ymax] '
       print *, '                 [-zmi zmin] '
       print *, '                 [-zma zmax] '
       print *, '                 [-lma lmax] '
       print *, '                 [-typ type] '
       print *, '                 [-fil filetype] '
       print *, '                 [-act action] '
       print *, 'ex: amr2map -inp output_00001 -out map.dat'// &
            &   ' -dir z -xmi 0.1 -xma 0.7 -lma 12'
       print *, ' '
       print *, ' axis : "x" = choose x-axis for direction of projection'
       print *, ' axis : "y" = choose y-axis for direction of projection'
       print *, ' axis : "z" = choose z-axis for direction of projection (default)'
       print *, ' '
       print *, ' type :-1 = cpu number'
       print *, ' type : 0 = ref. level (default)'
       print *, ' type : 1 = gas density'
       print *, ' type : 2 = X velocity'
       print *, ' type : 3 = Y velocity'
       print *, ' type : 4 = Z velocity'
       print *, ' type : 5 = gas pressure'
       print *, ' type : 6 = gas metallicity'
       print *, ' '
       print *, ' filetype : "bin" = binary output file (default)'
       print *, ' filetype : "ascii" = ascii output file'
       print *, ' '
       print *, ' action : 0 = mass-weighted average along line of sight (default)'
       print *, ' action : 1 = volume-weighted average along line of sight'
       print *, ' action : 2 = maximum along line of sight'
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
       case ('-lma')
          read (arg,*) lmax
       case ('-nx')
          read (arg,*) nx_sample
       case ('-ny')
          read (arg,*) ny_sample
       case ('-typ')
          read (arg,*) type
       case ('-fil')
          read (arg,*) filetype
       case ('-act')
          read (arg,*) action
       case default
          print '("unknown option ",a2," ignored")', opt
       end select
    end do

    return

  end subroutine read_params

end program amr2map
