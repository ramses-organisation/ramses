!===============================================================================
!
!  This module adds a random forcing module to RAMSES
!
!  It calculate a random force, confined to a shell in k-space, at
!  regular intervals in time.  We assume rho, P, and the sound
!  speed to be near unity.  The velocity amplitude should be of the
!  order ampl_turb.  The size (scale) of the driving motions is of the
!  order of mx*dx/k1.  Thus the turnover time is given by
!  t_turb*ampl_turb = mx*dx/k1; and the acceleration is of the order
!  ampl_turb/t_turb.
!
!  Note:  In order to restart properly, one would have to save the
!  previous value of iseed. Instead we trust the user to specify the same value
!  in the input file on consecutive runs
!
!  if you use this module in your scientific work please cite:
!
!    Paolo Padoan, Aake Nordlund, ApJ 526, 279 (1999)
!
!  @ARTICLE{1999ApJ...526..279P,
!     author = {{Padoan}, P. and {Nordlund}, {\AA}.},
!      title = "{A Super-Alfv{\'e}nic Model of Dark Clouds}",
!    journal = {\apj},
!     eprint = {astro-ph/9901288},
!   keywords = {ISM: CLOUDS, ISM: KINEMATICS AND DYNAMICS, ISM: MAGNETIC FIELDS, MAGNETOHYDRODYNAMICS: MHD, SHOCK WAVES, TURBULENCE, ISM: Clouds, ISM: Kinematics and Dynamics, ISM: Magnetic Fields, Magnetohydrodynamics: MHD, Shock Waves, Turbulence},
!       year = 1999,
!      month = nov,
!     volume = 526,
!      pages = {279-294},
!        doi = {10.1086/307956},
!     adsurl = {http://adsabs.harvard.edu/abs/1999ApJ...526..279P},
!    adsnote = {Provided by the SAO/NASA Astrophysics Data System}
!  }
!
!  Version history:
!  1997 Paolo padoan: original stagger-code version
!  2010 Troels Haugboelle: rewritten to support RAMSES, an AMR grid, and OpenMP parallelization
!
!===============================================================================
module forcing
  use amr_parameters, only : dp
  implicit none
  logical do_force, do_helmh
  integer kmax, iseed, nmodes
  real(dp) tprev, t_turb, t_turn, ampl_turb, k1, k2, pk, a_helmh
  real(dp) rhoav, pxav, pyav, pzav, aranx, arany, aranz, accel
  real(dp),    parameter  :: pi = 3.14159265358979323846264338328
  complex(dp), allocatable, dimension(:,:) :: ak_prev, ak_next
  complex(dp), allocatable, dimension(:)   :: kx, ky, kz, fxx, fyy, fzz
contains
! Extremely simple local random number generator
!===============================================================================
function ran1s(idum)
  INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
  real(dp) ran1s,AM,EPS,RNMX
  parameter (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836, &
  NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
  integer k,iy
  if (idum.le.0) then
     idum=max(-idum,1)
  end if
  k=idum/IQ
  idum=IA*(idum-k*IQ)-IR*k
  if (idum.lt.0) idum=idum+IM
  iy=idum
  ran1s=min(AM*iy,RNMX)
end function ran1s

! Init routine for forcing module
!===============================================================================
subroutine init_force
  use amr_commons,    only : myid
  use amr_parameters, only : ndim, boxlen
  implicit none
  namelist /force/ do_force, do_helmh, a_helmh, k1, k2, ampl_turb, t_turn, t_turb, pk, iseed
  integer :: jx, jy, jz

  ampl_turb = 0.1
  do_force = .false.
  do_helmh = .true.
  a_helmh=1.
  k1=1.
  k2=1.5
  pk=7./6.
  t_turn = 0.1
  t_turb = 0.
  iseed = 123456789

  rewind(1); read (1,NML=force); if (myid==1) write(*,force)                    ! namelist parameters

  if (.not. do_force) return

  if (ndim .ne. 3) then
     if (myid == 1) write(*,*) 'Forcing only works in 3 dimensions!'
     stop
  end if

  if (t_turb == 0.) t_turb = -1.5*t_turn

  kmax = ifix(real(k2,kind=4)+1.)

  ! count how many modes we have
  nmodes = 0
  do jx=-kmax,kmax
     do jy=-kmax,kmax
        do jz=-kmax,kmax
           if (float(jx**2+jy**2+jz**2) .ge. k1**2-1e-5 .and. &
               float(jx**2+jy**2+jz**2) .le. k2**2+1e-5) nmodes = nmodes + 1
        end do
     end do
  end do

  ! calculate average acceleration
  accel  = ampl_turb/t_turn/sqrt(real(nmodes,kind=dp)/8.0_dp)    ! rms=1./8. per comp.

  ! allocate mode arrays
  allocate (ak_prev(3,nmodes), ak_next(3,nmodes), &
            fxx(nmodes), fyy(nmodes), fzz(nmodes), &
             kx(nmodes),  ky(nmodes),  kz(nmodes))
  ak_prev=0.; ak_next=0.

  tprev = -1.

  ! init mode arrays
  nmodes = 0
  do jx=-kmax,kmax
     do jy=-kmax,kmax
        do jz=-kmax,kmax
           if (float(jx**2+jy**2+jz**2) .ge. k1**2-1e-5 .and. &
               float(jx**2+jy**2+jz**2) .le. k2**2+1e-5) then
              nmodes = nmodes + 1
              kx(nmodes) = cmplx (0., jx*2.*pi/boxlen)
              ky(nmodes) = cmplx (0., jy*2.*pi/boxlen)
              kz(nmodes) = cmplx (0., jz*2.*pi/boxlen)
           end if
        end do
     end do
  end do

end subroutine
!===============================================================================
end module

! Routine for updating the driving modes and calculating averages on root grid
!===============================================================================
subroutine update_random_forcing
  use amr_commons
  use hydro_commons
  use forcing
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer     :: ilevel,i,j,k,jx,jy,jz
  complex(dp) :: lfxx,lfyy,lfzz,lkx,lky,lkz,corr,expikr
  real(dp)    :: dx, x, y, z, w, fact, test, fpow, fk, average, f
  integer     :: ivar,idim,ind,ncache,igrid,iskip
  integer     :: info,nleaf,ngrid,nx_loc
  integer, dimension(1:nvector)       :: ind_grid,ind_cell
  real(dp),dimension(1:nvector,4)     :: uu
  real(dp),dimension(1:twotondim,1:3) :: xc
  double precision :: rhoavp, pxavp, pyavp, pzavp, aranxp, aranyp, aranzp
  double precision,dimension(7)       :: comm_buffin, comm_buffout
  logical     :: first_call = .true.

  if (first_call) then
     call init_force
     first_call = .false.
  end if

  if (.not. do_force) return

!-------------------------------------------------------------------------------
  if (t .ne. tprev) then    ! prevent repeat at same t
     tprev = t

     if (myid==1) then
        print *,'forceit: t,t_turb=',t,t_turb
        if (t_turb .le. t) print *,'new force'
     end if

!--------------------------------------------- new random numbers:
!  Loop until t_turb > t
     do while (t_turb .le. t)
        ak_prev = ak_next                                                        ! advances a turnover time
        t_turb = t_turb + t_turn
        if (myid==1) print *,'random_force:t=',myid,t,t_turb,nmodes,iseed
        ! To obtain a Kolmogorov slope of the driving alone, one should have amplitudes
        ! a(k) that drop with k^(-11./6.), to have a(k)^2*k^2 = k^(-5./3.).  This ASSUMES
        ! that the amplitudes are proportional to the driving.  But the energy drain may
        ! be rather inversely proportional to the turnover time, which goes as k^(2./3.)

        ! pk = 11./6.
        ! pk = 7./6.

        fpow = 0.
        j = 0
        do jx=-kmax,kmax
           lkx =  cmplx (0., jx*2.*pi/boxlen)
           do jy=-kmax,kmax
              lky =  cmplx (0., jy*2.*pi/boxlen)
              do jz=-kmax,kmax
                 lkz =  cmplx (0., jz*2.*pi/boxlen)
                 fk = sqrt(float(jx**2+jy**2+jz**2))
                 if (fk .ge. k1-1e-5 .and. fk .le. k2+1e-5) then
                    lfxx = exp(cmplx(0., 2.*pi*ran1s(iseed)))/fk**pk
                    lfyy = exp(cmplx(0., 2.*pi*ran1s(iseed)))/fk**pk
                    lfzz = exp(cmplx(0., 2.*pi*ran1s(iseed)))/fk**pk

                    !------------------
                    ! solenoidal field:
                    if (do_helmh) then
                       corr=(lkx*lfxx+lky*lfyy+lkz*lfzz)/(lkx*lkx+lky*lky+lkz*lkz+1e-20)

                       if (jx.ne.0) lfxx = lfxx - a_helmh*corr*lkx
                       if (jy.ne.0) lfyy = lfyy - a_helmh*corr*lky
                       if (jz.ne.0) lfzz = lfzz - a_helmh*corr*lkz
                    end if
                    !------------------
                    j = j + 1
                    ak_next(1,j)=accel*lfxx
                    ak_next(2,j)=accel*lfyy
                    ak_next(3,j)=accel*lfzz

                    fact=1.
                    if (jx.gt.0) fact=fact*0.5
                    if (jy.gt.0) fact=fact*0.5
                    if (jz.gt.0) fact=fact*0.5
                    fpow = fpow+fact*(abs(lfxx)**2+abs(lfyy)**2+abs(lfzz)**2)
                 end if
              end do
           end do
        end do

     end do ! while

  end if   ! t > t_prev

  !------------------------------------------------------------------
  !  Time interpolation of the force
  w = (t_turn-(t_turb-t))/t_turn
  w = 0.5*(1.-cos(w*pi))

  ! Precompute interpolated accelerations in k-space
  do j=1, nmodes
     fxx(j) = ak_prev(1,j) * (1.0 - w) + ak_next(1,j) * w
     fyy(j) = ak_prev(2,j) * (1.0 - w) + ak_next(2,j) * w
     fzz(j) = ak_prev(3,j) * (1.0 - w) + ak_next(3,j) * w
  end do

  if (myid==1 .and. verbose) &
     print *,'force: t-t_turb, t_turn, w:', t-t_turb, t_turn, w

! Force and velocity averages done on the root grid
!  uu = (\rho, \rho u, \rho v, \rho w, Etot, Bx, By, Bz)
!--------------------------------------------------

  ilevel = levelmin

  ! Mesh size at level ilevel in coarse cell units
  dx=0.5D0**ilevel

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     jz=(ind-1)/4
     jy=(ind-1-4*jz)/2
     jx=(ind-1-2*jy-4*jz)
     xc(ind,1)=(dble(jx)-0.5D0)*dx
     xc(ind,2)=(dble(jy)-0.5D0)*dx
     xc(ind,3)=(dble(jz)-0.5D0)*dx
  end do

  ! Zero averages
  rhoavp=0.
  pxavp=0.;  pyavp=0.; pzavp=0.
  aranxp=0.; aranyp=0.; aranzp=0.

  ncache=active(ilevel)%ngrid

  ! Loop over grids
  !!$omp parallel do private(igrid,i,ngrid,ind_grid,ind,iskip,ivar,ind_cell,uu,x,y,z,j,expikr) &
  !!$omp&            reduction(+:rhoavp,pxavp,pyavp,pzavp,aranxp,aranyp,aranzp) &
  !!$omp&            default(none) shared(ncache,active,ncoarse,ngridmax,uold,xg,xc,nmodes,kx,ky,kz,fxx,fyy,fzz,ilevel)
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do

     ! Loop over cells in grid
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=ind_grid(i)+iskip
        end do

        ! Gather hydro variables (rho, px, py, pz)
        do ivar=1,4
           do i=1,ngrid
              uu(i,ivar)=uold(ind_cell(i),ivar)
           end do
        end do

        ! Compute average density and velocities
        do i=1,ngrid
           rhoavp = rhoavp + uu(i,1)
           pxavp  = pxavp  + uu(i,2)
           pyavp  = pyavp  + uu(i,3)
           pzavp  = pzavp  + uu(i,4)
        end do

        ! Compute average acceleration * rho on grid
        do i=1,ngrid
          ! Gather cell center positions
          x = xg(ind_grid(i),1)+xc(ind,1)
          y = xg(ind_grid(i),2)+xc(ind,2)
          z = xg(ind_grid(i),3)+xc(ind,3)

          ! Calculate acceleration * rho
          do j=1,nmodes
            expikr = exp(kx(j)*x+ky(j)*y+kz(j)*z)
            aranxp = aranxp + real(fxx(j)*expikr) * uu(i,1)
            aranyp = aranyp + real(fyy(j)*expikr) * uu(i,1)
            aranzp = aranzp + real(fzz(j)*expikr) * uu(i,1)
          end do
        enddo

     end do
     ! End loop over cells

  end do
  ! End loop over grids

  ! Compute global quantities
#ifndef WITHOUTMPI
  comm_buffin(1)=rhoavp
  comm_buffin(2)=pxavp
  comm_buffin(3)=pyavp
  comm_buffin(4)=pzavp
  comm_buffin(5)=aranxp
  comm_buffin(6)=aranyp
  comm_buffin(7)=aranzp
  call MPI_ALLREDUCE(comm_buffin,comm_buffout,7,MPI_DOUBLE_PRECISION,MPI_SUM,&
       &mpi_comm_use,info)
  rhoav = comm_buffout(1)
  pxav  = comm_buffout(2)
  pyav  = comm_buffout(3)
  pzav  = comm_buffout(4)
  aranx = comm_buffout(5)
  arany = comm_buffout(6)
  aranz = comm_buffout(7)
#else
  rhoav = rhoavp
  pxav  = pxavp
  pyav  = pyavp
  pzav  = pzavp
  aranx = aranxp
  arany = aranyp
  aranz = aranzp
#endif

  aranx = (aranx + pxav / t_turn) / rhoav
  arany = (arany + pyav / t_turn) / rhoav
  aranz = (aranz + pzav / t_turn) / rhoav

  if (myid==1 .and. verbose) write(*,'(a,4e12.4)') 'rho_av, xav, yav, zav :', &
     (/ rhoav / (real(numbtot(1,levelmin),kind=dp)*twotondim),aranx,arany,aranz/)

end subroutine update_random_forcing

! Routine doing the actual driving on grid ilevel
!===============================================================================
subroutine do_random_force(ilevel,dt)
  use forcing
  use amr_commons
  use hydro_commons
  implicit none
  integer,  intent(in) :: ilevel
  real(dp), intent(in) :: dt
  integer     :: i, j, jx, jy, jz
  complex(dp) :: expikr
  real(dp)    :: dx, x, y, z, w
  integer     :: ivar, idim, ind, ncache, igrid, iskip, ngrid, nleaf
  integer, dimension(1:nvector)      :: ind_grid, ind_cell, ind_leaf
  real(dp),dimension(1:nvector)      :: rho, ein, ax, ay, az
  real(dp),dimension(1:twotondim,1:3):: xc

  if (.not. do_force) return

  !  Time interpolation of the force
  w = (t_turn-(t_turb-t))/t_turn
  w = 0.5*(1.-cos(w*pi))

  ! Precompute interpolated accelerations in k-space
  do j=1, nmodes
     fxx(j) = ak_prev(1,j) * (1.0 - w) + ak_next(1,j) * w
     fyy(j) = ak_prev(2,j) * (1.0 - w) + ak_next(2,j) * w
     fzz(j) = ak_prev(3,j) * (1.0 - w) + ak_next(3,j) * w
  end do

  ! Mesh size at level ilevel in units of coarse grid cells
 dx = 0.5D0**ilevel

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     jz=(ind-1)/4
     jy=(ind-1-4*jz)/2
     jx=(ind-1-2*jy-4*jz)
     xc(ind,1)=(dble(jx)-0.5D0)*dx
     xc(ind,2)=(dble(jy)-0.5D0)*dx
     xc(ind,3)=(dble(jz)-0.5D0)*dx
  end do

  ncache=active(ilevel)%ngrid

  ! Loop over grids
  !!$omp parallel do default(none) schedule(dynamic) private(igrid,ngrid, &
  !!$omp&            i,ind_grid,ind,iskip,nleaf,j,ind_cell,rho,ein,x,y,z, &
  !!$omp&            expikr,ax,ay,az) shared(ncache,active,ilevel,ncoarse,&
  !!$omp&            ngridmax,son,uold,xg,xc,nmodes,kx,ky,kz,fxx,fyy,fzz, &
  !!$omp&            dt,aranx,arany,aranz)
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do

     ! Loop over cells in grid
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        ! Gather leaf cells
        nleaf=0
        do i=1,ngrid
           j=ind_grid(i)+iskip
           if(son(j)==0)then
              nleaf=nleaf+1
              ind_cell(nleaf)=j
           end if
        end do

        ! Gather hydro variables (rho, internal energy)
        do i=1,nleaf
           rho(i)=uold(ind_cell(i),1)
           ein(i)=uold(ind_cell(i),5)-0.5*( &
                    uold(ind_cell(i),2)*uold(ind_cell(i),2) &
                  + uold(ind_cell(i),3)*uold(ind_cell(i),3) &
                  + uold(ind_cell(i),4)*uold(ind_cell(i),4) &
                                          )/uold(ind_cell(i),1)
        end do

        ! Compute average acceleration on grid
        do i=1,nleaf
           ! Gather cell center positions
           j = ind_cell(i)-iskip
           x = xg(j,1)+xc(ind,1)
           y = xg(j,2)+xc(ind,2)
           z = xg(j,3)+xc(ind,3)
           ! Calculate acceleration
           expikr = exp(kx(1)*x+ky(1)*y+kz(1)*z)
           ax(i) = real(fxx(1)*expikr)
           ay(i) = real(fyy(1)*expikr)
           az(i) = real(fzz(1)*expikr)
           do j=2,nmodes
              expikr = exp(kx(j)*x+ky(j)*y+kz(j)*z)
              ax(i) = ax(i) + real(fxx(j)*expikr)
              ay(i) = ay(i) + real(fyy(j)*expikr)
              az(i) = az(i) + real(fzz(j)*expikr)
           end do
           ax(i) = rho(i) * (ax(i) - aranx) * dt
           ay(i) = rho(i) * (ay(i) - arany) * dt
           az(i) = rho(i) * (az(i) - aranz) * dt
        end do

        ! Update momentum and total energy
        do i=1,nleaf
           uold(ind_cell(i),2)=uold(ind_cell(i),2) + ax(i)
           uold(ind_cell(i),3)=uold(ind_cell(i),3) + ay(i)
           uold(ind_cell(i),4)=uold(ind_cell(i),4) + az(i)
           uold(ind_cell(i),5) = ein(i) + 0.5*( &
                    uold(ind_cell(i),2)*uold(ind_cell(i),2) &
                  + uold(ind_cell(i),3)*uold(ind_cell(i),3) &
                  + uold(ind_cell(i),4)*uold(ind_cell(i),4) &
                                              )/uold(ind_cell(i),1)
        end do
     end do
     ! End loop over cells

  end do
  ! End loop over grids
end subroutine do_random_force
!===============================================================================
