!================================================================
! Jiang & Oh (2018) cosmic-ray test suite -- CR initial conditions patch.
!
! Patch override of cr/cr_condinit.f90 for the separated-CR module: it sets
! ONLY the CR conservative variables u(1:nn,1:ncrvars) (energy at iCRu=1, the
! ndim fluxes in the slots after it). The gas and magnetic field come from the
! region-based &INIT_PARAMS block in the namelist (region_condinit, unchanged),
! so no gas condinit override is needed for the region-based tests.
!
! The test is chosen at run time with the namelist string jiang_test (in
! &CR_PARAMS). jiang_test='' reproduces the default cr_condinit (smallcr energy
! floor, zero flux), so this patch is inert for non-test runs.
!
! Ported from ramses_cral patch/jiang_tests/condinit.f90 (the jiang_cr_init CR
! block), applying the central transformation: cral wrote the CR variables into
! the embedded u(:,icrU...) with icrU=nvar+4; here we write the separated CR
! buffer u(:,iCRu...) with iCRu=1. The CR advective flux F = 4/3 v_gas E_cr is
! zero for the static-gas region tests below (u_region=0), so the gas velocity
! is not needed in this CR-only buffer; flux is set to zero, identical to cral
! for these tests.
!
! Test '422' (colliding flow) takes its CR ENERGY from crmom_region per region
! geometry -- cral does this in region_condinit, which the separated CR module
! does not call -- so cr_region_condinit below replicates region_condinit's
! square-region geometry into the CR buffer, then the '422' branch sets the CR
! advective flux F = 4/3 v_gas E_cr from the namelist region velocity u_region
! (= the gas velocity q(:,2) cral reads), giving the identical CR state.
!
! Test 'TP_1D_shock' (two-pressure shock tube) is also region-based: like 422 the
! CR ENERGY (and here the zero flux) come from crmom_region per region geometry,
! but unlike 422 the CR flux is taken DIRECTLY from crmom_region(:,2)=0 (not from
! the gas velocity), so cr_region_condinit alone reproduces cral's region_condinit.
!
!   1D tests covered here: 411  411_triangular  413  414  421  422  424  TP_1D_shock
!================================================================
subroutine cr_condinit(x,u,dx,nn,ilevel)
  use amr_parameters
  use cr_parameters
  implicit none
  integer ::nn                              ! Number of cells
  integer:: ilevel                          ! Refinement level
  real(dp)::dx                              ! Cell size
  real(dp),dimension(1:nvector,1:ncrvars)::u ! CR conservative variables
  real(dp),dimension(1:nvector,1:ndim  )::x ! Cell center position in [0,boxlen]**ndim
  !----------------------------------------------------------------
  ! u(i,1:ncrvars) is the CR conservative vector for group g:
  !   energy at iCRu+(ndim+1)*(g-1), the ndim fluxes in the slots after it.
  !   (iCRu=1; CR indexing never references nvar.)
  !----------------------------------------------------------------
  integer::i,igrp,icrE
  real(dp),dimension(1:nvector)::tmp
  real(dp)::pi

  ! Default for every group: energy floor, zero flux. Test branches below
  ! overwrite the first group's energy (and, where applicable, flux).
  do i=1,nn
     do igrp=1,ncr
        icrE=iCRu+(ndim+1)*(igrp-1)
        u(i,icrE)=smallcr
        u(i,icrE+1:icrE+ndim)=0d0
     end do
  end do

  select case(trim(jiang_test))

  case('')
     ! No override: pure default (smallcr floor, zero flux).

  case('411','414')                          ! Gaussian CR-energy pulse
     ! cral: u(:,icrU)=exp(-40*(x-boxlen/2)^2); flux=4/3*v_gas*E_cr=0 (static gas).
     tmp(1:nn)=(x(1:nn,1)-boxlen*0.5d0)**2
     do i=1,nn
        u(i,iCRu)=exp(-40d0*tmp(i))
     end do

  case('411_triangular')                     ! Triangular CR-energy profile
     ! cral: u(:,icrU)=E0-slope*|x-boxlen/2| (E0=2,slope=1); the advective flux
     ! u(:,icrU+1)=4/3*v_gas*E_cr is zero here because the gas is static
     ! (u_region=0), so it stays at the default-zero set above. The streaming
     ! wave is driven entirely by the time-dependent analytic boundary
     ! (cr_boundana, bound_type=3). E0/slope match jiang_cr_init('411_triangular').
     tmp(1:nn)=(x(1:nn,1)-boxlen*0.5d0)**2
     do i=1,nn
        u(i,iCRu)=2d0-1d0*sqrt(tmp(i))
     end do

  case('421')                                ! Sinusoidal CR-energy profile
     pi=acos(-1d0)
     do i=1,nn
        u(i,iCRu)=20d0+10d0*sin(pi*(x(i,1)-boxlen*0.5d0))
     end do

  case('413','424')                          ! CR floor (CR injected at reflexive BC)
     ! Density step is part of the GAS state (region/condinit); here CR only.
     do i=1,nn
        u(i,iCRu)=1d-6
        u(i,iCRu+1:iCRu+ndim)=0d0
     end do

  case('422')                                ! CR energy from crmom_region; flux from gas velocity
     ! cral: region_condinit fills the CR energy E_cr=crmom_region(k,1) inside
     ! each region (here crmom_region(:,1)=50,50); then jiang_cr_init sets the
     ! CR advective flux u(:,icrU+1)=4/3 q(:,2) E_cr from the gas velocity q(:,2).
     ! The separated CR module never calls region_condinit, so we replicate its
     ! square-region geometry: cr_region_condinit fills E_cr (and any crmom flux),
     ! then we overwrite the x-flux with 4/3 u_region(k) E_cr per region. Inside
     ! each region q(:,2)=u_region(k) (region_condinit), so this matches cral.
     call cr_region_condinit(x,u,dx,nn)
     call cr_flux_from_region_velocity(x,u,dx,nn)

  case('TP_1D_shock')                        ! Two-pressure shock tube: CR energy+flux from crmom_region
     ! cral: a region-based test with NO condinit override -- jiang_cr_init is a
     ! no-op; the entire CR state comes from region_condinit reading crmom_region
     ! (here crmom_region(:,1)=3,1 for the energy and crmom_region(:,2)=0,0 for
     ! the x-flux), giving a CR-energy contact discontinuity (3|1) across x=5 with
     ! zero initial flux. The separated CR module never calls region_condinit, so
     ! cr_region_condinit replicates its square-region geometry into the CR buffer
     ! (filling BOTH the energy at iCRu and the flux at iCRu+1 from crmom_region),
     ! giving the identical CR state. No gas-velocity flux is needed: u_region=0,
     ! so cral's flux would be zero anyway, and crmom_region(:,2)=0 sets it here.
     call cr_region_condinit(x,u,dx,nn)

  case default
     write(*,*)'cr_condinit: unknown jiang_test = "'//trim(jiang_test)//'"'
     call clean_stop
  end select

end subroutine cr_condinit
!================================================================
!================================================================
!================================================================
!================================================================
subroutine cr_region_condinit(x,u,dx,nn)
  ! Fill the separated CR buffer u(1:nn,1:ncrvars) from the namelist
  ! crmom_region per region geometry. Faithful port of the CR-handling part
  ! of ramses_cral mhd/init_flow_fine.f90 region_condinit (the square/point
  ! region loop), keeping ONLY the CR slots: cral wrote q(i,icru+ivar-1) with
  ! icru=nvar+4; here we write u(i,ivar) with the separated layout (iCRu=1),
  ! ivar=1..ncrvars maps energy + ndim fluxes per group identically.
  ! Default (no region covers a cell) leaves the smallcr floor / zero flux
  ! that cr_condinit set before calling this routine; cral's default is 0,
  ! but every 422 cell is covered by a region so the difference never shows.
  use amr_parameters
  use cr_parameters
  implicit none
  integer ::nn
  real(dp)::dx
  real(dp),dimension(1:nvector,1:ncrvars)::u
  real(dp),dimension(1:nvector,1:ndim  )::x
  integer::i,k,ivar
  real(dp)::vol,r,xn,yn,zn,en

  ! Loop over initial conditions regions
  do k=1,nregion

     ! For "square" regions only:
     if(region_type(k) .eq. 'square')then
        en=exp_region(k)
        do i=1,nn
           xn=0.0d0; yn=0.0d0; zn=0.0d0
           xn=2.0d0*abs(x(i,1)-x_center(k))/length_x(k)
#if NDIM>1
           yn=2.0d0*abs(x(i,2)-y_center(k))/length_y(k)
#endif
#if NDIM>2
           zn=2.0d0*abs(x(i,3)-z_center(k))/length_z(k)
#endif
           if(exp_region(k)<10)then
              r=(xn**en+yn**en+zn**en)**(1.0/en)
           else
              r=max(xn,yn,zn)
           end if
           ! If cell lies within region, REPLACE CR variables by region values
           if(r<1.0)then
              do ivar=1,ncrvars
                 u(i,ivar)=crmom_region(k,ivar)
              end do
           end if
        end do
     end if

     ! For "point" regions only:
     if(region_type(k) .eq. 'point')then
        vol=dx**ndim
        do i=1,nn
           xn=1.0; yn=1.0; zn=1.0
           xn=max(1.0-abs(x(i,1)-x_center(k))/dx,0.0_dp)
#if NDIM>1
           yn=max(1.0-abs(x(i,2)-y_center(k))/dx,0.0_dp)
#endif
#if NDIM>2
           zn=max(1.0-abs(x(i,3)-z_center(k))/dx,0.0_dp)
#endif
           r=xn*yn*zn
           ! ADD region CR energy density (cral adds crmom_region(k,1)*r/vol)
           u(i,iCRu)=u(i,iCRu)+crmom_region(k,1)*r/vol
        end do
     end if
  end do

end subroutine cr_region_condinit
!================================================================
!================================================================
!================================================================
!================================================================
subroutine cr_flux_from_region_velocity(x,u,dx,nn)
  ! jiang_test='422': set the CR advective flux F = 4/3 v_gas E_cr from the gas
  ! velocity. cral does u(:,icrU+1)=4/3*q(:,2)*u(:,icrU) where q(:,2) is the gas
  ! x-velocity = u_region(k) inside region k (region_condinit). The separated CR
  ! buffer has no gas state, so we recover the same per-region velocity from the
  ! namelist region geometry (identical square/point test as region_condinit).
  ! Only the first CR group's x-flux is set, matching jiang_cr_init('422').
  use amr_parameters
  use cr_parameters
  use hydro_parameters, only: u_region
  implicit none
  integer ::nn
  real(dp)::dx
  real(dp),dimension(1:nvector,1:ncrvars)::u
  real(dp),dimension(1:nvector,1:ndim  )::x
  integer::i,k
  real(dp)::r,xn,yn,zn,en

  do k=1,nregion
     if(region_type(k) .ne. 'square')cycle
     en=exp_region(k)
     do i=1,nn
        xn=0.0d0; yn=0.0d0; zn=0.0d0
        xn=2.0d0*abs(x(i,1)-x_center(k))/length_x(k)
#if NDIM>1
        yn=2.0d0*abs(x(i,2)-y_center(k))/length_y(k)
#endif
#if NDIM>2
        zn=2.0d0*abs(x(i,3)-z_center(k))/length_z(k)
#endif
        if(exp_region(k)<10)then
           r=(xn**en+yn**en+zn**en)**(1.0/en)
        else
           r=max(xn,yn,zn)
        end if
        if(r<1.0)then
           u(i,iCRu+1)=4d0/3d0*u_region(k)*u(i,iCRu)
        end if
     end do
  end do

end subroutine cr_flux_from_region_velocity
