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
!   1D tests covered here: 411  411_triangular  413  414  421  424
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

  case default
     write(*,*)'cr_condinit: unknown jiang_test = "'//trim(jiang_test)//'"'
     call clean_stop
  end select

end subroutine cr_condinit
