!###############################################################
!###############################################################
subroutine crmom_step(ilevel)
  ! Two-moment cosmic-ray transport driver: subcycled explicit transport
  ! + implicit scattering source terms + cooling, ported from ramses_cral
  ! feat/CR_tests cr/cr_godunov_fine.f90.
  ! Phase 0: no-op scaffold. The Courant loop, cr_godunov_fine,
  ! add_cr_source_terms and cr_cooling_fine land in phases 1-4.
  use amr_commons
  use cr_parameters
  use cr_hydro_commons
  implicit none
  integer,intent(in)::ilevel

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

111 format('   Entering crmom_step for level ',I2)
end subroutine crmom_step
