!*************************************************************************
subroutine cr_init_hydro
  ! Allocate the cosmic-ray state arrays (separate from the hydro
  ! uold/unew) and initialise them to the smallcr floor.
  ! Restart read of CR snapshot data is added in phase 5; until then a
  ! restart re-initialises CR fields to the floor.
  use amr_commons
  use cr_parameters
  use cr_hydro_commons
  implicit none
  integer::ncell

  if(verbose)write(*,*)'Entering cr_init_hydro'

  ncell=ncoarse+twotondim*ngridmax
  allocate(cruold(1:ncell,1:ncrvars))
  allocate(crunew(1:ncell,1:ncrvars))
  cruold=smallcr ; crunew=smallcr

  if(verbose)write(*,*)'Allocate done for',ncrvars,'CR variables'

end subroutine cr_init_hydro
