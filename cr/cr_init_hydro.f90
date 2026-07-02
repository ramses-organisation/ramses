!*************************************************************************
subroutine cr_init_hydro
  ! Allocate the cosmic-ray state arrays (separate from the hydro
  ! uold/unew) and initialise them to zero.
  use amr_commons
  use cr_parameters
  use cr_hydro_commons
  implicit none
  integer::ncell

  if(verbose)write(*,*)'Entering cr_init_hydro'

  ncell=ncoarse+twotondim*ngridmax
  allocate(cruold(1:ncell,1:ncrvar))
  allocate(crunew(1:ncell,1:ncrvar))
  cruold=0d0 ; crunew=0d0

  if(verbose)write(*,*)'Allocate done for',ncrvar,'CR variables'

end subroutine cr_init_hydro
