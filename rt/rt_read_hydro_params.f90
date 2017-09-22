!************************************************************************
SUBROUTINE rt_read_hydro_params()

! 'Interface' for reading any additional parameters from .nml file.
! This routine can be overridden by any patch to include more parameters.
!------------------------------------------------------------------------
  implicit none
!------------------------------------------------------------------------
  call read_rt_params()

END SUBROUTINE rt_read_hydro_params




