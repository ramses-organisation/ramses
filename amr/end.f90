
subroutine clean_end
  !---------------------------
  ! Properly end the run. 
  !---------------------------
  ! use amr_commons
  ! use poisson_commons
  ! use pm_commons
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  character(LEN=80)::str

  call output_timer(.false., str)

#ifndef WITHOUTMPI
  call MPI_FINALIZE(info)
#endif

  call deallocate_amr
  call deallocate_pm
  call deallocate_poisson

  stop
end subroutine clean_end

subroutine clean_stop
  !-----------------------------------------------------
  ! This subroutine brings the program to a halt after
  ! an error.
  !-----------------------------------------------------
  use amr_commons
  use poisson_commons
  use pm_commons
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif

#ifndef WITHOUTMPI
  call MPI_ABORT(MPI_COMM_WORLD, 2, info)
#endif

  call deallocate_amr
  call deallocate_pm
  call deallocate_poisson

  stop
end subroutine clean_stop
