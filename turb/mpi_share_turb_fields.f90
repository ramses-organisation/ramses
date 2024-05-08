#if USE_TURB==1
! MPI_SHARE_TURB_FIELDS.F90
! A.McLeod - 24/11/2012
! Share turbulent acceleration fields to non-root tasks
! ============================================================================
#ifndef WITHOUTMPI
subroutine mpi_share_turb_fields(include_last)
  use turb_commons
  use mpi_mod
  implicit none
  logical, intent(in) :: include_last

  ! Set MPI_REAL_DP to the appropriate MPI type integer
#ifndef NPRE
  integer, parameter :: MPI_REAL_DP=MPI_REAL
#else
#if NPRE==4
  integer, parameter :: MPI_REAL_DP=MPI_REAL
#else
  integer, parameter :: MPI_REAL_DP=MPI_DOUBLE_PRECISION
#endif
#endif

   integer, parameter :: message_length=NDIM*TURB_GS**NDIM
                                               ! Length of flattened arrays
   integer            :: ierr                  ! MPI error variable

   ! Share afield_last and afield_next
   if (include_last) then
      call MPI_BCAST(afield_last, message_length, MPI_REAL_DP, 0, &
                     & MPI_COMM_WORLD, ierr)
      call MPI_BCAST(turb_last_time, 1, MPI_REAL_DP, 0, &
                     & MPI_COMM_WORLD, ierr)
   end if

   call MPI_BCAST(afield_next, message_length, MPI_REAL_DP, 0, &
                  & MPI_COMM_WORLD, ierr)
   call MPI_BCAST(turb_next_time, 1, MPI_REAL_DP, 0, &
                  & MPI_COMM_WORLD, ierr)

end subroutine mpi_share_turb_fields
#endif
#endif
