module mpi_mod
#ifndef WITHOUTMPI
#ifdef MPI_OLD
  include 'mpif.h'
#else
  use mpi
#endif
#endif
end module mpi_mod
