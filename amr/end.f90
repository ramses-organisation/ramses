
subroutine clean_end
  !---------------------------
  ! Properly end the run.
  !---------------------------
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

subroutine deallocate_amr

  use amr_commons
  implicit none
  integer :: ilevel

  ! allocations in read_params.f90
  if(allocated(remap_pscalar)) deallocate(remap_pscalar)


  ! allocations in init_amr.f90
  if(allocated(bound_key)) deallocate(bound_key)
  if(allocated(bound_key2)) deallocate(bound_key2)
  if(allocated(headl)) deallocate(headl)
  if(allocated(taill)) deallocate(taill)
  if(allocated(numbl)) deallocate(numbl)
  if(allocated(numbtot)) deallocate(numbtot)
  if(allocated(headb)) deallocate(headb)
  if(allocated(tailb)) deallocate(tailb)
  if(allocated(numbb)) deallocate(numbb)
  if(allocated(boundary)) deallocate(boundary)
  ! communicators
  if(allocated(active))then
     do ilevel=1,nlevelmax
        ! virtual_boundaries.f90
        ! TODO: cleaner solution (fortran 2003): s/pointer/allocatable/
        if(active(ilevel)%ngrid>0) deallocate(active(ilevel)%igrid)
     enddo
     deallocate(active)
  endif
  if(allocated(emission)) deallocate(emission)
#ifdef LIGHT_MPI_COMM
  if(allocated(emission_part)) deallocate(emission_part)
#endif
  if(allocated(reception)) deallocate(reception)
  !
  if(allocated(father)) deallocate(father)
  if(allocated(nbor)) deallocate(nbor)
  if(allocated(next)) deallocate(next)
  if(allocated(prev)) deallocate(prev)

  if(allocated(xg)) deallocate(xg)
  ! amr cell-based arrays
  if(allocated(flag1)) deallocate(flag1)
  if(allocated(flag2)) deallocate(flag2)
  if(allocated(son)) deallocate(son)
  ! mpi cell-based arrays
  if(allocated(cpu_map)) deallocate(cpu_map)
  if(allocated(cpu_map2)) deallocate(cpu_map2)
  if(allocated(hilbert_key)) deallocate(hilbert_key)

  ! allocations in init_time.f90
  if(allocated(aexp_frw)) deallocate(aexp_frw)
  if(allocated(hexp_frw)) deallocate(hexp_frw)
  if(allocated(tau_frw)) deallocate(tau_frw)
  if(allocated(t_frw)) deallocate(t_frw)

end subroutine deallocate_amr

subroutine deallocate_pm
  use pm_commons
  use amr_commons, only: pic
  implicit none

  if(pic)then
     if(allocated(headp)) deallocate(headp)
     if(allocated(tailp)) deallocate(tailp)
     if(allocated(numbp)) deallocate(numbp)

     ! init_part.f90 - in general, BIG deallocations
     if(allocated(idp)) deallocate(idp)
     if(allocated(nextp)) deallocate(nextp)
     if(allocated(prevp)) deallocate(prevp)
     if(allocated(levelp)) deallocate(levelp)
     if(allocated(mp)) deallocate(mp)
     if(allocated(vp)) deallocate(vp)
     if(allocated(xp)) deallocate(xp)
  endif

end subroutine deallocate_pm

subroutine deallocate_poisson
  use poisson_commons
  implicit none

  if(allocated(lookup_mg)) deallocate(lookup_mg)
  if(allocated(safe_mode)) deallocate(safe_mode)
  if(allocated(active_mg)) deallocate(active_mg)
  if(allocated(emission_mg)) deallocate(emission_mg)

  ! cell-centred variables
  if(allocated(rho)) deallocate(rho)
  if(allocated(phi)) deallocate(phi)
  if(allocated(phi_old)) deallocate(phi_old)
  if(allocated(f)) deallocate(f)

end subroutine deallocate_poisson
