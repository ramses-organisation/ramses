#if USE_TURB==1
! TURB_NEXT_FIELD.F90
! A.McLeod - 24/11/2012
! Generate next turbulent field
! ============================================================================
subroutine turb_next_field
   use turb_commons
   implicit none

   ! Set turbulent field times
   turb_last_time = turb_next_time
   turb_next_time = turb_last_time + turb_dt

   ! Copy current field, and add more turbulence
   turb_last = turb_next
   afield_last = afield_next
   call add_turbulence(turb_next, turb_dt)

   ! Subtract decay of turbulence
   call decay_turbulence(turb_last, turb_next)

   ! Fourier transform
#if NDIM==1
   call FFT_1D(turb_next(1,:,0,0), afield_next(1,:,0,0))
#elif NDIM==2
   call FFT_2D(turb_next(:,:,:,0), afield_next(:,:,:,0))
#else
   call FFT_3D(turb_next, afield_next)
#endif
   afield_next = afield_next * turb_norm * turb_rms

end subroutine turb_next_field
#endif
