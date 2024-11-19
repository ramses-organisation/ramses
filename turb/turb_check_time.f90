#if USE_TURB==1
subroutine turb_check_time
   use amr_commons
   use turb_commons
   implicit none

   real(kind=dp) :: turb_last_tfrac        ! Time fraction since last
   real(kind=dp) :: turb_next_tfrac        ! Time fraction until next

   if (turb_type == 3) then
      ! decaying turbulence
      afield_now = 0.0_dp
      fturb = 0.0_dp
      turb_next_time = huge(turb_next_time) / 10.0d0
   else if (turb_type == 1) then
      ! evolving forced turbulence
      do
         if (t >= turb_next_time) then
#ifndef WITHOUTMPI
            if (myid==1) then
               call turb_next_field
            else
               turb_last_time = turb_next_time
               afield_last = afield_next
            end if
            call mpi_share_turb_fields(.FALSE.)
#else
            call turb_next_field
#endif
         else
            exit
         end if
      end do

      turb_last_tfrac = real((t - turb_last_time) / turb_dt, dp)
      turb_next_tfrac = 1.0_dp - turb_last_tfrac

      afield_now = turb_last_tfrac*afield_last + turb_next_tfrac*afield_next
   end if

end subroutine turb_check_time
#endif
