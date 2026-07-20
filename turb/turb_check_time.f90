subroutine turb_check_time
   use amr_commons
   use turb_commons
   implicit none

   real(kind=dp) :: turb_tfrac

   select case (turb_type)
   case (3)
      ! decaying turbulence - the initial field set up by init_turb has already
      ! been applied as an initial velocity, so no forcing from here on
      afield_now = 0.0_dp
      fturb = 0.0_dp
      turb_next_time = huge(turb_next_time) / 10.0d0
   case (2)
      ! fixed forced turbulence - afield_now is the static field chosen by
      ! init_turb and is deliberately left untouched for the whole run
      continue
   case (1)
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

      ! Time fraction since last turbulence field evaluation
      turb_tfrac = real((t - turb_last_time) / turb_dt, dp)

      ! interpolate for current time between last and next turb field
      afield_now = (1.0_dp - turb_tfrac)*afield_last + turb_tfrac*afield_next
   end select

end subroutine turb_check_time
