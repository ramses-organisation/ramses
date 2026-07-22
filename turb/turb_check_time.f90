subroutine turb_check_time
   use amr_commons
   use turb_commons
   implicit none

   ! A non-evolving field is the static one chosen by init_turb and is
   ! deliberately left untouched for the whole run
   if (turb_evolving) then
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

      ! interpolate for current time between last and next turb field
      call turb_interpolate_now

      ! Optionally hold the injected rms exactly on turb_rms while the
      ! pattern keeps evolving
      if (turb_exact_rms) call turb_normalise_rms
   end if

end subroutine turb_check_time
