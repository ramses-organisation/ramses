!###########################################################
!###########################################################
!###########################################################
SUBROUTINE cr_courant_fine(ilevel)

   ! Cosmic-ray moment-transport contribution to the level timestep: refresh
   ! cr_vmax/Dcr_code, apply the cr_varvmax cap, clamp to c, reduce dtnew by the
   ! CR Courant condition. Must run at the tail of courant_fine, after the gas
   ! sweep has set dtnew(ilevel) (gas global-min) and cr_va_max that the cap reads.
   !-------------------------------------------------------------------------
   use amr_commons
   use cr_parameters, only: cr_vmax,cr_c_code,cr_nsubcycle,cr_varvmax, &
        & cr_varvmax_fudge,cr_varvmax_vdvs,cr_va_max, &
        & gamma_cr,Dcr_code,cr_streaming_diffusion,ncr_groups
   implicit none
   integer::ilevel
   real(dp)::dt_cr,dx,scale
   integer::nx_loc,igrp
   !-------------------------------------------------------------------------
   ! Mesh spacing at this level
   nx_loc=icoarse_max-icoarse_min+1
   scale=boxlen/dble(nx_loc)
   dx=0.5D0**ilevel*scale

   call update_cr_vmax_and_Dcr_code(cr_vmax(ilevel))

   ! Adaptive reduced light speed (cr_varvmax): keep cr_vmax at least
   ! cr_varvmax_fudge times larger than the gas signal speed dx/3/dt so the
   ! CR wave always outruns the moving gas. Inert when cr_varvmax=.false.
   if(cr_varvmax)then
      cr_vmax(ilevel) = max(cr_vmax(ilevel), dx/3d0/dtnew(ilevel) * cr_varvmax_fudge)
      if(cr_varvmax_vdvs)then
         if(cr_streaming_diffusion)then
            do igrp=1,ncr_groups
               cr_vmax(ilevel) = max(cr_vmax(ilevel), gamma_cr(igrp)*cr_va_max * cr_varvmax_fudge)
            end do
         endif
         do igrp=1,ncr_groups
            cr_vmax(ilevel) = max(cr_vmax(ilevel), Dcr_code(igrp)/dx * (gamma_cr(igrp)-1d0) * cr_varvmax_fudge)
         end do
      endif
   endif

   ! Finally, make sure vmax <= c
   if(cr_vmax(ilevel) .gt. cr_c_code) cr_vmax(ilevel)=cr_c_code

   call get_crmom_courant(dt_cr,ilevel)
   dtnew(ilevel) = MIN(dtnew(ilevel), dt_cr * cr_nsubcycle*0.99999d0)

END SUBROUTINE cr_courant_fine
