!###########################################################
!###########################################################
!###########################################################
SUBROUTINE cr_courant_fine(ilevel)

   ! Cosmic-ray moment-transport contribution to the level timestep.
   ! Extracted verbatim from mhd/courant_fine.f90 (the #ifdef CRPHYS block
   ! formerly at lines 234-268): refresh cr_vmax/Dcr_code in code units,
   ! apply the adaptive reduced-light-speed cap (cr_varvmax), clamp to c,
   ! and reduce dtnew(ilevel) by the CR Courant condition. Called from the
   ! tail of courant_fine so dtnew(ilevel) already holds the gas global-min
   ! that the cr_varvmax cap reads, and cr_va_max already holds
   ! the values cmpdt produced during the gas sweep. (This CR<->gas-Courant
   ! ordering is why CR is unlike the gas-independent RT timestep.)
   !-------------------------------------------------------------------------
   use amr_commons
   use cr_parameters, only: cr_vmax,c_cu,cr_nsubcycle,cr_varvmax, &
        & cr_varvmax_fudge,cr_varvmax_vdvs,cr_va_max, &
        & gamma_cr,Dcr_code,mom_streaming_diffusion,ncr
   implicit none
   integer::ilevel
   real(dp)::dt_cr,dx,scale
   integer::nx_loc,igrp
   !-------------------------------------------------------------------------
   ! Mesh spacing at this level (matches courant_fine:57-60)
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
         if(mom_streaming_diffusion)then
            do igrp=1,ncr
               cr_vmax(ilevel) = max(cr_vmax(ilevel), gamma_cr(igrp)*cr_va_max * cr_varvmax_fudge)
            end do
         endif
         do igrp=1,ncr
            cr_vmax(ilevel) = max(cr_vmax(ilevel), Dcr_code(igrp)/dx * (gamma_cr(igrp)-1d0) * cr_varvmax_fudge)
         end do
      endif
   endif

   ! Finally, make sure vmax <= c
   if(cr_vmax(ilevel) .gt. c_cu) cr_vmax(ilevel)=c_cu

   call get_crmom_courant(dt_cr,ilevel)
   dtnew(ilevel) = MIN(dtnew(ilevel), dt_cr * cr_nsubcycle*0.99999d0)

END SUBROUTINE cr_courant_fine
