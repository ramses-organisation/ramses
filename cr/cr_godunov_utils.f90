!###########################################################
!###########################################################
!###########################################################
SUBROUTINE get_crmom_courant(dt,ilevel)

   ! Determine the coarse CR timestep length set by the Courant condition
   !-------------------------------------------------------------------------
     use amr_parameters
     use cr_parameters, only: cr_vmax
     use hydro_parameters, only: courant_factor
     implicit none
     integer:: nx_loc,ilevel
     real(dp):: dt, scale, dx
   !-------------------------------------------------------------------------
     ! Mesh spacing at coarse level
     nx_loc=icoarse_max-icoarse_min+1
     scale=boxlen/dble(nx_loc)
     dx=0.5D0**ilevel*scale
     dt=courant_factor*dx/3d0/cr_vmax(ilevel)
END SUBROUTINE get_crmom_courant
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE update_cr_vmax_and_Dcr_code(crvmax)

! Update the maximum CR speed, in code units.
! This cannot be just a constant, since scale_v changes with time in
! cosmological simulations.
! Additionally update diffusion coefficient in code units.
!-------------------------------------------------------------------------
   use constants, only: c_cgs
   use cr_parameters
   implicit none
   real(dp),intent(out)::crvmax
   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_kappa
   integer::igrp
!-------------------------------------------------------------------------
    call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
    crvmax=c_cgs/scale_v * cr_c_fraction
    cr_c_code = c_cgs/scale_v

    scale_kappa = scale_l**2/scale_t
    DCRmax_code=DCRmax/scale_kappa
    do igrp=1,ncr_groups
       ! YD: I think the factor 3 should not appear there.
       ! It would be better to have it in cr_godunov_fine.f90
       ! in source term because sigma = 1/(3*kappa_E).
       ! Here that means that our definition of DCR_code
       ! corresponds to DCR_code=1/sigma. That could be
       ! a bit confusing when going through the code.
       DCR_code(igrp)=max(cr_smalld,DCR(igrp)/scale_kappa*3d0)
    end do

END SUBROUTINE update_cr_vmax_and_Dcr_code
