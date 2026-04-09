!##################################################################################################
!##################################################################################################
!##################################################################################################
!##################################################################################################

!  Function PLANCK_ANA:
!
!> Compute Planck average opacity.
!<
function planck_ana(dens,Tp,Tr,igroup,insink)

  use amr_commons
  use radiation_parameters
  use constants, only:pi

  implicit none

  integer ,intent(in)    :: igroup
  real(dp),intent(in)    :: dens,Tp,Tr
  logical, intent(in)    :: insink
  real(dp)               :: planck_ana,Tgd
  real(dp)               :: Tevap ! if sublimation_kuiper
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  Tgd = Tp
  if(sublimation_kuiper) then
     Tevap = 2000.0d0*dens**0.0195  !! Evaporation temperature
     if(Tp .gt. Tevap) Tgd = Tevap
  endif

  if(fit_semenov)then
     planck_ana = dens*0.5d0 * ( &
         (1.d0 - dtanh((Tgd - 145.d0)/2.d0)) * 0.001d0 * Tgd**((10.d0 - Tgd)/750.d0 + 2.d0) + &
         (1.d0 + dtanh((Tgd - 145.d0)/2.d0)) * &
         (2.25d0 - 1.8d0 * dtanh((Tgd - 1250.d0)/50.d0)))
  else 
     planck_ana = planck_params(1)*(dens**planck_params(2))*(Tgd**planck_params(3))
  endif  

  if(sublimation_kuiper) then
     !# RMR #### Sublimation of dust grains as in kuiper+10 ApJ ####
     !### The highest dust temperature is the evaporation temperature
     !### Opacities are taken at this temperature
     !### The input temperature is kept to compute the d/g ratio
     !# RMR ## Sublimation mimicked by a d/g ratio that decreases as a arctan function centered on Tevap ##
     planck_ana = planck_ana*(0.5d0 - 1./pi*atan(0.01d0*(Tp - Tevap) ) ) & !planck_ana contains the d/g ratio of 0.01
          +dens*0.01d0*(1.0d0-0.01d0*(0.5d0 - 1./pi*atan(0.01d0*(Tp - Tevap) ) )) !=> quasi full gas
  endif

  if (sinks_opt_thin .and. insink) planck_ana = min_optical_depth/(0.5D0**nlevelmax*boxlen*scale_l)

end function planck_ana

!##################################################################################################
!##################################################################################################
!##################################################################################################
!##################################################################################################

!  Function ROSSELAND_ANA:
!
!> Compute Rosseland mean opacity.
!<
function rosseland_ana(dens,Tp,Tr,igroup,insink)

  use amr_commons
  use radiation_parameters
  use constants, only:pi

  implicit none

  integer ,intent(in)    :: igroup
  real(dp),intent(in)    :: dens,Tp,Tr
  logical, intent(in)    :: insink
  real(dp)               :: rosseland_ana,Tgd
  real(dp)               :: Tevap ! if sublimation_kuiper
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  
  Tgd = Tp
  if(sublimation_kuiper) then
     Tevap = 2000.0d0*dens**0.0195  !! Evaporation temperature
     if(Tp .gt. Tevap) Tgd = Tevap
  endif  

  if(fit_semenov) then !! No dust grains above Tevap
     rosseland_ana = dens*0.5d0 * ( &
         (1.d0 - dtanh((Tgd - 145.d0)/2.d0)) * 0.00022d0 * Tgd**2 + &
         (1.d0 + dtanh((Tgd - 145.d0)/2.d0)) * &
         (1.51d0 - 1.5d0 * dtanh((Tgd - 1250.d0)/50.d0)))
  else
     rosseland_ana = rosseland_params(1)*(dens**rosseland_params(2))*(Tgd**rosseland_params(3))
  endif

  if(sublimation_kuiper) then !! No dust grains above Tevap
     !# RMR #### Sublimation of dust grains as in kuiper+10 ApJ ####
     !### The highest dust temperature is the evaporation temperature
     !### Opacities are taken at this temperature
     !### The input temperature is kept to compute the d/g ratio
     !# RMR ## Sublimation mimicked by a d/g ratio that decreases as a arctan function centered on Tevap ##
     rosseland_ana = rosseland_ana*(0.5d0 - 1./pi*atan(0.01d0*(Tp - Tevap) ) ) & !ross_ana contains the d/g ratio of 0.01
          +dens*0.01d0*(1.0d0-0.01d0*(0.5d0 - 1./pi*atan(0.01d0*(Tp - Tevap) ) )) !=> quasi full gas
  endif

  if (sinks_opt_thin .and. insink) rosseland_ana = min_optical_depth/(0.5D0**nlevelmax*boxlen*scale_l)

end function rosseland_ana