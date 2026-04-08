subroutine barotropic_eos_temperature(density, temperature)
   use amr_parameters
   implicit none
   !--------------------------------------------------------------
   ! This routine selects the chosen EOS and calculates the
   ! temperature T[in Kelvin]/mu from the density[in H/cc]
   !--------------------------------------------------------------
   real(dp), intent(in) ::density
   real(dp), intent(out)::temperature
   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v

   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

   SELECT CASE (barotropic_eos_form)
   CASE ('isothermal')
      temperature = T2_eos
   CASE ('polytrope')
      temperature = T2_eos*((density/scale_nH)/polytrope_rho_cu)**(polytrope_index-1.0d0)
   CASE ('double_polytrope')
      ! to convert n to rho: rho = density/scale_nH*scale_d
      temperature = T2_eos * (1 + ((density/scale_nH)/polytrope_rho_cu)**(polytrope_index-1.0d0))
   CASE ('custom')
      ! WRITE YOUR FAVORITE EOS HERE
      if((density/scale_nH)<polytrope_rho_cu)then
         temperature = T2_eos
      else
         temperature = T2_eos * ((density/scale_nH)/polytrope_rho_cu)**(polytrope_index-1.0d0)
      endif
   CASE DEFAULT
     write(*,*)'unknown barotropic eos form'
     call clean_stop
   END SELECT

end subroutine barotropic_eos_temperature
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pressure_eos(rho_temp,Enint_temp,Peos)
  use amr_parameters      ,only:dp
  use hydro_commons       ,only:gamma
  implicit none
  !--------------------------------------------------------------
  ! This routine computes the pressure from the density and 
  ! internal volumic energy. Inputs/output are in code units
  !--------------------------------------------------------------
  real(dp), intent(in) :: Enint_temp,rho_temp
  real(dp), intent(out):: Peos

  Peos = (gamma-1.d0)*Enint_temp

  return

end subroutine pressure_eos
!################################################################
!################################################################
!################################################################
!################################################################
subroutine temperature_eos(rho_temp,Enint_temp,Teos,ht)
  use amr_parameters      ,only:dp,mu_gas
  use hydro_commons       ,only:gamma
  use constants           ,only:kB,mH
  implicit none
  !--------------------------------------------------------------
  ! This routine computes the temperature from the density and 
  ! internal volumic energy. Inputs/output are in code units.
  !--------------------------------------------------------------
  real(dp), intent(in) :: Enint_temp,rho_temp
  integer , intent(out):: ht 
  real(dp), intent(out):: Teos
  real(dp)::rho,Enint
  real(dp)::scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  rho   = rho_temp*scale_d
  Enint = Enint_temp*scale_d*scale_v**2 

  Teos = Enint/(rho*kB/(mu_gas*mH*(gamma-1.0d0)))

  ht=1

  return

end subroutine temperature_eos
!################################################################
!################################################################
!################################################################
!################################################################
subroutine enerint_eos(rho_temp,temp_temp,Eeos)
  use amr_parameters      ,only:dp,mu_gas
  use hydro_commons       ,only:gamma
  use constants           ,only:kB,mH
  implicit none
  !--------------------------------------------------------------
  ! This routine computes the internal volumic energy from  
  ! the density and the temperature. Inputs/output are in code units.
  !--------------------------------------------------------------
  real(dp), intent(in) :: temp_temp,rho_temp
  real(dp), intent(out):: Eeos
  real(dp)::rho,temp
  real(dp)::scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  rho  = rho_temp * scale_d
  temp = temp_temp

  Eeos = rho*kB/(mu_gas*mH*(gamma-1.0))*temp/(scale_d*scale_v**2)

  return

end subroutine enerint_eos
!################################################################
!################################################################
!################################################################
!################################################################
subroutine soundspeed_eos(rho_temp,Enint_temp,Cseos)
  use amr_parameters      ,only:dp
  use hydro_commons       ,only:gamma
  implicit none
  !--------------------------------------------------------------
  ! This routine computes the sound speed from the internal volumic energy 
  ! and the temperature. Inputs/output are in code units.
  !--------------------------------------------------------------
  real(dp), intent(in) :: Enint_temp,rho_temp
  real(dp), intent(out):: Cseos

  Cseos = sqrt(gamma*(gamma-1.d0)*Enint_temp/rho_temp)

  return

end subroutine soundspeed_eos
!################################################################
!################################################################
!################################################################
!################################################################
function cmp_Cv_eos(rho,Enint)
  use amr_parameters      ,only:dp,mu_gas
  use hydro_commons       ,only:gamma
  use constants           ,only:kB,mH
  implicit none
  !--------------------------------------------------------------
  ! This function computes the Cv from the density and 
  ! internal volumic energy. Inputs/output are in code units.
  !--------------------------------------------------------------
  real(dp)   :: rho,Enint
  real(dp)   :: cmp_Cv_eos

  real(dp)::scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  
  cmp_Cv_eos = rho*kB/(mu_gas*mH*(gamma-1.0d0))/scale_v**2

end function cmp_Cv_eos
