subroutine barotropic_eos_temperature(density, temperature)
   use amr_parameters
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

!=====================================================================================================
!=====================================================================================================
!=====================================================================================================
!=====================================================================================================
subroutine pressure_eos(rho_temp,Enint_temp,Peos)
   use amr_parameters      ,only:dp,barotropic_eos
   use hydro_commons       ,only:gamma
   implicit none
   !--------------------------------------------------------------
   ! This routine computes the pressure from the density and 
   ! internal volumic energy. Inputs/output are in code units
   !--------------------------------------------------------------
   real(dp), intent(in) :: Enint_temp,rho_temp
   real(dp), intent(out):: Peos
   real(dp) :: Cseos
 
   if(barotropic_eos)then
       call soundspeed_eos(rho_temp,Enint_temp,Cseos)
       !Peos = rho_temp * Cseos**2 / gamma
       ! without gamma, otherwise jeans_refine is wrong!
       Peos = rho_temp * Cseos**2
   else
       Peos = (gamma-1d0)*Enint_temp
   endif
 
   return
 
 end subroutine pressure_eos
 !===========================================================================================
 !===========================================================================================
 !===========================================================================================
 !===========================================================================================
 subroutine temperature_eos(rho_temp,Enint_temp,Teos)
   use amr_parameters      ,only:dp,mu_gas,barotropic_eos
   use hydro_commons       ,only:gamma
   use constants
   implicit none
   !--------------------------------------------------------------
   ! This routine computes the temperature from the density and 
   ! internal volumic energy. Inputs/output are in code units.
   !--------------------------------------------------------------
   real(dp), intent(in) :: Enint_temp,rho_temp
   real(dp), intent(out):: Teos
   real(dp)::rho,Enint,temperature
 
   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
 
   if(barotropic_eos)then
      call barotropic_eos_temperature(rho_temp*scale_nH, temperature)
      Teos=temperature*mu_gas
   else
      rho   = rho_temp*scale_d
      Enint = Enint_temp*scale_d*scale_v**2
      Teos = Enint/(rho*kB/(mu_gas*mH*(gamma-1.0d0)))
   end if

   return
 
 end subroutine temperature_eos
 !==================================================================================
 !==================================================================================
 !==================================================================================
 !==================================================================================
 subroutine soundspeed_eos(rho_temp,Enint_temp,Cseos)
   use amr_parameters      ,only:dp,mu_gas,barotropic_eos
   use hydro_commons       ,only:gamma
   use constants
   implicit none
   !--------------------------------------------------------------
   ! This routine computes the sound speed from the internal volumic energy 
   ! and the temperature. Inputs/output are in code units.
   !--------------------------------------------------------------
   real(dp), intent(in) :: Enint_temp,rho_temp
   real(dp), intent(out):: Cseos
   real(dp):: Teos
    real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
    call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
 
   if(barotropic_eos)then
       call temperature_eos(rho_temp,Enint_temp,Teos)
       Cseos = sqrt(kB*Teos/(mu_gas*mH))/scale_v
       !Cseos = sqrt(gamma*kB*Teos/(mu_gas*mH))/scale_v
       !mu_gas cancels out
       ! without gamma, otherwise jeans_refine is wrong!
   else
       Cseos = sqrt(gamma*(gamma-1d0)*Enint_temp/rho_temp)
   endif
 
   return
 
 end subroutine soundspeed_eos