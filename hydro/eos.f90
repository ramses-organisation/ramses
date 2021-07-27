subroutine get_eos(rho, temperature)
   use amr_parameters      ,only:dp,T2_star,n_star,g_star,mu_gas
   use hydro_commons       ,only:gamma,eos_form
   real(dp)::rho, temperature
   real(dp)::barotrop1D
   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v

   !--------------------------------------------------------------
   ! This routine selects the chosen EOS
   ! Inputs/output are in H/cc units
   !--------------------------------------------------------------
   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

   SELECT CASE (eos_form)
   CASE ('isothermal')
      temperature = T2_star
   CASE ('polytrop')
      temperature = T2_star*(rho/n_star)**(g_star-1.0d0)
   CASE ('barotrop')
      temperature = barotrop1D(rho/scale_nH*scale_d)/mu_gas
   CASE DEFAULT
     write(*,*)'unknown eos form'
     call clean_stop
   END SELECT

end subroutine get_eos
!###########################################################
!###########################################################
!###########################################################
!###########################################################

!subroutine read_eos_table()
!end subroutine read_eos_table

!subroutine get_eos_value_from_table()
!end subroutine get_eos_value_from_table

!###########################################################
!###########################################################
!###########################################################
!###########################################################
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
  real(dp) :: Cseos

  !!if(barotrop)then
      call soundspeed_eos(rho_temp,Enint_temp,Cseos)
      !Peos = rho_temp * Cseos**2 / gamma
      ! without gamma, otherwise jeans_refine is wrong!
      Peos = rho_temp * Cseos**2
  !!else
  !!    Peos = (gamma-1d0)*Enint_temp
  !!endif
  !Peos = (gamma-1d0)*Enint_temp

  return

end subroutine pressure_eos
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine temperature_eos(rho_temp,Enint_temp,Teos)
  use amr_parameters      ,only:dp,mu_gas
  use hydro_commons       ,only:gamma
  use constants
  implicit none
  !--------------------------------------------------------------
  ! This routine computes the temperature from the density and 
  ! internal volumic energy. Inputs/output are in code units.
  !--------------------------------------------------------------
  real(dp), intent(in) :: Enint_temp,rho_temp
  real(dp), intent(out):: Teos
  real(dp)::rho,Enint,barotrop1D

   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  !!if(barotrop)then
     Teos=barotrop1D(rho_temp*scale_d)
  !!else
  !!   rho   = rho_temp*scale_d
  !!   Enint = Enint_temp*scale_d*scale_v**2

  !!   Teos = Enint/(rho*kB/(mu_gas*mH*(gamma-1.0d0)))

  !!end if
  !rho   = rho_temp*scale_d
  !Enint = Enint_temp*scale_d*scale_v**2 
  !
  !Teos = Enint/(rho*kB/(mu_gas*mH*(gamma-1.0d0)))
  !

  return

end subroutine temperature_eos
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine soundspeed_eos(rho_temp,Enint_temp,Cseos)
  use amr_parameters      ,only:dp,mu_gas
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

  !!if(barotrop)then
      call temperature_eos(rho_temp,Enint_temp,Teos)
      Cseos = sqrt(kB*Teos/(mu_gas*mH))/scale_v
      !Cseos = sqrt(gamma*kB*Teos/(mu_gas*mH))/scale_v
      !mu_gas cancels out
      ! without gamma, otherwise jeans_refine is wrong!
  !!else
  !!    Cseos = sqrt(gamma*(gamma-1d0)*Enint_temp/rho_temp)
  !!endif
  !Cseos = sqrt(gamma*(gamma-1d0)*Enint_temp/rho_temp)

  return

end subroutine soundspeed_eos
!###########################################################
!###########################################################
!###########################################################
!###########################################################
!double precision function barotrop1D(rhon)
function barotrop1D(rhon)
   use hydro_commons
   use amr_parameters!, only : n_star, T2_star, mu_gas
   implicit none
   !--------------------------------------------------------------
   ! Barotropic equation of state
   !--------------------------------------------------------------
   real(dp)::barotrop1D
   real(dp)::rhon,nH
   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

   ! Implement analytical EOS here
   ! rhon is the density in g/cm3
   !if(rhon<barotrop_knee)then
   !   barotrop1D = mu_gas * T2_star
   !else
   !   barotrop1D = mu_gas * T2_star * (rhon/barotrop_knee)**(barotrop_slope-1.0d0)
   !endif
   barotrop1D = mu_gas * T2_star * (1 + (rhon/barotrop_knee)**(barotrop_slope-1.0d0))

end function barotrop1D