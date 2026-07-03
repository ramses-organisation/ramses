subroutine barotropic_eos_temperature(nH, temperature)
   use amr_parameters
   implicit none
   !--------------------------------------------------------------
   ! This routine selects the chosen EOS and calculates the
   ! temperature T[in Kelvin]/mu from the density nH[in H/cc]
   !--------------------------------------------------------------
   real(dp), intent(in) ::nH
   real(dp), intent(out)::temperature
   real(dp)::factor1,factor2,factor3

   SELECT CASE (barotropic_eos_form)
   CASE ('isothermal')
      temperature = T2_eos
   CASE ('polytrope')
      temperature = T2_eos*(nH/polytrope_n(1))**(polytrope_index(1)-1.0d0)
   CASE ('double_polytrope')
      temperature = T2_eos * (1 + (nH/polytrope_n(1))**(polytrope_index(1)-1.0d0))
   CASE ('2nd_collapse')
      ! Machida & Inutsuka 2006, Marchand et al. 2016
      factor1 = sqrt(1 + (nH/polytrope_n(1))**(2*polytrope_index(1)))
      factor2 = (1 + (nH/polytrope_n(2)))**polytrope_index(2)
      factor3 = (1 + (nH/polytrope_n(3)))**polytrope_index(3)
      temperature = T2_eos * factor1 * factor2 * factor3
   CASE ('custom')
      ! WRITE YOUR FAVORITE EOS HERE
      if(nH<polytrope_n(1))then
         temperature = T2_eos
      else
         temperature = T2_eos * (nH/polytrope_n(1))**(polytrope_index(1)-1.0d0)
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
subroutine temperature_eos(rho_temp,Enint_temp,Teos)
  use amr_parameters      ,only:dp,mu_gas
  use hydro_commons       !,only:gamma
  use constants           ,only:kB,mH
  use cooling_module, only:barotropic_eos
  implicit none
  !--------------------------------------------------------------
  ! This routine computes the temperature from the density and
  ! internal volumic energy. Inputs/output are in code units.
  !--------------------------------------------------------------
  real(dp), intent(in) :: Enint_temp,rho_temp
  real(dp), intent(out):: Teos
  real(dp)::rho,Enint
  real(dp)::scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l

  integer::i_t,i_r,i
  real(dp)::logr,tt,uu,y1,y2,y3,y4
  real(dp):: le,lr
  real(dp):: dd1,dd2,de1,de2
  integer :: ir,ie
  real(dp):: xx,drho,dener

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  rho   = rho_temp*scale_d
  Enint = Enint_temp*scale_d*scale_v**2
  if (barotropic_eos) then
      call barotropic_eos_temperature(rho_temp, Teos)
      Teos = Teos*mu_gas
  else
      Teos = Enint/(rho*kB/(mu_gas*mH*(gamma-1.0d0)))
  endif

end subroutine temperature_eos
!################################################################
!################################################################
!################################################################
!################################################################
