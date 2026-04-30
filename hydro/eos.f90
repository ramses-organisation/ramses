subroutine barotropic_eos_temperature(nH, temperature)
   use amr_parameters
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
