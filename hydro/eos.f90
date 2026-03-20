subroutine barotropic_eos_temperature(nH, temperature)
   use amr_parameters
   !--------------------------------------------------------------
   ! This routine selects the chosen EOS and calculates the
   ! temperature T[in Kelvin]/mu from the density nH[in H/cc]
   !--------------------------------------------------------------
   real(dp), intent(in) ::nH
   real(dp), intent(out)::temperature
   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v

   SELECT CASE (barotropic_eos_form)
   CASE ('isothermal')
      temperature = T2_eos
   CASE ('polytrope')
      temperature = T2_eos*(nH/polytrope_n(1))**(polytrope_index(1)-1.0d0)
   CASE ('double_polytrope')
      ! to convert n to rho: rho = nH/scale_nH*scale_d
      temperature = T2_eos * (1 + (nH/polytrope_n(1))**(polytrope_index(1)-1.0d0))
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
