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
