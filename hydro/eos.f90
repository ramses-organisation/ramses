subroutine barotropic_eos_temperature(nH, temperature)
   use amr_parameters
   !--------------------------------------------------------------
   ! This routine selects the chosen EOS and calculates the
   ! temperature T[in Kelvin]/mu from the density nH[in H/cc]
   !--------------------------------------------------------------
   real(dp), intent(in) ::nH
   real(dp), intent(out)::temperature
   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v

   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

   SELECT CASE (barotropic_eos_form)
   CASE ('legacy')
      ! remark: not exactly the same for cosmo=.true. since n_star!=nISM in all cases
      temperature = T2_star*(nH/n_star)**(g_star-1.0d0)
   CASE ('isothermal')
      temperature = T2_eos
   CASE ('polytrope')
      temperature = T2_eos*((nH/scale_nH*scale_d)/polytrope_rho)**(polytrope_index-1.0d0)
   CASE ('double_polytrope')
      ! to convert n to rho: rho = nH/scale_nH*scale_d
      temperature = T2_eos * (1 + ((nH/scale_nH*scale_d)/polytrope_rho)**(polytrope_index-1.0d0))
   CASE ('custom')
      ! WRITE YOUR FAVORITE EOS HERE
      if(nH<polytrope_rho)then
         temperature = T2_eos
      else
         temperature = T2_eos * ((nH/scale_nH*scale_d)/polytrope_rho)**(polytrope_index-1.0d0)
      endif
   CASE DEFAULT
     write(*,*)'unknown barotropic eos form'
     call clean_stop
   END SELECT

end subroutine barotropic_eos_temperature