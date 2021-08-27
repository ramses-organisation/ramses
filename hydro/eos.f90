subroutine eos_temperature_from_density(nH, temperature)
   use amr_parameters
   !--------------------------------------------------------------
   ! This routine selects the chosen EOS and calculates the
   ! temperature T[in Kelvin]/mu from the density nH[in H/cc]
   !--------------------------------------------------------------
   real(dp), intent(in) ::nH
   real(dp), intent(out)::temperature
   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v

   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

   SELECT CASE (eos_form)
   CASE ('isothermal')
      temperature = T2_eos
   CASE ('polytrop')
      temperature = T2_eos*(nH/n_star)**(g_star-1.0d0)
   CASE ('barotrop')
      ! to convert n to rho: rho = nH/scale_nH*scale_d
      temperature = T2_eos * (1 + ((nH/scale_nH*scale_d)/barotrop_knee)**(barotrop_slope-1.0d0))
   CASE ('custom')
      ! WRITE YOUR FAVORITE EOS HERE
      if(nH<barotrop_knee)then
         temperature = T2_eos
      else
         temperature = T2_eos * ((nH/scale_nH*scale_d)/barotrop_knee)**(barotrop_slope-1.0d0)
      endif
   CASE DEFAULT
     write(*,*)'unknown eos form'
     call clean_stop
   END SELECT

end subroutine eos_temperature_from_density