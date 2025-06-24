module slope_types
   use amr_parameters, only:dp
   use const
   implicit none

contains

   !#######################################################
   pure function slope_minmod(dlft,drgt) result(slope)
      real(dp),intent(in)::dlft,drgt
      real(dp)::slope
      ! slope_type==1

      if((dlft*drgt)<=zero) then
         slope = zero
      else if(dlft>0) then
         slope = min(dlft,drgt)
      else
         slope = max(dlft,drgt)
      end if

   end function slope_minmod
   !#######################################################
   pure function slope_minmod_or_average(dlft,drgt,slope_type) result(slope)
      real(dp),intent(in)::dlft,drgt,slope_type
      real(dp)::slope
      real(dp)::dcen,dsgn,dlim

      dcen = half*(dlft+drgt)
      dsgn = sign(one, dcen)
      dlim = min(slope_type,two)*min(abs(dlft),abs(drgt))
      if((dlft*drgt)<=zero)dlim=zero
      slope = dsgn*min(dlim,abs(dcen))

   end function slope_minmod_or_average
   !#######################################################
   pure function slope_moncen(dlft,drgt) result(slope)
      real(dp),intent(in)::dlft,drgt
      real(dp)::slope
      ! slope_type==2
      real(dp)::dcen,dsgn,dlim

      dcen = half*(dlft+drgt)
      ! TC: what's the point of this?
      !     half and slopetype=2 are just going to cancel each other
      !     even if slopetype is an int
      dsgn = sign(one, dcen)
      dlim = 2*min(abs(dlft),abs(drgt))
      if((dlft*drgt)<=zero)dlim=zero
      slope = dsgn*min(dlim,abs(dcen))

   end function slope_moncen
   !#######################################################
   pure function slope_vanLeer(dlft,drgt) result(slope)
      real(dp),intent(in)::dlft,drgt
      real(dp)::slope
      ! slope_type==7

      if((dlft*drgt)<=zero) then
         slope=zero
      else
         slope=(2*dlft*drgt/(dlft+drgt))
      end if

   end function slope_vanLeer
   !#######################################################
   pure function slope_vanLeer_bis(dlft,drgt) result(slope)
      real(dp),intent(in)::dlft,drgt
      real(dp)::slope
      ! slope_type==8
      ! generalized moncen/minmod parameterisation (van Leer 1979)
      real(dp)::dcen,dsgn,dlim
      real(dp),parameter::slope_theta=1.5d0

      dcen = half*(dlft+drgt)
      dsgn = sign(one, dcen)
      dlim = min(slope_theta*abs(dlft),slope_theta*abs(drgt))
      if((dlft*drgt)<=zero)dlim=zero
      slope = dsgn*min(dlim,abs(dcen))

   end function slope_vanLeer_bis
   !#######################################################

end module slope_types
