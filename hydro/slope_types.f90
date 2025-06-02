module slope_types
   use amr_parameters, only:dp
   use const
   implicit none

contains

   !#######################################################
   pure function slope_minmod(dlft,drgt) result(slope)
      real(dp),intent(in)::dlft,drgt
      real(dp)::slope
 
      if((dlft*drgt)<=zero) then
         slope = zero
      else if(dlft>0) then
         slope = min(dlft,drgt)
      else
         slope = max(dlft,drgt)
      end if
 
   end function slope_minmod
   !#######################################################
   pure function slope_moncen(dlft,drgt) result(slope)
      real(dp),intent(in)::dlft,drgt
      real(dp)::slope
      real(dp)::dcen,dsgn,dlim
      integer,parameter::slope_type=2
 
      dcen = half*(dlft+drgt)/slope_type
      ! TC: what's the point of this? 
      !     half and slopetype=2 are just going to cancel each other
      !     even if slopetype is an int
      dsgn = sign(one, dcen)
      dlim = min(abs(dlft),abs(drgt))
      if((dlft*drgt)<=zero)dlim=zero
      slope = dsgn*min(dlim,abs(dcen))
 
   end function slope_moncen
   !#######################################################
   pure function slope_vanLeer(dlft,drgt) result(slope)
      real(dp),intent(in)::dlft,drgt
      real(dp)::slope
 
      if((dlft*drgt)<=zero) then
         slope=zero
      else
         slope=(2*dlft*drgt/(dlft+drgt))
      end if
 
   end function slope_vanLeer
   !#######################################################


end module slope_types