module slope_types
   use amr_parameters, only:dp
   use const,  only:zero
   implicit none
   public::slope_minmod

contains

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
 

end module slope_types