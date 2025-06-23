module slope_types
   use amr_parameters, only:dp
   use const
   implicit none

contains

   ! HELPER FUNCTIONS

   !#######################################################
   pure function gather_local_values(q,l,i,j,k,n) result(qloc)
      use amr_parameters,   only:dp,nvector,ndim
      use hydro_parameters, only:iu1,iu2,ju1,ju2,ku1,ku2,nvar
      implicit none
      real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),intent(in)::q
      integer,intent(in)::l,i,j,k,n
      real(dp),dimension(0:2*ndim)::qloc
      ! store the center value at index 0
      integer,parameter::icen=0
      ! indices of left/right, bottom/top, back/front cells in q_neighbors array (min and plus)
      integer,parameter::im=1,ip=2,jm=3,jp=4,km=5,kp=6

      ! Gather values at center cell and its neighbors
      qloc(icen) = q(l,i,j,k,n)
      qloc(im)   = q(l,i-1,j,k,n)
      qloc(ip)   = q(l,i+1,j,k,n)
#if NDIM>1
      qloc(jm)   = q(l,i,j-1,k,n)
      qloc(jp)   = q(l,i,j+1,k,n)
#endif
#if NDIM>2
      qloc(km)   = q(l,i,j,k-1,n)
      qloc(kp)   = q(l,i,j,k+1,n)
#endif
   end function gather_local_values
   !#######################################################

   ! SLOPE TYPES

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
   pure function slope_moncen(dlft,drgt) result(slope)
      real(dp),intent(in)::dlft,drgt
      real(dp)::slope
      ! slope_type==2
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