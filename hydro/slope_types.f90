module slope_types
   use amr_parameters, only:dp,nvector,ndim
   use const
   implicit none

contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! USLOPE SUBROUTINE FOR EACH SLOPE TYPE !!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine calc_uslope_minmod_average(q,dq,i,j,k,ngrid,slope_type_real)
      use hydro_parameters
      implicit none
      real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),intent(in)::q
      real(dp),dimension(1:nvector,1:ndim),intent(out)::dq
      integer,intent(in)::ngrid
      integer,intent(in)::i, j, k
      real(dp),intent(in)::slope_type_real
      !-----------------------------------------
      ! 
      !-----------------------------------------
      integer::l
      real(dp)::dlft, drgt, qcen

      !DIR$ IVDEP
      !DIR$ SIMD
      do l = 1, ngrid
         qcen = q(l,i,j,k)

         ! slopes in first coordinate direction
         dlft = qcen - q(l,i-1,j,k)
         drgt = q(l,i+1,j,k) - qcen
#if NDIM==1 || NDIM==2
         dq(l,1) = slope_minmod_or_average(dlft,drgt,slope_type_real)
#else
         dq(l,1) = slope_minmod(dlft,drgt)
#endif
#if NDIM>1
         ! slopes in second coordinate direction
         dlft = qcen - q(l,i,j-1,k)
         drgt = q(l,i,j+1,k) - qcen
#if NDIM==2
         dq(l,2) = slope_minmod_or_average(dlft,drgt,slope_type_real)
#else
         dq(l,2) = slope_minmod(dlft,drgt)
#endif
#endif
#if NDIM>2
         ! slopes in third coordinate direction
         dlft = qcen - q(l,i,j,k-1)
         drgt = q(l,i,j,k+1) - qcen
         dq(l,3) = slope_minmod(dlft,drgt)
#endif
      end do

   end subroutine calc_uslope_minmod_average
   !#######################################################
   subroutine calc_uslope_vanLeer_bis(q,dq,i,j,k,ngrid)
      use hydro_parameters
      implicit none
      real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),intent(in)::q
      real(dp),dimension(1:nvector,1:ndim),intent(out)::dq
      integer,intent(in)::ngrid
      integer,intent(in)::i, j, k
      !-----------------------------------------
      ! Apply van Leer (bis) slope
      !-----------------------------------------
      integer::l
      real(dp)::dlft, drgt, qcen

      !DIR$ IVDEP
      !DIR$ SIMD
      do l = 1, ngrid
         ! Gather values at center cell and its neighbors
         qcen = q(l,i,j,k)

         ! slopes in first coordinate direction
         dlft = qcen - q(l,i-1,j,k)
         drgt = q(l,i+1,j,k) - qcen
         dq(l,1) = slope_vanLeer_bis(dlft,drgt)
#if NDIM>1
         ! slopes in second coordinate direction
         dlft = qcen - q(l,i,j-1,k)
         drgt = q(l,i,j+1,k) - qcen
         dq(l,2) = slope_vanLeer_bis(dlft,drgt)
#endif
#if NDIM>2
         ! slopes in third coordinate direction
         dlft = qcen - q(l,i,j,k-1)
         drgt = q(l,i,j,k+1) - qcen
         dq(l,3) = slope_vanLeer_bis(dlft,drgt)
#endif
      end do

   end subroutine calc_uslope_vanLeer_bis

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! SLOPE TYPE EQUATIONS !!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !#######################################################
   !DIR$ ATTRIBUTES FORCEINLINE :: slope_minmod
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
   !DIR$ ATTRIBUTES FORCEINLINE :: slope_minmod_or_average
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
   !DIR$ ATTRIBUTES FORCEINLINE :: slope_moncen
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
   !DIR$ ATTRIBUTES FORCEINLINE :: slope_vanLeer
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
   !DIR$ ATTRIBUTES FORCEINLINE :: slope_vanLeer_bis
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
