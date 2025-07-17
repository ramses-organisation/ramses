module slope_types
   use amr_parameters, only:dp,nvector,ndim
   use hydro_parameters, only:nvar,iu1,iu2,ju1,ju2,ku1,ku2
   use const
   implicit none

contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! USLOPE SUBROUTINE FOR EACH SLOPE TYPE !!!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   pure subroutine calc_uslope_minmod_average(q,dq,i,j,k,n,ngrid,slope_type_real)
      implicit none
      real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),intent(in)::q
      real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim),intent(out)::dq
      integer,intent(in)::ngrid
      integer,intent(in)::i, j, k, n
      real(dp),intent(in)::slope_type_real
      !-----------------------------------------------------
      ! Calculate minmod or average slope for each direction
      !-----------------------------------------------------
      integer::l
      real(dp)::dlft, drgt, qcen

      do l = 1, ngrid
         qcen = q(l,i,j,k,n)

         ! slopes in first coordinate direction
         dlft = qcen - q(l,i-1,j,k,n)
         drgt = q(l,i+1,j,k,n) - qcen
#if NDIM==1 || NDIM==2 || defined(SOLVERmhd)
         dq(l,i,j,k,n,1) = slope_minmod_or_average(dlft,drgt,slope_type_real)
#else
         dq(l,i,j,k,n,1) = slope_minmod(dlft,drgt)
#endif
#if NDIM>1
         ! slopes in second coordinate direction
         dlft = qcen - q(l,i,j-1,k,n)
         drgt = q(l,i,j+1,k,n) - qcen
#if NDIM==2 || defined(SOLVERmhd)
         dq(l,i,j,k,n,2) = slope_minmod_or_average(dlft,drgt,slope_type_real)
#else
         dq(l,i,j,k,n,2) = slope_minmod(dlft,drgt)
#endif
#endif
#if NDIM>2
         ! slopes in third coordinate direction
         dlft = qcen - q(l,i,j,k-1,n)
         drgt = q(l,i,j,k+1,n) - qcen
#if defined(SOLVERmhd)
         dq(l,i,j,k,n,3) = slope_minmod_or_average(dlft,drgt,slope_type_real)
#else
         dq(l,i,j,k,n,3) = slope_minmod(dlft,drgt)
#endif
#endif
      end do

   end subroutine calc_uslope_minmod_average
   !#######################################################
#if NDIM==3 && defined(SOLVERhydro)
   pure subroutine calc_uslope_moncen(q,dq,i,j,k,n,ngrid)
      implicit none
      real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),intent(in)::q
      real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim),intent(out)::dq
      integer,intent(in)::ngrid
      integer,intent(in)::i, j, k, n
      !------------------------------------------
      ! Calculate moncen slope for each direction
      !------------------------------------------
      integer::l
      real(dp)::dlft, drgt, qcen

      do l = 1, ngrid
         ! Gather values at center cell and its neighbors
         qcen = q(l,i,j,k,n)

         ! slopes in first coordinate direction
         dlft = qcen - q(l,i-1,j,k,n)
         drgt = q(l,i+1,j,k,n) - qcen
         dq(l,i,j,k,n,1) = slope_moncen(dlft,drgt)

         ! slopes in second coordinate direction
         dlft = qcen - q(l,i,j-1,k,n)
         drgt = q(l,i,j+1,k,n) - qcen
         dq(l,i,j,k,n,2) = slope_moncen(dlft,drgt)

         ! slopes in third coordinate direction
         dlft = qcen - q(l,i,j,k-1,n)
         drgt = q(l,i,j,k+1,n) - qcen
         dq(l,i,j,k,n,3) = slope_moncen(dlft,drgt)
      end do

   end subroutine calc_uslope_moncen
#endif
   !#######################################################
#if NDIM==2
   pure subroutine calc_uslope_positivity_preserving(q,dq,i,j,k,n,ngrid)
      implicit none
      real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),intent(in)::q
      real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim),intent(out)::dq
      integer,intent(in)::ngrid
      integer,intent(in)::i, j, k, n
      !--------------------------------------------------------------------------
      ! Calculate positivity preserving 2D unsplit slope slope for each direction
      !--------------------------------------------------------------------------
      integer::l
      real(dp)::dlft,drgt,qcen,dlim
      real(dp)::dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr
      real(dp)::vmin,vmax,dfx,dfy,dff

      do l = 1, ngrid
         ! Gather values at center cell and its neighbors
         qcen = q(l,i,j,k,n)

         dfll = q(l,i-1,j-1,k,n) - qcen
         dflm = q(l,i-1,j  ,k,n) - qcen
         dflr = q(l,i-1,j+1,k,n) - qcen
         dfml = q(l,i  ,j-1,k,n) - qcen
         dfmm = 0
         dfmr = q(l,i  ,j+1,k,n) - qcen
         dfrl = q(l,i+1,j-1,k,n) - qcen
         dfrm = q(l,i+1,j  ,k,n) - qcen
         dfrr = q(l,i+1,j+1,k,n) - qcen

         vmin = min(dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr)
         vmax = max(dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr)

         dfx  = half*(q(l,i+1,j,k,n)-q(l,i-1,j,k,n))
         dfy  = half*(q(l,i,j+1,k,n)-q(l,i,j-1,k,n))
         dff  = half*(abs(dfx)+abs(dfy))

         if(dff>zero)then
            dlim = min(one,min(abs(vmin),abs(vmax))/dff)
         else
            dlim = one
         endif

         dq(l,i,j,k,n,1) = dlim*dfx
         dq(l,i,j,k,n,2) = dlim*dfy
      end do

   end subroutine calc_uslope_positivity_preserving
#elif NDIM==3
   pure subroutine calc_uslope_positivity_preserving(q,dq,i,j,k,n,ngrid)
      implicit none
      real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),intent(in)::q
      real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim),intent(out)::dq
      integer,intent(in)::ngrid
      integer,intent(in)::i, j, k, n
      !--------------------------------------------------------------------------
      ! Calculate positivity preserving 3D unsplit slope slope for each direction
      !--------------------------------------------------------------------------
      integer::l
      real(dp)::dlft,drgt,qcen,dlim
      real(dp)::dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl
      real(dp)::dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm
      real(dp)::dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr
      real(dp)::vmin,vmax,dfx,dfy,dfz,dff

      do l = 1, ngrid
         qcen = q(l,i,j,k,n)

         dflll = q(l,i-1,j-1,k-1,n) - qcen
         dflml = q(l,i-1,j  ,k-1,n) - qcen
         dflrl = q(l,i-1,j+1,k-1,n) - qcen
         dfmll = q(l,i  ,j-1,k-1,n) - qcen
         dfmml = q(l,i  ,j  ,k-1,n) - qcen
         dfmrl = q(l,i  ,j+1,k-1,n) - qcen
         dfrll = q(l,i+1,j-1,k-1,n) - qcen
         dfrml = q(l,i+1,j  ,k-1,n) - qcen
         dfrrl = q(l,i+1,j+1,k-1,n) - qcen

         dfllm = q(l,i-1,j-1,k  ,n) - qcen
         dflmm = q(l,i-1,j  ,k  ,n) - qcen
         dflrm = q(l,i-1,j+1,k  ,n) - qcen
         dfmlm = q(l,i  ,j-1,k  ,n) - qcen
         dfmmm = 0
         dfmrm = q(l,i  ,j+1,k  ,n) - qcen
         dfrlm = q(l,i+1,j-1,k  ,n) - qcen
         dfrmm = q(l,i+1,j  ,k  ,n) - qcen
         dfrrm = q(l,i+1,j+1,k  ,n) - qcen

         dfllr = q(l,i-1,j-1,k+1,n) - qcen
         dflmr = q(l,i-1,j  ,k+1,n) - qcen
         dflrr = q(l,i-1,j+1,k+1,n) - qcen
         dfmlr = q(l,i  ,j-1,k+1,n) - qcen
         dfmmr = q(l,i  ,j  ,k+1,n) - qcen
         dfmrr = q(l,i  ,j+1,k+1,n) - qcen
         dfrlr = q(l,i+1,j-1,k+1,n) - qcen
         dfrmr = q(l,i+1,j  ,k+1,n) - qcen
         dfrrr = q(l,i+1,j+1,k+1,n) - qcen

         vmin = min(dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl, &
               &     dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm, &
               &     dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr)
         vmax = max(dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl, &
               &     dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm, &
               &     dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr)

         dfx  = half*(q(l,i+1,j,k,n) - q(l,i-1,j,k,n))
         dfy  = half*(q(l,i,j+1,k,n) - q(l,i,j-1,k,n))
         dfz  = half*(q(l,i,j,k+1,n) - q(l,i,j,k-1,n))
         dff  = half*(abs(dfx)+abs(dfy)+abs(dfz))

         if(dff>zero)then
            dlim = min(one,min(abs(vmin),abs(vmax))/dff)
         else
            dlim = one
         endif

         dq(l,i,j,k,n,1) = dlim*dfx
         dq(l,i,j,k,n,2) = dlim*dfy
         dq(l,i,j,k,n,3) = dlim*dfz
      end do

   end subroutine calc_uslope_positivity_preserving
#endif
   !#######################################################
#if NDIM==1
   pure subroutine calc_uslope_superbee(q,dq,i,j,k,n,ngrid,dtdx)
      implicit none
      real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),intent(in)::q
      real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim),intent(out)::dq
      integer,intent(in)::ngrid
      integer,intent(in)::i, j, k, n
      real(dp),intent(in)::dtdx
      !--------------------------------------------
      ! Calculate superbee slope for each direction
      !--------------------------------------------
      integer::l
      real(dp)::dlft, drgt, qcen,dcen,dsgn,dlim

      do l = 1, ngrid
         ! Gather values at center cell and its neighbors
         qcen = q(l,i,j,k,n)
         dcen = q(l,i,j,k,2)*dtdx
         ! slopes in first coordinate direction
         dlft = two/(one+dcen)*(qcen - q(l,i-1,j,k,n))
         drgt = two/(one-dcen)*(q(l,i+1,j,k,n) - qcen)
         !dcen = half*(q(l,i+1,j,k,n)-q(l,i-1,j,k,n))
         dsgn = sign(one, dlft)
         dlim = min(abs(dlft),abs(drgt))
         if((dlft*drgt)<=zero)dlim=zero
         dq(l,i,j,k,n,1) = dsgn*dlim !min(dlim,abs(dcen))
      end do

   end subroutine calc_uslope_superbee
   !#######################################################
   pure subroutine calc_uslope_ultrabee(q,dq,i,j,k,n,ngrid,dtdx)
      implicit none
      real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),intent(in)::q
      real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim),intent(out)::dq
      integer,intent(in)::ngrid
      integer,intent(in)::i, j, k, n
      real(dp),intent(in)::dtdx
      !--------------------------------------------
      ! Calculate ultrabee slope for each direction
      !--------------------------------------------
      integer::l
      real(dp)::dlft, drgt, qcen,dcen,dsgn,dlim

      do l = 1, ngrid
         ! Gather values at center cell and its neighbors
         qcen = q(l,i,j,k,n)
         dcen = q(l,i,j,k,2)*dtdx
         if(dcen>=0)then
            dlft = two/(zero+dcen+1d-10)*(qcen - q(l,i-1,j,k,n))
            drgt = two/(one -dcen      )*(q(l,i+1,j,k,n) - qcen)
         else
            dlft = two/(one +dcen      )*(qcen - q(l,i-1,j,k,n))
            drgt = two/(zero-dcen+1d-10)*(q(l,i+1,j,k,n) - qcen)
         endif
         dsgn = sign(one, dlft)
         dlim = min(abs(dlft),abs(drgt))
         !dcen = half*(q(l,i+1,j,k,n)-q(l,i-1,j,k,n))
         if((dlft*drgt)<=zero)dlim=zero
         dq(l,i,j,k,n,1) = dsgn*dlim !min(dlim,abs(dcen))
      end do

   end subroutine calc_uslope_ultrabee
   !#######################################################
   pure subroutine calc_uslope_unstable(q,dq,i,j,k,n,ngrid)
      implicit none
      real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),intent(in)::q
      real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim),intent(out)::dq
      integer,intent(in)::ngrid
      integer,intent(in)::i, j, k, n
      !--------------------------------------------
      ! Calculate unstable slope for each direction
      !--------------------------------------------
      integer::l
      real(dp)::dlft, drgt, qcen

      do l = 1, ngrid
         ! Gather values at center cell and its neighbors
         qcen = q(l,i,j,k,n)
         ! slopes in first coordinate direction
         dlft = qcen - q(l,i-1,j,k,n)
         drgt = q(l,i+1,j,k,n) - qcen
         dq(l,i,j,k,n,1) = 0.5d0*(dlft+drgt)
      end do

   end subroutine calc_uslope_unstable
#endif
   !#######################################################
   pure subroutine calc_uslope_vanLeer(q,dq,i,j,k,n,ngrid)
      implicit none
      real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),intent(in)::q
      real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim),intent(out)::dq
      integer,intent(in)::ngrid
      integer,intent(in)::i, j, k, n
      !--------------------------------------------
      ! Calculate van Leer slope for each direction
      !--------------------------------------------
      integer::l
      real(dp)::dlft, drgt, qcen

      do l = 1, ngrid
         ! Gather values at center cell and its neighbors
         qcen = q(l,i,j,k,n)

         ! slopes in first coordinate direction
         dlft = qcen - q(l,i-1,j,k,n)
         drgt = q(l,i+1,j,k,n) - qcen
         dq(l,i,j,k,n,1) = slope_vanLeer(dlft,drgt)
#if NDIM>1
         ! slopes in second coordinate direction
         dlft = qcen - q(l,i,j-1,k,n)
         drgt = q(l,i,j+1,k,n) - qcen
         dq(l,i,j,k,n,2) = slope_vanLeer(dlft,drgt)
#endif
#if NDIM>2
         ! slopes in third coordinate direction
         dlft = qcen - q(l,i,j,k-1,n)
         drgt = q(l,i,j,k+1,n) - qcen
         dq(l,i,j,k,n,3) = slope_vanLeer(dlft,drgt)
#endif
      end do

   end subroutine calc_uslope_vanLeer
   !#######################################################
   pure subroutine calc_uslope_vanLeer_bis(q,dq,i,j,k,n,ngrid)
      implicit none
      real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),intent(in)::q
      real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim),intent(out)::dq
      integer,intent(in)::ngrid
      integer,intent(in)::i, j, k, n
      !--------------------------------------------------
      ! Calculate van Leer (bis) slope for each direction
      !--------------------------------------------------
      integer::l
      real(dp)::dlft, drgt, qcen

      do l = 1, ngrid
         ! Gather values at center cell and its neighbors
         qcen = q(l,i,j,k,n)

         ! slopes in first coordinate direction
         dlft = qcen - q(l,i-1,j,k,n)
         drgt = q(l,i+1,j,k,n) - qcen
         dq(l,i,j,k,n,1) = slope_vanLeer_bis(dlft,drgt)
#if NDIM>1
         ! slopes in second coordinate direction
         dlft = qcen - q(l,i,j-1,k,n)
         drgt = q(l,i,j+1,k,n) - qcen
         dq(l,i,j,k,n,2) = slope_vanLeer_bis(dlft,drgt)
#endif
#if NDIM>2
         ! slopes in third coordinate direction
         dlft = qcen - q(l,i,j,k-1,n)
         drgt = q(l,i,j,k+1,n) - qcen
         dq(l,i,j,k,n,3) = slope_vanLeer_bis(dlft,drgt)
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
