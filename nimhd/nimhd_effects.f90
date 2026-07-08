! by Jacques Masson, Benoit Commercon and Neil Vaytet
! refactored by Tine Colman
!
! This file contains the core routines for the non-ideal MHD (NIMHD)
! effects: ambipolar diffusion and Ohmic dissipation. They are called
! from the MHD Godunov update (mhd/umuscl.f90) and add the non-ideal
! contributions to the electromotive fields (EMFs) that evolve the
! magnetic field, together with the associated non-ideal energy fluxes.
!
! Overview of the routines in this file:
!   - compute_bemf / compute_bmagij / compute_bmagijbis : interpolate the
!     magnetic field to the various staggered locations needed below.
!   - computejb2   : build the current density J at the EMF/face locations
!                    and the ideal-MHD flux building blocks (fluxmd, fluxad)
!                    shared by both effects.
!   - computdifmag : Ohmic dissipation EMF (emfohmdiss) and energy flux,
!                    using the fixed resistivity eta = etaMD.
!   - computambip  : ambipolar diffusion EMF (emfambdiff) and energy flux,
!                    using the fixed coefficient beta = 1/(gammaAD*rho).
!   - helpers      : finite-difference derivatives (computd{x,y,z}bis)
!                    and cross products (crossprod*).

!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine compute_bemf(u,q,ngrid,bemfx,bemfy,bemfz)

   USE amr_parameters
   use hydro_commons
   IMPLICIT NONE
   !-------------------------------------------
   ! compute magnetic field at location of EMF
   !-------------------------------------------
   ! inputs
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::u
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q
   integer::ngrid
   ! output
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::bemfx,bemfy,bemfz
   ! local variables
   integer ::i, j, k, l

   bemfx=0d0
   bemfy=0d0
   bemfz=0d0

   !!!!!!!!!!!!!!!!!!
   ! EMF x
   !!!!!!!!!!!!!!!!!!

   do k=min(1,ku1+1),ku2
      do j=min(1,ju1+1),ju2
         do i=iu1,iu2
            do l=1,ngrid
               bemfx(l,i,j,k,1)=0.25d0*( q(l,i,j,k,6)+q(l,i,j-1,k,6)+q(l,i,j,k-1,6)+q(l,i,j-1,k-1,6) )
            end do
         end do
      end do
   end do

   do k=min(1,ku1+1),ku2
      do j=ju1,ju2
         do i=iu1,iu2
            do l=1,ngrid
               bemfx(l,i,j,k,2)=0.5d0*( u(l,i,j,k,7)+u(l,i,j,k-1,7) )
            end do
         end do
      end do
   end do

   do k=ku1,ku2
      do j=min(1,ju1+1),ju2
         do i=iu1,iu2
            do l=1,ngrid
               bemfx(l,i,j,k,3)=0.5d0*(u(l,i,j,k,8)+u(l,i,j-1,k,8))
            end do
         end do
      end do
   end do

   !!!!!!!!!!!!!!!!!!
   ! EMF y
   !!!!!!!!!!!!!!!!!!

   do k=min(1,ku1+1),ku2
      do j=ju1,ju2
         do i=iu1,iu2
            do l=1,ngrid
               bemfy(l,i,j,k,1)=0.5d0*(u(l,i,j,k,6)+u(l,i,j,k-1,6))
            end do
         end do
      end do
   end do

   do k=min(1,ku1+1),ku2
      do j=ju1,ju2
         do i=min(1,iu1+1),iu2
            do l=1,ngrid
               bemfy(l,i,j,k,2)=0.25d0*(q(l,i,j,k,7)+q(l,i-1,j,k,7)+q(l,i,j,k-1,7)+q(l,i-1,j,k-1,7))
            end do
         end do
      end do
   end do

   do k=ku1,ku2
      do j=ju1,ju2
         do i=min(1,iu1+1),iu2
            do l=1,ngrid
               bemfy(l,i,j,k,3)=0.5d0*(u(l,i-1,j,k,8)+u(l,i,j,k,8))
            end do
         end do
      end do
   end do

   !!!!!!!!!!!!!!!!!!
   ! EMF z
   !!!!!!!!!!!!!!!!!!

   do k=ku1,ku2
      do j=min(1,ju1+1),ju2
         do i=iu1,iu2
            do l=1,ngrid
               bemfz(l,i,j,k,1)=0.5d0*(u(l,i,j,k,6)+u(l,i,j-1,k,6))
            end do
         end do
      end do
   end do

   do k=ku1,ku2
      do j=ju1,ju2
         do i=min(1,iu1+1),iu2
            do l=1,ngrid
               bemfz(l,i,j,k,2)=0.5d0*(u(l,i,j,k,7)+u(l,i-1,j,k,7))
            end do
         end do
      end do
   end do

   do k=ku1,ku2
      do j=min(1,ju1+1),ju2
         do i=min(1,iu1+1),iu2
            do l=1,ngrid
               bemfz(l,i,j,k,3)=0.25d0*(q(l,i,j,k,8)+q(l,i-1,j,k,8)+q(l,i,j-1,k,8)+q(l,i-1,j-1,k,8))
            end do
         end do
      end do
   end do

end subroutine compute_bemf
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine compute_bmagij(u,q,ngrid,bmagij)
   USE amr_parameters
   use hydro_commons
   IMPLICIT NONE
   !-----------------------------------------------------------------
   ! bmagij is the value of the magnetic field Bi where Bj
   ! is naturally defined; Ex bmagij(l,i,j,k,1,2) is Bx at i,j-1/2,k
   ! and we can write it Bx,y
   !-----------------------------------------------------------------
   ! inputs
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::u
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q
   integer::ngrid
   ! output
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:3)::bmagij
   ! declare local variables
   INTEGER ::i, j, k, l, m

  bmagij=0d0

   do k=ku1,ku2
      do j=ju1,ju2
         do i=iu1,iu2
            do l=1,ngrid
               do m=1,3
                  !! m+5 mandatory cf Bx=uin(l,i,j,k,6)
                  bmagij(l,i,j,k,m,m)=u(l,i,j,k,m+5)
               end do
            end do
         end do
      end do
   end do

   ! case Bx,y
   do k=ku1,ku2
      do j=min(1,ju1+1),ju2
         do i=iu1,max(1,iu2-1)
            do l=1,ngrid
               bmagij(l,i,j,k,1,2)=0.5d0*(q(l,i,j,k,6)+q(l,i,j-1,k,6))
            end do
         end do
      end do
   end do

   ! case Bx,z
   do k=min(1,ku1+1),ku2
      do j=ju1,ju2
         do i=iu1,max(1,iu2-1)
            do l=1,ngrid
               bmagij(l,i,j,k,1,3)=0.5d0*(q(l,i,j,k,6)+q(l,i,j,k-1,6))
            end do
         end do
      end do
   end do

   ! case By,x
   do k=ku1,ku2
      do j=ju1,max(1,ju2-1)
         do i=min(1,iu1+1),iu2
            do l=1,ngrid
               bmagij(l,i,j,k,2,1)=0.5d0*(q(l,i,j,k,7)+q(l,i-1,j,k,7))
            end do
         end do
      end do
   end do

   ! case By,z
   do k=min(1,ku1+1),ku2
      do j=ju1,max(1,ju2-1)
         do i=iu1,iu2
            do l=1,ngrid
               bmagij(l,i,j,k,2,3)=0.5d0*(q(l,i,j,k,7)+q(l,i,j,k-1,7))
            end do
         end do
      end do
   end do

   ! case Bz,x
   do k=ku1,max(1,ku2-1)
      do j=ju1,ju2
         do i=min(1,iu1+1),iu2
            do l=1,ngrid
               bmagij(l,i,j,k,3,1)=0.5d0*(q(l,i,j,k,8)+q(l,i-1,j,k,8))
            end do
         end do
      end do
   end do

   ! case Bz,y
   do k=ku1,max(1,ku2-1)
      do j=min(1,ju1+1),ju2
         do i=iu1,iu2
            do l=1,ngrid
               bmagij(l,i,j,k,3,2)=0.5d0*(q(l,i,j,k,8)+q(l,i,j-1,k,8))
            end do
         end do
      end do
   end do

end subroutine compute_bmagij
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine compute_bmagijbis(u,ngrid,bmagijbis)
   use amr_parameters
   use hydro_commons
   IMPLICIT NONE
   !-----------------------------------------------------------------
   ! bmagijbis(l,i,j,k,n) is the value of the magnetic field component
   ! Bn at i-1/2,j-1/2,k-1/2
   !-----------------------------------------------------------------
   ! inputs
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::u
   integer::ngrid
   ! output
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::bmagijbis
   ! declare local variables
   INTEGER ::i, j, k, l

   bmagijbis=0d0

   ! case Bx for Lorentz force EMF
   do k=min(1,ku1+1),ku2
      do j=min(1,ju1+1),ju2
         do i=iu1,iu2
            do l=1,ngrid
               bmagijbis(l,i,j,k,1)=0.25d0*(u(l,i,j,k,6)+u(l,i,j-1,k,6)+u(l,i,j,k-1,6)+u(l,i,j-1,k-1,6))
            end do
         end do
      end do
   end do

   ! case By for Lorentz force EMF
   do k=min(1,ku1+1),ku2
      do j=ju1,ju2
         do i=min(1,iu1+1),iu2
            do l=1,ngrid
               bmagijbis(l,i,j,k,2)=0.25d0*(u(l,i,j,k,7)+u(l,i-1,j,k,7)+u(l,i,j,k-1,7)+u(l,i-1,j,k-1,7))
            end do
         end do
      end do
   end do

   ! case Bz for Lorentz force EMF
   do k=ku1,ku2
      do j=min(1,ju1+1),ju2
         do i=min(1,iu1+1),iu2
            do l=1,ngrid
               bmagijbis(l,i,j,k,3)=0.25d0*(u(l,i,j,k,8)+u(l,i-1,j,k,8)+u(l,i,j-1,k,8)+u(l,i-1,j-1,k,8))
            end do
         end do
      end do
   end do

end subroutine compute_bmagijbis
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine computejb2(u,q,ngrid,dx,dy,dz,dt,bemfx,bemfy,bemfz,jemfx,jemfy,jemfz,bmagij,fluxmd,fluxad)

   USE amr_parameters
   use hydro_commons
   use nimhd_parameters
   IMPLICIT NONE
   !-----------------------------------------------------------------
   ! Build the geometric quantities shared by both non-ideal effects.
   ! Interpolates the magnetic field to the EMF edges (bemfx/y/z), the
   ! face locations (bmagij) and computes the current density J at the
   ! EMF edges (jemfx/y/z) and on the cell faces (jface). From J and B
   ! it forms the ideal-MHD energy-flux building blocks:
   !   fluxmd = J x B          (used by Ohmic dissipation)
   !   fluxad = (J x B) x B x B (used by ambipolar diffusion, if active)
   ! These are later scaled by the resistivities in computdifmag and
   ! computambip. Must be called before those two routines.
   !-----------------------------------------------------------------
   ! inputs
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::u
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q
   INTEGER ::ngrid
   REAL(dp)::dx,dy,dz,dt

   ! outputs
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::bemfx,bemfy,bemfz
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::jemfx,jemfy,jemfz
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:3)::bmagij
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::fluxmd,fluxad

   ! declare local variables
   INTEGER ::i, j, k, l, m
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::bmagijbis
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:3)::jface
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::bcenter
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:3)::fluxbis,fluxter,fluxquat
   real(dp)::bsquare
   real(dp)::computdxbis,computdybis,computdzbis

   fluxmd=0d0
   fluxad=0d0
   bmagij=0d0

   bemfx=0d0
   bemfy=0d0
   bemfz=0d0
   jemfx=0d0
   jemfy=0d0
   jemfz=0d0

   ! magnetic field at center of cells
   do k=ku1,ku2
      do j=ju1,ju2
         do i=iu1,iu2
            do l=1,ngrid
               bcenter(l,i,j,k,1)=q(l,i,j,k,6)
               bcenter(l,i,j,k,2)=q(l,i,j,k,7)
               bcenter(l,i,j,k,3)=q(l,i,j,k,8)
            end do
         end do
      end do
   end do

   call compute_bemf(u,q,ngrid,bemfx,bemfy,bemfz)

   call compute_bmagij(u,q,ngrid,bmagij)

   call compute_bmagijbis(u,ngrid,bmagijbis)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! computation of the component of j where EMFs are located
   ! jemfx(l,i,j,k,n) is the component Jn at i,j-1/2,k-1/2
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   do k=min(1,ku1+1),max(1,ku2-1)
      do j=min(1,ju1+1),max(1,ju2-1)
         do i=min(1,iu1+1),max(1,iu2-1)
            do l=1,ngrid
               jemfx(l,i,j,k,1)=(u(l,i,j,k,8)-u(l,i,j-1,k,8))/dy-(u(l,i,j,k,7)-u(l,i,j,k-1,7))/dz
               jemfx(l,i,j,k,2)=(bmagij(l,i,j,k,1,2)-bmagij(l,i,j,k-1,1,2))/dz- (bmagijbis(l,i+1,j,k,3)-bmagijbis(l,i,j,k,3))/dx
               jemfx(l,i,j,k,3)=(bmagijbis(l,i+1,j,k,2) -bmagijbis(l,i,j,k,2))/dx- (bmagij(l,i,j,k,1,3)-bmagij(l,i,j-1,k,1,3))/dy

               jemfy(l,i,j,k,1)=(bmagijbis(l,i,j+1,k,3)-bmagijbis(l,i,j,k,3))/dy-(bmagij(l,i,j,k,2,1) - bmagij(l,i,j,k-1,2,1) )/dz
               jemfy(l,i,j,k,2)=(u(l,i,j,k,6)-u(l,i,j,k-1,6))/dz-(u(l,i,j,k,8)-u(l,i-1,j,k,8))/dx
               jemfy(l,i,j,k,3)=(bmagij(l,i,j,k,2,3)-bmagij(l,i-1,j,k,2,3))/dx-(bmagijbis(l,i,j+1,k,1)-bmagijbis(l,i,j,k,1))/dy

               jemfz(l,i,j,k,1)=(bmagij(l,i,j,k,3,1) -bmagij(l,i,j-1,k,3,1))/dy-(bmagijbis(l,i,j,k+1,2)-bmagijbis(l,i,j,k,2))/dz
               jemfz(l,i,j,k,2)=( bmagijbis(l,i,j,k+1,1)-bmagijbis(l,i,j,k,1))/dz-(bmagij(l,i,j,k,3,2)-bmagij(l,i-1,j,k,3,2))/dx
               jemfz(l,i,j,k,3)=(u(l,i,j,k,7)-u(l,i-1,j,k,7))/dx-(u(l,i,j,k,6)-u(l,i,j-1,k,6))/dy
            end do
         end do
      end do
   end do

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! computation of the component of j at center of cell
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   ! computation of current on faces

   ! face at i-1/2,j,k

   do k=min(1,ku1+1),max(1,ku2-1)
      do j=min(1,ju1+1),max(1,ju2-1)
         do i=min(1,iu1+1),iu2
            do l=1,ngrid
               jface(l,i,j,k,1,1)=computdybis(bemfz,3,l,i,j,k,dy)-computdzbis(bemfy,2,l,i,j,k,dz)
            end do
         end do
      end do
   end do

   do k=min(1,ku1+1),max(1,ku2-1)
      do j=ju1,ju2
         do i=min(1,iu1+1),iu2
            do l=1,ngrid
               jface(l,i,j,k,2,1)=computdzbis(bemfy,1,l,i,j,k,dz)-computdxbis(bcenter,3,l,i-1,j,k,dx)
            end do
         end do
      end do
   end do

   do k=ku1,ku2
      do j=min(1,ju1+1),max(1,ju2-1)
         do i=min(1,iu1+1),iu2
            do l=1,ngrid
               jface(l,i,j,k,3,1)=computdxbis(bcenter,2,l,i-1,j,k,dx)-computdybis(bemfz,1,l,i,j,k,dy)
            end do
         end do
      end do
   end do

   ! face at i,j-1/2,k

   do k=min(1,ku1+1),max(1,ku2-1)
      do j=min(1,ju1+1),ju2
         do i=iu1,iu2
            do l=1,ngrid
               jface(l,i,j,k,1,2)=computdybis(bcenter,3,l,i,j-1,k,dy)-computdzbis(bemfx,2,l,i,j,k,dz)
            end do
         end do
      end do
   end do

   do k=min(1,ku1+1),max(1,ku2-1)
      do j=min(1,ju1+1),ju2
         do i=min(1,iu1+1),max(1,iu2-1)
            do l=1,ngrid
               jface(l,i,j,k,2,2)=computdzbis(bemfx,1,l,i,j,k,dz)-computdxbis(bemfz,3,l,i,j,k,dx)
            end do
         end do
      end do
   end do

   do k=ku1,ku2
      do j=min(1,ju1+1),ju2
         do i=min(1,iu1+1),max(1,iu2-1)
            do l=1,ngrid
               jface(l,i,j,k,3,2)=computdxbis(bemfz,2,l,i,j,k,dx)-computdybis(bcenter,1,l,i,j-1,k,dy)
            end do
         end do
      end do
   end do

   ! face at i,j,k-1/2

   do k=min(1,ku1+1),ku2
      do j=min(1,ju1+1),max(1,ju2-1)
         do i=iu1,iu2
            do l=1,ngrid
               jface(l,i,j,k,1,3)=computdybis(bemfx,3,l,i,j,k,dy)-computdzbis(bcenter,2,l,i,j,k-1,dz)
            end do
         end do
      end do
   end do

   do k=min(1,ku1+1),ku2
      do j=ju1,ju2
         do i=min(1,iu1+1),max(1,iu2-1)
            do l=1,ngrid
               jface(l,i,j,k,2,3)=computdzbis(bcenter,1,l,i,j,k-1,dz)-computdxbis(bemfy,3,l,i,j,k,dx)
            end do
         end do
      end do
   end do

   do k=min(1,ku1+1),ku2
      do j=min(1,ju1+1),max(1,ju2-1)
         do i=min(1,iu1+1),max(1,iu2-1)
            do l=1,ngrid
               jface(l,i,j,k,3,3)=computdxbis(bemfy,2,l,i,j,k,dx)-computdybis(bemfx,1,l,i,j,k,dx)
            end do
         end do
      end do
   end do


   do k=min(1,ku1+1),max(1,ku2-1)
      do j=min(1,ju1+1),max(1,ju2-1)
         do i=min(1,iu1+1),max(1,iu2-1)
            do l = 1, ngrid
               call crossprodbis(jface,bmagij,fluxbis,l,i,j,k)
               fluxmd(l,i,j,k,1)=fluxbis(l,i,j,k,1,1)
               fluxmd(l,i,j,k,2)=fluxbis(l,i,j,k,2,2)
               fluxmd(l,i,j,k,3)=fluxbis(l,i,j,k,3,3)
            end do
         end do
      end do
   end do

   if(nambipolar) then
      do k=min(1,ku1+1),max(1,ku2-1)
         do j=min(1,ju1+1),max(1,ju2-1)
            do i=min(1,iu1+1),max(1,iu2-1)
               do l = 1, ngrid
                  call crossprodbis(fluxbis,bmagij,fluxter,l,i,j,k)
                  call crossprodbis(fluxter,bmagij,fluxquat,l,i,j,k)
                  fluxad(l,i,j,k,1)=fluxquat(l,i,j,k,1,1)
                  fluxad(l,i,j,k,2)=fluxquat(l,i,j,k,2,2)
                  fluxad(l,i,j,k,3)=fluxquat(l,i,j,k,3,3)
               end do
            end do
         end do
      end do
   endif

end subroutine computejb2
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine computdifmag(u,ngrid,dx,dy,dz,dt,bemfx,bemfy,bemfz,jemfx,jemfy,jemfz,bmagij,fluxmd,emfohmdiss,fluxohm)

   use amr_parameters
   use hydro_commons
   use nimhd_parameters
   implicit none
   !-----------------------------------------------------------------
   ! Ohmic dissipation. Computes the Ohmic contribution to the EMF,
   !   emfohmdiss = -eta * J   (from dB/dt = -curl(eta*J)),
   ! evaluated at the EMF edges using the fixed resistivity eta = etaMD.
   ! If nimhdheating_in_flux is set, it also builds the associated Ohmic
   ! energy flux (fluxohm = eta * fluxmd) on the faces.
   ! Expects jemf* and fluxmd from computejb2.
   !-----------------------------------------------------------------
   ! inputs
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::u
   integer ::ngrid
   real(dp)::dx,dy,dz,dt
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::bemfx,bemfy,bemfz
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::jemfx,jemfy,jemfz
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:3)::bmagij
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::fluxmd

   ! outputs
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3):: emfohmdiss,fluxohm

   ! local variables
   integer ::i,j,k,l,h
   emfohmdiss = 0d0
   fluxohm = 0d0

   do k=min(1,ku1+1),max(1,ku2-1)
      do j=min(1,ju1+1),max(1,ju2-1)
         do i=min(1,iu1+1),max(1,iu2-1)

            do l=1,ngrid
               ! WARNING dB/dt=-curl(eta*J)
               emfohmdiss(l,i,j,k,1)=-etaMD*jemfx(l,i,j,k,1)
               emfohmdiss(l,i,j,k,2)=-etaMD*jemfy(l,i,j,k,2)
               emfohmdiss(l,i,j,k,3)=-etaMD*jemfz(l,i,j,k,3)
            end do

            if(nimhdheating_in_flux) then
               do h = 1,3
                  do l=1,ngrid
                     fluxohm(l,i,j,k,h)=etaMD*fluxmd(l,i,j,k,h)
                  enddo
               enddo
            endif

         end do
      end do
   end do

end subroutine computdifmag
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine computambip(u,ngrid,dx,dy,dz,dt,bemfx,bemfy,bemfz,jemfx,jemfy,jemfz,bmagij,fluxad,emfambdiff,fluxambdiff)

   use amr_commons
   use amr_parameters
   use hydro_commons
   use nimhd_parameters
   use const
   implicit none
   !-----------------------------------------------------------------
   ! Ambipolar diffusion. Computes the ambipolar contribution to the
   ! EMF from the Lorentz force F = J x B:
   !   emfambdiff = beta * (F x B)   (from dB/dt = curl(beta*(JxB)xB)),
   ! evaluated at the EMF edges, with the fixed coefficient
   ! beta = 1/(gammaAD*rho). If nimhdheating_in_flux is set, it also
   ! builds the ambipolar energy flux (fluxambdiff = -beta * fluxad) on
   ! the faces. Expects bemf*, jemf* and fluxad from computejb2.
   !-----------------------------------------------------------------
   ! inputs
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::u
   integer ::ngrid
   real(dp)::dx,dy,dz,dt
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::bemfx,bemfy,bemfz
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::jemfx,jemfy,jemfz
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::fluxad
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:3)::bmagij

   ! outputs
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::emfambdiff
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::fluxambdiff
   ! declare local variables
   INTEGER ::i, j, k, l
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::florentzx,florentzy,florentzz

   real(dp)::v1x,v1y,v1z,v2x,v2y,v2z
   real(dp)::rhox,rhoy,rhoz,rhofx,rhofy,rhofz
   real(dp)::crossprodx,crossprody,crossprodz

   emfambdiff=0d0
   fluxambdiff=0d0

   florentzx=0d0
   florentzy=0d0
   florentzz=0d0

   ! compute Loretz force
   do k=min(1,ku1+1),max(1,ku2-1)
      do j=min(1,ju1+1),max(1,ju2-1)
         do i=min(1,iu1+1),max(1,iu2-1)
            do l = 1, ngrid
               call crossprod(jemfx,bemfx,florentzx,l,i,j,k)
               call crossprod(jemfy,bemfy,florentzy,l,i,j,k)
               call crossprod(jemfz,bemfz,florentzz,l,i,j,k)
            end do
         end do
      end do
   end do


   !dtlim=dt!*coefalfven
   !dt est deja dtnew, qui a été choisi comme le dt normal (avec la condition de courant) ou le dt normal seuillé si le dtAD est trop faible(bricolo)

   do k=min(1,ku1+1),max(1,ku2-1)
      do j=min(1,ju1+1),max(1,ju2-1)
         do i=min(1,iu1+1),max(1,iu2-1)

            do l = 1, ngrid

               rhox=0.25d0*(u(l,i,j,k,1)+u(l,i,j-1,k,1)+u(l,i,j,k-1,1)+u(l,i,j-1,k-1,1))
               rhoy=0.25d0*(u(l,i,j,k,1)+u(l,i-1,j,k,1)+u(l,i,j,k-1,1)+u(l,i-1,j,k-1,1))
               rhoz=0.25d0*(u(l,i,j,k,1)+u(l,i-1,j,k,1)+u(l,i,j-1,k,1)+u(l,i-1,j-1,k,1))

               ! EMF x

               v1x=florentzx(l,i,j,k,1)
               v1y=florentzx(l,i,j,k,2)
               v1z=florentzx(l,i,j,k,3)
               v2x=bemfx(l,i,j,k,1)
               v2y=bemfx(l,i,j,k,2)
               v2z=bemfx(l,i,j,k,3)
               emfambdiff(l,i,j,k,1)=crossprodx(v1x,v1y,v1z,v2x,v2y,v2z)/(gammaAD*rhox)

               ! EMF y

               v1x=florentzy(l,i,j,k,1)
               v1y=florentzy(l,i,j,k,2)
               v1z=florentzy(l,i,j,k,3)
               v2x=bemfy(l,i,j,k,1)
               v2y=bemfy(l,i,j,k,2)
               v2z=bemfy(l,i,j,k,3)
               emfambdiff(l,i,j,k,2)=crossprody(v1x,v1y,v1z,v2x,v2y,v2z)/(gammaAD*rhoy)

               ! EMF z

               v1x=florentzz(l,i,j,k,1)
               v1y=florentzz(l,i,j,k,2)
               v1z=florentzz(l,i,j,k,3)
               v2x=bemfz(l,i,j,k,1)
               v2y=bemfz(l,i,j,k,2)
               v2z=bemfz(l,i,j,k,3)
               emfambdiff(l,i,j,k,3)=crossprodz(v1x,v1y,v1z,v2x,v2y,v2z)/(gammaAD*rhoz)

               if(nimhdheating_in_flux) then
                  ! energy flux on faces
                  rhofx=0.5d0*(u(l,i,j,k,1)+u(l,i-1,j,k,1))
                  rhofy=0.5d0*(u(l,i,j,k,1)+u(l,i,j-1,k,1))
                  rhofz=0.5d0*(u(l,i,j,k,1)+u(l,i,j,k-1,1))

                  fluxambdiff(l,i,j,k,1)=-fluxad(l,i,j,k,1)/(gammaAD*rhofx)
                  fluxambdiff(l,i,j,k,2)=-fluxad(l,i,j,k,2)/(gammaAD*rhofy)
                  fluxambdiff(l,i,j,k,3)=-fluxad(l,i,j,k,3)/(gammaAD*rhofz)
               endif
            end do
         end do
      end do
   end do

end subroutine computambip
!###########################################################
!###########################################################
!###########################################################
!###########################################################

! VECTOR HELPER FUNCTIONS
! computd{x,y,z}bis : first-order forward finite difference of component
!   n2 of the vector field vec along the x, y or z direction, at cell dx.
! crossprod / crossprodbis : cross product of two vector (rank-1) or
!   tensor (rank-2, done column by column) fields at cell (l,i,j,k).
! crossprod{x,y,z} : the three scalar components of a single cross product.

!###########################################################
! forward x-derivative of component n2 of vec, divided by dx
double precision function computdxbis(vec,n2,l,i,j,k,dx)

   use amr_parameters,only:dp,nvector
   use hydro_parameters,only:iu1,iu2,ju1,ju2,ku1,ku2
   implicit none
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::vec
   real(dp)::dx
   integer::n2,l,i,j,k

   computdxbis = (vec(l,i+1,j,k,n2) - vec(l,i,j,k,n2)) / dx

end function computdxbis

! forward y-derivative of component n2 of vec, divided by dx
double precision  function computdybis(vec,n2,l,i,j,k,dx)

   use amr_parameters,only:dp,nvector
   use hydro_parameters,only:iu1,iu2,ju1,ju2,ku1,ku2
   implicit none
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::vec
   real(dp)::dx
   integer::n2,l,i,j,k

   computdybis=(vec(l,i,j+1,k,n2) - vec(l,i,j,k,n2)) / dx

end function computdybis

! forward z-derivative of component n2 of vec, divided by dx
double precision  function computdzbis(vec,n2,l,i,j,k,dx)

   use amr_parameters,only:dp,nvector
   use hydro_parameters,only:iu1,iu2,ju1,ju2,ku1,ku2
   implicit none
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::vec
   real(dp)::dx
   integer::n2,l,i,j,k

   computdzbis = (vec(l,i,j,k+1,n2) - vec(l,i,j,k,n2)) / dx

end function computdzbis
!###########################################################
! column-wise cross product of two rank-2 (3x3) fields:
! v1crossv2(:,n) = vec1(:,n) x vec2(:,n) for each column n=1,3
subroutine crossprodbis(vec1,vec2,v1crossv2,l,i,j,k)

   use amr_parameters,only:dp,nvector
   use hydro_parameters,only:iu1,iu2,ju1,ju2,ku1,ku2
   implicit none
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:3)::vec1,vec2,v1crossv2
   integer ::l,i,j,k

   real(dp)::v1x,v1y,v1z,v2x,v2y,v2z,crossprodx,crossprody,crossprodz

   integer::n

   do n=1,3
      v1x=vec1(l,i,j,k,1,n)
      v1y=vec1(l,i,j,k,2,n)
      v1z=vec1(l,i,j,k,3,n)

      v2x=vec2(l,i,j,k,1,n)
      v2y=vec2(l,i,j,k,2,n)
      v2z=vec2(l,i,j,k,3,n)

      v1crossv2(l,i,j,k,1,n)=crossprodx(v1x,v1y,v1z,v2x,v2y,v2z)
      v1crossv2(l,i,j,k,2,n)=crossprody(v1x,v1y,v1z,v2x,v2y,v2z)
      v1crossv2(l,i,j,k,3,n)=crossprodz(v1x,v1y,v1z,v2x,v2y,v2z)
   end do

end subroutine crossprodbis

! cross product of two vector fields at cell (l,i,j,k): v1crossv2 = vec1 x vec2
subroutine crossprod(vec1,vec2,v1crossv2,l,i,j,k)

   use amr_parameters,only:dp,nvector
   use hydro_parameters,only:iu1,iu2,ju1,ju2,ku1,ku2
   implicit none
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::vec1,vec2,v1crossv2
   integer ::l,i,j,k

   real(dp)::v1x,v1y,v1z,v2x,v2y,v2z,crossprodx,crossprody,crossprodz

   v1x=vec1(l,i,j,k,1)
   v1y=vec1(l,i,j,k,2)
   v1z=vec1(l,i,j,k,3)

   v2x=vec2(l,i,j,k,1)
   v2y=vec2(l,i,j,k,2)
   v2z=vec2(l,i,j,k,3)

   v1crossv2(l,i,j,k,1)=crossprodx(v1x,v1y,v1z,v2x,v2y,v2z)
   v1crossv2(l,i,j,k,2)=crossprody(v1x,v1y,v1z,v2x,v2y,v2z)
   v1crossv2(l,i,j,k,3)=crossprodz(v1x,v1y,v1z,v2x,v2y,v2z)

end subroutine crossprod
!###########################################################
double precision function  crossprodx(v1x,v1y,v1z,v2x,v2y,v2z)

   ! x component of a cross product
   use amr_parameters,only:dp
   implicit none
   real(dp)::v1x,v1y,v1z,v2x,v2y,v2z

   crossprodx=v1y*v2z-v1z*v2y

end function crossprodx

double precision function crossprody(v1x,v1y,v1z,v2x,v2y,v2z)

   ! y component of a cross product
   use amr_parameters,only:dp
   implicit none
   real(dp)::v1x,v1y,v1z,v2x,v2y,v2z

   crossprody=v1z*v2x-v1x*v2z

end function crossprody

double precision function crossprodz(v1x,v1y,v1z,v2x,v2y,v2z)

   ! z component of a cross product
   use amr_parameters,only:dp
   implicit none
   real(dp)::v1x,v1y,v1z,v2x,v2y,v2z

   crossprodz=v1x*v2y-v2x*v1y

end function crossprodz
!###########################################################
!###########################################################
!###########################################################
