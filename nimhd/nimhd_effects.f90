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
   implicit none
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3),intent(in)::u
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),intent(in)::q
   integer,intent(in)::ngrid
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3),intent(out)::bemfx,bemfy,bemfz
   !-------------------------------------------
   ! Interpolates the magnetic field at location of EMF (edges)
   !-------------------------------------------
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
   implicit none
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3),intent(in)::u
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),intent(in)::q
   integer,intent(in)::ngrid
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:3),intent(out)::bmagij
   !-----------------------------------------------------------------
   ! bmagij is the value of the magnetic field Bi where Bj
   ! is naturally defined; Ex bmagij(l,i,j,k,1,2) is Bx at i,j-1/2,k
   ! and we can write it Bx,y
   !-----------------------------------------------------------------
   integer ::i, j, k, l, m

   bmagij=0d0

   ! Diagonal: Bx x, By y, Bz z
   do k=ku1,ku2
      do j=ju1,ju2
         do i=iu1,iu2
            do l=1,ngrid
               bmagij(l,i,j,k,1,1)=u(l,i,j,k,6)
               bmagij(l,i,j,k,2,2)=u(l,i,j,k,7)
               bmagij(l,i,j,k,3,3)=u(l,i,j,k,8)
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
   implicit none
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3),intent(in)::u
   integer,intent(in)::ngrid
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3),intent(out)::bmagijbis
   !-----------------------------------------------------------------
   ! bmagijbis(l,i,j,k,n) is the value of the magnetic field component
   ! Bn at i-1/2,j-1/2,k-1/2
   !-----------------------------------------------------------------
   integer ::i, j, k, l

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
   use amr_parameters
   use hydro_commons
   use nimhd_parameters
   implicit none
   ! inputs
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3),intent(in)::u
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),intent(in)::q
   integer,intent(in)::ngrid
   real(dp),intent(in)::dx,dy,dz,dt
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3),intent(in)::bemfx,bemfy,bemfz
   ! outputs
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3),intent(out)::jemfx,jemfy,jemfz
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:3),intent(out)::bmagij
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3),intent(out)::fluxmd,fluxad
   !-----------------------------------------------------------------
   ! Build the geometric quantities shared by both non-ideal effects.
   ! the
   ! face locations (bmagij) and computes the current density J at the
   ! EMF edges (jemfx/y/z) and on the cell faces (jface). From J and B
   ! it forms the ideal-MHD energy-flux building blocks:
   !   fluxmd = J x B          (used by Ohmic dissipation)
   !   fluxad = (J x B) x B x B (used by ambipolar diffusion, if active)
   ! These are later scaled by the resistivities in computdifmag and
   ! computambip. Must be called before those two routines.
   !-----------------------------------------------------------------
   integer ::i, j, k, l, m, n
   real(dp)::v1x,v1y,v1z,v2x,v2y,v2z
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::bmagijbis
   real(dp),dimension(1:nvector,1:3,1:3)::jface
   real(dp),dimension(1:nvector,1:3,1:3)::fluxbis
   real(dp)::bsquare
   real(dp)::computdxbis,computdybis,computdzbis  !forward derivatives

   fluxmd=0d0
   fluxad=0d0
   jemfx=0d0
   jemfy=0d0
   jemfz=0d0

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

   do k=min(1,ku1+1),max(1,ku2-1)
      do j=min(1,ju1+1),max(1,ju2-1)
         do i=min(1,iu1+1),max(1,iu2-1)

            ! computation of current on faces
            ! (q contains the magnetic field at center of cells)
            do l = 1, ngrid
               ! face at i-1/2,j,k
               computdybis = (bemfz(l,i,j+1,k  ,3) - bemfz(l,i,j,k,3)) / dy
               computdzbis = (bemfy(l,i,j  ,k+1,2) - bemfy(l,i,j,k,2)) / dz
               jface(l,1,1) = computdybis - computdzbis

               computdzbis = (bemfy(l,i,j,k+1,1) - bemfy(l,i  ,j,k,1)) / dz
               computdxbis = (    q(l,i,j,k  ,8) -     q(l,i-1,j,k,8)) / dx
               jface(l,2,1) = computdzbis - computdxbis

               computdxbis = (    q(l,i,j  ,k,7) -     q(l,i-1,j,k,7)) / dx
               computdybis = (bemfz(l,i,j+1,k,1) - bemfz(l,i  ,j,k,1)) / dy
               jface(l,3,1) = computdxbis - computdybis

               ! face at i,j-1/2,k
               computdybis = (    q(l,i,j,k  ,8) -     q(l,i,j-1,k,8)) / dy
               computdzbis = (bemfx(l,i,j,k+1,2) - bemfx(l,i,j  ,k,2)) / dz
               jface(l,1,2) = computdybis - computdzbis

               computdzbis = (bemfx(l,i  ,j,k+1,1) - bemfx(l,i,j,k,1)) / dz
               computdxbis = (bemfz(l,i+1,j,k  ,3) - bemfz(l,i,j,k,3)) / dx
               jface(l,2,2) = computdzbis - computdxbis

               computdxbis = (bemfz(l,i+1,j,k,2) - bemfz(l,i,j  ,k,2)) / dx
               computdybis = (    q(l,i  ,j,k,6) -     q(l,i,j-1,k,6)) / dy
               jface(l,3,2) = computdxbis - computdybis

               ! face at i,j,k-1/2
               computdybis = (bemfx(l,i,j+1,k,3) - bemfx(l,i,j,k  ,3)) / dy
               computdzbis = (    q(l,i,j  ,k,7) -     q(l,i,j,k-1,7)) / dz
               jface(l,1,3) = computdybis - computdzbis

               computdzbis = (    q(l,i  ,j,k,6) -     q(l,i,j,k-1,6)) / dz
               computdxbis = (bemfy(l,i+1,j,k,3) - bemfy(l,i,j,k  ,3)) / dx
               jface(l,2,3) = computdzbis - computdxbis

               computdxbis = (bemfy(l,i+1,j  ,k,2) - bemfy(l,i,j,k,2)) / dx
               computdybis = (bemfx(l,i  ,j+1,k,1) - bemfx(l,i,j,k,1)) / dx
               jface(l,3,3) = computdxbis - computdybis
            end do 

            do l = 1, ngrid
               do n=1,3
                  v1x=jface(l,1,n)
                  v1y=jface(l,2,n)
                  v1z=jface(l,3,n)

                  v2x=bmagij(l,i,j,k,1,n)
                  v2y=bmagij(l,i,j,k,2,n)
                  v2z=bmagij(l,i,j,k,3,n)

                  fluxbis(l,1,n)=v1y*v2z-v1z*v2y
                  fluxbis(l,2,n)=v1z*v2x-v1x*v2z
                  fluxbis(l,3,n)=v1x*v2y-v2x*v1y
               end do

               fluxmd(l,i,j,k,1)=fluxbis(l,1,1)
               fluxmd(l,i,j,k,2)=fluxbis(l,2,2)
               fluxmd(l,i,j,k,3)=fluxbis(l,3,3)
            end do

            if(nambipolar) then
               ! Compute fluxad from crossproduct ((J x B) x B) x B
               do l = 1, ngrid

                  ! x-component
                  v1x=fluxbis(l,1,1)
                  v2x=bmagij(l,i,j,k,1,1)
                  v2y=bmagij(l,i,j,k,2,1)
                  v2z=bmagij(l,i,j,k,3,1)
                  v1y=fluxbis(l,3,1)*v2x-v1x*v2z
                  v1z=v1x*v2y-v2x*fluxbis(l,2,1)
                  fluxad(l,i,j,k,1)=v1y*v2z-v1z*v2y

                  ! y-component
                  v1y=fluxbis(l,2,2)
                  v2x=bmagij(l,i,j,k,1,2)
                  v2y=bmagij(l,i,j,k,2,2)
                  v2z=bmagij(l,i,j,k,3,2)
                  v1x=v1y*v2z-fluxbis(l,3,2)*v2y
                  v1z=fluxbis(l,1,2)*v2y-v2x*v1y
                  fluxad(l,i,j,k,2)=v1z*v2x-v1x*v2z

                  ! z-component
                  v1z=fluxbis(l,3,3)
                  v2x=bmagij(l,i,j,k,1,3)
                  v2y=bmagij(l,i,j,k,2,3)
                  v2z=bmagij(l,i,j,k,3,3)
                  v1x=fluxbis(l,2,3)*v2z-v1z*v2y
                  v1y=v1z*v2x-fluxbis(l,1,3)*v2z
                  fluxad(l,i,j,k,3)=v1x*v2y-v2x*v1y

               end do
            endif
         end do
      end do
   end do

end subroutine computejb2
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine computdifmag(ngrid,jemfx,jemfy,jemfz,emfohmdiss)

   use amr_parameters
   use hydro_commons
   use nimhd_parameters
   implicit none
   ! inputs
   integer,intent(in)::ngrid
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3),intent(in)::jemfx,jemfy,jemfz
   ! outputs
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3),intent(out):: emfohmdiss
   !-----------------------------------------------------------------
   ! Ohmic dissipation. Computes the Ohmic contribution to the EMF,
   !   emfohmdiss = -eta * J   (from dB/dt = -curl(eta*J)),
   ! evaluated at the EMF edges using the fixed resistivity eta = etaMD.
   ! Expects jemf* from computejb2.
   !-----------------------------------------------------------------
   integer ::i,j,k,l,h

   emfohmdiss = 0d0

   do k=min(1,ku1+1),max(1,ku2-1)
      do j=min(1,ju1+1),max(1,ju2-1)
         do i=min(1,iu1+1),max(1,iu2-1)

            do l=1,ngrid
               ! WARNING dB/dt=-curl(eta*J)
               emfohmdiss(l,i,j,k,1)=-etaMD*jemfx(l,i,j,k,1)
               emfohmdiss(l,i,j,k,2)=-etaMD*jemfy(l,i,j,k,2)
               emfohmdiss(l,i,j,k,3)=-etaMD*jemfz(l,i,j,k,3)
            end do

         end do
      end do
   end do

end subroutine computdifmag
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine compute_heating_difmag(ngrid,fluxmd,fluxohm)
   use amr_parameters
   use hydro_commons
   use nimhd_parameters
   implicit none
   ! inputs
   integer,intent(in)::ngrid
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3),intent(in)::fluxmd
   ! outputs
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3),intent(out)::fluxohm
   !-----------------------------------------------------------------
   ! Compute the Ohmic energy flux (fluxohm = eta * fluxmd) on the faces.
   ! Expects fluxmd from computejb2.
   ! Resistivity is assumed to be constant: eta = etaMD
   !-----------------------------------------------------------------
   integer ::i,j,k,l,h

   fluxohm = 0d0

   do k=min(1,ku1+1),max(1,ku2-1)
      do j=min(1,ju1+1),max(1,ju2-1)
         do i=min(1,iu1+1),max(1,iu2-1)
            do l=1,ngrid
               fluxohm(l,i,j,k,1)=etaMD*fluxmd(l,i,j,k,1)
               fluxohm(l,i,j,k,2)=etaMD*fluxmd(l,i,j,k,2)
               fluxohm(l,i,j,k,3)=etaMD*fluxmd(l,i,j,k,3)
            enddo
         end do
      end do
   end do

end subroutine compute_heating_difmag
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine computambip(u,ngrid,bemfx,bemfy,bemfz,jemfx,jemfy,jemfz,emfambdiff)
   use amr_commons
   use amr_parameters
   use hydro_commons
   use nimhd_parameters
   use const
   implicit none
   ! inputs
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3),intent(in)::u
   integer,intent(in)::ngrid
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3),intent(in)::bemfx,bemfy,bemfz
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3),intent(in)::jemfx,jemfy,jemfz
   ! outputs
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3),intent(out)::emfambdiff
   !-----------------------------------------------------------------
   ! Ambipolar diffusion. Computes the ambipolar contribution to the
   ! EMF from the Lorentz force F = J x B:
   !   emfambdiff = beta * (F x B)   (from dB/dt = curl(beta*(JxB)xB)),
   ! evaluated at the EMF edges, with the fixed coefficient
   ! beta = 1/(gammaAD*rho).
   ! Expects bemf*, jemf* from computejb2.
   !-----------------------------------------------------------------
   integer ::i, j, k, l
   real(dp)::jx,jy,jz,bx,by,bz,fx,fy,fz,beta_x,beta_y,beta_z
   real(dp)::rhox,rhoy,rhoz

   emfambdiff=0d0

   ! Compute (J x B) x B
   do k=min(1,ku1+1),max(1,ku2-1)
      do j=min(1,ju1+1),max(1,ju2-1)
         do i=min(1,iu1+1),max(1,iu2-1)

            do l = 1, ngrid

               rhox=0.25d0*(u(l,i,j,k,1)+u(l,i,j-1,k,1)+u(l,i,j,k-1,1)+u(l,i,j-1,k-1,1))
               rhoy=0.25d0*(u(l,i,j,k,1)+u(l,i-1,j,k,1)+u(l,i,j,k-1,1)+u(l,i-1,j,k-1,1))
               rhoz=0.25d0*(u(l,i,j,k,1)+u(l,i-1,j,k,1)+u(l,i,j-1,k,1)+u(l,i-1,j-1,k,1))

               ! EMF x
               jx=jemfx(l,i,j,k,1)
               jy=jemfx(l,i,j,k,2)
               jz=jemfx(l,i,j,k,3)
               bx=bemfx(l,i,j,k,1)
               by=bemfx(l,i,j,k,2)
               bz=bemfx(l,i,j,k,3)
               fy=jz*bx-jx*bz !florentzx(l,2)
               fz=jx*by-bx*jy !florentzx(l,3)
               beta_x=1d0/(gammaAD*rhox)
               emfambdiff(l,i,j,k,1) = (fy*bz-fz*by) * beta_x

               ! EMF y
               jx=jemfy(l,i,j,k,1)
               jy=jemfy(l,i,j,k,2)
               jz=jemfy(l,i,j,k,3)
               bx=bemfy(l,i,j,k,1)
               by=bemfy(l,i,j,k,2)
               bz=bemfy(l,i,j,k,3)
               fx=jy*bz-jz*by !florentzy(l,1)
               fz=jx*by-bx*jy !florentzy(l,3)
               beta_y=1d0/(gammaAD*rhoy)
               emfambdiff(l,i,j,k,2)=(fz*bx-fx*bz) * beta_y

               ! EMF z
               jx=jemfz(l,i,j,k,1)
               jy=jemfz(l,i,j,k,2)
               jz=jemfz(l,i,j,k,3)
               bx=bemfz(l,i,j,k,1)
               by=bemfz(l,i,j,k,2)
               bz=bemfz(l,i,j,k,3)
               fx=jy*bz-jz*by !florentzz(l,1)
               fy=jz*bx-jx*bz !florentzz(l,2)
               beta_z=1d0/(gammaAD*rhoz)
               emfambdiff(l,i,j,k,3)=(fx*by-bx*fy) * beta_z

            end do

         end do
      end do
   end do

end subroutine computambip
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine compute_heating_ambip(u,ngrid,fluxad,fluxambdiff)
   use amr_commons
   use amr_parameters
   use hydro_commons
   use nimhd_parameters
   use const
   implicit none
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3),intent(in)::u
   integer,intent(in)::ngrid
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3),intent(in)::fluxad
   ! outputs
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3),intent(out)::fluxambdiff
   !-----------------------------------------------------------------
   ! Calculate the ambipolar energy flux
   !   fluxambdiff = -beta * fluxad on
   ! the faces.
   ! Assume fixed coefficient
   !   beta = 1/(gammaAD*rho)
   !-----------------------------------------------------------------
   integer ::i, j, k, l
   real(dp)::beta_x,beta_y,beta_z
   real(dp)::rhofx,rhofy,rhofz

   fluxambdiff=0d0

   ! Compute (J x B) x B
   do k=min(1,ku1+1),max(1,ku2-1)
      do j=min(1,ju1+1),max(1,ju2-1)
         do i=min(1,iu1+1),max(1,iu2-1)

            do l = 1, ngrid
               ! energy flux on faces
               rhofx=0.5d0*(u(l,i,j,k,1)+u(l,i-1,j,k,1))
               rhofy=0.5d0*(u(l,i,j,k,1)+u(l,i,j-1,k,1))
               rhofz=0.5d0*(u(l,i,j,k,1)+u(l,i,j,k-1,1))

               beta_x=1d0/(gammaAD*rhofx)
               beta_y=1d0/(gammaAD*rhofy)
               beta_z=1d0/(gammaAD*rhofz)

               fluxambdiff(l,i,j,k,1)=-fluxad(l,i,j,k,1) * beta_x
               fluxambdiff(l,i,j,k,2)=-fluxad(l,i,j,k,2) * beta_y
               fluxambdiff(l,i,j,k,3)=-fluxad(l,i,j,k,3) * beta_z
            end do

         end do
      end do
   end do

end subroutine compute_heating_ambip
!###########################################################
!###########################################################
!###########################################################
!###########################################################

