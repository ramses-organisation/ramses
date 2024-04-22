!  by Jacques Masson, Benoit Commercon and Neil Vaytet

!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine computejb2(u,q,ngrid,dx,dy,dz,dt,bemfx,bemfy,bemfz,jemfx,jemfy,jemfz,bmagij,florentzx,florentzy,florentzz,fluxmd,fluxh,fluxad)

   USE amr_parameters
   use hydro_commons
   use nimhd_parameters
   IMPLICIT NONE

   ! inputs
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::u 
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q 

   INTEGER ::ngrid
   REAL(dp)::dx,dy,dz,dt

   ! outputs
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::bemfx,bemfy,bemfz
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::jemfx,jemfy,jemfz
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::florentzx,florentzy,florentzz
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::fluxmd,fluxh,fluxad
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:3)::bmagij

   ! declare local variables
   INTEGER ::i, j, k, l, m, n 

   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::bmagijbis
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::flxmagxx,flxmagxy,flxmagxz,flxmagyx,flxmagyy,flxmagyz,flxmagzx,flxmagzy,flxmagzz
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:3)::jface
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::bcenter
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:3)::fluxbis,fluxter,fluxquat
   real(dp)::b12x,b12y,b12z,emag,bsquare
   real(dp)::computdivbisx,computdivbisy,computdivbisz
   real(dp)::computdxbis,computdybis,computdzbis

   ! magnetic field at center of cells
   do k=ku1,ku2
      do j=ju1,ju2
         do i=iu1,iu2
            do l=1,ngrid
               bcenter(l,i,j,k,nxx)=q(l,i,j,k,6)
               bcenter(l,i,j,k,nyy)=q(l,i,j,k,7)
               bcenter(l,i,j,k,nzz)=q(l,i,j,k,8)
            end do
         end do
      end do
   end do

   !!!!!!!!!!!!!!!!!!
   ! EMF x
   !!!!!!!!!!!!!!!!!!

   ! magnetic field at location of EMF

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

   ! magnetic field at location of EMF

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

   ! magnetic field at location of EMF

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

   ! bmagij is the value of the magnetic field Bi where Bj 
   ! is naturally defined; Ex bmagij(l,i,j,k,1,2) is Bx at i,j-1/2,k
   ! and we can write it Bx,y

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

   !!!!!!!!!!!!!!!!!!
   ! bmagijbis(l,i,j,k,n) is the value of the magnetic field component
   ! Bn at i-1/2,j-1/2,k-1/2
   !!!!!!!!!!!!!!!!!!

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

   if(nambipolar) then
   ! EMF x
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
   endif

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
                  fluxh(l,i,j,k,1)=fluxter(l,i,j,k,1,1)
                  fluxh(l,i,j,k,2)=fluxter(l,i,j,k,2,2)
                  fluxh(l,i,j,k,3)=fluxter(l,i,j,k,3,3)
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
   real(dp)::rhox,rhoy,rhoz,epsx,epsy,epsz,bsquarex,bsquarey,bsquarez
   real(dp)::tcellx,tcelly,tcellz,etaod2x,etaod2y,etaod2z
   real(dp)::ionisrate
   real(dp)::rhof,bsqf,epsf,tcellf
   real(dp)::etaohmdiss,etaod2,etaohmdissbricolo
   integer , dimension(1:3) :: index_i,index_j,index_k

   index_i = (/1,0,0/)
   index_j = (/0,1,0/)
   index_k = (/0,0,1/)

   do k=min(1,ku1+1),max(1,ku2-1)
      do j=min(1,ju1+1),max(1,ju2-1)
         do i=min(1,iu1+1),max(1,iu2-1)

            do l=1,ngrid

               rhox=0.25d0*(u(l,i,j,k,   1)+u(l,i  ,j-1,k,   1)+u(l,i,j  ,k-1,   1)+u(l,i  ,j-1,k-1,   1))
               rhoy=0.25d0*(u(l,i,j,k,   1)+u(l,i-1,j  ,k,   1)+u(l,i,j  ,k-1,   1)+u(l,i-1,j  ,k-1,   1))
               rhoz=0.25d0*(u(l,i,j,k,   1)+u(l,i-1,j  ,k,   1)+u(l,i,j-1,k  ,   1)+u(l,i-1,j-1,k  ,   1))

               epsx=0.25d0*(u(l,i,j,k,5)+u(l,i  ,j-1,k,5)+u(l,i,j  ,k-1,5)+u(l,i  ,j-1,k-1,5))
               epsy=0.25d0*(u(l,i,j,k,5)+u(l,i-1,j  ,k,5)+u(l,i,j  ,k-1,5)+u(l,i-1,j  ,k-1,5))
               epsz=0.25d0*(u(l,i,j,k,5)+u(l,i-1,j  ,k,5)+u(l,i,j-1,k  ,5)+u(l,i-1,j-1,k  ,5))

               bsquarex=bemfx(l,i,j,k,1)**2+bemfx(l,i,j,k,2)**2+bemfx(l,i,j,k,3)**2
               bsquarey=bemfy(l,i,j,k,1)**2+bemfy(l,i,j,k,2)**2+bemfy(l,i,j,k,3)**2
               bsquarez=bemfz(l,i,j,k,1)**2+bemfz(l,i,j,k,2)**2+bemfz(l,i,j,k,3)**2

               if(ntestDADM.eq.1)then
                     tcellx=1.0d0
                     tcelly=1.0d0
                     tcellz=1.0d0
               else
                     call ideal_gas_temperature(rhox, epsx, tcellx)
                     call ideal_gas_temperature(rhoy, epsy, tcelly)
                     call ideal_gas_temperature(rhoz, epsz, tcellz)
               endif

               ionisrate=default_ionisrate 
               etaod2x=etaohmdissbricolo(rhox,bsquarex,tcellx,dt,dx,ionisrate)
               etaod2y=etaohmdissbricolo(rhoy,bsquarey,tcelly,dt,dx,ionisrate)
               etaod2z=etaohmdissbricolo(rhoz,bsquarez,tcellz,dt,dx,ionisrate)
               ! TC: shouldn't dy and dz be used here?  

               ! WARNING dB/dt=-curl(eta*J)
               emfohmdiss(l,i,j,k,nxx)=-etaod2x*jemfx(l,i,j,k,1)
               emfohmdiss(l,i,j,k,nyy)=-etaod2y*jemfy(l,i,j,k,2)
               emfohmdiss(l,i,j,k,nzz)=-etaod2z*jemfz(l,i,j,k,3)

               do h = 1,3
                  rhof=0.5d0*(u(l,i,j,k,1)+u(l,i-index_i(h),j-index_j(h),k-index_k(h),1))
                  epsf=0.5d0*(u(l,i,j,k,5)+u(l,i-index_i(h),j-index_j(h),k-index_k(h),5))
                  bsqf=bmagij(l,i,j,k,1,h)**2+bmagij(l,i,j,k,2,h)**2+bmagij(l,i,j,k,3,h)**2

                  ! Compute gas temperature in cgs
                  if(ntestDADM.eq.1)then
                        tcellf=1.0d0
                  else 
                        call ideal_gas_temperature(rhof, epsf, tcellf)
                  endif
                     
                  etaod2=etaohmdiss(rhof,bsqf,tcellf,ionisrate)
                  fluxohm(l,i,j,k,h)=etaod2*fluxmd(l,i,j,k,h)
               enddo

            end do
         end do
      end do
   end do
  
end subroutine computdifmag
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine computambip(u,ngrid,dx,dy,dz,dt,bemfx,bemfy,bemfz,florentzx,florentzy,florentzz,fluxad,bmagij,emfambdiff,fluxambdiff)

   use amr_commons
   use amr_parameters
   use hydro_commons
   use nimhd_parameters
   use const
   implicit none

   ! inputs
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::u 
   integer ::ngrid
   real(dp)::dx,dy,dz,dt
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::bemfx,bemfy,bemfz
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::florentzx,florentzy,florentzz
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::fluxad
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:3)::bmagij

   ! outputs
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::emfambdiff
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::fluxambdiff

   ! declare local variables
   INTEGER ::i, j, k, l

   real(dp)::v1x,v1y,v1z,v2x,v2y,v2z
   real(dp)::rhox,rhoy,rhoz,rhofx,rhofy,rhofz
   real(dp)::bsquarex,bsquarey,bsquarez,bsquare
   real(dp)::bsquarexx,bsquareyy,bsquarezz
   real(dp)::betaad2,betaadbricolo
   real(dp)::rhocell,bcell,tcell,ionisrate
   real(dp)::crossprodx,crossprody,crossprodz

   ! do NOT change value below Variation of betaad
   ! to avoid too small time step allowed
   !ntest=0

   !dtlim=dt!*coefalfven
   !dt est deja dtnew, qui a été choisi comme le dt normal (avec la condition de courant) ou le dt normal seuillé si le dtAD est trop faible(bricolo)

   do k=min(1,ku1+1),max(1,ku2-1)
      do j=min(1,ju1+1),max(1,ju2-1)
         do i=min(1,iu1+1),max(1,iu2-1)
           
            do l = 1, ngrid

               rhox=0.25d0*(u(l,i,j,k,1)+u(l,i,j-1,k,1)+u(l,i,j,k-1,1)+u(l,i,j-1,k-1,1))
               rhoy=0.25d0*(u(l,i,j,k,1)+u(l,i-1,j,k,1)+u(l,i,j,k-1,1)+u(l,i-1,j,k-1,1))
               rhoz=0.25d0*(u(l,i,j,k,1)+u(l,i-1,j,k,1)+u(l,i,j-1,k,1)+u(l,i-1,j-1,k,1))

               rhofx=0.5d0*(u(l,i,j,k,1)+u(l,i-1,j,k,1))
               rhofy=0.5d0*(u(l,i,j,k,1)+u(l,i,j-1,k,1))
               rhofz=0.5d0*(u(l,i,j,k,1)+u(l,i,j,k-1,1))

               rhocell = min(rhox,rhoy,rhoz,rhofx,rhofy,rhofz)

               ! Compute gas temperature in cgs
               if(ntestDADM.eq.1) then
                  tcell=1.0d0
               else
                  call ideal_gas_temperature(u(l,i,j,k,1), u(l,i,j,k,5), tcell)
               end if

               bsquarex=bemfx(l,i,j,k,1)**2+bemfx(l,i,j,k,2)**2+bemfx(l,i,j,k,3)**2
               bsquarey=bemfy(l,i,j,k,1)**2+bemfy(l,i,j,k,2)**2+bemfy(l,i,j,k,3)**2
               bsquarez=bemfz(l,i,j,k,1)**2+bemfz(l,i,j,k,2)**2+bemfz(l,i,j,k,3)**2

               bsquarexx=bmagij(l,i,j,k,1,1)**2+bmagij(l,i,j,k,2,1)**2+bmagij(l,i,j,k,3,1)**2
               bsquareyy=bmagij(l,i,j,k,1,2)**2+bmagij(l,i,j,k,2,2)**2+bmagij(l,i,j,k,3,2)**2
               bsquarezz=bmagij(l,i,j,k,1,3)**2+bmagij(l,i,j,k,2,3)**2+bmagij(l,i,j,k,3,3)**2

               bcell = max(bsquarex,bsquarey,bsquarez,bsquarexx,bsquareyy,bsquarezz)

               ionisrate=default_ionisrate

               ! EMF x

               v1x=florentzx(l,i,j,k,1)
               v1y=florentzx(l,i,j,k,2)
               v1z=florentzx(l,i,j,k,3)
               v2x=bemfx(l,i,j,k,1)
               v2y=bemfx(l,i,j,k,2)
               v2z=bemfx(l,i,j,k,3)
               emfambdiff(l,i,j,k,1)=crossprodx(v1x,v1y,v1z,v2x,v2y,v2z)

               rhox=0.25d0*(u(l,i,j,k,1)+u(l,i,j-1,k,1)+u(l,i,j,k-1,1)+u(l,i,j-1,k-1,1))
               betaad2=betaadbricolo(rhocell,rhox,dt,bcell,bcell,dx,tcell,ionisrate)
               emfambdiff(l,i,j,k,1)=emfambdiff(l,i,j,k,1)*betaad2 

               ! EMF y

               v1x=florentzy(l,i,j,k,1)
               v1y=florentzy(l,i,j,k,2)
               v1z=florentzy(l,i,j,k,3)
               v2x=bemfy(l,i,j,k,1)
               v2y=bemfy(l,i,j,k,2)
               v2z=bemfy(l,i,j,k,3)
               emfambdiff(l,i,j,k,2)=crossprody(v1x,v1y,v1z,v2x,v2y,v2z)

               rhoy=0.25d0*(u(l,i,j,k,1)+u(l,i-1,j,k,1)+u(l,i,j,k-1,1)+u(l,i-1,j,k-1,1))
               betaad2=betaadbricolo(rhocell,rhoy,dt,bcell,bcell,dx,tcell,ionisrate)
               emfambdiff(l,i,j,k,2)=emfambdiff(l,i,j,k,2)*betaad2            
                    
               ! EMF z

               v1x=florentzz(l,i,j,k,1)
               v1y=florentzz(l,i,j,k,2)
               v1z=florentzz(l,i,j,k,3)
               v2x=bemfz(l,i,j,k,1)
               v2y=bemfz(l,i,j,k,2)
               v2z=bemfz(l,i,j,k,3)
               emfambdiff(l,i,j,k,3)=crossprodz(v1x,v1y,v1z,v2x,v2y,v2z)

               rhoz=0.25d0*(u(l,i,j,k,1)+u(l,i-1,j,k,1)+u(l,i,j-1,k,1)+u(l,i-1,j-1,k,1))
               betaad2=betaadbricolo(rhocell,rhoz,dt,bcell,bcell,dx,tcell,ionisrate)
               emfambdiff(l,i,j,k,3)=emfambdiff(l,i,j,k,3)*betaad2

               ! energy flux on faces

               rhofx=0.5d0*(u(l,i,j,k,1)+u(l,i-1,j,k,1))
               betaad2=betaadbricolo(rhocell,rhofx,dt,bcell,bcell,dx,tcell,ionisrate)
               fluxambdiff(l,i,j,k,1)=-betaad2*fluxad(l,i,j,k,1)

               rhofy=0.5d0*(u(l,i,j,k,1)+u(l,i,j-1,k,1))
               betaad2=betaadbricolo(rhocell,rhofy,dt,bcell,bcell,dx,tcell,ionisrate) !TC:dy?
               fluxambdiff(l,i,j,k,2)=-betaad2*fluxad(l,i,j,k,2)

               rhofz=0.5d0*(u(l,i,j,k,1)+u(l,i,j,k-1,1))
               betaad2=betaadbricolo(rhocell,rhofz,dt,bcell,bcell,dx,tcell,ionisrate) !TC: dz?
               fluxambdiff(l,i,j,k,3)=-betaad2*fluxad(l,i,j,k,3)

            end do
         end do
      end do
   end do

end subroutine computambip
!###########################################################
!###########################################################
!###########################################################
!###########################################################

! VECTOR FUNCTION

!###########################################################
double precision function computdxbis(vec,n2,l,i,j,k,dx)

   use amr_parameters,only:dp,nvector
   use hydro_parameters,only:iu1,iu2,ju1,ju2,ku1,ku2
   implicit none 
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::vec
   real(dp)::dx
   integer::n2,l,i,j,k

   computdxbis = (vec(l,i+1,j,k,n2) - vec(l,i,j,k,n2)) / dx

end function computdxbis

double precision  function computdybis(vec,n2,l,i,j,k,dx)

   use amr_parameters,only:dp,nvector
   use hydro_parameters,only:iu1,iu2,ju1,ju2,ku1,ku2
   implicit none 
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::vec
   real(dp)::dx
   integer::n2,l,i,j,k

   computdybis=(vec(l,i,j+1,k,n2) - vec(l,i,j,k,n2)) / dx

end function computdybis

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

! NIMHD COEFFICIENTS

!###########################################################
!###########################################################
!###########################################################
double precision function gammaadbis(rhon,BBcell,BBcellold,temper,ionisrate)

   use hydro_parameters
   use amr_parameters,only:mu_gas
   use nimhd_parameters
   use constants
   implicit none

   real(dp)::rhon,rhoH,n_H_max,BBcell,temper,BBcellold
   real(dp)::eta_AD_chimie,ionisrate

   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

   ! function which computes the coefficient gamma which
   ! appears in ambipolar diffusion dB/dt=1/(gamma*rhoi*rhon)curl*(j*B)*B)+...
   ! see Duffin & Pudritz 2008, astro-ph 08/10/08 eq (6)
   ! WARNING no mu_0 needed here

   n_H_max = 2.5d+17 !H/cc
   ! TC: Why this value?

   ! C shock Duffin et Pudritz
   ! gammaadbis in CGS
   !gammaadbis=gammaAD

   ! TC: what is xmolaire? Why is it constant to 2? Because molecular hydrogen?
   !rhoH=rhon*xmolaire*H2_fraction*scale_d/(mu_gas*mH) ! convert in H/cc
   rhoH=rhon*2.0d0*H2_fraction*scale_d/(mu_gas*mH) ! convert in H/cc
   rhoH = max(rhoH, n_H_max)

   gammaadbis=eta_AD_chimie(rhoH,BBcell,BBcellold,temper,ionisrate)

   gammaadbis=gammaadbis*scale_t*scale_d ! in code units

end function gammaadbis
!###########################################################
!###########################################################
!###########################################################
! TC: eta_AD_chimie and eta_ohm_chimie contain a lot of duplicate code
double precision function eta_AD_chimie(rhon,BBcell,BBcellold,temper,ionisrate)

   use hydro_commons
   use constants
   use nimhd_commons
   use nimhd_parameters
   implicit none

   real(dp):: sigO,sigH,sigP,densionbis,BBcgs, bbcell,BBcellold
   real(dp)::inp,ll,rhon,ii,temper,lb,j_dp,xx
   integer :: i,j,k,ib
   real(dp)::ionisrate

   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

   ! TC: extrapolate from table[density,temperature,magnetic field]
   if(use_x2d==1)then
      call sig_x2d(rhon,temper,BBcellold,sigO,sigH,sigP) 
      !inp=rhon/xmolaire/H2_fraction     ! inp is neutrals.cc, to fit densionbis
      inp=rhon/2d0/H2_fraction     ! inp is neutrals.cc, to fit densionbis
      eta_AD_chimie=(sigO/(sigO**2+sigH**2)-1d0/sigP)   ! resistivity in s

      BBcgs=sqrt(BBcell*(4d0*pi*scale_d*(scale_v)**2))
      eta_AD_chimie=BBcgs**2/(eta_AD_chimie*densionbis(inp)*inp*scale_d*scale_d*c_cgs**2)  ! need B in G, output is gammaad in cgs

   ! TC: extrapolate from table[density,temperature,ionisation rate,magnetic field]
   else if(use_x3d==1)then
      ll=(1d0+(log10(rhon)-log10(nminchimie))/dnchimie)
      j=floor(ll)
      ii = max(0d0, (log10(temper)-log10(tminchimie))/dtchimie)
      i=floor(1d0+ii)
      xx=(1d0+(log10(ionisrate)-log10(ximinchimie))/dxichimie)
      k=floor(xx)
      BBcgs=sqrt(BBcellold*(4d0*pi*scale_d*(scale_v)**2))
      lb=(1d0+(log10(BBcgs)-log10(bminchimie))/dbchimie)
      ib=floor(lb)

      call sig_x3d(ll,ii,xx,j,i,k,lb,ib,sigO,sigH,sigP) 
      inp=rhon/2d0/H2_fraction     ! inp is neutrals.cc, to fit densionbis
      eta_AD_chimie=(sigO/(sigO**2+sigH**2)-1d0/sigP)   ! resistivity in s

      BBcgs=sqrt(BBcell*(4d0*pi*scale_d*(scale_v)**2))

      eta_AD_chimie=BBcgs**2/(eta_AD_chimie*densionbis(inp)*inp*scale_d*scale_d*c_cgs**2)  ! need B in G, output is gammaad in cgs

   endif

end function eta_AD_chimie
!###########################################################
!###########################################################
!###########################################################
double precision function eta_ohm_chimie(rhon,BBcell,temper,ionisrate)

   use hydro_commons
   use constants
   use nimhd_commons
   use nimhd_parameters
   implicit none

   real(dp) :: inp,ll,ii,lb,rhon,BBcell
   real(dp) :: temper,sigO,sigH,sigP,BBcgs
   real(dp) :: j_dp,ionisrate,xx
   integer  :: j,i,ib,k

   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

   ! TC: extrapolate from table[density,temperature,magnetic field] ?
   if(use_x2d==1)then
      call sig_x2d(rhon,temper,BBcell,sigO,sigH,sigP) 
      eta_ohm_chimie = (1d0 / sigP) * c_cgs * c_cgs / (4.0_dp*pi)
      eta_ohm_chimie = max(eta_ohm_chimie * (1.0d0-tanh(rhon/1.0d15)), 1d-36)

   else if(use_x3d==1)then
      ll=(1d0+(log10(rhon)-log10(nminchimie))/dnchimie)
      j=floor(ll)
      ii=(1d0+(log10(temper)-log10(tminchimie))/dtchimie)
      i=floor(ii)   
      xx=(1d0+(log10(ionisrate)-log10(ximinchimie))/dxichimie)
      k=floor(xx)
      BBcgs=sqrt(BBcell*(4d0*pi*scale_d*(scale_v)**2))
      lb=(1d0+(log10(BBcgs)-log10(bminchimie))/dbchimie)
      ib=floor(lb)

      call sig_x3d(ll,ii,xx,j,i,k,lb,ib,sigO,sigH,sigP)
      eta_ohm_chimie = (1d0 / sigP) * c_cgs * c_cgs / (4.0_dp*pi)
      ! TC: why don't we do the adhoc thing here?
   endif

   ! Ad-hoc modification to ensure that the ohmic resistivity falls to zero when the density exceeds 1.0e15
   ! when alkali metals are ionized.
   !eta_ohm_chimie = eta_ohm_chimie * (1.0d0-tanh(rhon/1.0d15))
   ! eta_ohm_chimie = max(eta_ohm_chimie * (1.0d0-tanh(rhon/1.0d15)), 1d-36)

end function eta_ohm_chimie
!###########################################################
!###########################################################
!###########################################################
double precision function densionbis(rhon)

   use nimhd_parameters, only : coefionis,default_ionisrate,ntestDADM,dp

   implicit none 
   real(dp)::rhon
   real(dp)::xn, rhoncgs

   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

   ! density of neutral in g/cm3  
   rhoncgs=rhon*scale_d

   ! function which computes the density in g/cm3 of ions 
   ! see Duffin & Pudritz 2008, astro-ph 08/10/08 eq (14)

   ! Mellon & Li 2009 (?) or Hennebelle & Teyssier 2007
   ! WARNING 3d-16 si in cgs
   densionbis=coefionis*sqrt(rhoncgs*default_ionisrate/1.0d-17)

   !!!!!!!!!!!!!!!!!!!!!!!!!
   ! densionbis in USER UNITS !!!
   !!!!!!!!!!!!!!!!!!!!!!!!!

   ! transformation coefionis in user units
   !densionbis=densionbis/sqrt(scale_d)

   ! back in user units
   densionbis=densionbis/scale_d

   ! test C shock Duffin et Pudritz
   if(ntestDADM.eq.1) then
      densionbis=1d0
   endif

end function densionbis
!###########################################################
!###########################################################
!###########################################################
double precision function etaohmdiss(rhon,BBcell,temper,ionisrate)

   use hydro_commons
   use amr_parameters,only:mu_gas
   use nimhd_parameters
   use constants

   implicit none 
   real(dp) ::rhon,xpressure,rhoH,rhotemp,BBcell
   real(dp)::gammaadbis,densionbis
   real(dp)::xionisation,temper,scale_p,xpcgs,rhocgs,xnbcgs,n_H_max
   real(dp)::eta_ohm_chimie,ionisrate

   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

   if(ntestDADM.eq.0) then
      ! function which computes the coefficient eta which
      ! appears in ohmic dissipation dB/dt=-curl(eta*curl(B))+...
      ! see Machida, Inutsuka, Matsumoto, ApJ, 670,1198-1213, 2007

      n_H_max = 2.5d+17

      ! convert to CGS

      ! scale_p = scale_d*(scale_v**2.)
      ! xpcgs=xpressure*scale_p
      ! rhocgs=rhon*scale_d
      ! ! nb per cm3
      ! xnbcgs=rhocgs/xmneutre
      ! ! temperature in cgs
      ! temper=xpcgs*xmolaire/(rhocgs*rperfectgaz)
      ! !write(*,*)'temper',temper
      ! 
      ! ! degree of ionisation
      ! ! Machida et al 2007 
      ! xionisation=5.7d-4/(xnbcgs)
      ! ! Shu 1987 27.7 p 363 and p 361 m_n=2.33 m_i=29
      ! !xionisation=2.33d0/29d0*densionbis(rhon)/rhon
      ! 
      ! ! Machida et al 2007 : etaMD=740
      ! !etaohmdiss=etaMD*sqrt(temper/10d0)*(1d0-tanh(xnbcgs/1d15))/xionisation
      ! ! Dapp & Basu 2010
      ! ! etaohmdiss=etaMD*1.3d18*(xnbcgs/1d12)*sqrt(temper/10d0)*(1d0-tanh(xnbcgs/1d15))
      ! ! back to user units
      ! !print*, etaohmdiss
      ! !stop

      !rhoH=rhon*xmolaire*H2_fraction*scale_d/(mu_gas*mH) ! convert in H/cc
      rhoH=rhon*2.0d0*H2_fraction*scale_d/(mu_gas*mH) ! convert in H/cc

      rhotemp = MAX(rhoH,rho_threshold)

      if(rhotemp < n_H_max)then
         etaohmdiss=eta_ohm_chimie(rhotemp,BBcell,temper,ionisrate)
      else
         etaohmdiss=eta_ohm_chimie(n_H_max,BBcell,temper,ionisrate)
      endif
      
      etaohmdiss=etaohmdiss*scale_t/(scale_l)**2

   elseif(ntestDADM.eq.1) then
      ! test Alfven Lessaffre
      !etaohmdiss=2d-2

      ! test heat diffusion
      !etaohmdiss=1d0

      ! test oblique shock
      !etaohmdiss=0.15d0

      etaohmdiss=etaMD
   endif
 
end function etaohmdiss
!###########################################################
!###########################################################
!###########################################################
double precision function etaohmdissbricolo(rhon,BBcell,temper,dt,dx,ionisrate)

   use hydro_commons
   use amr_parameters,only:mu_gas
   use nimhd_parameters
   use constants

   implicit none 
   real(dp) ::rhon,xpressure,rhoH,rhotemp,BBcell
   real(dp)::gammaadbis,densionbis,ionisrate
   real(dp)::xionisation,temper,scale_p,xpcgs,rhocgs,xnbcgs,n_H_max
   real(dp)::eta_ohm_chimie,dx,dt,xx,dtt

   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

   ! function which computes the coefficient eta which
   ! appears in ohmic dissipation dB/dt=-curl(eta*curl(B))+...
   ! see Machida, Inutsuka, Matsumoto, ApJ, 670,1198-1213, 2007

   if(ntestDADM.eq.0) then

      n_H_max = 2.5d+17

      ! convert to CGS

      ! scale_p = scale_d*(scale_v**2.)
      ! xpcgs=xpressure*scale_p
      ! rhocgs=rhon*scale_d
      ! ! nb per cm3
      ! xnbcgs=rhocgs/xmneutre
      ! ! temperature in cgs
      ! temper=xpcgs*xmolaire/(rhocgs*rperfectgaz)
      ! !write(*,*)'temper',temper
      ! 
      ! ! degree of ionisation
      ! ! Machida et al 2007 
      ! xionisation=5.7d-4/(xnbcgs)
      ! ! Shu 1987 27.7 p 363 and p 361 m_n=2.33 m_i=29
      ! !xionisation=2.33d0/29d0*densionbis(rhon)/rhon
      ! 
      ! ! Machida et al 2007 : etaMD=740
      ! !etaohmdiss=etaMD*sqrt(temper/10d0)*(1d0-tanh(xnbcgs/1d15))/xionisation
      ! ! Dapp & Basu 2010
      ! ! etaohmdiss=etaMD*1.3d18*(xnbcgs/1d12)*sqrt(temper/10d0)*(1d0-tanh(xnbcgs/1d15))
      ! ! back to user units
      ! !print*, etaohmdiss
      ! !stop

      !rhoH=rhon*xmolaire*H2_fraction*scale_d/(mu_gas*mH) ! convert in H/cc
      rhoH=rhon*2.0d0*H2_fraction*scale_d/(mu_gas*mH) ! convert in H/cc

      rhotemp = MAX(rhoH,rho_threshold)

      if(rhotemp < n_H_max)then
         etaohmdissbricolo=eta_ohm_chimie(rhotemp,BBcell,temper,ionisrate)
      else
         etaohmdissbricolo=eta_ohm_chimie(n_H_max,BBcell,temper,ionisrate)
      endif

      etaohmdissbricolo=etaohmdissbricolo*scale_t/(scale_l)**2

      ! TC: this subroutine is exactly the same as etaohmdiss except for the following part...
      ! robbery to avoid too small time step
      if(nminitimestep) then
         if(dt.ne.0d0) then
            xx=etaohmdissbricolo
            if(xx.ne.0d0) then
            dtt=coefohm*dx*dx/xx   !dtohm pour la cellule
            else
               dtt=1d39
            endif
            if (dtt.le.dt) then
               etaohmdissbricolo=coefohm*dx*dx/(dt)
            endif
         endif
      endif

   elseif(ntestDADM.eq.1) then
      ! test Alfven Lessaffre
      !etaohmdiss=2d-2

      ! test heat diffusion
      !etaohmdiss=1d0

      ! test oblique shock
      !etaohmdiss=0.15d0

      etaohmdissbricolo=etaMD
   endif

end function etaohmdissbricolo
!###########################################################
!###########################################################
!###########################################################
double precision function betaad(rhon,bsquare,temper,ionisrate)

   use hydro_parameters
   use nimhd_parameters
   implicit none

   ! rhon in code units
   real(dp)::rhon,rhotemp,bsquare,temper
   real(dp)::gammaadbis,densionbis,ionisrate
   real(dp)::xx

   ! function which computes the coefficient beta which
   ! appears in ambipolar diffusion dB/dt=curl(gamma(j*B)*B)+...
   ! see Duffin & Pudritz 2008, astro-ph 08/10/08 eq (5)
   ! WARNING no mu_0 needed here because F_Lorentz used

   if(ntestDADM.eq.0) then

      ! Warning gammaadbis and densionbis already in code units
      rhotemp = MAX(rhon,rho_threshold)
      ! TC: rho_threshold is a units disaster waiting to happen

      xx=gammaadbis(rhotemp,bsquare,bsquare,temper,ionisrate)*densionbis(rhotemp)*rhotemp
      if(xx.ne.0d0) then
         betaad=1d0/xx 
      else
         betaad=1d39
         if(rhotemp < 1.0d+14)then
         !TC: hard coded value in code units (not good)
            write(*,*)'WARNING gammaadbis(rhotemp,bsquare,temper,ionisrate)*densionbis(rhon)*rhon equal 0',gammaadbis(rhotemp,bsquare,bsquare,temper,ionisrate),densionbis(rhotemp),rhotemp
         endif
      endif

   elseif(ntestDADM.eq.1) then
      ! test Barenblatt
      !betaad=1d0
      ! test C shock
      betaad=1d0/(gammaAD*rhon)
   endif

   !rhon, gammaadbis(rhon) and densionbis(rhon) already in user units
   !!betaad=betaad/scale_d

end function betaad
!###########################################################
!###########################################################
!###########################################################
double precision function betaadbricolo(rhocelln,rhon,dt,bsquare,bsquareold,dx,temper,ionisrate)

   use hydro_parameters
   use amr_commons
   use cooling_module
   use nimhd_parameters
   use constants

   implicit none
   real(dp) ::rhocelln,rhon,betaadbricolotemp,dt,bsquare,bsquareold,dx,temper
   real(dp)::gammaadbis,densionbis,rhotemp,rhotemp_cell,ionisrate
   real(dp)::xx,dtt,bbcgs

   ! function which computes the coefficient beta which
   ! appears in ambipolar diffusion dB/dt=curl(gamma(j*B)*B)+...
   ! see Duffin & Pudritz 2008, astro-ph 08/10/08 eq (5)
   ! WARNING no mu_0 needed here because F_Lorentz used

   if(ntestDADM.eq.0) then

      rhotemp = MAX(rhon,rho_threshold)
      rhotemp_cell = MAX(rhocelln,rho_threshold)

      xx=gammaadbis(rhotemp_cell,bsquare,bsquareold,temper,ionisrate)*densionbis(rhotemp_cell)*rhotemp_cell  ! dans la cellule
      ! gammaadbis and densionbis already in user units

      if(xx.ne.0d0) then
         betaadbricolo=1d0/xx 
      else
         betaadbricolo=1d39
         if(rhotemp < 1.0d+14)then
            write(*,*)'WARNING gammaadbis(rhocelln,bsquare,bsquareold,temper,ionisrate)*densionbis(rhocelln)*rhocelln in the cell equals 0',gammaadbis(rhotemp_cell,bsquare,bsquareold,temper,ionisrate),densionbis(rhocelln),rhocelln,bsquare,bsquareold,temper
         endif
      endif

      !xx=gammaadbis(rhotemp,bsquare,bsquareold,temper)*densionbis(rhon)*rhon   ! a l'interface : cote ou coin selon les cas. A utiliser si l'on est pas dans un cas seuille
      xx=gammaadbis(rhotemp,bsquare,bsquareold,temper,ionisrate)*densionbis(rhotemp)*rhotemp  

      ! a l'interface : cote ou coin selon les cas. A utiliser si l'on est pas dans un cas seuille

      if(xx.ne.0d0) then
         betaadbricolotemp=1d0/xx 
      else
         betaadbricolotemp=1d39
         if(rhotemp < 1.0d+14)then
            write(*,*)'WARNING gammaadbis(rhon,bsquare,bsquareold,temper,ionisrate)*densionbis(rhon)*rhon at the interface equals 0',gammaadbis(rhotemp,bsquare,bsquareold,temper,ionisrate),densionbis(rhon),rhon
         endif
      endif

      ! robbery to avoid too small time step
      if(nminitimestep) then
         if(dt.ne.0d0) then
            xx=bsquare*betaadbricolo
            if(xx.ne.0d0) then
               dtt=coefad*dx*dx/xx   !dtAD pour la cellule
            else
               dtt=1d39
            endif
            if (dtt.le.dt) then   ! on compare bien dtAD calcule pour la cellule (rhocelln) avec le temps de la simu
               betaadbricolo=coefad*dx*dx/(dt*bsquare)
         !      write(*,*) 'la où ça seuille rho et B valent : ', rhocelln, bsquare
               !ici dtlim est le dt le plus petit : normal, ou seuillé si besoin est.
            else
               betaadbricolo=betaadbricolotemp  ! le betaad normal calcule avec rho a l'interface
            endif
         endif
      endif

   elseif(ntestDADM.eq.1) then
      ! test Barenblatt
      !betaadbricolo=1d0
      ! test C shock
      betaadbricolo=1d0/(gammaAD*rhon)
   endif

end function betaadbricolo
!###########################################################
!###########################################################
!###########################################################