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
   real(dp)::rhof,bsqf,epsf,tcellf
   real(dp)::etaod2,etaohmdiss
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

               call ideal_gas_temperature(rhox, epsx, tcellx)
               call ideal_gas_temperature(rhoy, epsy, tcelly)
               call ideal_gas_temperature(rhoz, epsz, tcellz)

               etaod2x=etaohmdiss(rhox,bsquarex,tcellx,dt,dx,.true.)
               etaod2y=etaohmdiss(rhoy,bsquarey,tcelly,dt,dx,.true.)
               etaod2z=etaohmdiss(rhoz,bsquarez,tcellz,dt,dx,.true.)
               ! TC: shouldn't dy and dz be used here in principle? (in practice they are the same)  

               ! WARNING dB/dt=-curl(eta*J)
               emfohmdiss(l,i,j,k,nxx)=-etaod2x*jemfx(l,i,j,k,1)
               emfohmdiss(l,i,j,k,nyy)=-etaod2y*jemfy(l,i,j,k,2)
               emfohmdiss(l,i,j,k,nzz)=-etaod2z*jemfz(l,i,j,k,3)

               do h = 1,3
                  rhof=0.5d0*(u(l,i,j,k,1)+u(l,i-index_i(h),j-index_j(h),k-index_k(h),1))
                  epsf=0.5d0*(u(l,i,j,k,5)+u(l,i-index_i(h),j-index_j(h),k-index_k(h),5))
                  bsqf=bmagij(l,i,j,k,1,h)**2+bmagij(l,i,j,k,2,h)**2+bmagij(l,i,j,k,3,h)**2

                  ! Compute gas temperature in cgs
                  call ideal_gas_temperature(rhof, epsf, tcellf)
                     
                  etaod2=etaohmdiss(rhof,bsqf,tcellf,0d0,0d0,.false.)
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
   real(dp)::rhocell,bcell,tcell
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
               call ideal_gas_temperature(u(l,i,j,k,1), u(l,i,j,k,5), tcell)

               bsquarex=bemfx(l,i,j,k,1)**2+bemfx(l,i,j,k,2)**2+bemfx(l,i,j,k,3)**2
               bsquarey=bemfy(l,i,j,k,1)**2+bemfy(l,i,j,k,2)**2+bemfy(l,i,j,k,3)**2
               bsquarez=bemfz(l,i,j,k,1)**2+bemfz(l,i,j,k,2)**2+bemfz(l,i,j,k,3)**2

               bsquarexx=bmagij(l,i,j,k,1,1)**2+bmagij(l,i,j,k,2,1)**2+bmagij(l,i,j,k,3,1)**2
               bsquareyy=bmagij(l,i,j,k,1,2)**2+bmagij(l,i,j,k,2,2)**2+bmagij(l,i,j,k,3,2)**2
               bsquarezz=bmagij(l,i,j,k,1,3)**2+bmagij(l,i,j,k,2,3)**2+bmagij(l,i,j,k,3,3)**2

               bcell = max(bsquarex,bsquarey,bsquarez,bsquarexx,bsquareyy,bsquarezz)

               ! EMF x

               v1x=florentzx(l,i,j,k,1)
               v1y=florentzx(l,i,j,k,2)
               v1z=florentzx(l,i,j,k,3)
               v2x=bemfx(l,i,j,k,1)
               v2y=bemfx(l,i,j,k,2)
               v2z=bemfx(l,i,j,k,3)
               emfambdiff(l,i,j,k,1)=crossprodx(v1x,v1y,v1z,v2x,v2y,v2z)

               rhox=0.25d0*(u(l,i,j,k,1)+u(l,i,j-1,k,1)+u(l,i,j,k-1,1)+u(l,i,j-1,k-1,1))
               betaad2=betaadbricolo(rhocell,rhox,dt,bcell,bcell,dx,tcell,.true.)
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
               betaad2=betaadbricolo(rhocell,rhoy,dt,bcell,bcell,dx,tcell,.true.)
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
               betaad2=betaadbricolo(rhocell,rhoz,dt,bcell,bcell,dx,tcell,.true.)
               emfambdiff(l,i,j,k,3)=emfambdiff(l,i,j,k,3)*betaad2

               ! energy flux on faces

               rhofx=0.5d0*(u(l,i,j,k,1)+u(l,i-1,j,k,1))
               betaad2=betaadbricolo(rhocell,rhofx,dt,bcell,bcell,dx,tcell,.true.)
               fluxambdiff(l,i,j,k,1)=-betaad2*fluxad(l,i,j,k,1)

               rhofy=0.5d0*(u(l,i,j,k,1)+u(l,i,j-1,k,1))
               betaad2=betaadbricolo(rhocell,rhofy,dt,bcell,bcell,dx,tcell,.true.) !TC:dy?
               fluxambdiff(l,i,j,k,2)=-betaad2*fluxad(l,i,j,k,2)

               rhofz=0.5d0*(u(l,i,j,k,1)+u(l,i,j,k-1,1))
               betaad2=betaadbricolo(rhocell,rhofz,dt,bcell,bcell,dx,tcell,.true.) !TC: dz?
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
double precision function betaadbricolo(rhocelln,rhon,dt,bsquare,bsquareold,dx,temper,limit)

   use hydro_parameters
   use amr_commons
   use cooling_module
   use nimhd_parameters
   use constants

   implicit none
   real(dp) ::rhocelln,rhon,betaadbricolotemp,dt,bsquare,bsquareold,dx,temper
   real(dp)::gammaadbis,densionbis,rhotemp,rhotemp_cell
   real(dp)::xx,dtt,bbcgs
   logical::limit

   ! function which computes the coefficient beta which
   ! appears in ambipolar diffusion dB/dt=curl(gamma(j*B)*B)+...
   ! see Duffin & Pudritz 2008, astro-ph 08/10/08 eq (5)
   ! WARNING no mu_0 needed here because F_Lorentz used

   if(resistivity_method==0)then
      ! fixed resistivity
      betaadbricolo=1d0/(gammaAD*rhon)
   elseif(resistivity_method==1)then
      ! *** put your analytic resistivity here ***
      !analytical model resitivity(rho,T), Shu?
      gammaAD = 1
      betaadbricolo=1d0/(gammaAD*rhon)
   else
      ! table
      rhotemp = MAX(rhon,rho_threshold)
      rhotemp_cell = MAX(rhocelln,rho_threshold)

      xx=gammaadbis(rhotemp_cell,bsquare,bsquareold,temper)*densionbis(rhotemp_cell)*rhotemp_cell  ! dans la cellule
      ! gammaadbis and densionbis already in user units

      if(xx.ne.0d0) then
         betaadbricolo=1d0/xx 
      else
         betaadbricolo=1d39
         if(rhotemp < 1.0d+14)then
            write(*,*)'WARNING gammaadbis(rhocelln,bsquare,bsquareold,temper)*densionbis(rhocelln)*rhocelln in the cell equals 0',gammaadbis(rhotemp_cell,bsquare,bsquareold,temper),densionbis(rhocelln),rhocelln,bsquare,bsquareold,temper
         endif
      endif

      !xx=gammaadbis(rhotemp,bsquare,bsquareold,temper)*densionbis(rhon)*rhon   ! a l'interface : cote ou coin selon les cas. A utiliser si l'on est pas dans un cas seuille
      xx=gammaadbis(rhotemp,bsquare,bsquareold,temper)*densionbis(rhotemp)*rhotemp  

      ! a l'interface : cote ou coin selon les cas. A utiliser si l'on est pas dans un cas seuille

      if(xx.ne.0d0) then
         betaadbricolotemp=1d0/xx 
      else
         betaadbricolotemp=1d39
         if(rhotemp < 1.0d+14)then
            write(*,*)'WARNING gammaadbis(rhon,bsquare,bsquareold,temper)*densionbis(rhon)*rhon at the interface equals 0',gammaadbis(rhotemp,bsquare,bsquareold,temper),densionbis(rhon),rhon
         endif
      endif

      ! if the timestep has been limited, the resistivity needs to be adjusted
      if(limit.and.nminitimestep) then
         if(dt.ne.0d0) then
            ! recalculate the ambipolar diffusion timestep for the current cell
            xx=bsquare*betaadbricolo
            if(xx.ne.0d0) then
               dtt=coefad*dx*dx/xx
            else
               dtt=1d39
            endif
            ! check whether it is smaller than the global timestep that has been determined
            if (dtt.le.dt) then
               ! if so, adjust the resistivity to match the timestep
               betaadbricolo=coefad*dx*dx/(dt*bsquare)
            else
               betaadbricolo=betaadbricolotemp
            endif
         endif
      endif
   endif

end function betaadbricolo
!###########################################################
!###########################################################
!###########################################################
double precision function betaad(rhon,bsquare,temper)

   use hydro_parameters
   use nimhd_parameters
   implicit none

   ! rhon in code units
   real(dp)::rhon,rhotemp,bsquare,temper
   real(dp)::gammaadbis,densionbis
   real(dp)::xx

   ! function which computes the coefficient beta which
   ! appears in ambipolar diffusion dB/dt=curl(gamma(j*B)*B)+...
   ! see Duffin & Pudritz 2008, astro-ph 08/10/08 eq (5)
   ! WARNING no mu_0 needed here because F_Lorentz used

   if(resistivity_method==0)then
      ! fixed resistivity
      betaad=1d0/(gammaAD*rhon)
   elseif(resistivity_method==1)then
      ! *** put your analytic resistivity here ***
      !analytical model resitivity(rho,T), Shu?
      gammaAD = 1
      betaad=1d0/(gammaAD*rhon)
   else
      ! table

      ! Warning gammaadbis and densionbis already in code units
      rhotemp = MAX(rhon,rho_threshold)
      ! TC: rho_threshold is a units disaster waiting to happen

      xx=gammaadbis(rhotemp,bsquare,bsquare,temper)*densionbis(rhotemp)*rhotemp
      if(xx.ne.0d0) then
         betaad=1d0/xx 
      else
         betaad=1d39
         if(rhotemp < 1.0d+14)then
         !TC: hard coded value in code units (not good)
            write(*,*)'WARNING gammaadbis(rhotemp,bsquare,temper)*densionbis(rhon)*rhon equal 0',gammaadbis(rhotemp,bsquare,bsquare,temper),densionbis(rhotemp),rhotemp
         endif
      endif
   endif

end function betaad
!###########################################################
!###########################################################
!###########################################################
double precision function gammaadbis(rhon,BBcell,BBcellold,temper)

   use hydro_parameters
   use amr_parameters,only:mu_gas
   use nimhd_parameters
   use resistivity_table
   use constants
   implicit none

   real(dp)::rhon,rhoH,n_H_max,BBcell,temper,BBcellold
   real(dp)::eta_AD_chimie
   real(dp) :: n_H_max=2.5d+17     ! cm**-3

   real(dp):: sigO,sigH,sigP,densionbis,BBcgs
   real(dp)::inp

   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

   ! function which computes the coefficient gamma which
   ! appears in ambipolar diffusion dB/dt=1/(gamma*rhoi*rhon)curl*(j*B)*B)+...
   ! see Duffin & Pudritz 2008, astro-ph 08/10/08 eq (6)
   ! WARNING no mu_0 needed here

   rhoH=rhon*2.0d0*H2_fraction*scale_d/(mu_gas*mH) ! convert in H/cc
   rhoH = max(rhoH, n_H_max)

   ! TC: extrapolate from table[density,temperature,magnetic field]
   call interpolate_table(rhoH,temper,BBcellold,sigO,sigH,sigP) 
   !inp=rhoH/xmolaire/H2_fraction     ! inp is neutrals.cc, to fit densionbis
   inp=rhoH/2d0/H2_fraction     ! inp is neutrals.cc, to fit densionbis
   eta_AD_chimie=(sigO/(sigO**2+sigH**2)-1d0/sigP)   ! resistivity in s

   BBcgs=sqrt(BBcell*(4d0*pi*scale_d*(scale_v)**2))
   eta_AD_chimie=BBcgs**2/(eta_AD_chimie*densionbis(inp)*inp*scale_d*scale_d*c_cgs**2)  ! need B in G, output is gammaad in cgs

   gammaadbis=eta_AD_chimie*scale_t*scale_d ! in code units

end function gammaadbis
!###########################################################
!###########################################################
!###########################################################
double precision function densionbis(rhon)

   use amr_parameters, only : dp

   implicit none 
   real(dp)::rhon
   real(dp)::xn, rhoncgs

   ! Mellon & Li 2009 (?) or Hennebelle & Teyssier 2007
   real(dp):: coefionis=3d-16 !in cgs !remove ! coefionis*sqrt(n_H)=n_i , empirical value from Shu book 2, p. 363
   real(dp):: default_ionisrate=1d-17

   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

   ! density of neutral in g/cm3  
   rhoncgs=rhon*scale_d

   ! function which computes the density in g/cm3 of ions 
   ! see Duffin & Pudritz 2008, astro-ph 08/10/08 eq (14)

   densionbis=coefionis*sqrt(rhoncgs*default_ionisrate/1.0d-17)

   ! back in code units
   densionbis=densionbis/scale_d

end function densionbis
!###########################################################
!###########################################################
!###########################################################
double precision function etaohmdiss(rhon,BBcell,temper,dt,dx,limit)

   use amr_parameters,    only:dp,mu_gas
   use nimhd_parameters
   use resistivity_table, only:interpolate_table
   use constants,         only:c_cgs,pi,mH
   !-----------------------------------------------------------
   ! Function which computes the coefficient eta which appears
   ! in ohmic dissipation dB/dt=-curl(eta*curl(B))+...
   ! See Machida, Inutsuka, Matsumoto, ApJ, 670,1198-1213, 2007
   !-----------------------------------------------------------
   implicit none 
   real(dp) :: rhon,BBcell,temper  ! input cell variables
   real(dp) :: dx,dt               ! input cell size and simulation time step
   logical  :: limit               ! take into account limitation of timestep or not
   real(dp) :: rhoH,n_H_max,BBcgs
   real(dp) :: n_H_max=2.5d+17     ! cm**-3
   real(dp) :: sigO,sigH,sigP,eta_ohm_chimie
   real(dp) :: dtt
   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v

   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

   if(resistivity_method==0) then ! fixed value
      etaohmdiss=etaMD

   elseif(resistivity_method==1) then ! analytical formula
      ! TODO
      etaohmdiss=etaMD

   else ! table
      ! convert to CGS
      rhoH=rhon*2.0d0*H2_fraction*scale_d/(mu_gas*mH) ! convert in H/cc
      rhoH = MAX(rhoH,rho_threshold)
      rhoH = MIN(rhoH,n_H_max)

      ! extrapolate from table[density,temperature,magnetic field]
      call interpolate_table(rhoH,temper,BBcell,sigO,sigH,sigP) 
      eta_ohm_chimie = (1d0 / sigP) * c_cgs * c_cgs / (4.0_dp*pi)

      ! Ad-hoc modification to ensure that the ohmic resistivity falls to zero when the density exceeds 1.0e15
      ! when alkali metals are ionized.
      eta_ohm_chimie = max(eta_ohm_chimie * (1.0d0-tanh(rhoH/1.0d15)), 1d-36)

      ! convert to code units
      etaohmdiss=etaohmdiss*scale_t/(scale_l)**2

      ! if the timestep was limited in courant fine, we need to adjust the resistivity to make things consistent.
      if(limit.and.nminitimestep) then
         if(dt.ne.0d0) then
            if(etaohmdiss.ne.0d0) then
               ! recalculate the ohmic timestep for the cell
               dtt=coefohm*dx*dx/etaohmdiss
            else
               dtt=1d39
            endif
            if (dtt.le.dt) then
               ! if it is smaller than the global timestep, we need to adjust the resistivity
               etaohmdiss=coefohm*dx*dx/dt
            endif
         endif
      endif
   endif

end function etaohmdiss
!###########################################################
!###########################################################
!###########################################################