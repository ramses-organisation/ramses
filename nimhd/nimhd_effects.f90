!  by Jacques Masson, Benoit Commercon and Neil Vaytet

!###########################################################
!###########################################################
!###########################################################
!###########################################################
! to remove! viscosity will be a separate update
! modif cmm
subroutine computevisco(q,ngrid,dx,dy,dz,dt,fvisco)

  USE amr_parameters
  use hydro_commons
  USE const
  IMPLICIT NONE

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q 
  INTEGER ::ngrid
  REAL(dp)::dx,dy,dz,dt
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:3) :: fvisco
  real(dp) :: muvisco

! declare local variables
  INTEGER :: i, j, k, l
  real(dp) :: rhox,rhoy,rhoz

 do k=min(1,ku1+1),ku2
     do j=min(1,ju1+1),ju2       
        do i=min(1,iu1+1),iu2
           do l=1,ngrid
              rhox=0.5d0*(q(l,i,j,k,1)+q(l,i-1,j,k,1))
              rhoy=0.5d0*(q(l,i,j,k,1)+q(l,i,j-1,k,1))
              rhoz=0.5d0*(q(l,i,j,k,1)+q(l,i,j,k-1,1))

! WARNING Flux F defined as dU/dt+dF/dx=0 
              fvisco(l,i,j,k,1,1)=-muvisco(rhox)*(q(l,i,j,k,2)-q(l,i-1,j,k,2))/dx
              fvisco(l,i,j,k,1,2)=-muvisco(rhoy)*(q(l,i,j,k,2)-q(l,i,j-1,k,2))/dy
              fvisco(l,i,j,k,1,3)=-muvisco(rhoz)*(q(l,i,j,k,2)-q(l,i,j,k-1,2))/dz
              fvisco(l,i,j,k,2,1)=-muvisco(rhox)*(q(l,i,j,k,3)-q(l,i-1,j,k,3))/dx
              fvisco(l,i,j,k,2,2)=-muvisco(rhoy)*(q(l,i,j,k,3)-q(l,i,j-1,k,3))/dy
              fvisco(l,i,j,k,2,3)=-muvisco(rhoz)*(q(l,i,j,k,3)-q(l,i,j,k-1,3))/dz
              fvisco(l,i,j,k,3,1)=-muvisco(rhox)*(q(l,i,j,k,4)-q(l,i-1,j,k,4))/dx
              fvisco(l,i,j,k,3,2)=-muvisco(rhoy)*(q(l,i,j,k,4)-q(l,i,j-1,k,4))/dy
              fvisco(l,i,j,k,3,3)=-muvisco(rhoz)*(q(l,i,j,k,4)-q(l,i,j,k-1,4))/dz

           end do
        end do
     end do
  end do

end subroutine computevisco
! fin modif cmm

!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine computejb(u,q,ngrid,dx,dy,dz,dt,bemfx,bemfy,bemfz,jemfx,jemfy,jemfz,bmagij,florentzx,florentzy,florentzz,fluxmd,fluxh,fluxad)
! comment Tine: Not used

  USE amr_parameters
  use hydro_commons
  USE const
  IMPLICIT NONE

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::u 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q 
  INTEGER ::ngrid
  REAL(dp)::dx,dy,dz,dt

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::bemfx,bemfy,bemfz
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::jemfx,jemfy,jemfz
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::florentzx,florentzy,florentzz
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::fluxmd,fluxh,fluxad
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:3)::bmagij

  ! declare local variables
  INTEGER ::i, j, k, l, m, n 

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::bmagijbis
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::flxmagxx,flxmagxy,flxmagxz,flxmagyx,flxmagyy,flxmagyz,flxmagzx,flxmagzy,flxmagzz
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:3)::jface,fluxbis,fluxter,fluxquat
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::bcenter
  real(dp)::v1x,v1y,v1z,v2x,v2y,v2z
  real(dp)::b12x,b12y,b12z,emag,bsquare
  real(dp)::computdivbisx,computdivbisy,computdivbisz
  real(dp)::computdxbis,computdybis,computdzbis
  real(dp)::crossprodx,crossprody,crossprodz

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

  if((nambipolar.eq.1).or.(nhall.eq.1)) then

     ! Fx,x
     do k=min(1,ku1+1),ku2
        do j=min(1,ju1+1),ju2       
           do i=min(1,iu1+1),iu2
              do l=1,ngrid
                 b12x=bmagijbis(l,i,j,k,1)
                 b12y=bmagijbis(l,i,j,k,2)
                 b12z=bmagijbis(l,i,j,k,3)
                 emag=0.5d0*(b12x*b12x+b12y*b12y+b12z*b12z)
                 flxmagxx(l,i,j,k,1)=b12x*b12x-emag
              end do
           end do
        end do
     end do
   
     ! Fx,y
     do k=min(1,ku1+1),ku2
        do j=ju1,max(1,ju2-1)
           do i=iu1,max(1,iu2-1)
              do l=1,ngrid
                 b12x=bmagij(l,i,j,k,1,3)
                 b12y=bmagij(l,i,j,k,2,3)
                 flxmagxx(l,i,j,k,2)=b12x*b12y
              end do
           end do
        end do
     end do
   
     ! Fx,z
     do k=ku1,max(1,ku2-1)
        do j=min(1,ju1+1),ju2
           do i=iu1,max(1,iu2-1)
              do l=1,ngrid
                 b12x=bmagij(l,i,j,k,1,2)
                 b12z=bmagij(l,i,j,k,3,2)
                 flxmagxx(l,i,j,k,3)=b12x*b12z
              end do
           end do
        end do
     end do
   
     ! Fy,x
     do k=min(1,ku1+1),ku2
        do j=min(1,ju1+1),ju2       
           do i=min(1,iu1+1),iu2
              do l=1,ngrid
                 b12x=bmagijbis(l,i,j,k,1)
                 b12y=bmagijbis(l,i,j,k,2)
                 flxmagyx(l,i,j,k,1)=b12y*b12x
              end do
           end do
        end do
     end do
   
     ! Fy,y
     do k=min(1,ku1+1),ku2
        do j=ju1,max(1,ju2-1)
           do i=iu1,max(1,iu2-1)
              do l=1,ngrid
                 b12x=bmagij(l,i,j,k,1,3)
                 b12y=bmagij(l,i,j,k,2,3)
                 b12z=bmagij(l,i,j,k,3,3)
                 emag=0.5d0*(b12x*b12x+b12y*b12y+b12z*b12z)
                 flxmagyx(l,i,j,k,2)=b12y*b12y-emag
              end do
           end do
        end do
     end do
   
     ! Fy,z
     do k=ku1,max(1,ku2-1)
        do j=min(1,ju1+1),ju2
           do i=iu1,max(1,iu2-1)
              do l=1,ngrid
                 b12y=bmagij(l,i,j,k,2,2)
                 b12z=bmagij(l,i,j,k,3,2)
                 flxmagyx(l,i,j,k,3)=b12y*b12z
              end do
           end do
        end do
     end do
   
     ! Fz,x
     do k=min(1,ku1+1),ku2
        do j=min(1,ju1+1),ju2       
           do i=min(1,iu1+1),iu2
              do l=1,ngrid
                 b12x=bmagijbis(l,i,j,k,1)
                 b12z=bmagijbis(l,i,j,k,3)
                 flxmagzx(l,i,j,k,1)=b12z*b12x
              end do
           end do
        end do
     end do
   
     ! Fz,y
     do k=min(1,ku1+1),ku2
        do j=ju1,max(1,ju2-1)
           do i=iu1,max(1,iu2-1)
              do l=1,ngrid
                 b12y=bmagij(l,i,j,k,2,3)
                 b12z=bmagij(l,i,j,k,3,3)               
                 flxmagzx(l,i,j,k,2)=b12z*b12y
              end do
           end do
        end do
     end do
   
     ! Fz,z
     do k=ku1,max(1,ku2-1)
        do j=min(1,ju1+1),ju2
           do i=iu1,max(1,iu2-1)
              do l=1,ngrid
                 b12x=bmagij(l,i,j,k,1,2) 
                 b12y=bmagij(l,i,j,k,2,2)
                 b12z=bmagij(l,i,j,k,3,2)
                 emag=0.5d0*(b12x*b12x+b12y*b12y+b12z*b12z)
                 flxmagzx(l,i,j,k,3)=b12z*b12z-emag
              end do
           end do
        end do
     end do
   
     !!!!!!!!!!!!!!!!!!!!!!!!!!!
   
     ! Fx,x
     do k=min(1,ku1+1),ku2
        do j= ju1,max(1,ju2-1)     
           do i=iu1,max(1,iu2-1)
              do l=1,ngrid
                 b12x=bmagij(l,i,j,k,1,3)
                 b12y=bmagij(l,i,j,k,2,3)
                 b12z=bmagij(l,i,j,k,3,3)
                 emag=0.5d0*(b12x*b12x+b12y*b12y+b12z*b12z)
                 flxmagxy(l,i,j,k,1)=b12x*b12x-emag
              end do
           end do
        end do
     end do
   
     ! Fx,y
     do k=min(1,ku1+1),ku2
        do j=min(1,ju1+1),ju2       
           do i=min(1,iu1+1),iu2
              do l=1,ngrid
                 b12x=bmagijbis(l,i,j,k,1)
                 b12y=bmagijbis(l,i,j,k,2)
                 flxmagxy(l,i,j,k,2)=b12x*b12y
              end do
           end do
        end do
     end do
   
     ! Fx,z
     do k=ku1,max(1,ku2-1)
       do j=ju1,max(1,ju2-1)
          do i=min(1,iu1+1),iu2
              do l=1,ngrid
                 b12x=bmagij(l,i,j,k,1,1)
                 b12z=bmagij(l,i,j,k,3,1)
                 flxmagxy(l,i,j,k,3)=b12x*b12z
              end do
           end do
        end do
     end do
   
     ! Fy,x
     do k=min(1,ku1+1),ku2
        do j= ju1,max(1,ju2-1)     
           do i=iu1,max(1,iu2-1)
              do l=1,ngrid
                 b12x=bmagij(l,i,j,k,1,3)
                 b12y=bmagij(l,i,j,k,2,3)
                 flxmagyy(l,i,j,k,1)=b12y*b12x
              end do
           end do
        end do
     end do
   
     ! Fy,y
     do k=min(1,ku1+1),ku2
        do j=min(1,ju1+1),ju2       
           do i=min(1,iu1+1),iu2
              do l=1,ngrid
                 b12x=bmagijbis(l,i,j,k,1)
                 b12y=bmagijbis(l,i,j,k,2)
                 b12z=bmagijbis(l,i,j,k,3)
                 emag=0.5d0*(b12x*b12x+b12y*b12y+b12z*b12z)
                 flxmagyy(l,i,j,k,2)=b12y*b12y-emag
              end do
           end do
        end do
     end do
   
     ! Fy,z
     do k=ku1,max(1,ku2-1)
       do j=ju1,max(1,ju2-1)
          do i=min(1,iu1+1),iu2
              do l=1,ngrid
                 b12y=bmagij(l,i,j,k,2,1)
                 b12z=bmagij(l,i,j,k,3,1)
                 flxmagyy(l,i,j,k,3)=b12y*b12z
              end do
           end do
        end do
     end do
   
     ! Fz,x
     do k=min(1,ku1+1),ku2
        do j= ju1,max(1,ju2-1)     
           do i=iu1,max(1,iu2-1)
              do l=1,ngrid
                 b12x=bmagij(l,i,j,k,1,3)
                 b12z=bmagij(l,i,j,k,3,3)
                 flxmagzy(l,i,j,k,1)=b12z*b12x
              end do
           end do
        end do
     end do
   
     ! Fz,y
     do k=min(1,ku1+1),ku2
        do j=min(1,ju1+1),ju2       
           do i=min(1,iu1+1),iu2
              do l=1,ngrid
                 b12y=bmagijbis(l,i,j,k,2)
                 b12z=bmagijbis(l,i,j,k,3)
                 flxmagzy(l,i,j,k,2)=b12z*b12y
              end do
           end do
        end do
     end do
   
     ! Fz,z
     do k=ku1,max(1,ku2-1)
       do j=ju1,max(1,ju2-1)
          do i=min(1,iu1+1),iu2
              do l=1,ngrid
                 b12x=bmagij(l,i,j,k,1,1)
                 b12y=bmagij(l,i,j,k,2,1)
                 b12z=bmagij(l,i,j,k,3,1)
                 emag=0.5d0*(b12x*b12x+b12y*b12y+b12z*b12z)
                 flxmagzy(l,i,j,k,3)=b12z*b12z-emag
              end do
           end do
        end do
     end do
   
     ! Fx,x
     do k=ku1,max(1,ku2-1)
        do j=min(1,ju1+1),ju2     
           do i=iu1,max(1,iu2-1)
              do l=1,ngrid
                 b12x=bmagij(l,i,j,k,1,2)
                 b12y=bmagij(l,i,j,k,2,2)
                 b12z=bmagij(l,i,j,k,3,2)
                 emag=0.5d0*(b12x*b12x+b12y*b12y+b12z*b12z)
                 flxmagxz(l,i,j,k,1)=b12x*b12x-emag
              end do
           end do
        end do
     end do
   
     ! Fx,y
     do k=ku1,max(1,ku2-1)
        do j=ju1,max(1,ju2-1)      
           do i=min(1,iu1+1),iu2 
              do l=1,ngrid
                 b12x=bmagij(l,i,j,k,1,1)
                 b12y=bmagij(l,i,j,k,2,1)
                 flxmagxz(l,i,j,k,2)=b12x*b12y
              end do
           end do
        end do
     end do
   
     ! Fx,z
     do k=min(1,ku1+1),ku2
        do j=min(1,ju1+1),ju2       
           do i=min(1,iu1+1),iu2
              do l=1,ngrid
                 b12x=bmagijbis(l,i,j,k,1)
                 b12z=bmagijbis(l,i,j,k,3)
                 flxmagxz(l,i,j,k,3)=b12x*b12z
              end do
           end do
        end do
     end do
   
     ! Fy,x
     do k=ku1,max(1,ku2-1)
        do j=min(1,ju1+1),ju2     
           do i=iu1,max(1,iu2-1)
              do l=1,ngrid
                 b12x=bmagij(l,i,j,k,1,2)
                 b12y=bmagij(l,i,j,k,2,2)
                 flxmagyz(l,i,j,k,1)=b12y*b12x
              end do
           end do
        end do
     end do
   
     ! Fy,y
     do k=ku1,max(1,ku2-1)
        do j=ju1,max(1,ju2-1)      
           do i=min(1,iu1+1),iu2 
              do l=1,ngrid
                 b12x=bmagij(l,i,j,k,1,1)
                 b12y=bmagij(l,i,j,k,2,1)
                 b12z=bmagij(l,i,j,k,3,1)
                 emag=0.5d0*(b12x*b12x+b12y*b12y+b12z*b12z)
                 flxmagyz(l,i,j,k,2)=b12y*b12y-emag
              end do
           end do
        end do
     end do
   
     ! Fy,z
     do k=min(1,ku1+1),ku2
        do j=min(1,ju1+1),ju2       
           do i=min(1,iu1+1),iu2
              do l=1,ngrid
                 b12y=bmagijbis(l,i,j,k,2)
                 b12z=bmagijbis(l,i,j,k,3)
                 flxmagyz(l,i,j,k,3)=b12y*b12z
              end do
           end do
        end do
     end do
   
     ! Fz,x
     do k=ku1,max(1,ku2-1)
        do j=min(1,ju1+1),ju2     
           do i=iu1,max(1,iu2-1)
              do l=1,ngrid
                 b12x=bmagij(l,i,j,k,1,2)
                 b12z=bmagij(l,i,j,k,3,2)
                 flxmagzz(l,i,j,k,1)=b12z*b12x
              end do
           end do
        end do
     end do
   
     ! Fz,y
     do k=ku1,max(1,ku2-1)
        do j=ju1,max(1,ju2-1)      
           do i=min(1,iu1+1),iu2 
              do l=1,ngrid
                 b12y=bmagij(l,i,j,k,2,1)
                 b12z=bmagij(l,i,j,k,3,1)
                 flxmagzz(l,i,j,k,2)=b12z*b12y
              end do
           end do
        end do
     end do
   
     ! Fz,z
     do k=min(1,ku1+1),ku2
        do j=min(1,ju1+1),ju2       
           do i=min(1,iu1+1),iu2
              do l=1,ngrid
                 b12x=bmagijbis(l,i,j,k,1)
                 b12y=bmagijbis(l,i,j,k,2)
                 b12z=bmagijbis(l,i,j,k,3)
                 emag=0.5d0*(b12x*b12x+b12y*b12y+b12z*b12z)
                 flxmagzz(l,i,j,k,3)=b12z*b12z-emag
              end do
           end do
        end do
     end do
   
     !!!!!!!!!!!!!!!!!!!!!!!
     do k=min(1,ku1+1),max(1,ku2-1)
        do j=min(1,ju1+1),max(1,ju2-1)
           do i=min(1,iu1+1),max(1,iu2-1)
              do l = 1, ngrid
                 florentzx(l,i,j,k,1)=computdivbisx(flxmagxx,l,i,j,k,dx,dy,dz)
                 florentzx(l,i,j,k,2)=computdivbisx(flxmagyx,l,i,j,k,dx,dy,dz)
                 florentzx(l,i,j,k,3)=computdivbisx(flxmagzx,l,i,j,k,dx,dy,dz)
   
                 florentzy(l,i,j,k,1)=computdivbisy(flxmagxy,l,i,j,k,dx,dy,dz)
                 florentzy(l,i,j,k,2)=computdivbisy(flxmagyy,l,i,j,k,dx,dy,dz)
                 florentzy(l,i,j,k,3)=computdivbisy(flxmagzy,l,i,j,k,dx,dy,dz)
   
                 florentzz(l,i,j,k,1)=computdivbisz(flxmagxz,l,i,j,k,dx,dy,dz)
                 florentzz(l,i,j,k,2)=computdivbisz(flxmagyz,l,i,j,k,dx,dy,dz)
                 florentzz(l,i,j,k,3)=computdivbisz(flxmagzz,l,i,j,k,dx,dy,dz)
              end do
           end do
        end do
     end do
   
   ! end if((nambipolar.eq.1).or.(nhall.eq.1)) then
   endif

  ! computation of current on faces

  if((nambipolar.eq.1).or.(nhall.eq.1).or.(nmagdiffu.eq.1).or.(nmagdiffu2.eq.1)) then

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
   ! end if((nambipolar.eq.1).or.(nhall.eq.1).or.(nmagdiffu.eq.1)) then
  endif

  if((nambipolar.eq.1).or.(nhall.eq.1)) then
     do k=min(1,ku1+1),max(1,ku2-1)
        do j=min(1,ju1+1),max(1,ju2-1)
           do i=min(1,iu1+1),max(1,iu2-1)
              do l = 1, ngrid
                 call crossprodbis(fluxbis,bmagij,fluxter,l,i,j,k)
                 fluxh(l,i,j,k,1)=fluxter(l,i,j,k,1,1)
                 fluxh(l,i,j,k,2)=fluxter(l,i,j,k,2,2)
                 fluxh(l,i,j,k,3)=fluxter(l,i,j,k,3,3)
              end do
           end do
        end do
     end do
  endif

  if(nambipolar.eq.1)then
     do k=min(1,ku1+1),max(1,ku2-1)
        do j=min(1,ju1+1),max(1,ju2-1)
           do i=min(1,iu1+1),max(1,iu2-1)
              do l = 1, ngrid
                 call crossprodbis(fluxter,bmagij,fluxquat,l,i,j,k)
                 fluxad(l,i,j,k,1)=fluxquat(l,i,j,k,1,1)
                 fluxad(l,i,j,k,2)=fluxquat(l,i,j,k,2,2)
                 fluxad(l,i,j,k,3)=fluxquat(l,i,j,k,3,3)
              end do
           end do
        end do
     end do
  endif

end subroutine computejb
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine computejb2(u,q,ngrid,dx,dy,dz,dt,bemfx,bemfy,bemfz,jemfx,jemfy,jemfz,bmagij,florentzx,florentzy,florentzz,fluxmd,fluxh,fluxad,jcell)

  USE amr_parameters
  use hydro_commons
  USE const
  IMPLICIT NONE

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::u 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q 

  INTEGER ::ngrid
  REAL(dp)::dx,dy,dz,dt

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::bemfx,bemfy,bemfz
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::jemfx,jemfy,jemfz
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::florentzx,florentzy,florentzz
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::fluxmd,fluxh,fluxad
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:3)::bmagij
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::jcell
 
  ! declare local variables
  INTEGER ::i, j, k, l, m, n 

  real(dp)::computdx,computdy,computdz
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
!
! EMF x
!
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
!
! EMF y
!
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
!
! EMF z
!
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
!
! bmagijbis(l,i,j,k,n) is the value of the magnetic field component
! Bn at i-1/2,j-1/2,k-1/2
!
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
!
! computation of the component of j where EMFs are located
! jemfx(l,i,j,k,n) is the component Jn at i,j-1/2,k-1/2
!
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
!
! computation of the component of j at center of cell
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do k=min(1,ku1+1),max(1,ku2-1)
     do j=min(1,ju1+1),max(1,ju2-1)
        do i=min(1,iu1+1),max(1,iu2-1)
           do l=1,ngrid
              jcell(l,i,j,k,1)=computdy(bmagij,nzz,nyy,l,i,j,k,dy)-computdz(bmagij,nyy,nzz,l,i,j,k,dy)
              jcell(l,i,j,k,2)=computdz(bmagij,nxx,nzz,l,i,j,k,dy)-computdx(bmagij,nzz,nxx,l,i,j,k,dy)
              jcell(l,i,j,k,3)=computdx(bmagij,nyy,nxx,l,i,j,k,dy)-computdy(bmagij,nxx,nyy,l,i,j,k,dy)
           end do
        end do
     end do
  end do

  if((nambipolar.eq.1).or.(nhall.eq.1).or.(nambipolar2.eq.1)) then
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

  if((nambipolar.eq.1).or.(nhall.eq.1).or.(nmagdiffu.eq.1).or.(nambipolar2.eq.1).or.(nmagdiffu2.eq.1)) then

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

  endif

  if((nambipolar.eq.1).or.(nhall.eq.1).or.(nambipolar2==1)) then

     do k=min(1,ku1+1),max(1,ku2-1)
        do j=min(1,ju1+1),max(1,ju2-1)
           do i=min(1,iu1+1),max(1,iu2-1)
              do l = 1, ngrid
                 call crossprodbis(fluxbis,bmagij,fluxter,l,i,j,k)
                 fluxh(l,i,j,k,1)=fluxter(l,i,j,k,1,1)
                 fluxh(l,i,j,k,2)=fluxter(l,i,j,k,2,2)
                 fluxh(l,i,j,k,3)=fluxter(l,i,j,k,3,3)
              end do
           end do
        end do
     end do

  endif

  if((nambipolar.eq.1).or.(nambipolar2==1))then

     do k=min(1,ku1+1),max(1,ku2-1)
        do j=min(1,ju1+1),max(1,ju2-1)
           do i=min(1,iu1+1),max(1,iu2-1)
              do l = 1, ngrid
                 call crossprodbis(fluxter,bmagij,fluxquat,l,i,j,k)
                 fluxad(l,i,j,k,1)=fluxquat(l,i,j,k,1,1)
                 fluxad(l,i,j,k,2)=fluxquat(l,i,j,k,2,2)
                 fluxad(l,i,j,k,3)=fluxquat(l,i,j,k,3,3)
              end do
           end do
        end do
     end do

  endif

  do k=min(1,ku1+1),max(1,ku2-1)
     do j=min(1,ju1+1),max(1,ju2-1)
        do i=min(1,iu1+1),max(1,iu2-1)   
           do l = 1, ngrid
              call crossprodbis(fluxter,bmagij,fluxquat,l,i,j,k)
              fluxad(l,i,j,k,1)=fluxquat(l,i,j,k,1,1)
              fluxad(l,i,j,k,2)=fluxquat(l,i,j,k,2,2)
              fluxad(l,i,j,k,3)=fluxquat(l,i,j,k,3,3)
           end do
        end do
     end do
  end do

end subroutine computejb2
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine computdifmag(u,q,ngrid,dx,dy,dz,dt,bemfx,bemfy,bemfz,jemfx,jemfy,jemfz,bmagij,fluxmd,emfohmdiss,fluxohm,jcentersquare)

  USE amr_parameters
  use hydro_commons
  USE const
  use constants
  IMPLICIT NONE

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::u 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q 

  INTEGER ::ngrid
  REAL(dp)::dx,dy,dz,dt

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::bemfx,bemfy,bemfz
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::jemfx,jemfy,jemfz
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:3)::bmagij
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::fluxmd

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3):: emfohmdiss,fluxohm 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::jcentersquare

  ! declare local variables
  INTEGER ::i, j, k, l, m, n,h

! WARNING following quantities defined with three components even
! if ndim<3 !
real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::jcenter
real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:3)::jface
real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::jemf


real(dp)::computdx,computdy,computdz,tcell,rhocell,bcell,ionisrate
real(dp)::crossprodx,crossprody,crossprodz, Cv

real(dp)::etaohmdiss,etaod2,tcellx,tcelly,tcellz,bsquarex,bsquarey,bsquarez,etaohmdissbricolo,dtlim
real(dp)::pressurex,pressurey,pressurez,rhox,rhoy,rhoz,epsx,epsy,epsz
real(dp)::etaod2x,etaod2y,etaod2z,rhof,pf,bsqf,epsf,tcellf,barotrop1D

integer , dimension(1:3) :: index_i,index_j,index_k

   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

index_i = (/1,0,0/)
index_j = (/0,1,0/)
index_k = (/0,0,1/)

dtlim = dt !neil

do k=min(1,ku1+1),max(1,ku2-1)
     do j=min(1,ju1+1),max(1,ju2-1)
        do i=min(1,iu1+1),max(1,iu2-1)
           
           do l=1,ngrid

              jemf(l,i,j,k,1)=jemfx(l,i,j,k,1)
              jemf(l,i,j,k,2)=jemfy(l,i,j,k,2)
              jemf(l,i,j,k,3)=jemfz(l,i,j,k,3)


              rhox=0.25d0*(u(l,i,j,k,   1)+u(l,i  ,j-1,k,   1)+u(l,i,j  ,k-1,   1)+u(l,i  ,j-1,k-1,   1))
              rhoy=0.25d0*(u(l,i,j,k,   1)+u(l,i-1,j  ,k,   1)+u(l,i,j  ,k-1,   1)+u(l,i-1,j  ,k-1,   1))
              rhoz=0.25d0*(u(l,i,j,k,   1)+u(l,i-1,j  ,k,   1)+u(l,i,j-1,k  ,   1)+u(l,i-1,j-1,k  ,   1))
              epsx=0.25d0*(u(l,i,j,k,nvar)+u(l,i  ,j-1,k,nvar)+u(l,i,j  ,k-1,nvar)+u(l,i  ,j-1,k-1,nvar))
              epsy=0.25d0*(u(l,i,j,k,nvar)+u(l,i-1,j  ,k,nvar)+u(l,i,j  ,k-1,nvar)+u(l,i-1,j  ,k-1,nvar))
              epsz=0.25d0*(u(l,i,j,k,nvar)+u(l,i-1,j  ,k,nvar)+u(l,i,j-1,k  ,nvar)+u(l,i-1,j-1,k  ,nvar))
               if(nmagdiffu .eq.1)then
                 bsquarex=bemfx(l,i,j,k,1)**2+bemfx(l,i,j,k,2)**2+bemfx(l,i,j,k,3)**2
                 bsquarey=bemfy(l,i,j,k,1)**2+bemfy(l,i,j,k,2)**2+bemfy(l,i,j,k,3)**2
                 bsquarez=bemfz(l,i,j,k,1)**2+bemfz(l,i,j,k,2)**2+bemfz(l,i,j,k,3)**2
               else if(nmagdiffu2 .eq.1)then
                  bsquarex=u(l,i,j,k,   2)
                  bsquarey=u(l,i,j,k,   2)
                  bsquarez=u(l,i,j,k,   2)
                  epsx=u(l,i,j,k,nvar)
                  epsy=u(l,i,j,k,nvar)
                  epsz=u(l,i,j,k,nvar)
                  rhox=u(l,i,j,k,1)
                  rhoy=u(l,i,j,k,1)
                  rhoz=u(l,i,j,k,1)
                  if(epsx .ne.u(l,i,j,k,3))then
                     ! Attention, on est sur les boundary du domaine, divu et enew ne sont pas connus....
                     bsquarex=bemfx(l,i,j,k,1)**2+bemfx(l,i,j,k,2)**2+bemfx(l,i,j,k,3)**2
                     bsquarey=bemfy(l,i,j,k,1)**2+bemfy(l,i,j,k,2)**2+bemfy(l,i,j,k,3)**2
                     bsquarez=bemfz(l,i,j,k,1)**2+bemfz(l,i,j,k,2)**2+bemfz(l,i,j,k,3)**2
                  end if
               end if

               if(ntestDADM.eq.1)then
                  tcellx=1.0d0
                  tcelly=1.0d0
                  tcellz=1.0d0
               else
!                  print*,'x',rhox,epsx,u(l,i,j,k,2),bemfx(l,i,j,k,1)**2+bemfx(l,i,j,k,2)**2+bemfx(l,i,j,k,3)**2
!                  if(epsy* scale_d*scale_v**2  .lt. 1d-16) print*,'y',rhoy,epsy
!                  if(epsz* scale_d*scale_v**2  .lt. 1d-16) print*,'z',rhoz,epsz
!                  print*,rhox,epsx,rhoy,epsy,rhoz,epsz
                  call temperature_eos(rhox,epsx,tcellx)
                  call temperature_eos(rhoy,epsy,tcelly)
                  call temperature_eos(rhoz,epsz,tcellz)
!!$                  tcelly=10.
!!$                  tcellx=10.
!!$                  tcellz=10.
               endif
              ionisrate=default_ionisrate 
!               etaod2x=etaohmdiss(rhox,bsquarex,tcellx)
!               etaod2y=etaohmdiss(rhoy,bsquarey,tcelly)
!               etaod2z=etaohmdiss(rhoz,bsquarez,tcellz)
              etaod2x=etaohmdissbricolo(rhox,bsquarex,tcellx,dtlim,dx,ionisrate)
              etaod2y=etaohmdissbricolo(rhoy,bsquarey,tcelly,dtlim,dx,ionisrate)
              etaod2z=etaohmdissbricolo(rhoz,bsquarez,tcellz,dtlim,dx,ionisrate)
              
! WARNING dB/dt=-curl(eta*J)
              emfohmdiss(l,i,j,k,nxx)=-etaod2x*jemf(l,i,j,k,1)
              emfohmdiss(l,i,j,k,nyy)=-etaod2y*jemf(l,i,j,k,2)
              emfohmdiss(l,i,j,k,nzz)=-etaod2z*jemf(l,i,j,k,3)

! !!!!!!!!!!!!!!!!!!!!!!!
! !
! ! compute j at center of cells
! !
! ! mandatory for non isotherm case
! 
! ! bmagij is the value of the magnetic field Bi where Bj 
! ! is naturally defined; Ex bmagij(l,i,j,k,1,2) is Bx at i,j-1/2,k
! ! and we can write it Bx,y

              jcenter(l,i,j,k,1)=computdy(bmagij,nzz,nyy,l,i,j,k,dy)-computdz(bmagij,nyy,nzz,l,i,j,k,dy)
              jcenter(l,i,j,k,2)=computdz(bmagij,nxx,nzz,l,i,j,k,dy)-computdx(bmagij,nzz,nxx,l,i,j,k,dy)
              jcenter(l,i,j,k,3)=computdx(bmagij,nyy,nxx,l,i,j,k,dy)-computdy(bmagij,nxx,nyy,l,i,j,k,dy)

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !
! ! computation of j at faces
! ! jface is the value of the current where Bj 
! ! is naturally defined; Ex jface(l,i,j,k,1,2) is Jx at i,j-1/2,k
! ! and we can write it Jx,y
              if(nmagdiffu2 .eq.0)then
                 do h = 1,3
                    
                    rhof=0.5d0*(u(l,i,j,k,   1)+u(l,i-index_i(h),j-index_j(h),k-index_k(h),   1))
!                 epsf=u(l,i,j,k,3)
                    epsf=0.5d0*(u(l,i,j,k,nvar)+u(l,i-index_i(h),j-index_j(h),k-index_k(h),nvar))
                    bsqf=bmagij(l,i,j,k,1,h)**2+bmagij(l,i,j,k,2,h)**2+bmagij(l,i,j,k,3,h)**2
                    
                 ! Compute gas temperature in cgs

                    if(barotrop)then
                       tcellf=barotrop1D(rhof*scale_d)
                    elseif(ntestDADM.eq.1)then
                       tcellf=1.0d0
                    else 
                       call temperature_eos(rhof,epsf,tcellf)
                    endif
                    
                    etaod2=etaohmdiss(rhof,bsqf,tcellf,ionisrate)
                    fluxohm(l,i,j,k,h)=etaod2*fluxmd(l,i,j,k,h)
                    
                    !               rhof=0.5d0*(u(l,i,j,k,1)+u(l,i-1,j,k,1))
                    !               pf=0.5d0*(q(l,i,j,k,5)+q(l,i-1,j,k,5))
!               etaod2=etaohmdiss(rhof,pf)
!               fluxohm(l,i,j,k,1)=etaod2*fluxmd(l,i,j,k,1)
                    
                 enddo
              end if

!            end do
!         end do
!      end do
!   end do
! 
! 
!   do k=min(1,ku1+1),max(1,ku2-1) 
! !     do j=min(1,ju1+1),ju2 
!      do j=min(1,ju1+1),max(1,ju2-1)
!         do i=min(1,iu1+1),max(1,iu2-1) 
!            
!            do l=1,ngrid
! 
!               rhof=0.5d0*(u(l,i,j,k,1)+u(l,i,j-1,k,1))
!               pf=0.5d0*(q(l,i,j,k,5)+q(l,i,j-1,k,5))
!               etaod2=etaohmdiss(rhof,pf)
!               fluxohm(l,i,j,k,2)=etaod2*fluxmd(l,i,j,k,2)
! 
!            end do
!         end do
!      end do
!   end do
! 
! 
!  do k=min(1,ku1+1),max(1,ku2-1) 
! !     do j=min(1,ju1+1),ju2  
!     do j=min(1,ju1+1),max(1,ju2-1)
!         do i=min(1,iu1+1),max(1,iu2-1) 
!            
!            do l=1,ngrid
! 
!               rhof=0.5d0*(u(l,i,j,k,1)+u(l,i,j,k-1,1))
!               pf=0.5d0*(q(l,i,j,k,5)+q(l,i,j,k-1,5))
!               etaod2=etaohmdiss(rhof,pf)
!               fluxohm(l,i,j,k,3)=etaod2*fluxmd(l,i,j,k,3)
! 
!            end do
!         end do
!      end do
!   end do


! compute contribution to energy flux +eta*I*B


!            do l=1,ngrid
              if(nmagdiffu2 .eq. 1)then
                 jcentersquare(l,i,j,k)=jcenter(l,i,j,k,1)*jcenter(l,i,j,k,1)+jcenter(l,i,j,k,2)*jcenter(l,i,j,k,2)+jcenter(l,i,j,k,3)*jcenter(l,i,j,k,3)
                 
                 rhocell = u(l,i,j,k,1)
                 bcell   = u(l,i,j,k,2)
                 if(u(l,i,j,k,nvar) .ne.u(l,i,j,k,3))then
                    ! Attention, on est sur les boundary du domaine, divu et enew ne sont pas connus....
                    bcell=(0.5*(u(l,i,j,k,6)+u(l,i,j,k,nvar+1)))**2 + (0.5*(u(l,i,j,k,7)+u(l,i,j,k,nvar+2)))**2 +(0.5*(u(l,i,j,k,8)+u(l,i,j,k,nvar+3)))**2
                 end if

                 if(ntestDADM.eq.1)then
                    tcell=1.0d0
                 else 
                    call temperature_eos(rhocell,u(l,i,j,k,nvar),tcell)
!                    if(nmagdiffu2.eq.1)call temperature_eos(rhocell,u(l,i,j,k,3),tcell)
                    end if
                    
                    jcentersquare(l,i,j,k) = jcentersquare(l,i,j,k)*etaohmdiss(rhocell,bcell,tcell,ionisrate)*dt
                    
                 end if
              end do
        end do
     end do
  end do
  
end subroutine computdifmag
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE  computambip(u,q,ngrid,dx,dy,dz,dt,bemfx,bemfy,bemfz,florentzx,florentzy,florentzz,fluxad,bmagij,emfambdiff,fluxambdiff,jxbsquare)

  use amr_commons
  USE amr_parameters
  use hydro_commons
  USE const
  IMPLICIT NONE
  
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::u 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q 
  
  INTEGER ::ngrid
  REAL(dp)::dx,dy,dz,dt
  
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::bemfx,bemfy,bemfz
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::rhocellmin,bsquaremax
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::florentzx,florentzy,florentzz
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::fluxad
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:3)::bmagij
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::emfambdiff
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::fluxambdiff

! declare local variables
  INTEGER ::i, j, k, l, m, n, ntest,ic,ivar
  real(dp)::computdx,computdy,computdz

  real(dp)::v1x,v1y,v1z,v2x,v2y,v2z
  real(dp)::rhofx,rhofy,rhofz
  real(dp)::bsquarex,bsquarey,bsquarez,bsquare
  real(dp)::bsquarexx,bsquareyy,bsquarezz
  real(dp)::betaad2,betaadbricolo,betaad
  real(dp)::rhox,rhoy,rhoz,rhocell,bcell,bcellold,tcell,ionisrate
  real(dp)::dtlim,Cv,eps
  real(dp)::crossprodx,crossprody,crossprodz

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::florentz

  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::jxbsquare
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::jcenter
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::jxb

   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

!modif pour voir les lieux du seuil
!  real(dp),dimension(1:3)             :: skip_loc
!  real(dp),dimension(1:twotondim,1:3) :: xc
!  integer                             :: nx_loc
!  real(dp)                            :: scale


!nx_loc = (icoarse_max -icoarse_min+1)
!scale = dble(nx_loc)/boxlen
!print*, dx, 0.5d0**11, 0.5d0**7/scale
!
!do  ind=1,twotondim
!    iz=(ind-1)/4
!    iy=



! do NOT change value below Variation of betaad
! to avoid too small time step allowed
  ntest=0

  dtlim=dt!*coefalfven
!dt est deja dtnew, qui a été choisi comme le dt normal (avec la condition de courant) ou le dt normal seuillé si le dtAD est trop faible(bricolo)

jxb=0.0d0

jxbsquare=0.0d0
jcenter=0.0d0

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
              
              rhocellmin(l,i,j,k)=min(rhox,rhoy,rhoz,rhofx,rhofy,rhofz)

              rhocell = u(l,i,j,k,1)

              ! Compute gas temperature in cgs
              if(ntestDADM.eq.1) then
                 tcell=1.0d0
              else
                 call temperature_eos(u(l,i,j,k,1),u(l,i,j,k,nvar),tcell)
                 if(nambipolar2.eq.1)call temperature_eos(u(l,i,j,k,1),u(l,i,j,k,3),tcell)
              end if
             

              bsquarex=bemfx(l,i,j,k,1)**2+bemfx(l,i,j,k,2)**2+bemfx(l,i,j,k,3)**2
              bsquarey=bemfy(l,i,j,k,1)**2+bemfy(l,i,j,k,2)**2+bemfy(l,i,j,k,3)**2
              bsquarez=bemfz(l,i,j,k,1)**2+bemfz(l,i,j,k,2)**2+bemfz(l,i,j,k,3)**2


              bsquarexx=bmagij(l,i,j,k,1,1)**2+bmagij(l,i,j,k,2,1)**2+bmagij(l,i,j,k,3,1)**2
              bsquareyy=bmagij(l,i,j,k,1,2)**2+bmagij(l,i,j,k,2,2)**2+bmagij(l,i,j,k,3,2)**2
              bsquarezz=bmagij(l,i,j,k,1,3)**2+bmagij(l,i,j,k,2,3)**2+bmagij(l,i,j,k,3,3)**2

              bsquaremax(l,i,j,k)=max(bsquarex,bsquarey,bsquarez,bsquarexx,bsquareyy,bsquarezz)
                 
! EMF x
  
              v1x=florentzx(l,i,j,k,1)
              v1y=florentzx(l,i,j,k,2)
              v1z=florentzx(l,i,j,k,3)
              v2x=bemfx(l,i,j,k,1)
              v2y=bemfx(l,i,j,k,2)
              v2z=bemfx(l,i,j,k,3)
             
              emfambdiff(l,i,j,k,1)=crossprodx(v1x,v1y,v1z,v2x,v2y,v2z)
              rhox=0.25d0*(u(l,i,j,k,1)+u(l,i,j-1,k,1)+u(l,i,j,k-1,1)+u(l,i,j-1,k-1,1))
              bcell = bsquaremax(l,i,j,k)
              bcellold=bcell
              if(nambipolar2.eq.1)then
!                 bcell=v2x*v2x+v2y*v2y+v2z*v2z
                 bcellold=u(l,i,j,k,2)
!!$                 rhox=rhocell
!!$                 rhoy=rhocell
!!$                 rhoz=rhocell
              end if

              if(nambipolar2 .eq. 0)rhocell = rhocellmin(l,i,j,k)
! alfven time alone maybe not correct
!             betaad2=betaadbricolo(rhox,dtlim,bsquare,dx,ntest)
! comparison with hydro+idealMHD
!!$              rhox= u(l,i,j,k,3)
!!$              rhoy= u(l,i,j,k,3)
!!$              rhoz= u(l,i,j,k,3)
              ionisrate=default_ionisrate
              betaad2=betaadbricolo(rhocell,rhox,dtlim,bcell,bcellold,dx,ntest,tcell,ionisrate)
!              betaad2=betaadbricolo(rhocell,rhox,dtlim,bcellold,bcellold,dx,ntest,tcell)

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
              bcell = bsquaremax(l,i,j,k)
              bcellold=bcell
              if(nambipolar2.eq.1)then
!                 bcell=v2x*v2x+v2y*v2y+v2z*v2z
                 bcellold=u(l,i,j,k,2)
              end if

              if(nambipolar2 .eq. 0)rhocell = rhocellmin(l,i,j,k)
! alfven time alone maybe not correct
!             betaad2=betaadbricolo(rhoy,dtlim,bsquare,dx,ntest)
! comparison with hydro+idealMHD 

             betaad2=betaadbricolo(rhocell,rhoy,dtlim,bcell,bcellold,dx,ntest,tcell,ionisrate)
!             betaad2=betaadbricolo(rhocell,rhoy,dtlim,bcellold,bcellold,dx,ntest,tcell)

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
              bcell = bsquaremax(l,i,j,k)
              bcellold=bcell
              if(nambipolar2.eq.1)then
 !                bcell=v2x*v2x+v2y*v2y+v2z*v2z
                 bcellold=u(l,i,j,k,2)
              end if
              if(nambipolar2 .eq. 0) rhocell = rhocellmin(l,i,j,k)
             
             betaad2=betaadbricolo(rhocell,rhoz,dtlim,bcell,bcellold,dx,ntest,tcell,ionisrate)
!             betaad2=betaadbricolo(rhocell,rhoz,dtlim,bcellold,bcellold,dx,ntest,tcell)

             emfambdiff(l,i,j,k,3)=emfambdiff(l,i,j,k,3)*betaad2

! energy flux on faces

              v2x=bmagij(l,i,j,k,1,1)
              v2y=bmagij(l,i,j,k,2,1)
              v2z=bmagij(l,i,j,k,3,1)

             bcell = bsquaremax(l,i,j,k)
             if(nambipolar2.eq.1)then
                 bcell=v2x*v2x+v2y*v2y+v2z*v2z
              end if


              rhocell = rhocellmin(l,i,j,k)
              rhofx=0.5d0*(u(l,i,j,k,1)+u(l,i-1,j,k,1))
              betaad2=betaadbricolo(rhocell,rhofx,dtlim,bcell,bcell,dx,ntest,tcell,ionisrate)
              fluxambdiff(l,i,j,k,1)=-betaad2*fluxad(l,i,j,k,1)

              v2x=bmagij(l,i,j,k,1,2)
              v2y=bmagij(l,i,j,k,2,2)
              v2z=bmagij(l,i,j,k,3,2)

              rhofy=0.5d0*(u(l,i,j,k,1)+u(l,i,j-1,k,1))
              bcell = bsquaremax(l,i,j,k)
              if(nambipolar2.eq.1)then
!                 bcell=0.5d0*(u(l,i,j,k,2)+u(l,i,j-1,k,2))
                 bcell=v2x*v2x+v2y*v2y+v2z*v2z
              end if

              betaad2=betaadbricolo(rhocell,rhofy,dtlim,bcell,bcell,dx,ntest,tcell,ionisrate)
              fluxambdiff(l,i,j,k,2)=-betaad2*fluxad(l,i,j,k,2)

              v2x=bmagij(l,i,j,k,1,3)
              v2y=bmagij(l,i,j,k,2,3)
              v2z=bmagij(l,i,j,k,3,3)
!              bsquare=v2x*v2x+v2y*v2y+v2z*v2z
              rhofz=0.5d0*(u(l,i,j,k,1)+u(l,i,j,k-1,1))
              bcell = bsquaremax(l,i,j,k)
              if(nambipolar2.eq.1)then
!                 bcell=0.5d0*(u(l,i,j,k,2)+u(l,i,j,k-1,2))
                 bcell=v2x*v2x+v2y*v2y+v2z*v2z
              end if

              betaad2=betaadbricolo(rhocell,rhofz,dtlim,bcell,bcell,dx,ntest,tcell,ionisrate)
              fluxambdiff(l,i,j,k,3)=-betaad2*fluxad(l,i,j,k,3)

              v2x=u(l,i,j,k,6)
              v2y=u(l,i,j,k,7)
              v2z=u(l,i,j,k,8)
             !              bsquare=v2x*v2x+v2y*v2y+v2z*v2z
              bcellold=bcell
             if(nambipolar2.eq.1)then
                 bcellold=u(l,i,j,k,2)
                 bcell=v2x*v2x+v2y*v2y+v2z*v2z
              end if

              jcenter(l,i,j,k,1)=computdy(bmagij,nzz,nyy,l,i,j,k,dy)-computdz(bmagij,nyy,nzz,l,i,j,k,dy)
              jcenter(l,i,j,k,2)=computdz(bmagij,nxx,nzz,l,i,j,k,dy)-computdx(bmagij,nzz,nxx,l,i,j,k,dy)
              jcenter(l,i,j,k,3)=computdx(bmagij,nyy,nxx,l,i,j,k,dy)-computdy(bmagij,nxx,nyy,l,i,j,k,dy)

              call crossprod(jcenter,u(:,:,:,:,6:8),jxb,l,i,j,k)

              jxbsquare(l,i,j,k)=(jxb(l,i,j,k,1)*jxb(l,i,j,k,1)+jxb(l,i,j,k,2)*jxb(l,i,j,k,2)+jxb(l,i,j,k,3)*jxb(l,i,j,k,3))*&
              & betaad(u(l,i,j,k,1),bcell,tcell,ionisrate)*dtlim



           end do
        end do
     end do
  end do

end SUBROUTINE computambip
!########################################################################################
!########################################################################################
!########################################################################################
!########################################################################################
#if HALL==1
subroutine computvhall(q,dx,dy,dz,ngrid,bpred,rppred,vhall)
  USE amr_parameters
  use hydro_commons
  USE const
  implicit none
  integer :: l,i,j,k,kk,m
  integer :: imm1,jmm1,kmm1
  integer :: ngrid
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,3,3)::bpred
  REAL(dp),DIMENSION(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,2)::rppred
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q 
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:4,1:3)::vhall
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::bijbis
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::jx,jy,jz
  real(dp) :: dx,dy,dz
  real(dp) :: rhoedge,nrjedge,tedge,Bedge,bedge2
  real(dp) :: eta,ionisrate
  real(dp) :: eta_hall_chimie

!!!!!!!!!!!!!!!!!!
!
! bijbis(l,i,j,k,n) is the value of the magnetic field component
! Bn at i-1/2,j-1/2,k-1/2
!
!!!!!!!!!!!!!!!!!!

do k=min(1,ku1+1),ku2
     do j=min(1,ju1+1),ju2
        do i=iu1,iu2
           do l=1,ngrid
              bijbis(l,i,j,k,1)=0.25d0*(bpred(l,i,j,k,1,1)+bpred(l,i,j-1,k,1,1)+bpred(l,i,j,k-1,1,1)+bpred(l,i,j-1,k-1,1,1))
           end do
        end do
     end do
  end do

  ! case By for Lorentz force EMF 
  do k=min(1,ku1+1),ku2
     do j=ju1,ju2
        do i=min(1,iu1+1),iu2
           do l=1,ngrid
              bijbis(l,i,j,k,2)=0.25d0*(bpred(l,i,j,k,2,2)+bpred(l,i-1,j,k,2,2)+bpred(l,i,j,k-1,2,2)+bpred(l,i-1,j,k-1,2,2))
           end do
        end do
     end do
  end do
 
  ! case Bz for Lorentz force EMF 
  do k=ku1,ku2
     do j=min(1,ju1+1),ju2
        do i=min(1,iu1+1),iu2
           
           do l=1,ngrid

              bijbis(l,i,j,k,3)=0.25d0*(bpred(l,i,j,k,3,3)+bpred(l,i,j-1,k,3,3)+bpred(l,i-1,j,k,3,3)+bpred(l,i-1,j-1,k,3,3))
 
           end do
        end do
     end do
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! computation of the component of j where EMFs are located
! jemfx(l,i,j,k,n) is the component Jn at i,j-1/2,k-1/2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do k=min(1,ku1+1),max(1,ku2-1)
     do j=min(1,ju1+1),max(1,ju2-1)
        do i=min(1,iu1+1),max(1,iu2-1)
           
           do l=1,ngrid

              jx(l,i,j,k,1)=(bpred(l,i,j,k,3,3)-bpred(l,i,j-1,k,3,3))/dy-(bpred(l,i,j,k,2,2)-bpred(l,i,j,k-1,2,2))/dz 
              jx(l,i,j,k,2)=(bpred(l,i,j,k,2,1)-bpred(l,i,j,k-1,2,1))/dz- (bijbis(l,i+1,j,k,3)-bijbis(l,i,j,k,3))/dx
              jx(l,i,j,k,3)=(bijbis(l,i+1,j,k,2) -bijbis(l,i,j,k,2))/dx- (bpred(l,i,j,k,3,1)-bpred(l,i,j-1,k,3,1))/dy

              jy(l,i,j,k,1)=(bijbis(l,i,j+1,k,3)-bijbis(l,i,j,k,3))/dy-(bpred(l,i,j,k,1,2) - bpred(l,i,j,k-1,1,2) )/dz
              jy(l,i,j,k,2)=(bpred(l,i,j,k,1,1)-bpred(l,i,j,k-1,1,1))/dz-(bpred(l,i,j,k,3,3)-bpred(l,i-1,j,k,3,3))/dx
              jy(l,i,j,k,3)=(bpred(l,i,j,k,3,2)-bpred(l,i-1,j,k,3,2))/dx-(bijbis(l,i,j+1,k,1)-bijbis(l,i,j,k,1))/dy

              jz(l,i,j,k,1)=(bpred(l,i,j,k,1,3) -bpred(l,i,j-1,k,1,3))/dy-(bijbis(l,i,j,k+1,2)-bijbis(l,i,j,k,2))/dz
              jz(l,i,j,k,2)=( bijbis(l,i,j,k+1,1)-bijbis(l,i,j,k,1))/dz-(bpred(l,i,j,k,2,3)-bpred(l,i-1,j,k,2,3))/dx
              jz(l,i,j,k,3)=(bpred(l,i,j,k,2,2)-bpred(l,i-1,j,k,2,2))/dx-(bpred(l,i,j,k,1,1)-bpred(l,i,j-1,k,1,1))/dy

        end do
     end do
  end do
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! vhall(l,i,j,k,1,n) is nth componement of vhall at i,j-1/2,k-1/2
do k=min(1,ku1+1),max(1,ku2-1)
  do j=min(1,ju1+1),max(1,ju2-1)
    do i=min(1,iu1+1),max(1,iu2-1)
        
      do l=1,ngrid
        if(nhall==1) then
          do kk=1,3
            if(kk==1) then
              imm1=0
              jmm1=1
              kmm1=1
              bedge=dsqrt(sum((0.25*(bpred(l,i,j,k,2,1:3)+bpred(l,i,j,k,3,1:3)+bpred(l,i,j-1,k,3,1:3)+bpred(l,i,j,k-1,2,1:3)))**2d0))
            else if (kk==2) then
              imm1=1
              jmm1=0
              kmm1=1
              bedge=dsqrt(sum((0.25*(bpred(l,i,j,k,1,1:3)+bpred(l,i,j,k,3,1:3)+bpred(l,i-1,j,k,3,1:3)+bpred(l,i,j,k-1,1,1:3)))**2d0))
            else
              imm1=1
              jmm1=1
              kmm1=0
              bedge=dsqrt(sum((0.25*(bpred(l,i,j,k,2,1:3)+bpred(l,i,j,k,1,1:3)+bpred(l,i,j-1,k,1,1:3)+bpred(l,i-1,j,k,2,1:3)))**2d0))
            end if

            ! rppred (l,i,j,k,1) is the density in the cell, and rppred(l,i,j,k,2) its pressure
            rhoedge=0.25*(rppred(l,i,j,k,1)+rppred(l,i-imm1,j-jmm1,k,1)+rppred(l,i-imm1,j,k-kmm1,1)+rppred(l,i,j-jmm1,k-kmm1,1))
            nrjedge=0.25*(rppred(l,i,j,k,2)+rppred(l,i-imm1,j-jmm1,k,2)+rppred(l,i-imm1,j,k-kmm1,2)+rppred(l,i,j-jmm1,k-kmm1,2))

            ! We take non-updated ionisation rate !!
            ionisrate=1d-17
            bedge2=bedge*bedge
            call temperature_eos(rhoedge,nrjedge/(gamma-1),tedge)
            eta=eta_hall_chimie(rhoedge,tedge,ionisrate,bedge2)/bedge

            if (kk==1)  vhall(l,i,j,k,1,1:3)=jx(l,i,j,k,1:3)*eta
            if (kk==2)  vhall(l,i,j,k,2,1:3)=jy(l,i,j,k,1:3)*eta
            if (kk==3)  vhall(l,i,j,k,3,1:3)=jz(l,i,j,k,1:3)*eta
            vhall(l,i,j,k,4,kk)=ionisrate
           end do
         else
           vhall(l,i,j,k,:,:)=0d0
         end if

       end do
     end do
  end do
end do

end subroutine computvhall
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
!###########################################################

! fonctions de produits vectoriels et coef nimhd

!###########################################################
double precision  function computdx(vec,n2,n3,l,i,j,k,dx)

   use hydro_parameters
   implicit none 
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:3)::vec
   real(dp)::dx
   integer::n2,n3,l,i,j,k

   computdx=(vec(l,i+1,j,k,n2,n3)-vec(l,i,j,k,n2,n3))/dx

end function computdx
!###########################################################
!###########################################################
!###########################################################
double precision  function computdy(vec,n2,n3,l,i,j,k,dx)

   use hydro_parameters
   implicit none 
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:3)::vec
   real(dp)::dx
   integer::n2,n3,l,i,j,k

   computdy=(vec(l,i,j+1,k,n2,n3)-vec(l,i,j,k,n2,n3))/dx

end function computdy
!###########################################################
!###########################################################
!###########################################################
double precision  function computdz(vec,n2,n3,l,i,j,k,dx)

   use hydro_parameters
   implicit none 
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:3)::vec
   real(dp)::dx
   integer::n2,n3,l,i,j,k

   computdz=(vec(l,i,j,k+1,n2,n3)-vec(l,i,j,k,n2,n3))/dx

end function computdz
!###########################################################
!###########################################################
!###########################################################
double precision  function computdxbis(vec,n2,l,i,j,k,dx)

   use hydro_parameters
   implicit none 
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::vec
   real(dp)::dx
   integer::n2,l,i,j,k

   computdxbis=(vec(l,i+1,j,k,n2)-vec(l,i,j,k,n2))/dx

end function computdxbis
!###########################################################
!###########################################################
!###########################################################
double precision  function computdybis(vec,n2,l,i,j,k,dx)

   use hydro_parameters
   implicit none 
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::vec
   real(dp)::dx
   integer::n2,l,i,j,k

   computdybis=(vec(l,i,j+1,k,n2)-vec(l,i,j,k,n2))/dx

end function computdybis
!###########################################################
!###########################################################
!###########################################################
double precision  function computdzbis(vec,n2,l,i,j,k,dx)

   use hydro_parameters
   implicit none
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::vec
   real(dp)::dx
   integer::n2,l,i,j,k

   computdzbis=(vec(l,i,j,k+1,n2)-vec(l,i,j,k,n2))/dx

end function computdzbis
!###########################################################
!###########################################################
!###########################################################
double precision function gammaadbis(rhon,BBcell,BBcellold,temper,ionisrate)

   use hydro_parameters
   use amr_parameters,only:mu_gas
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

   n_H_max = 2.5d+17

   ! C shock Duffin et Pudritz
   ! gammaadbis in CGS
   !gammaadbis=gammaAD

   !rhoH=rhon*xmolaire*H2_fraction*scale_d/(mu_gas*mH) ! convert in H/cc
   rhoH=rhon*2.0d0*H2_fraction*scale_d/(mu_gas*mH) ! convert in H/cc

   if(rhoH < n_H_max)then
      gammaadbis=eta_AD_chimie(rhoH,BBcell,BBcellold,temper,ionisrate)
   else
      gammaadbis=eta_AD_chimie(n_H_max,BBcell,BBcellold,temper,ionisrate)
   endif

   gammaadbis=gammaadbis*scale_t*scale_d ! in code units

   ! test
   !gammaadbis=gammaAD

end function gammaadbis
!###########################################################
!###########################################################
!###########################################################
subroutine sig_x2d(ll,ii,j,k,lb,ib,sigO,sigH,sigP,bsquare)

   use amr_parameters,    only : dp
   use nimhd_commons
   use amr_commons, only : myid
   implicit none

   integer, intent(in)             :: j,k,ib
   real(dp)                        :: B,nH,temper,sigav
   real(dp)                        :: j_dp,k_dp,b_dp
   real(dp), dimension(nvarchimie) :: x
   real(dp), intent(in)            :: ll,ii,lb,bsquare
   real(dp), intent(out)           :: sigO,sigH,sigP
   integer                         :: i,kk

   j_dp = real(j,dp)
   kk=min(k,tchimie-1)
   k_dp = real(kk,dp)
   b_dp = real(ib,dp)

   x(1:3)=(1d0-(ll-j_dp))*(1d0-(ii-k_dp))*(1d0-(lb-b_dp))*(resistivite_chimie(1:3,j,kk,ib,1))+&
            &((ll-j_dp))*(1d0-(ii-k_dp))*(1d0-(lb-b_dp))*(resistivite_chimie(1:3,j+1,kk,ib,1))+&
            &(1d0-(ll-j_dp))*((ii-k_dp))*(1d0-(lb-b_dp))*(resistivite_chimie(1:3,j,kk+1,ib,1))+&
                  &((ll-j_dp))*((ii-k_dp))*(1d0-(lb-b_dp))*(resistivite_chimie(1:3,j+1,kk+1,ib,1))+&
               (1d0-(ll-j_dp))*(1d0-(ii-k_dp))*(lb-b_dp)*(resistivite_chimie(1:3,j,kk,ib+1,1))+&
                     &((ll-j_dp))*(1d0-(ii-k_dp))*(lb-b_dp)*(resistivite_chimie(1:3,j+1,kk,ib+1,1))+&
                     &(1d0-(ll-j_dp))*((ii-k_dp))*(lb-b_dp)*(resistivite_chimie(1:3,j,kk+1,ib+1,1))+&
                        &((ll-j_dp))*((ii-k_dp))*(lb-b_dp)*(resistivite_chimie(1:3,j+1,kk+1,ib+1,1))
               
   sigP= 10.0d0**x(1)
   sigO= 10.0d0**x(2)

   ! modification since x(3) can be negative we simply use the sign of the leftmost
   ! point. If there is a sign inversion, we set it to zero.
   ! If you are using Hall resisitvities, this could be improved by using a linear
   ! interpolation instead of log.
   sigH=(10.0d0**x(3))*resistivite_chimie(0,j,kk,ib,1)
   sigav = sum(resistivite_chimie(0,j:j+1,kk:kk+1,ib:ib+1,1)) / 8.0d0
   if(sigav .ne. resistivite_chimie(0,j,kk,ib,1))then
      sigH = 0.0_dp
      !if(myid==1) write(*,*)'Sign inversion in Hall resistivity'
   endif

   return

end subroutine sig_x2d
!###########################################################
!###########################################################
!###########################################################
subroutine sig_x3d(ll,ii,xx,j,k,xi,lb,ib,sigO,sigH,sigP,bsquare)

   use amr_parameters,    only : dp
   use nimhd_commons
   use amr_commons, only : myid
   implicit none

   integer, intent(in)             :: j,k,xi,ib
   real(dp)                        :: B,nH,temper,sigav
   real(dp)                        :: j_dp,k_dp,xi_dp,b_dp
   real(dp), dimension(0:3)        :: x
   real(dp), intent(in)            :: ll,ii,xx,lb,bsquare
   real(dp), intent(out)           :: sigO,sigH,sigP
   integer                         :: i,kk

   j_dp = real(j,dp)
   kk=min(k,tchimie-1)
   k_dp = real(kk,dp)
   xi_dp=real(xi,dp)
   b_dp = real(ib,dp)

   x(0:3)=(1d0-(ll-j_dp))*(1d0-(ii-k_dp))*(1d0-(xx-xi_dp))*(1d0-(lb-b_dp))*(resistivite_chimie(0:3,j,kk,xi,ib))+&
            &((ll-j_dp))*(1d0-(ii-k_dp))*(1d0-(xx-xi_dp))*(1d0-(lb-b_dp))*(resistivite_chimie(0:3,j+1,kk,xi,ib))+&
            &(1d0-(ll-j_dp))*((ii-k_dp))*(1d0-(xx-xi_dp))*(1d0-(lb-b_dp))*(resistivite_chimie(0:3,j,kk+1,xi,ib))+&
                  &((ll-j_dp))*((ii-k_dp))*(1d0-(xx-xi_dp))*(1d0-(lb-b_dp))*(resistivite_chimie(0:3,j+1,kk+1,xi,ib))+&
            &(1d0-(ll-j_dp))*(1d0-(ii-k_dp))*((xx-xi_dp))*(1d0-(lb-b_dp))*(resistivite_chimie(0:3,j,kk,xi+1,ib))+&
                  &((ll-j_dp))*(1d0-(ii-k_dp))*((xx-xi_dp))*(1d0-(lb-b_dp))*(resistivite_chimie(0:3,j+1,kk,xi+1,ib))+&
                  &(1d0-(ll-j_dp))*((ii-k_dp))*((xx-xi_dp))*(1d0-(lb-b_dp))*(resistivite_chimie(0:3,j,kk+1,xi+1,ib))+&
                        &((ll-j_dp))*((ii-k_dp))*((xx-xi_dp))*(1d0-(lb-b_dp))*(resistivite_chimie(0:3,j+1,kk+1,xi+1,ib))+&
               (1d0-(ll-j_dp))*(1d0-(ii-k_dp))*(1d0-(xx-xi_dp))*((lb-b_dp))*(resistivite_chimie(0:3,j,kk,xi,ib+1))+&
                  &((ll-j_dp))*(1d0-(ii-k_dp))*(1d0-(xx-xi_dp))*((lb-b_dp))*(resistivite_chimie(0:3,j+1,kk,xi,ib+1))+&
                  &(1d0-(ll-j_dp))*((ii-k_dp))*(1d0-(xx-xi_dp))*((lb-b_dp))*(resistivite_chimie(0:3,j,kk+1,xi,ib+1))+&
                        &((ll-j_dp))*((ii-k_dp))*(1d0-(xx-xi_dp))*((lb-b_dp))*(resistivite_chimie(0:3,j+1,kk+1,xi,ib+1))+&
                  &(1d0-(ll-j_dp))*(1d0-(ii-k_dp))*((xx-xi_dp))*((lb-b_dp))*(resistivite_chimie(0:3,j,kk,xi+1,ib+1))+&
                        &((ll-j_dp))*(1d0-(ii-k_dp))*((xx-xi_dp))*((lb-b_dp))*(resistivite_chimie(0:3,j+1,kk,xi+1,ib+1))+&
                        &(1d0-(ll-j_dp))*((ii-k_dp))*((xx-xi_dp))*((lb-b_dp))*(resistivite_chimie(0:3,j,kk+1,xi+1,ib+1))+&
                           &((ll-j_dp))*((ii-k_dp))*((xx-xi_dp))*((lb-b_dp))*(resistivite_chimie(0:3,j+1,kk+1,xi+1,ib+1))

               
   sigP= 10.0d0**x(1)
   sigO= 10.0d0**x(2)

   ! modification since x(3) can be negative we simply use the sign of the leftmost
   ! point. If there is a sign inversion, we set it to zero.
   ! If you are using Hall resisitvities, this could be improved by using a linear
   ! interpolation instead of log.
   sigH=(10.0d0**x(3))*sign(1d0,x(0))
   sigav = sum(resistivite_chimie(0,j:j+1,kk:kk+1,xi:xi+1,ib:ib+1)) / 16.0d0
   if(sigav .ne. resistivite_chimie(0,j,kk,xi,ib))then
      sigH = 0.0_dp
      !if(myid==1) write(*,*)'Sign inversion in Hall resistivity'
   endif

   return

end subroutine sig_x3d
!###########################################################
!###########################################################
!###########################################################
double precision function eta_AD_chimie(rhon,BBcell,BBcellold,temper,ionisrate)

   use hydro_commons
   use constants
   use nimhd_commons
   implicit none

   real(dp)     :: sigO,sigH,sigP,densionbis,BBcgs, bbcell,BBcellold
   real(dp)::inp,ll,rhon,ii,temper,lb,j_dp,xx
   integer :: i,j,k,ib
   logical :: notfound
   real(dp)::ionisrate

   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

   if(use_res==1)then
      inp=rhon
      ll=(1d0+(log10(inp)-log10(nminchimie))/dnchimie)
      j=floor(ll)
      j_dp=real(j,dp)
      !ll=(1d0+(log10(inp)-log10(300d0))/(17d0/35d0))
      !j=dble(floor(ll))
      eta_AD_chimie=(ll-j_dp)*log10(resistivite_chimie_res(6,j+1))+(1d0-(ll-j_dp))*log10(resistivite_chimie_res(6,j))
      eta_AD_chimie=10**eta_AD_chimie
      !print*, rhon,temper,eta_AD_chimie
      !print*, eta_AD_chimie, inp,(1.43d-7*sqrt(inp))**2
      !stop
   else if(use_x2d==1)then
      inp=rhon
      ll=(1d0+(log10(inp)-log10(nminchimie))/dnchimie)
      j=floor(ll)
      inp=temper
      ii=(1d0+(log10(inp)-log10(tminchimie))/dtchimie)
      ii=max(ii,1.0d0)
      !ii=(1d0+(log10(inp)-log10(5d0))/(3d0/50d0))
      i=floor(ii)
      BBcgs=sqrt(BBcellold*(4d0*pi*scale_d*(scale_v)**2))
      !!$   bbcgs=1.43d-7*sqrt(rhon/2d0/H2_fraction)
      !!$   print*, bbcgs, sqrt(BBcellold*(4d0*pi*scale_d*(scale_v)**2)),rhon
      inp=BBcgs
      lb=(1d0+(log10(inp)-log10(bminchimie))/dbchimie)
      ib=floor(lb)

      call sig_x2d(ll,ii,j,i,lb,ib,sigO,sigH,sigP,BBcgs) 
      !inp=rhon/xmolaire/H2_fraction     ! inp is neutrals.cc, to fit densionbis
      inp=rhon/2d0/H2_fraction     ! inp is neutrals.cc, to fit densionbis
      eta_AD_chimie=(sigO/(sigO**2+sigH**2)-1d0/sigP)   ! resistivity in s
      !print*,   eta_AD_chimie,inp*xmolaire*H2_fraction,BBcgs,inp,densionbis(inp),scale_d
      !print*, eta_AD_chimie,inp*xmolaire*H2_fraction,BBcgs

      BBcgs=sqrt(BBcell*(4d0*pi*scale_d*(scale_v)**2))
      !!$   (eta_AD_chimie = max(eta_AD_chimie * (1.0d0-tanh(ii/(dble(tchimie))), 1d-36)

      eta_AD_chimie=BBcgs**2/(eta_AD_chimie*densionbis(inp)*inp*scale_d*scale_d*c_cgs**2)  ! need B in G, output is gammaad in cgs

   else if(use_x3d==1)then
      ll=(1d0+(log10(rhon)-log10(nminchimie))/dnchimie)
      j=floor(ll)
      !ii=(1d0+(log10(temper)-log10(tminchimie))/dtchimie)
      !    ii=max(ii,1.0d0)
      !i=floor(ii)
      ii = max(0d0, (log10(temper)-log10(tminchimie))/dtchimie)
      i=floor(1d0+ii)

      xx=(1d0+(log10(ionisrate)-log10(ximinchimie))/dxichimie)
      k=floor(xx)

      BBcgs=sqrt(BBcellold*(4d0*pi*scale_d*(scale_v)**2))
      lb=(1d0+(log10(BBcgs)-log10(bminchimie))/dbchimie)
      ib=floor(lb)

      call sig_x3d(ll,ii,xx,j,i,k,lb,ib,sigO,sigH,sigP,BBcgs) 
      inp=rhon/2d0/H2_fraction     ! inp is neutrals.cc, to fit densionbis
      eta_AD_chimie=(sigO/(sigO**2+sigH**2)-1d0/sigP)   ! resistivity in s

      BBcgs=sqrt(BBcell*(4d0*pi*scale_d*(scale_v)**2))

      eta_AD_chimie=BBcgs**2/(eta_AD_chimie*densionbis(inp)*inp*scale_d*scale_d*c_cgs**2)  ! need B in G, output is gammaad in cgs

   endif

   !!print*, inp, ll, j,resistivite_chimie(1,1),resistivite_chimie(1,2),resistivite_chimie(6,1),resistivite_chimie(1,35)
   !eta_AD_chimie=(ll-j)*log10(resistivite_chimie(6,j+1))+(1d0-(ll-j))*log10(resistivite_chimie(6,j))
   !!print*, eta_AD_chimie
   !eta_AD_chimie=10**eta_AD_chimie

   ! Ad-hoc modification to ensure that the ambipolar resistivity falls to zero when the density exceeds 5.0e13
   !eta_AD_chimie = eta_AD_chimie * (1.0d0-tanh(rhon/5.0d13))

end function eta_AD_chimie
!###########################################################
!###########################################################
!###########################################################
double precision function eta_ohm_chimie(rhon,BBcell,temper,ionisrate)

   use hydro_commons
   use constants
   use nimhd_commons
   implicit none

   real(dp) :: inp,ll,ii,lb,rhon,BBcell
   real(dp) :: temper,sigO,sigH,sigP,BBcgs
   real(dp) :: j_dp,ionisrate,xx
   integer  :: j,i,ib,k

   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

   if(use_res==1)then
      inp=rhon
      ll=(1d0+(log10(inp)-log10(nminchimie))/dnchimie)
      j=floor(ll)
      j_dp=real(j,dp)
      eta_ohm_chimie=(ll-j_dp)*log10(resistivite_chimie_res(7,j+1))+(1d0-(ll-j_dp))*log10(resistivite_chimie_res(7,j))
      eta_ohm_chimie=10.0d0**eta_ohm_chimie
      eta_ohm_chimie = max(eta_ohm_chimie * (1.0d0-tanh(rhon/1.0d15)), 1d-36)
   else if(use_x2d==1)then
      inp=rhon
      ll=(1d0+(log10(inp)-log10(nminchimie))/dnchimie)
      j=floor(ll)
      inp=temper
      ii=(1d0+(log10(inp)-log10(tminchimie))/dtchimie)
      ii=max(ii,1.0d0)
      i=floor(ii)
      BBcgs=sqrt(BBcell*(4d0*pi*scale_d*(scale_v)**2))
      inp=BBcgs
      lb=(1d0+(log10(inp)-log10(bminchimie))/dbchimie)
      ib=floor(lb)
      call sig_x2d(ll,ii,j,i,lb,ib,sigO,sigH,sigP,BBcgs)
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
      call sig_x3d(ll,ii,xx,j,i,k,lb,ib,sigO,sigH,sigP,BBcgs)
      eta_ohm_chimie = (1d0 / sigP) * c_cgs * c_cgs / (4.0_dp*pi)
   endif

   ! Ad-hoc modification to ensure that the ohmic resistivity falls to zero when the density exceeds 1.0e15
   ! when alkali metals are ionized.
   !eta_ohm_chimie = eta_ohm_chimie * (1.0d0-tanh(rhon/1.0d15))
   ! eta_ohm_chimie = max(eta_ohm_chimie * (1.0d0-tanh(rhon/1.0d15)), 1d-36)

end function eta_ohm_chimie
!###########################################################
!###########################################################
!###########################################################
#if HALL==1
double precision function eta_hall_chimie(rhon,temper,ionisrate,BBcell)

   use hydro_commons
   use constants
   use nimhd_commons,ONLY:dtchimie,dnchimie,nminchimie,tminchimie,ximinchimie,&
                     &dbchimie,bminchimie,nchimie,tchimie,xichimie,dxichimie,&
                     &bchimie
   implicit none

   real(dp) :: inp,ll,ii,lb,rhon,BBcell
   real(dp) :: temper,sigO,sigH,sigP,BBcgs
   real(dp) :: j_dp,ionisrate,xx
   integer  :: j,i,ib,k

   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

   if(ntestDADM==1) then
      eta_hall_chimie=rhall
   else
      if(use_x2d==1)then
         inp=rhon
         ll=(1d0+(log10(inp)-log10(nminchimie))/dnchimie)
         j=floor(ll)
         inp=temper
         ii=(1d0+(log10(inp)-log10(tminchimie))/dtchimie)
         ii=max(ii,1.0d0)
         i=floor(ii)
         BBcgs=sqrt(BBcell*(4d0*pi*scale_d*(scale_v)**2))
         inp=BBcgs
         lb=(1d0+(log10(inp)-log10(bminchimie))/dbchimie)
         ib=floor(lb)
         call sig_x2d(ll,ii,j,i,lb,ib,sigO,sigH,sigP,BBcgs)
         eta_hall_chimie=sigH/(sigO**2+sigH**2) 
         eta_hall_chimie=eta_hall_chimie * c_cgs * c_cgs / (4.0_dp*pi)
         eta_hall_chimie=eta_hall_chimie*scale_t/(scale_l*scale_l)
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
         call sig_x3d(ll,ii,xx,j,i,k,lb,ib,sigO,sigH,sigP,BBcgs)
         eta_hall_chimie=sigH/(sigO**2+sigH**2) 
         eta_hall_chimie=eta_hall_chimie * c_cgs * c_cgs / (4.0_dp*pi)
         eta_hall_chimie=eta_hall_chimie*scale_t/(scale_l*scale_l)
      endif
   endif
end function eta_hall_chimie
#endif
!###########################################################
!###########################################################
!###########################################################
double precision function densionbis(rhon)

   use hydro_parameters, only : coefionis,default_ionisrate,ntestDADM,rhoi0,dp

   implicit none 
   real(dp)::rhon
   real(dp)::xn, rhoncgs

   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

   ! density of neutral in g/cm3  
   rhoncgs=rhon*scale_d

   ! function which computes the density in g/cm3 of ions 
   ! see Duffin & Pudritz 2008, astro-ph 08/10/08 eq (14)

   ! density of neutral in number per cm3
   !xn=rhoncgs/xmneutre

   ! density of ions in g/cm3 
   !densionbis=densionbis*xmion


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
      densionbis=rhoi0
   endif

end function densionbis
!###########################################################
!###########################################################
!###########################################################
double precision function computdxvx(vec,l,i,j,k,dx,dy,dz)

   use hydro_parameters

   implicit none
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::vec
   integer l,i,j,k
   real(dp)::dx,dy,dz

   computdxvx=(vec(l,i+1,j,k,nxx)-vec(l,i,j,k,nxx))/dx 

end function computdxvx
!###########################################################
!###########################################################
!###########################################################
double precision function computdyvy(vec,l,i,j,k,dx,dy,dz)

   use hydro_parameters

   implicit none
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::vec
   integer l,i,j,k
   real(dp)::dx,dy,dz

   computdyvy=(vec(l,i,j+1,k,nyy)-vec(l,i,j,k,nyy))/dy 

end function computdyvy
!###########################################################
!###########################################################
!###########################################################
double precision function computdzvz(vec,l,i,j,k,dx,dy,dz)

   use hydro_parameters

   implicit none
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::vec
   integer l,i,j,k
   real(dp)::dx,dy,dz

   computdzvz=(vec(l,i,j,k+1,nzz)-vec(l,i,j,k,nzz))/dz

end function computdzvz
!###########################################################
!###########################################################
!###########################################################
double precision function computdiv(vec,l,i,j,k,dx,dy,dz)

   use hydro_parameters

   implicit none
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::vec
   integer l,i,j,k
   real(dp)::dx,dy,dz

   if(ndim.eq.1) then
      computdiv=(vec(l,i+1,j,k,nxx)-vec(l,i,j,k,nxx))/dx 
   endif
   if(ndim.eq.2) then
      computdiv=(vec(l,i+1,j,k,nxx)-vec(l,i,j,k,nxx))/dx + (vec(l,i,j+1,k,nyy)-vec(l,i,j,k,nyy))/dy 
   endif
   if(ndim.eq.3) then
      computdiv=(vec(l,i+1,j,k,nxx)-vec(l,i,j,k,nxx))/dx + (vec(l,i,j+1,k,nyy)-vec(l,i,j,k,nyy))/dy + (vec(l,i,j,k+1,nzz)-vec(l,i,j,k,nzz))/dz
   endif

end function computdiv
!###########################################################
!###########################################################
!###########################################################
double precision function computdivbisx(vec,l,i,j,k,dx,dy,dz)

   use hydro_parameters

   implicit none
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::vec
   integer l,i,j,k
   real(dp)::dx,dy,dz
   real(dp)::computdxvx,computdyvy,computdzvz

   if(ndim.eq.1) then
      computdivbisx=computdxvx(vec,l,i,j,k,dx,dy,dz)
   endif
   if(ndim.eq.2) then
      computdivbisx=computdxvx(vec,l,i,j,k,dx,dy,dz)+computdyvy(vec,l,i,j-1,k,dx,dy,dz)
   endif
   if(ndim.eq.3) then
      computdivbisx=computdxvx(vec,l,i,j,k,dx,dy,dz)+computdyvy(vec,l,i,j-1,k,dx,dy,dz)+computdzvz(vec,l,i,j,k-1,dx,dy,dz)
   endif

end function computdivbisx
!###########################################################
!###########################################################
!###########################################################
double precision function computdivbisy(vec,l,i,j,k,dx,dy,dz)

   use hydro_parameters

   implicit none
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::vec
   integer l,i,j,k
   real(dp)::dx,dy,dz
   real(dp)::computdxvx,computdyvy,computdzvz

   if(ndim.eq.1) then
      computdivbisy=computdxvx(vec,l,i-1,j,k,dx,dy,dz)
   endif
   if(ndim.eq.2) then
      computdivbisy=computdxvx(vec,l,i-1,j,k,dx,dy,dz)+computdyvy(vec,l,i,j,k,dx,dy,dz)
   endif
   if(ndim.eq.3) then
      computdivbisy=computdxvx(vec,l,i-1,j,k,dx,dy,dz)+computdyvy(vec,l,i,j,k,dx,dy,dz)+computdzvz(vec,l,i,j,k-1,dx,dy,dz)
   endif

end function computdivbisy
!###########################################################
!###########################################################
!###########################################################
double precision function computdivbisz(vec,l,i,j,k,dx,dy,dz)

   use hydro_parameters

   implicit none
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::vec
   integer l,i,j,k
   real(dp)::dx,dy,dz
   real(dp)::computdxvx,computdyvy,computdzvz

   if(ndim.eq.1) then
      computdivbisz=computdxvx(vec,l,i-1,j,k,dx,dy,dz)
   endif
   if(ndim.eq.2) then
      computdivbisz=computdxvx(vec,l,i-1,j,k,dx,dy,dz)+computdyvy(vec,l,i,j-1,k,dx,dy,dz)
   endif
   if(ndim.eq.3) then
      computdivbisz=computdxvx(vec,l,i-1,j,k,dx,dy,dz)+computdyvy(vec,l,i,j-1,k,dx,dy,dz)+computdzvz(vec,l,i,j,k,dx,dy,dz)
   endif

end function computdivbisz
!###########################################################
!###########################################################
!###########################################################
subroutine crossprodbis(vec1,vec2,v1crossv2,l,i,j,k)

   use hydro_parameters

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
!###########################################################
!###########################################################
!###########################################################
subroutine crossprod(vec1,vec2,v1crossv2,l,i,j,k)

   use hydro_parameters

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
!###########################################################
!###########################################################
double precision function  crossprodx(v1x,v1y,v1z,v2x,v2y,v2z)

   ! function which gives the x component of a cross product of two
   ! vectors of coordinates v1x,v1y,v1z,v2x,v2y,v2z
   use hydro_parameters
   implicit none
   real(dp)::v1x,v1y,v1z,v2x,v2y,v2z

   crossprodx=v1y*v2z-v1z*v2y

end function crossprodx
!###########################################################
!###########################################################
!###########################################################
double precision function crossprody(v1x,v1y,v1z,v2x,v2y,v2z)

   ! function which gives the y component of a cross product of two
   ! vectors of coordinates v1x,v1y,v1z,v2x,v2y,v2z
   use hydro_parameters
   implicit none
   real(dp)::v1x,v1y,v1z,v2x,v2y,v2z

   crossprody=v1z*v2x-v1x*v2z

end function crossprody
!###########################################################
!###########################################################
!###########################################################
double precision function crossprodz(v1x,v1y,v1z,v2x,v2y,v2z)

   ! function which gives the z component of a cross product of two
   ! vectors of coordinates v1x,v1y,v1z,v2x,v2y,v2z
   use hydro_parameters
   implicit none
   real(dp)::v1x,v1y,v1z,v2x,v2y,v2z

   crossprodz=v1x*v2y-v2x*v1y

end function crossprodz
!###########################################################
!###########################################################
!###########################################################
double precision function muvisco(rhon)

   use hydro_parameters
   implicit none 
   real(dp) ::rhon

   muvisco=visco
   if(ntestDADM.eq.1) then
      muvisco=visco
   endif

end function muvisco
!###########################################################
!###########################################################
!###########################################################
double precision function etaohmdiss(rhon,BBcell,temper,ionisrate)

   use hydro_commons
   use amr_parameters,only:mu_gas
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
double precision function etaohmdissbricolo(rhon,BBcell,temper,dtlim,dx,ionisrate)

   use hydro_commons
   use amr_parameters,only:mu_gas
   use constants

   implicit none 
   real(dp) ::rhon,xpressure,rhoH,rhotemp,BBcell
   real(dp)::gammaadbis,densionbis,ionisrate
   real(dp)::xionisation,temper,scale_p,xpcgs,rhocgs,xnbcgs,n_H_max
   real(dp)::eta_ohm_chimie,dx,dtlim,xx,dtt

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
         etaohmdissbricolo=eta_ohm_chimie(rhotemp,BBcell,temper,ionisrate)
      else
         etaohmdissbricolo=eta_ohm_chimie(n_H_max,BBcell,temper,ionisrate)
      endif

      etaohmdissbricolo=etaohmdissbricolo*scale_t/(scale_l)**2

      ! robbery to avoid too small time step
      if(nminitimestep.eq.1 .and. nmagdiffu2.eq.0) then
         if(dtlim.ne.0d0) then
            xx=etaohmdissbricolo
            if(xx.ne.0d0) then
            dtt=coefohm*dx*dx/xx   !dtohm pour la cellule
         !    if (myid ==1) print*, dtt,bsquare,betaadbricolo,betaadbricolotemp
            else
               dtt=1d39
            endif
            if (dtt.le.dtlim) then
               etaohmdissbricolo=coefohm*dx*dx/(dtlim)
            endif
         endif
      endif
      
   !   etaohmdissbricolo=etaohmdissbricolo*scale_t/(scale_l)**2

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
   implicit none

   real(dp) ::rhon,rhotemp,bsquare,temper
   real(dp)::gammaadbis,densionbis,ionisrate
   real(dp)::xx

   if(ntestDADM.eq.0) then
      ! function which computes the coefficient beta which
      ! appears in ambipolar diffusion dB/dt=curl(gamma(j*B)*B)+...
      ! see Duffin & Pudritz 2008, astro-ph 08/10/08 eq (5)
      ! WARNING no mu_0 needed here because F_Lorentz used

      ! Warning gammaadbis and densionbis already in user units
      ! but NOT rhon/xmneutre

      !betaad=1.4d0/(gammaadbis(rhon)*densionbis(rhon)*rhon/xmneutre )
      ! no xmneutre for Duffin and Pudritz

      rhotemp = MAX(rhon,rho_threshold)

      !xx=gammaadbis(rhotemp,bsquare,temper)*densionbis(rhon)*rhon
      xx=gammaadbis(rhotemp,bsquare,bsquare,temper,ionisrate)*densionbis(rhotemp)*rhotemp
      !xx=gammaadbis(rhotemp)*densionbis(rhotemp)*rhotemp
      if(xx.ne.0d0) then
         betaad=1d0/xx 
      else
         betaad=1d39
         if(rhotemp < 1.0d+14)then
            write(*,*)'WARNING gammaadbis(rhotemp,bsquare,temper,ionisrate)*densionbis(rhon)*rhon equal 0',gammaadbis(rhotemp,bsquare,bsquare,temper,ionisrate),densionbis(rhotemp),rhotemp
         endif
      endif

      ! Barenblatt
      !betaad=1d0

   elseif(ntestDADM.eq.1) then
      ! test Barenblatt
      !betaadbricolo=1d0
      ! test C shock
         betaad=1d0/(gammaAD*rhoi0*rhon)
      !betaadbricolo=0d0
   endif

   !rhon, gammaadbis(rhon) and densionbis(rhon) already in user units
   !!betaad=betaad/scale_d

end function betaad
!###########################################################
!###########################################################
!###########################################################
double precision function betaadbricolo(rhocelln,rhon,dtlim,bsquare,bsquareold,dx,ntest,temper,ionisrate)

   use hydro_parameters
   use amr_commons
   use cooling_module
   use constants

   implicit none
   integer :: ntest
   real(dp) ::rhocelln,rhon,betaadbricolotemp,dtlim,bsquare,bsquareold,dx,temper
   real(dp)::gammaadbis,densionbis,rhotemp,rhotemp_cell,ionisrate
   real(dp)::xx,dtt,bbcgs

   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

   if(ntestDADM.eq.0) then
      ! function which computes the coefficient beta which
      ! appears in ambipolar diffusion dB/dt=curl(gamma(j*B)*B)+...
      ! see Duffin & Pudritz 2008, astro-ph 08/10/08 eq (5)
      ! WARNING no mu_0 needed here because F_Lorentz used

      ! Warning gammaadbis and densionbis already in user units
      ! but NOT rhon/xmneutre

      !betaad=1.4d0/(gammaadbis(rhon)*densionbis(rhon)*rhon/xmneutre )
      ! no xmneutre for Duffin and Pudritz

      rhotemp = MAX(rhon,rho_threshold)
      rhotemp_cell = MAX(rhocelln,rho_threshold)

      !xx=gammaadbis(rhotemp_cell,bsquare,temper)*densionbis(rhocelln)*rhocelln  ! dans la cellule
      xx=gammaadbis(rhotemp_cell,bsquare,bsquareold,temper,ionisrate)*densionbis(rhotemp_cell)*rhotemp_cell  ! dans la cellule

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
      if(nminitimestep.eq.1 .and. nambipolar2.eq.0) then
         if(dtlim.ne.0d0) then
            xx=bsquare*betaadbricolo
            if(xx.ne.0d0) then
            dtt=coefad*dx*dx/xx   !dtAD pour la cellule
            else
               dtt=1d39
            endif
            if (dtt.le.dtlim) then   ! on compare bien dtAD calcule pour la cellule (rhocelln) avec le temps de la simu
               betaadbricolo=coefad*dx*dx/(dtlim*bsquare)
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
         betaadbricolo=1d0/(gammaAD*rhoi0*rhon)
      !betaadbricolo=0d0
   endif

end function betaadbricolo
!###########################################################
!###########################################################
!###########################################################