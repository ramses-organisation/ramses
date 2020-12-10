program rungekutta
implicit none
integer :: i,imax
integer ::ib=1
integer :: recentre=1
integer, parameter :: n=10001
double precision,dimension(n) :: vx,vy,vz,by,bz,rho
double precision,dimension(7) :: state_l,state_r
double precision :: Bx=1d0
double precision :: eta_o,eta_h,eta_ad
double precision :: Co=1.d-9
double precision :: Ch=-2d-2
double precision :: Cad=3.5d-3
double precision :: Q,cs
double precision :: B2,Kx,Ky,Kz
double precision,dimension(2) :: rho0,P0,vx0,vy0,vz0,By0,Bz0
double precision :: Bycur,Bzcur,rhocur,vxcur
double precision :: f11,f12,f21,f22,f31,f32,f41,f42
double precision :: h,drho
character(60) :: infile


!ib=1 -> on part de la gauche
!ib=2 -> on part de la droite

open(20,file='nimhd-cshock-hall-ana.dat',status='replace')

state_l = (/1.7942, 0.017942, -0.9759, -0.6561, 0.0, 1.74885, 0.0/)
state_r = (/1.0, 0.01, -1.751, 0.0, 0.0, 0.6, 0.0/)

rho0(1)=state_l(1) 
P0(1  )=state_l(2) 
vx0(1) =state_l(3) 
vy0(1) =state_l(4) 
vz0(1) =state_l(5) 
By0(1) =state_l(6) 
Bz0(1) =state_l(7) 

rho0(2)=state_r(1) 
P0(2  )=state_r(2) 
vx0(2) =state_r(3) 
vy0(2) =state_r(4) 
vz0(2) =state_r(5) 
By0(2) =state_r(6) 
Bz0(2) =state_r(7) 

Q = rho0(ib)*vx0(ib)

cs = 0.1

Kx=Q*vx0(ib)+Q/vx0(ib)*cs*cs+0.5d0*(Bx**2+By0(ib)**2+Bz0(ib)**2)
Ky=Q*vy0(ib)-Bx*By0(ib)
Kz=Q*vz0(ib)-Bx*Bz0(ib)

By(1)=By0(ib)-1e-4
Bz(1)=Bz0(ib)
rho(1)=rho0(ib)
vx(1)=vx0(ib)
vy(1)=vy0(ib)
vz(1)=vz0(ib)

h = 1d0/dble(n-1)
i=1

drho=0d0
imax=1

do i=2,n

  Bycur=By(i-1)
  Bzcur=Bz(i-1)
  rhocur=rho(i-1)


  call ff(Bycur,Bzcur,f11,f12,Co,Ch,Cad,ib)
  Bycur=By(i-1)+f11*h/2d0
  Bzcur=Bz(i-1)+f12*h/2d0
  call ff(Bycur,Bzcur,f21,f22,Co,Ch,Cad,ib)
  Bycur=By(i-1)+f21*h/2d0
  Bzcur=Bz(i-1)+f22*h/2d0
  call ff(Bycur,Bzcur,f31,f32,Co,Ch,Cad,ib)
  Bycur=By(i-1)+f31*h
  Bzcur=Bz(i-1)+f32*h
  call ff(Bycur,Bzcur,f41,f42,Co,Ch,Cad,ib)

  By(i)=By(i-1)+h/6d0*(f11+2*f21+2*f31+f41)
  Bz(i)=Bz(i-1)+h/6d0*(f12+2*f22+2*f32+f42)
  vx(i)=1/(2d0*Q) * (Kx-(Bx**2+By(i)**2+Bz(i)**2)/2d0+dsqrt((Kx-(Bx**2+By(i)**2+Bz(i)**2)/2d0)**2-4*cs**2*Q**2))
  vy(i)=(Ky+Bx*By(i))/Q
  vz(i)=(Kz+Bx*Bz(i))/Q
  rho(i)=Q/vx(i)  
   
  if(recentre==1) then
    if(abs(rho(i)-rho(i-1))>drho) then 
      drho=abs(rho(i)-rho(i-1))
      imax=i
    end if
  end if

  if(recentre==0) write(20,'(8(e15.7,2X))') dble(i-1)*h,rho(i),vx(i),vy(i),vz(i),By(i),Bz(i)
end do

if(recentre==1) then
  write(20,'(8(e15.7,2X))') 0d0,rho(1),vx(1),vy(1),vz(1),By(1),Bz(1)
  do i=1,n
    write(20,'(8(e15.7,2X))') 0.6+dble(i-imax)*h,rho(i),vx(i),vy(i),vz(i),By(i),Bz(i)
  end do
  write(20,'(8(e15.7,2X))') 1d0,rho(n),vx(n),vy(n),vz(n),By(n),Bz(n)
end if





end program





subroutine ff(By,Bz, f1,f2,Co,Ch,Cad,ib)
implicit none
double precision, intent(in) :: By,Bz
double precision,intent(in) :: Co,Ch,Cad
integer,intent(in) :: ib
double precision,intent(out) :: f1,f2
double precision,dimension(7) :: state_l,state_r
double precision :: Bx=1d0
double precision :: eta_o,eta_h,eta_ad
double precision :: Q,cs,rho
double precision :: B2,Kx,Ky,Kz,vx,vy,vz
double precision,dimension(2) :: rho0,P0,vx0,vy0,vz0,By0,Bz0
double precision :: M1,M2,R11,R12,R21,R22

state_l = (/1.7942, 0.017942, -0.9759, -0.6561, 0.0, 1.74885, 0.0/)
state_r = (/1.0, 0.01, -1.751, 0.0, 0.0, 0.6, 0.0/)

rho0(1)=state_l(1) 
P0(1  )=state_l(2) 
vx0(1) =state_l(3) 
vy0(1) =state_l(4) 
vz0(1) =state_l(5) 
By0(1) =state_l(6) 
Bz0(1) =state_l(7) 

rho0(2)=state_r(1) 
P0(2  )=state_r(2) 
vx0(2) =state_r(3) 
vy0(2) =state_r(4) 
vz0(2) =state_r(5) 
By0(2) =state_r(6) 
Bz0(2) =state_r(7) 

Q = rho0(1)*vx0(1)

cs = 0.1

Kx=Q*vx0(ib)+Q/vx0(ib)*cs*cs+0.5d0*(Bx**2+By0(ib)**2+Bz0(ib)**2)
Ky=Q*vy0(ib)-Bx*By0(ib)
Kz=Q*vz0(ib)-Bx*Bz0(ib)



B2=Bx**2+By**2+Bz**2
vx=1d0/(2d0*Q)*(Kx-B2/2d0 + dsqrt((Kx-B2/2d0)**2 - 4*cs*cs*Q*Q))
vy=(Ky+Bx*By)/Q
vz=(Kz+Bx*Bz)/Q
rho=Q/vx


eta_o=Co
eta_h=Ch*dsqrt(B2)
eta_ad=Cad*B2/rho

M1=vx*By-vx0(ib)*By0(ib)+vy0(ib)*Bx-vy*Bx
M2=vx*Bz-vx0(ib)*Bz0(ib)+vz0(ib)*Bx-vz*Bx
R11=(eta_o+eta_ad*(1d0-Bz**2/B2))
R12=(eta_h*Bx/dsqrt(B2)+eta_ad*By*Bz/B2)
R21=(-eta_h*Bx/dsqrt(B2)+eta_ad*By*Bz/B2)
R22=(eta_o+eta_ad*(1d0-By**2/B2))

f1=(M1*R22-M2*R12)/(R11*R22-R21*R12)
f2=(M1*R21-M2*R11)/(-R11*R22+R21*R12)

!f2=-1d0*f2


end subroutine ff
