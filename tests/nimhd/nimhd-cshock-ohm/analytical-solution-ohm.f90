program anal
implicit none

double precision              :: x                   ! position
double precision              :: xoffset             ! position decalee
double precision              :: epsi                ! pour le départ
double precision              :: xmil                ! position du dirac initial
double precision              :: box                 ! taille de la boite
double precision              :: dx                  ! pas
integer                       :: nx                  ! nombre de points (moins 1)
double precision              :: rho                 ! densité 
double precision              :: b,bnew,bold,boldnew              ! By
double precision              :: vx,rac1,rac2        ! Vitesse des neutres suivant x
double precision              :: vy                  ! Vitesse des neutres suivant y
double precision,parameter    :: rho0 = 0.4          ! densite initiale 
double precision, parameter   :: B0 = 1.             ! champ initial  
double precision,parameter    :: gammaAD = 75.       ! 
double precision,parameter    :: gammagp = 1.666667   ! 
double precision,parameter    :: eta = 0.1          ! 
double precision,parameter    :: cs = 0.5            ! vitesse du son 
double precision,parameter    :: vx0 = 3.            !  
double precision,parameter    :: vy0 = 0.            !  
double precision,parameter    :: By0 = sqrt(2.)/2.   !  
double precision,parameter    :: Bx0 = sqrt(2.)/2.   !  
integer                       :: i
double precision              :: K1,K2,K3,K4,K5,K6,etak7,etot   ! constantes
double precision, parameter :: M  = 10.            ! mach number
double precision, parameter :: b_0 = sqrt(2.)/2.   ! sin theta
double precision, parameter :: p_0 = 0.4           ! p_0
double precision            :: s                   ! deflection du champ
double precision            :: r_n,r_i,r           ! vs/vnx, vs/vix, 1-r_n/r_i
double precision            :: p                   ! pression

open(10,file='ressteady.txt')

xmil = 0.5   ! a changer avec le dirac à t=0
nx =  1000
box = 1.
dx = box/dble(nx)

epsi = 0.16d-4
!write(*,*) "valeur de epsilon ?"
!read(*,*) epsi
!!!!!!!! unités rationnelles de merdeeee

!read(*,*) xoffset 
xoffset= 0. 
x=0.

K1=rho0*vx0
K2=Bx0
K3=rho0*vx0**2+0.5*By0**2+p_0
K4=rho0*vx0*vy0-Bx0*By0
K5=(p_0/(gammagp-1.)+0.5*rho0*(vx0**2+vy0**2)+0.5*(Bx0**2+By0**2)+&
   &p_0+0.5*(Bx0**2+By0**2))*vx0-vx0*Bx0**2-vy0*Bx0*By0
K6=vx0*By0-vy0*Bx0
etot=p_0/(gammagp-1.)+0.5*rho0*vx0**2+0.5*(Bx0**2+By0**2)
print*, 'les konstantes',K1,K2,K3,K4,K5,K6,etot


! Etat initial
b=b_0+epsi
p=p_0
vx=vx0
vy=vy0
rho=rho0



do  i=1,nx

    vy   = (K2*b+K4)/K1
    
    call vnx(gammagp,b,K1,K2,K3,K4,K5,K6,vy,Bx0,etot,rac1,rac2)
    vx=rac2
    !if(vx>=1.715.or.vx<1.4) then
        !vx=rac2
    !else
        !vx=rac1
    !end if
    p = K3-0.5*b**2-K1*vx
    rho = K1/vx
    if (i==1) then
        bold=b
        b=b+dx/eta*(vx*b-vy*Bx0-K6)
    else
        b=b+dx/eta*(vx*b-vy*Bx0-K6)
    end if
    x=i*dx
    !call vnx(gammagp,b,K3,K1,K6,K5,K4,vy,Bx0,vx)
    !print*, 'vx :',vx, 'vy :', vy, 'by :', b,'rho :', rho, 'p :', p
    !write(*,*)x+xoffset,b*B0,p*rho0*vx0**2,vnx,vny,rho                ! x,By,P,vnx,vny,rho
    write(10,*)x+xoffset,rho,vx,vy,Bx0,b,p 
    !write(*,*)x+xoffset,rho,vx,vy,Bx0,b,p 
end do
    write(*,*)'x+xoffset,rho,vx,vy,Bx0,b,p '
    write(*,*)x+xoffset,rho,vx,vy,Bx0,b,p 

close(10)

end program anal


subroutine vnx(gammagp,b,K1,K2,K3,K4,K5,K6,vy,Bx0,etot,rac1,rac2)
implicit none


double precision              :: K1,K2,K3,K4,K5,K6,gammagp,b,vy,Bx0,etot
double precision              :: A,BB,C,delt,rac1,rac2

A=1.
BB=gammagp*2./((gammagp+1.)*K1)*(0.5*b**2-K3)
C=2.*(1.-gammagp)*(b*K6+0.5*K1*vy**2-K5)/(K1*(gammagp+1.))
 

delt=BB**2-4.*A*C
!print*, 'delt', delt
if  (delt <= 0.) then 
    print*, 'probleeeeemeee'
    stop
end if
rac1=(-BB-sqrt(delt))/2.
rac2=(-BB+sqrt(delt))/2.
print*,'racines',rac1,rac2

end subroutine vnx
