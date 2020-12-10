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
double precision              :: b,bnew              ! By
double precision              :: vnx                 ! Vitesse des neutres suivant x
double precision              :: vny                 ! Vitesse des neutres suivant y
double precision,parameter    :: rho0 = 0.5          ! densite initiale 
double precision,parameter    :: rhoi0 = 1.          ! densite ionique initiale 
double precision, parameter   :: B0 = 2.             ! champ initial  
double precision,parameter    :: gammaAD = 75.       ! 
double precision,parameter    :: gammagp = 1.66667   ! 
double precision,parameter    :: cs = 0.5            ! vitesse du son 
double precision,parameter    :: vx0 = 5.            !  
double precision,parameter    :: vy0 = 0.            !  
double precision,parameter    :: By0 = sqrt(2.)            !  
double precision,parameter    :: Bx0 = sqrt(2.)            !  
integer                       :: i
double precision              :: A                 ! alfven mach number
double precision, parameter :: M  = 10.            ! mach number
double precision, parameter :: b_0 = sqrt(2.)/2.   ! sin theta
double precision, parameter :: p_0 = 0.125/0.5/5.**2  ! p_0
double precision            :: s                   ! deflection du champ
double precision            :: r_n,r_i,r           ! vs/vnx, vs/vix, 1-r_n/r_i
double precision            :: p                   ! pression

open(10,file='res.txt')

xmil = 0.5   ! a changer avec le dirac à t=0
nx =  1000
box = 1.
dx = box/dble(nx)

epsi = 4.6d-7
!write(*,*) "valeur de epsilon ?"
!read(*,*) epsi
!!!!!!!! unités rationnelles de merdeeee
!A = (4.*cos(-1.)*rho0)**0.5*5./B0
A = (rho0)**0.5*5./B0
print*, "A :", A

!read(*,*) xoffset 
xoffset= -0.03 ! 0.0 


! Etat initial
b=b_0+epsi
p=p_0

call calcul(p,p_0,b_0,b,A,s,r_n,r_i,r)

do  i=1,nx

    bnew = b+dx*gammaAD*rhoi0*A**2*r/(vx0*b)
    p = dx*gammaAD*rhoi0*r/vx0*(1./r_n+gammagp/(gammagp-1.)*p-(s+b_0)/b)*(gammagp-1.)*r_n/(1.-gammagp*r_n*p)+p
    !p = dx*gammaAD*rhoi0*r/vx0*(gammagp*p)/(1.-gammagp*p*r_n)+p
    b=bnew
    
    call calcul(p,p_0,b_0,b,A,s,r_n,r_i,r)
   
    ! autres variables
    vnx = vx0/r_n
    vny = s*vx0/b_0
    rho = rho0*vx0/vnx
    
    x = i*dx
    !write(*,*)x+xoffset,b*B0,p*rho0*vx0**2,vnx,vny,rho                ! x,By,P,vnx,vny,rho
    write(10,*)x+xoffset,b*B0,p*rho0*vx0**2,vnx,vny,rho 
end do

close(10)

end program anal

!!!!!!!!!!!! calcul de b
subroutine calcul(p,p_0,b_0,b,A,s,r_n,r_i,r)
implicit none

double precision            :: b_0                 ! sin theta
double precision            :: A                   ! nombres de Mach
double precision            :: b                   ! champ normalisé
double precision            :: s                   ! deflection du champ
double precision            :: r_n,r_i,r           ! vs/vnx, vs/vix, 1-r_n/r_i
double precision            :: p,p_0               ! pressions

s = (b-b_0)*b_0**2/A**2
r_n = 1./(1.-(p-p_0)-(b**2-b_0**2)/(2.*A**2))
r_i = r_n*( (b**2+b_0**2)/(b*r_n*(s+b_0)+b_0**2) )
r = (1.-r_n/r_i)


end subroutine calcul
