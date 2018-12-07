!================================================================
!================================================================
!================================================================
!================================================================

subroutine condinit(x,u,dx,nn)
  use amr_parameters
  use hydro_parameters
  implicit none

  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:ndim+1): d.u,d.v,d.w and U(i,ndim+2): E.
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2): vx, Q(i,3):vy , Q(i,4): vz ,Q(i,5): p
  ! If nvar >= 5, remaining variables are treated as passive
  ! scalars in the hydro solver.
   ! U(:,:) and Q(:,:) are in user units.
  !================================================================
  integer::ivar
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables
  integer ::i
  real(dp)::h,lor,tau
  real(dp)::rl,rr,ul,vl,ur,vr,pl,pr

  ! Call built-in initial condition generator
  call region_condinit(x,q,dx,nn)






  do i=1,nn
     ! Convert primitive to conservative variables
     ! specific enthalpy

     h=1.0d0+gamma/(gamma-1.0d0)*q(i,5)/q(i,1)
     if (eos.eq. 'TM') then
        tau=q(i,5)/q(i,1)
        h=5d0/2d0*tau+3d0/2d0*sqrt(tau**2+4d0/9d0)
     endif
     ! Lorentz factor
     lor=(1.-(q(i,2)**2+q(i,3)**2+q(i,4)**2))**(-1d0/2d0)

     ! proper density  -> density in lab frame
     u(i,1)= q(i,1)*lor
     ! proper 3-velocity -> momentum in lab frame
     u(i,2)=lor**2*q(i,1)*h*q(i,2)
     u(i,3)=lor**2*q(i,1)*h*q(i,3)
     u(i,4)=lor**2*q(i,1)*h*q(i,4)
     ! isotropic gas pressure -> total fluid energy
     u(i,5)=lor**2*q(i,1)*h-q(i,5)
    !passive scalars
     do ivar=6,nvar
        u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)*lor
     end do

  end do


end subroutine condinit
