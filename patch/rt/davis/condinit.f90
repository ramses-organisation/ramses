!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(x,u,dx,nn)
  use amr_parameters
  use poisson_parameters
  use hydro_parameters
  use pm_commons
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
  ! Q(i,1): d, Q(i,2:ndim+1):u,v,w and Q(i,ndim+2): P.
  ! If nvar >= ndim+3, remaining variables are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
  integer::ivar,i
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables
  real(dp)::rho_star, scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(kind=8)::hl,hu,RandNum

  ! Call built-in initial condition generator
  call region_condinit(x,q,dx,nn)

  ! Add here, if you wish, some user-defined initial conditions
  ! BEGIN KRUMHOLZ/DAVIS LEVITATION EXPERIMENT PATCH
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  rho_star=1d5 ! Max density in code units (unit_d = 1d-5 * rho_star)
  do i=1,nn
     hl = x(i,2)-dx/2d0
     hu = x(i,2)+dx/2d0
     q(i,1) = rho_star / dx * ( exp(-hl) - exp(-hu) )
     q(i,1) = max(q(i,1),1d-10*rho_star)
  end do

  ! Introduce sinusoidal fluctuations:
  q(1:nn,1) = q(1:nn,1) * (1d0 + 0.25 * sin(4d0*3.14*x(1:nn,1)/0.5/boxlen))
  ! Add randomness:
  do i=1,nn
     call ranf(localseed,RandNum)
     q(i,1) = q(i,1) * (1+0.5*(0.5+RandNum))
print*,(1+0.5*(0.5+RandNum))
  end do
  ! Set temperature to 82 K everywhere:
  q(1:nn,4) = 82./scale_T2 / 2.33 * q(1:nn,1) ! Const T
  ! END KRUMHOLZ/DAVIS LEVITATION EXPERIMENT PATCH

  if(metal) q(1:nn,ndim+3)=z_ave*0.02

  ! Convert primitive to conservative variables
  ! density -> density
  u(1:nn,1)=q(1:nn,1)
  ! velocity -> momentum
  u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
#if NDIM>1
  u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
#endif
#if NDIM>2
  u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
#endif
  ! kinetic energy
  u(1:nn,ndim+2)=0.0d0
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,2)**2
#if NDIM>1
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,3)**2
#endif
#if NDIM>2
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,4)**2
#endif
  ! thermal pressure -> total fluid energy
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+q(1:nn,ndim+2)/(gamma-1.0d0)
#if NENER>0
  ! radiative pressure -> radiative energy
  ! radiative energy -> total fluid energy
  do ivar=1,nener
     u(1:nn,ndim+2+ivar)=q(1:nn,ndim+2+ivar)/(gamma_rad(ivar)-1.0d0)
     u(1:nn,ndim+2)=u(1:nn,ndim+2)+u(1:nn,ndim+2+ivar)
  enddo
#endif
#if NVAR>NDIM+2+NENER
  ! passive scalars
  do ivar=ndim+3+nener,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do
#endif

end subroutine condinit
