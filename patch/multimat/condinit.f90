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
  ! Q(i,1): d, Q(i,2:ndim+1):u,v,w and Q(i,ndim+2): P.
  ! If nvar >= ndim+3, remaining variables are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
  integer::ivar,imat
  real(dp),dimension(1:nvector,1:npri),save::q   ! Primitive variables
  real(dp),dimension(1:nvector,1:nmat),save::f,g,kappa_mat ! Volume fraction and densities
  real(dp),dimension(1:nvector),save::ekin,dtot,eint,cs,kappa_hat

  ! Call built-in initial condition generator
  call region_condinit(x,q,f,g,dx,nn)

  ! Add here, if you wish, some user-defined initial conditions
  ! ........

  ! Convert primitive to conservative variables

  ! call inverse eos routine (f,g,d,p) -> (e,c)
  call eosinv(f,g,q,eint,cs,kappa_mat,kappa_hat,nn)

  ! density -> density
  u(1:nn,1)=q(1:nn,1)
  ! velocity -> momentum
  u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
  ekin(1:nn)=0.5*q(1:nn,2)**2
#if NDIM>1
  u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
  ekin(1:nn)=ekin(1:nn)+0.5*q(1:nn,3)**2
#endif
#if NDIM>2
  u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
  ekin(1:nn)=ekin(1:nn)+0.5*q(1:nn,4)**2
#endif
  ! total energy (E = e + 0.5 rho u**2)
  u(1:nn,npri)=eint(1:nn)+q(1:nn,1)*ekin(1:nn)
  ! volume fraction -> volume fraction
  do imat=1,nmat
     u(1:nn,imat+npri)=f(1:nn,imat)
  end do
  ! fluid densities -> fluid densities
  do imat=1,nmat
     u(1:nn,imat+npri+nmat)=g(1:nn,imat)
  end do

end subroutine condinit
