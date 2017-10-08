!================================================================
!================================================================
!================================================================
!================================================================
subroutine rt_condinit(x,u,dx,nn)
  use amr_parameters
  use rt_parameters                                                    !RT
  implicit none
  integer ::nn                              ! Number of cells
  real(dp)::dx                              ! Cell size
  real(dp),dimension(1:nvector,1:nrtvar)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim  )::x ! Cell center position.
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2 &        !DAVIS
       ,scale_Np,scale_Fp                                              !DAVIS
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative RT variable vector. Conventions are here:
  ! U(i,1): N, U(i,2:ndim+1): Fx,Fy,Fy per group.
  ! U(:,:) are in user units.
  !================================================================

  ! Call built-in initial condition generator
  call rt_region_condinit(x,u,dx,nn)

  ! Add here, if you wish, some user-defined initial conditions
  ! For Davis/Krumholz experiment: Initialize to F_y=F_star
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)        !DAVIS
  call rt_units(scale_Np, scale_Fp)                                    !DAVIS
  u(1:nn,1)=1.03d4/scale_Fp/rt_c                                       !DAVIS
  u(1:nn,2)=0d0                                                        !DAVIS
  u(1:nn,3)=1.03d4/scale_Fp                                            !DAVIS
  ! ........

end subroutine rt_condinit
