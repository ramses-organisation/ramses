!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(x,u,dx,pert,nn)
  use amr_parameters
  use hydro_parameters
  use poisson_parameters
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar+3)::u ! Conservative variables
  real(dp)::pert                          ! perturbations on/off
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:4): d.u,d.v,d.w, U(i,5): E, U(i,6:8): Bleft,
  ! U(i,nvar+1:nvar+3): Bright
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:4):u,v,w, Q(i,5): P, Q(i,6:8): Bleft,
  ! Q(i,nvar+1:nvar+3): Bright
  ! If nvar > 8, remaining variables (9:nvar) are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
#if NENER>0
  integer::irad
#endif
#if NVAR>8+NENER
  integer::ivar
#endif
  integer::i,id,iu,iv,iw,ip,is                     ! Indices Euler
  integer::iAL,iBL,iCL,iAR,iBR,iCR                 ! Indices Magnetic field
  real(dp)::p1,p2,p3,rho1,rho2,rho3,x1,x2,g,T0,drho ! Variables
  real(dp)::A1,A2,A3,B1,B2,B3,C1,C2,C3             ! Magnetic fields                   
  real(dp)::gammainit1,gammainit2,gammainit3       ! Polytropic indices
  real(dp),dimension(1:nvector,1:nvar+3),save::q   ! Primitive variables

  ! Initialize 
  id=1; iu=2; iv=3; iw=4; ip=5; is=9
  iAL=6; iBL=7; iCL=8
  iAR=nvar+1; iBR=nvar+2; iCR=nvar+3
  x1=x_center(1); x2=x_center(2)
  rho1=d_region(1); rho2=d_region(2); rho3=d_region(3)
  p1=p_region(1); p2=p_region(2); p3=p_region(3)
  A1=A_region(1); A2=A_region(2); A3=A_region(3)
  B1=B_region(1); B2=B_region(2); B3=B_region(3)
  C1=C_region(1); C2=C_region(2); C3=C_region(3)
  gammainit1=gamma_region(1)
  gammainit2=gamma_region(2)
  gammainit3=gamma_region(3)
  g = abs(gravity_params(1))

  ! Call built-in initial condition generator
  call region_condinit(x,q,dx,nn)

  do i=1,nn

    ! Bottom stable zone
    if(x(i,1) .le. x1)then
      T0 = 1 - ((gammainit1-1.0d0)/gammainit1)*g*(rho1/p1)*(x(i,1))
      q(i,id)=rho1*((T0)**(1.0/(gammainit1-1.0d0)))
      q(i,ip)=p1*((T0)**(gammainit1/(gammainit1-1.0d0)))

      q(i,iAL)=A1
      q(i,iBL)=B1
      q(i,iCL)=C1
      q(i,iAR)=A1
      q(i,iBR)=B1
      q(i,iCR)=C1
    
    ! Convective zone
    else if ((x(i,1) .gt. x1) .and. (x(i,1) .lt. x2)) then
      T0 = 1 - ((gammainit2-1.0d0)/gammainit2)*g*(rho2/p2)*((x(i,1)-x1))
      drho = 0.0d0
      if ((pert_r .gt. 0.) .and. (pert_dx .gt. 0.)) then
        if (x(i,1) .lt. x1+pert_dx) then
          !! produce density perturbation in small layer of convection zone!
          call random_number(drho)
          drho = pert_r*pert*2.0*(drho-0.5)*10.0**(-2.0)
        end if 
      end if
      q(i,id)=rho2*((T0)**(1.0/(gammainit2-1.0d0)))*(1.0-drho)
      q(i,ip)=p2*((T0)**(gammainit2/(gammainit2-1.0d0)))

      q(i,iAL)=A2
      q(i,iBL)=B2
      q(i,iCL)=C2
      q(i,iAR)=A2
      q(i,iBR)=B2
      q(i,iCR)=C2

    ! Top stable zone
    else 
      T0 = 1 - ((gammainit3-1.0d0)/gammainit3)*g*(rho3/p3)*((x(i,1)-x2))
      q(i,id)=rho3*((T0)**(1.0/(gammainit3-1.0d0)))
      q(i,ip)=p3*((T0)**(gammainit3/(gammainit3-1.0d0)))

      q(i,iAL)=A3
      q(i,iBL)=B3
      q(i,iCL)=C3
      q(i,iAR)=A3
      q(i,iBR)=B3
      q(i,iCR)=C3
    endif
    
    !! u,v,w
    q(i,iu)=0.0d0
    q(i,iv)=0.0d0
    q(i,iw)=0.0d0

#if NVAR>8
    ! Set entropy as a passive scalar
    q(i,is)=q(i,ip)/q(i,id)**gamma
#endif
  end do

  ! Convert primitive to conservative variables
  ! density -> density
  u(1:nn,1)=q(1:nn,1)
  ! velocity -> momentum
  u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
  u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
  u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
  ! kinetic energy
  u(1:nn,5)=0.0d0
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,2)**2
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,3)**2
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,4)**2
  ! pressure -> total fluid energy
  u(1:nn,5)=u(1:nn,5)+q(1:nn,5)/(gamma-1.0d0)
  ! magnetic energy -> total fluid energy
  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,6)+q(1:nn,nvar+1))**2
  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,7)+q(1:nn,nvar+2))**2
  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,8)+q(1:nn,nvar+3))**2
  u(1:nn,6:8)=q(1:nn,6:8)
  u(1:nn,nvar+1:nvar+3)=q(1:nn,nvar+1:nvar+3)
#if NENER>0
  ! radiative pressure -> radiative energy
  ! radiative energy -> total fluid energy
  do irad=1,nener
     u(1:nn,8+irad)=q(1:nn,8+irad)/(gamma_rad(irad)-1.0d0)
     u(1:nn,5)=u(1:nn,5)+u(1:nn,8+irad)
  enddo
#endif
#if NVAR>8+NENER
  ! passive scalars
  do ivar=9+nener,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do
#endif

end subroutine condinit
!================================================================
!================================================================
!================================================================
!================================================================
subroutine eneana(x,e,dx,t,ncell)
  use amr_parameters
  use hydro_parameters

  implicit none
  integer ::ncell                         ! Size of input arrays
  real(dp)::dx                            ! Cell size
  real(dp)::t                             ! Current time
  real(dp),dimension(1:nvector)::e        ! Energy 
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! Heat bottom of convection zone with 
  ! e = 2*10**10 erg/g/s * rho0       for HeFlash
  ! e = 2.489*10.0**12 erg/g/s * rho0 for Model S
  ! where rho0 is the density at the bottom of convection zone
  !
  ! We heat the first layer of 
  !   2000 km for Model S
  !   500 km  for HeFlash
  ! of the convection zone
  !================================================================
  integer :: i, idim
  real(dp):: rho0, e0, x1, x2, dxq, pi

  pi = 4.d0*ATAN(1.d0)

  x1 = x_center(1)
  x2 = x_center(2)
  rho0 = d_region(2)

  e0 = heating_r
  dxq = heating_dx
  ! HeFlash
  !e0 = 2.0e-6 ! erg/g/s (Normalized 10**16)
  !dxq = 0.5d0
  ! Model S
  !e0 = 2.489e-4 ! erg/g/s (Normalized 10**16)
  ! dxq = 2.0d0
  
  ! Initialize
  do i=1,ncell
    e(i) = 0.0d0
  end do 

  ! Heating loop
  do i=1,ncell
    if ((x(i,1) .gt. x1) .and. (x(i,1) .lt. x1+dxq)) then
      ! heating
      !e(i) = e0*(1.0 + cos(2.0*pi*(x(i,1)-x1-dxq/2.0)/dxq))/dxq
      e(i) = e0*rho0 ! erg/s/cm^3
      !e(i) = 0.d0
    else if ((x(i,1) .gt. x2-dxq) .and. (x(i,1) .lt. x2)) then
      ! cooling
      !e(i) = e0*(-1.0 - cos(2.0*pi*(x(i,1)-x2+dxq/2.0)/dxq))/dxq
      e(i) = (-e0)*rho0 ! erg/s/cm^3
      !e(i) = 0.d0
    end if
  end do 

end subroutine eneana
!================================================================
!================================================================
!================================================================
!================================================================
subroutine spongelayers(x,u,req,peq,t,ncell)

  use amr_parameters
  use hydro_parameters

  implicit none
  integer ::ncell                             ! Size of input arrays
  real(dp)::t                                 ! Current time
  real(dp),dimension(1:nvector,1:ndim)::x     ! Cell center position.
  real(dp),dimension(1:nvector,1:nvar+3)::u   ! Conservative variables 
  real(dp),dimension(1:nvector)::req,peq      ! Equilibrium profiles
  !================================================================
  integer :: i, irad, ie
  real(dp):: x1, x2, crho, cp, cv
  
  ! c: factor of "damping"
  crho=0.1d0
  cv=0.1d0
  cp=0.1d0
  ie = 5
  x1 = x_center(1); x2 = x_center(2)
  
!   ! Sponge loop 
!   do i=1,ncell
!     if ((x(i,1) .lt. x1-0.5) .or. (x(i,1) .gt. x2+0.5)) then
       
!       ! Damp the internal energy
!       u(i,ie) = u(i,ie)*(1.0d0-cp) + cp*peq(i)/(gamma-1.0d0)
!       ! damp density
!       u(i,1) = u(i,1)*(1.0d0 - crho) + crho*req(i)
!       ! damp momentum components
!       u(i,2) = u(i,2)*(1.0d0-cv)
!       u(i,3) = u(i,3)*(1.0d0-cv)
!       u(i,4) = u(i,4)*(1.0d0-cv)

!     ! else
! !       if (t .lt. 500.0) then
! !         u(i,2) = u(i,2)*(t/500.0)
! ! #if NDIM>1
! !         u(i,3) = u(i,3)*(t/500.0)
! ! #endif
! ! #if NDIM>2
! !         u(i,4) = u(i,4)*(t/500.0)
! ! #endif
! !       end if 
!     end if
!   end do 

end subroutine spongelayers
!================================================================
!================================================================
!================================================================
!================================================================
subroutine velana(x,v,dx,t,ncell)
  use amr_parameters
  use hydro_parameters
  implicit none
  integer ::ncell                         ! Size of input arrays
  real(dp)::dx                            ! Cell size
  real(dp)::t                             ! Current time
  real(dp),dimension(1:nvector,1:3)::v    ! Velocity field
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine computes the user defined velocity fields.
  ! x(i,1:ndim) are cell center position in [0,boxlen] (user units).
  ! v(i,1:3) is the imposed 3-velocity in user units.
  !================================================================
  integer::i
  real(dp)::xx,yy=0.,zz=0.,vx,vy,vz,aa,twopi
!!$  real(dp)::rr,tt,omega
  
  ! Add here, if you wish, some user-defined initial conditions
  aa=1.0+0.*t
  twopi=2d0*ACOS(-1d0)
  do i=1,ncell

     xx=x(i,1)+0.*dx
#if NDIM > 1
     yy=x(i,2)
#endif
#if NDIM > 2
     zz=x(i,3)
#endif
     ! ABC
     vx=aa*(cos(twopi*yy)+sin(twopi*zz))
     vy=aa*(sin(twopi*xx)+cos(twopi*zz))
     vz=aa*(cos(twopi*xx)+sin(twopi*yy))

!!$     ! 1D advection test
!!$     vx=1.0_dp
!!$     vy=0.0_dp
!!$     vz=0.0_dp

!!$     ! Ponomarenko
!!$     xx=xx-boxlen/2.0
!!$     yy=yy-boxlen/2.0
!!$     rr=sqrt(xx**2+yy**2)
!!$     if(yy>0)then
!!$        tt=acos(xx/rr)
!!$     else
!!$        tt=-acos(xx/rr)+twopi
!!$     endif
!!$     if(rr<1.0)then
!!$        omega=0.609711
!!$        vz=0.792624
!!$     else
!!$        omega=0.0
!!$        vz=0.0
!!$     endif
!!$     vx=-sin(tt)*rr*omega
!!$     vy=+cos(tt)*rr*omega

     v(i,1)=vx
#if NDIM > 1
     v(i,2)=vy
#endif
#if NDIM > 2
     v(i,3)=vz
#endif
  end do


end subroutine velana
