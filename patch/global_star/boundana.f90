!############################################################
!############################################################
!############################################################
!############################################################
subroutine boundana(x,u,dx,ibound,nn)
  use amr_parameters
  use hydro_parameters
  use poisson_parameters
  implicit none
  integer ::ibound                        ! Index of boundary region
  integer ::nn                            ! Number of active cells
  real(dp)::dx                            ! Cell size
#ifdef SOLVERmhd
  real(dp),dimension(1:nvector,1:nvar+3)::u ! Conservative variables
#else
  real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
#endif
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates boundary conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:4): d.u,d.v,d.w, U(i,5): E,
  ! U(i,6:8): Bleft, U(i,nvar+1:nvar+3): Bright
  ! U is in user units.
  ! ibound is the index of the boundary region defined in the namelist.
  !================================================================

#if NENER>0
  integer::irad
#endif

  integer::i,ivar,id,iu,iv,iw,ip,is                     ! Indices Euler
  integer::iAL,iBL,iCL,iAR,iBR,iCR                 ! Indices Magnetic field
  real(dp)::p1,p2,p3,rho1,rho2,rho3,x1,x2,g,T0,drho! Variables
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

#ifdef SOLVERmhd
  do ivar=1,nvar+3
#else
  do ivar=1,nvar
#endif
     do i=1,nn
        u(i,ivar)=boundary_var(ibound,ivar)
     end do
  end do

  ! User defined boundary conditions

  do i=1,nn
    !! rho, P
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
      q(i,id)=rho2*((T0)**(1.0/(gammainit2-1.0d0)))
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

end subroutine boundana
