subroutine FirstAndSecondBorisKick(nn,dt,dtarr,ctm,ts,b,u,v)
  use amr_parameters
  use hydro_parameters
  implicit none
  integer ::kick ! kick number
  integer ::nn ! number of cells
  real(dp) ::dt ! timestep
  real(dp) ::ctm ! charge-to-mass ratio
  real(dp) ::ts ! stopping time
  real(dp),dimension(1:nvector,1:ndim) ::b ! magnetic field components
  real(dp),dimension(1:nvector,1:ndim) ::u ! fluid velocity
  real(dp),dimension(1:nvector,1:ndim) ::v ! grain velocity
  real(dp),dimension(1:nvector,1:ndim),save ::vo ! grain velocity "new"
  ! real(dp),dimension(1:nvector,1:nvar+3),save ::vo ! velocity output
  !real(dp),dimension(1:nvector,1:nvar+3),save ::q   ! Primitive variables
  integer ::i ! Just an index
  !if (kick==1) then
  do i=1,nn
    vo(i,1) = v(i,1) + (2*ctm*dt*(-(b(i,2)*b(i,2)*ctm*dt*v(i,1))&
              + b(i,2)*(b(i,1)*ctm*dt*v(i,2)&
              - 2*v(i,3)) + b(i,3)*(-(b(i,3)*ctm*dt*v(i,1)) + 2*v(i,2)&
              + b(i,1)*ctm*dt*v(i,3))))/(4 +&
              (b(i,1)*b(i,1) + b(i,2)*b(i,2) + b(i,3)*b(i,3))*ctm*ctm*dt*dt)
    vo(i,2) = v(i,2) + (2*ctm*dt*(-(b(i,3)*b(i,3)*ctm*dt*v(i,2)) &
              + b(i,1)*(b(i,2)*ctm*dt*v(i,1)&
              - b(i,1)*ctm*dt*v(i,2) + 2*v(i,3)) + b(i,3)*(-2*v(i,1)&
              + b(i,2)*ctm*dt*v(i,3))))/(4&
              + (b(i,1)*b(i,1) + b(i,2)*b(i,2) + b(i,3)*b(i,3))*ctm*ctm*dt*dt)
    vo(i,3) = v(i,3) + (2*ctm*dt*(2*b(i,2)*v(i,1) &
              + b(i,1)*b(i,3)*ctm*dt*v(i,1) - 2*b(i,1)*v(i,2)&
              + b(i,2)*b(i,3)*ctm*dt*v(i,2) - (b(i,1)*b(i,1)&
              + b(i,2)*b(i,2))*ctm*dt*v(i,3)))/(4 +&
              (b(i,1)*b(i,1) + b(i,2)*b(i,2) + b(i,3)*b(i,3))*ctm*ctm*dt*dt)
  end do
  v(1:nvector,1:ndim)=vo(1:nvector,1:ndim)
  !else
  do i=1,nn
    vo(i,1) = (v(i,1) - 0.5*dt*(ctm*(u(i,2)*b(i,3)-u(i,3)*b(i,2))&
              -u(i,1)/ts))/(1.+0.5*dt/ts)
    vo(i,2) = (v(i,2) - 0.5*dt*(ctm*(u(i,3)*b(i,1)-u(i,1)*b(i,3))&
              -u(i,2)/ts))/(1.+0.5*dt/ts)
    vo(i,3) = (v(i,3) - 0.5*dt*(ctm*(u(i,1)*b(i,2)-u(i,2)*b(i,1))&
              -u(i,3)/ts))/(1.+0.5*dt/ts)
  end do
  !end if
  ! do i=1,nn
  !   write(*,*)kick,nn,dt,ctm,ts
  ! end do
  v(1:nvector,1:ndim)=vo(1:nvector,1:ndim)
end subroutine FirstAndSecondBorisKick

subroutine ThirdBorisKick(nn,dtarr,ctm,ts,b,u,v)
  use amr_parameters
  use hydro_parameters
  implicit none
  integer ::kick ! kick number
  integer ::nn ! number of cells
  real(dp) ::dt ! timestep
  real(dp) ::ctm ! charge-to-mass ratio
  real(dp) ::ts ! stopping time
  real(dp),dimension(1:nvector,1:ndim) ::b ! magnetic field components
  real(dp),dimension(1:nvector,1:ndim) ::u ! fluid velocity
  real(dp),dimension(1:nvector,1:ndim) ::v ! grain velocity
  real(dp),dimension(1:nn)::dtarr
  real(dp),dimension(1:nvector,1:ndim),save ::vo ! grain velocity "new"
  ! real(dp),dimension(1:nvector,1:nvar+3),save ::vo ! velocity output
  !real(dp),dimension(1:nvector,1:nvar+3),save ::q   ! Primitive variables
  integer ::i ! Just an index
  !if (kick==1) then
    !do i=1,nn
    !  vo(i,1) = v(i,1) + (2*ctm*dtarr(i)*(-(b(i,2)*b(i,2)*ctm*dtarr(i)*v(i,1))&
    !            + b(i,2)*(b(i,1)*ctm*dtarr(i)*v(i,2)&
    !            - 2*v(i,3)) + b(i,3)*(-(b(i,3)*ctm*dtarr(i)*v(i,1)) + 2*v(i,2)&
    !            + b(i,1)*ctm*dtarr(i)*v(i,3))))/(4 +&
    !            (b(i,1)*b(i,1) + b(i,2)*b(i,2) + b(i,3)*b(i,3))*ctm*ctm*dtarr(i)*dtarr(i))
    !  vo(i,2) = v(i,2) + (2*ctm*dtarr(i)*(-(b(i,3)*b(i,3)*ctm*dtarr(i)*v(i,2)) &
    !            + b(i,1)*(b(i,2)*ctm*dtarr(i)*v(i,1)&
    !            - b(i,1)*ctm*dtarr(i)*v(i,2) + 2*v(i,3)) + b(i,3)*(-2*v(i,1)&
    !            + b(i,2)*ctm*dtarr(i)*v(i,3))))/(4&
    !            + (b(i,1)*b(i,1) + b(i,2)*b(i,2) + b(i,3)*b(i,3))*ctm*ctm*dtarr(i)*dtarr(i))
    !  vo(i,3) = v(i,3) + (2*ctm*dtarr(i)*(2*b(i,2)*v(i,1) &
    !            + b(i,1)*b(i,3)*ctm*dtarr(i)*v(i,1) - 2*b(i,1)*v(i,2)&
    !            + b(i,2)*b(i,3)*ctm*dtarr(i)*v(i,2) - (b(i,1)*b(i,1)&
    !            + b(i,2)*b(i,2))*ctm*dtarr(i)*v(i,3)))/(4 +&
    !            (b(i,1)*b(i,1) + b(i,2)*b(i,2) + b(i,3)*b(i,3))*ctm*ctm*dtarr(i)*dtarr(i))
    !end do
  !else
  do i=1,nn
    vo(i,1) = (v(i,1) - 0.5*dtarr(i)*(ctm*(u(i,2)*b(i,3)-u(i,3)*b(i,2))&
              -u(i,1)/ts))/(1.+0.5*dtarr(i)/ts)
    vo(i,2) = (v(i,2) - 0.5*dtarr(i)*(ctm*(u(i,3)*b(i,1)-u(i,1)*b(i,3))&
              -u(i,2)/ts))/(1.+0.5*dtarr(i)/ts)
    vo(i,3) = (v(i,3) - 0.5*dtarr(i)*(ctm*(u(i,1)*b(i,2)-u(i,2)*b(i,1))&
              -u(i,3)/ts))/(1.+0.5*dtarr(i)/ts)
  end do
  !end if
  ! do i=1,nn
  !   write(*,*)kick,nn,dt,ctm,ts
  ! end do
  v(1:nvector,1:ndim)=vo(1:nvector,1:ndim)
end subroutine ThirdBorisKick

! This is a subroutine to compute the boris timestep.
subroutine BorisKick(kick,nn,dt,ctm,ts,b,u,v)
  use amr_parameters
  use hydro_parameters
  implicit none
  integer ::kick ! kick number
  integer ::nn ! number of cells
  real(dp) ::dt ! timestep
  real(dp) ::ctm ! charge-to-mass ratio
  real(dp) ::ts ! stopping time
  real(dp),dimension(1:nvector,1:ndim) ::b ! magnetic field components
  real(dp),dimension(1:nvector,1:ndim) ::u ! fluid velocity
  real(dp),dimension(1:nvector,1:ndim) ::v ! grain velocity
  real(dp),dimension(1:nvector,1:ndim),save ::vo ! grain velocity "new"
  ! real(dp),dimension(1:nvector,1:nvar+3),save ::vo ! velocity output
  !real(dp),dimension(1:nvector,1:nvar+3),save ::q   ! Primitive variables
  integer ::i ! Just an index
  if (kick==1) then
    do i=1,nn
      vo(i,1) = v(i,1) + (2*ctm*dt*(-(b(i,2)*b(i,2)*ctm*dt*v(i,1))&
                + b(i,2)*(b(i,1)*ctm*dt*v(i,2)&
                - 2*v(i,3)) + b(i,3)*(-(b(i,3)*ctm*dt*v(i,1)) + 2*v(i,2)&
                + b(i,1)*ctm*dt*v(i,3))))/(4 +&
                (b(i,1)*b(i,1) + b(i,2)*b(i,2) + b(i,3)*b(i,3))*ctm*ctm*dt*dt)
      vo(i,2) = v(i,2) + (2*ctm*dt*(-(b(i,3)*b(i,3)*ctm*dt*v(i,2)) &
                + b(i,1)*(b(i,2)*ctm*dt*v(i,1)&
                - b(i,1)*ctm*dt*v(i,2) + 2*v(i,3)) + b(i,3)*(-2*v(i,1)&
                + b(i,2)*ctm*dt*v(i,3))))/(4&
                + (b(i,1)*b(i,1) + b(i,2)*b(i,2) + b(i,3)*b(i,3))*ctm*ctm*dt*dt)
      vo(i,3) = v(i,3) + (2*ctm*dt*(2*b(i,2)*v(i,1) &
                + b(i,1)*b(i,3)*ctm*dt*v(i,1) - 2*b(i,1)*v(i,2)&
                + b(i,2)*b(i,3)*ctm*dt*v(i,2) - (b(i,1)*b(i,1)&
                + b(i,2)*b(i,2))*ctm*dt*v(i,3)))/(4 +&
                (b(i,1)*b(i,1) + b(i,2)*b(i,2) + b(i,3)*b(i,3))*ctm*ctm*dt*dt)
    end do
  else
    do i=1,nn
      vo(i,1) = (v(i,1) - 0.5*dt*(ctm*(u(i,2)*b(i,3)-u(i,3)*b(i,2))&
                -u(i,1)/ts))/(1.+0.5*dt/ts)
      vo(i,2) = (v(i,2) - 0.5*dt*(ctm*(u(i,3)*b(i,1)-u(i,1)*b(i,3))&
                -u(i,2)/ts))/(1.+0.5*dt/ts)
      vo(i,3) = (v(i,3) - 0.5*dt*(ctm*(u(i,1)*b(i,2)-u(i,2)*b(i,1))&
                -u(i,3)/ts))/(1.+0.5*dt/ts)
    end do
  end if
  ! do i=1,nn
  !   write(*,*)kick,nn,dt,ctm,ts
  ! end do
  v(1:nvector,1:ndim)=vo(1:nvector,1:ndim)
end subroutine BorisKick
