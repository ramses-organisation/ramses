program main

  use amr_parameters
  use hydro_parameters

  real(dp),dimension(1:nvector,1:ndim)::x
  real(dp),dimension(1:nvector,1:ndim)::v
  real(dp),dimension(1:nvector,1:ndim)::u
  real(dp),dimension(1:nvector,1:ndim)::b
  integer::n=10

  
  ! Intialize coordinates at random
  do i=1,n
     call random_number(x(i,1))
     call random_number(x(i,2))
     call random_number(x(i,3))
  end do
  
  ! Initialize particle velocity with gas velocity
  call velocity_field(x,u,n)
  v(1:n,1:ndim)=u(1:n,1:ndim)

  ! Write initial conditions to screen
  do i=1,n
     write(*,*)x(i,1),x(i,2),x(i,3),v(i,1),v(i,2),v(i,3)
  end do

  
  call magnetic_field(x,b,n)
  

end program main

subroutine magnetic_field(x,b,n)
  use amr_parameters
  use hydro_parameters
  integer::n
  real(dp),dimension(1:nvector,1:ndim)::x
  real(dp),dimension(1:nvector,1:ndim)::b

  integer::i  
  do i=1,n
     b(i,1)=cos(x(i,1))
     b(i,2)=cos(x(i,2))
     b(i,3)=cos(x(i,3))
  end do

end subroutine magnetic_field

subroutine velocity_field(x,v,n)
  use amr_parameters
  use hydro_parameters
  integer::n
  real(dp),dimension(1:nvector,1:ndim)::x
  real(dp),dimension(1:nvector,1:ndim)::v

  integer::i  
  do i=1,n
     v(i,1)=0.0d0
     v(i,2)=0.0d0
     v(i,3)=1.0d0
  end do

end subroutine velocity_field




