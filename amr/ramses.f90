program ramses
  implicit none

  ! Read run parameters
  call read_params

#ifndef CRAY
  ! Set signal handler
  call set_signal_handler
#endif

  ! Start time integration
  call adaptive_loop

end program ramses

#ifndef CRAY
! sets the hook to catch signal 10, doesn't work with CRAY
subroutine set_signal_handler
  implicit none
  external output_signal
  integer::istatus=-1
#ifdef NOSYSTEM
  integer::jsigact,jhandle
#endif

#ifndef NOSYSTEM
  call SIGNAL(10,output_signal,istatus)
#else
  call PXFSTRUCTCREATE("sigaction",jsigact,istatus)
  call PXFGETSUBHANDLE(output_signal,jhandle,istatus)
  call PXFINTSET(jsigact,"sa_handler",jhandle,istatus)
  call PXFSIGACTION(10,jsigact,0,istatus)
#endif
end subroutine set_signal_handler

! signal handler subroutine
subroutine output_signal
  use amr_commons
  implicit none

  if (myid==1) write (*,*) 'SIGNAL 10: Output will be written to disk during next main step.'

  ! output will be written to disk at next main step
  output_now = .true.

end subroutine output_signal
#endif
