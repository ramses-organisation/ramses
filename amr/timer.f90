!=======================================================================
real(kind=8) function wallclock()
  use mpi_mod
  implicit none
#ifdef WITHOUTMPI
  integer,      save :: tstart
  integer            :: tcur
  integer            :: count_rate
#else
  real(kind=8), save :: tstart
  real(kind=8)       :: tcur
#endif
  logical,      save :: first_call=.true.
  real(kind=8), save :: norm, offset=0
  !---------------------------------------------------------------------
  if (first_call) then
#ifdef WITHOUTMPI
     call system_clock(count=tstart, count_rate=count_rate)
     norm=1d0/count_rate
#else
     norm = 1d0
     tstart = MPI_Wtime()
#endif
     first_call=.false.
  end if
#ifdef WITHOUTMPI
  call system_clock(count=tcur)
#else
  tcur = MPI_Wtime()
#endif
  wallclock = (tcur-tstart)*norm + offset
  if (wallclock < 0.) then
     offset = offset + 24d0*3600d0
     wallclock = wallclock + 24d0*3600d0
  end if
end function wallclock
!=======================================================================
module timer_m
  implicit none
  integer,            parameter         :: mtimer=200                    ! max nr of timers
  real(kind=8),       dimension(mtimer) :: start, time
  integer                               :: ntimer=0, itimer
  character(len=72), dimension(mtimer)  :: labels
contains
!-----------------------------------------------------------------------
subroutine findit (label)
  implicit none
  character(len=*) label
  do itimer=1,ntimer
     if (trim(label) == trim(labels(itimer))) return
  end do
  ntimer = ntimer+1
  itimer = ntimer
  labels(itimer) = label
  time(itimer) = 0
end subroutine
end module
!=======================================================================
subroutine timer (label, cmd)
  use timer_m
  implicit none
  character(len=*)::label,cmd
  real(kind=8)::wallclock,current
!-----------------------------------------------------------------------
  current = wallclock()                                                 ! current time
  if (itimer > 0) then                                                  ! if timer is active ..
     time(itimer) = time(itimer) + current - start(itimer)              ! add to it
  end if
  call findit (label)                                                   ! locate timer slot
  if (cmd == 'start') then                                              ! start command
     start(itimer) = current                                            ! register start time
  else if (cmd == 'stop') then                                          ! stop command
     itimer = 0                                                         ! turn off timer
  end if
end subroutine
!=======================================================================
subroutine output_timer(write_file, filename)
  use amr_parameters
  use amr_commons
  use timer_m
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  real(kind=8) :: gtotal, avtime, rmstime
  real(kind=8), dimension(ncpu) :: vtime
  integer,      dimension(ncpu) :: all_ntimer
  logical,      dimension(ncpu) :: gprint_timer
  integer      :: imn, imx, mpi_err, icpu
  logical      :: print_timer
#endif
  real(kind=8) :: total
  integer      :: i
  logical      :: id_is_one, write_file
  integer      :: ilun=11
  character(LEN=80)::filename, fileloc !Optional for writing timing info
!-----------------------------------------------------------------------
  id_is_one = myid == 1
  total = 1e-9
  if (.not. write_file) ilun=6 ! 6 = std output
  if (id_is_one .and. write_file) then
     fileloc=TRIM(filename) ! Open file for timing info
     open(unit=ilun,file=fileloc,form='formatted')
  endif

  if (id_is_one .and. ncpu==1) write (ilun,'(/a,i7,a)') '     seconds         %    STEP (rank=',myid,')'
  do i = 1,ntimer
     total = total + time(i)
  end do
  if (ncpu==1) then
     do i = 1,ntimer
        if (id_is_one .and. time(i)/total > 0.001) write (ilun,'(f12.3,4x,f6.1,4x,a24)') &
          time(i), 100.*time(i)/total,labels(i)
     end do
     if (id_is_one) write (ilun,'(f12.3,4x,f6.1,4x,a)') total, 100., 'TOTAL'
  end if
#ifndef WITHOUTMPI
  if (ncpu > 1) then
     ! Check that timers are consistent across ranks
     call MPI_BARRIER(MPI_COMM_WORLD,mpi_err)
     call MPI_GATHER(ntimer,1,MPI_INTEGER,all_ntimer,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
     if (id_is_one) then
        if (maxval(all_ntimer) .ne. minval(all_ntimer)) then
           write (ilun,*)
           write (ilun,*) '--------------------------------------------------------------------'
           write (ilun,*) 'Error: Inconsistent number of timers on each rank. Min, max nr:', minval(all_ntimer), maxval(all_ntimer)
           write (ilun,*) 'Timing summary below can be misleading'
           write (ilun,*) 'Labels of timer on rank==1 :'
           write (ilun,*) '--------------------------------------------------------------------'
           do i=1,ntimer
              write(ilun,'(i3,1x,a)') i, labels(i)
           enddo
        endif
        ! Find first occurence of a rank with a different number of timers -- if it exists
        gprint_timer=.false.
        do icpu=1,ncpu
           if (all_ntimer(icpu) .ne. ntimer) then
              gprint_timer(icpu) = .true.
              exit
           endif
        enddo
        if (any(gprint_timer)) call sleep(1) ! Make sure that master rank finished, before we print from other rank.
     endif
     call MPI_SCATTER(gprint_timer,1,MPI_LOGICAL,print_timer,1,MPI_LOGICAL,0,MPI_COMM_WORLD,mpi_err)
     if (print_timer) then
        write (ilun,*)
        write (ilun,*) 'Labels of timer on rank==',myid
        write (ilun,*) '--------------------------------------------------------------------'
        do i=1,ntimer
           write(ilun,'(i3,1x,a)') i, labels(i)
        enddo
        write (ilun,*)
     endif

     call MPI_BARRIER(MPI_COMM_WORLD,mpi_err)
     call MPI_ALLREDUCE(total,gtotal,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpi_err)
     gtotal = gtotal / ncpu

     if (id_is_one) write (ilun,*) '--------------------------------------------------------------------'
     if (id_is_one) write (ilun,'(/a)') '     minimum       average       maximum' // &
                  '  standard dev        std/av       %   rmn   rmx  TIMER'
     do i = 1,ntimer
        call MPI_GATHER(real(time(i),kind=8),1,MPI_REAL8,vtime,1,MPI_REAL8,0,MPI_COMM_WORLD,mpi_err)
        if (id_is_one) then
           if (maxval(vtime)/gtotal > 0.001) then
              avtime  = sum(vtime) / ncpu ! average time used
              imn     = minloc(vtime,1)
              imx     = maxloc(vtime,1)
              rmstime = sqrt(sum((vtime - avtime)**2)/ncpu)
              write (ilun,'(5(f12.3,2x),f6.1,2x,2i4,4x,a24)') &
                 vtime(imn), avtime, vtime(imx), rmstime, rmstime/avtime, 100.*avtime/gtotal, imn, imx, labels(i)
           endif
        endif
     end do
     if (id_is_one) write (ilun,'(f12.3,4x,f6.1,4x,a)') total, 100., 'TOTAL'
  endif
#endif
  if (id_is_one) close(ilun)
end subroutine
!=======================================================================
subroutine reset_timer
   use timer_m
   use mpi_mod
   implicit none

!-----------------------------------------------------------------------
   do itimer = 1,ntimer
      time(itimer)=0
   end do
end subroutine
!=======================================================================
