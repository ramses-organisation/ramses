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
subroutine update_time(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use cooling_module
#if USE_TURB==1
  use turb_commons
#endif
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  real(kind=8)::ttend
  real(kind=8),save::ttstart=0
#endif
  integer::ilevel

  real(dp)::dt,econs,mcons
#ifdef SOLVERmhd
  real(dp)::sqrt_aexp_prev
#endif
#if USE_TURB==1
  real(kind=dp) :: cur_turb_rms
#endif
  integer::i,itest
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Added for new output file with energy data.
  character(LEN=80)::filename,fileloc
  character(LEN=5)::nchar
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Local constants
  dt=dtnew(ilevel)
  itest=0

#ifndef WITHOUTMPI
  if(myid==1)then
     if(ttstart.eq.0.0)ttstart=MPI_WTIME()
  endif
#endif

  !-------------------------------------------------------------
  ! At this point, IF nstep_coarse has JUST changed, all levels
  ! are synchronised, and all new refinements have been done.
  !-------------------------------------------------------------
  if(nstep_coarse .ne. nstep_coarse_old)then

     !--------------------------
     ! Check mass conservation
     !--------------------------
     if(mass_tot_0==0.0D0)then
        mass_tot_0=mass_tot
        mcons=0.0D0
     else
        mcons=(mass_tot-mass_tot_0)/mass_tot_0
     end if

     !----------------------------
     ! Check energy conservation
     !----------------------------
     if(epot_tot_old.ne.0)then
        epot_tot_int=epot_tot_int + &
             & 0.5D0*(epot_tot_old+epot_tot)*log(aexp/aexp_old)
     end if
     epot_tot_old=epot_tot
     aexp_old=aexp
     if(einit==0.0D0)then
        einit=epot_tot+ekin_tot  ! initial total energy
        econs=0.0D0
     else
        econs=(ekin_tot+epot_tot-epot_tot_int-einit) / &
             &(-(epot_tot-epot_tot_int-einit)+ekin_tot)
     end if

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Added
     ! filename='energy.dat'
     ! if(myid==1)then
     !   call title(myid,nchar)
     !   fileloc=TRIM(filename)//TRIM(nchar)
     !   open(25+myid, file = fileloc, status = 'unknown', access = 'append')
     ! endif
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     if(mod(nstep_coarse,ncontrol)==0.or.output_done)then
        if(myid==1)then

           !-------------------------------
           ! Output AMR structure to screen
           !-------------------------------
           write(*,*)'Mesh structure'
           do i=1,nlevelmax
              if(numbtot(1,i)>0)write(*,999)i,numbtot(1:4,i)
           end do

           !----------------------------------------------
           ! Output mass and energy conservation to screen
           !----------------------------------------------
           if(cooling.or.pressure_fix)then
              write(*,778)nstep_coarse,mcons,econs,epot_tot,ekin_tot,eint_tot
           else
              write(*,777)nstep_coarse,mcons,econs,epot_tot,ekin_tot
           end if
#ifdef SOLVERmhd
           write(*,'(" emag=",ES9.2)') emag_tot
#endif
           if(pic)then
              write(*,888)nstep,t,dt,aexp,&
                   & real(100.0D0*dble(used_mem_tot)/dble(ngridmax+1)),&
                   & real(100.0D0*dble(npartmax-numbp_free_tot)/dble(npartmax+1))
           else
              write(*,888)nstep,t,dt,aexp,&
                   & real(100.0D0*dble(used_mem_tot)/dble(ngridmax+1))
           endif
           itest=1
        end if
        output_done=.false.
     end if

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Added
     ! if(myid==1)then
     !   write(25+myid,*)t,emag_tot,ekin_tot-eint_tot-emag_tot
     ! endif
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     !---------------
     ! Exit program
     !---------------
     if(t>=tout(noutput).or.aexp>=aout(noutput).or. &
          & nstep_coarse>=nstepmax)then
        if(myid==1)then
           write(*,*)'Run completed'
#ifndef WITHOUTMPI
           ttend=MPI_WTIME()
           write(*,*)'Total elapsed time:',ttend-ttstart
#endif
        endif
        call clean_end
     end if

  end if
  nstep_coarse_old=nstep_coarse

  !----------------------------
  ! Output controls to screen
  !----------------------------
  if(mod(nstep,ncontrol)==0)then
     if(myid==1.and.itest==0)then
        if(pic)then
           write(*,888)nstep,t,dt,aexp,&
                & real(100.0D0*dble(used_mem_tot)/dble(ngridmax+1)),&
                & real(100.0D0*dble(npartmax-numbp_free_tot)/dble(npartmax+1))
        else
           write(*,888)nstep,t,dt,aexp,&
                & real(100.0D0*dble(used_mem_tot)/dble(ngridmax+1))
        endif
     end if
  end if

  !------------------------
  ! Update time variables
  !------------------------
  t=t+dt
  nstep=nstep+1
  if(cosmo)then
#ifdef SOLVERmhd
     ! Keep for magnetic field expansion
     sqrt_aexp_prev = SQRT(aexp)
#endif
     ! Find neighboring times
     i=1
     do while(tau_frw(i)>t.and.i<n_frw)
        i=i+1
     end do
     ! Interpolate expansion factor
     aexp = aexp_frw(i  )*(t-tau_frw(i-1))/(tau_frw(i  )-tau_frw(i-1))+ &
          & aexp_frw(i-1)*(t-tau_frw(i  ))/(tau_frw(i-1)-tau_frw(i  ))
     hexp = hexp_frw(i  )*(t-tau_frw(i-1))/(tau_frw(i  )-tau_frw(i-1))+ &
          & hexp_frw(i-1)*(t-tau_frw(i  ))/(tau_frw(i-1)-tau_frw(i  ))
     texp =    t_frw(i  )*(t-tau_frw(i-1))/(tau_frw(i  )-tau_frw(i-1))+ &
          &    t_frw(i-1)*(t-tau_frw(i  ))/(tau_frw(i-1)-tau_frw(i  ))

#ifdef SOLVERmhd
     do i=1,ilevel
       call update_cosmomag(i,SQRT(aexp)/sqrt_aexp_prev)
     end do
#endif
  else
     aexp = 1
     hexp = 0
     texp = t
  end if

#if USE_TURB==1
  if (turb) then
     call turb_check_time
     if (myid==1) then
        call current_turb_rms(cur_turb_rms)
        write (6,*) ' Current turbulent rms: ', cur_turb_rms
     end if
  end if
#endif

777 format(' Main step=',i7,' mcons=',1pe9.2,' econs=',1pe9.2, &
         & ' epot=',1pe9.2,' ekin=',1pe9.2)
778 format(' Main step=',i7,' mcons=',1pe9.2,' econs=',1pe9.2, &
         & ' epot=',1pe9.2,' ekin=',1pe9.2,' eint=',1pe9.2)
888 format(' Fine step=',i7,' t=',1pe12.5,' dt=',1pe10.3, &
         & ' a=',1pe10.3,' mem=',0pF4.1,'% ',0pF4.1,'%')
999 format(' Level ',I2,' has ',I10,' grids (',3(I8,','),')')

end subroutine update_time

!------------------------------------------------------------------------
SUBROUTINE getProperTime(tau,tproper)
! Calculate proper time tproper corresponding to conformal time tau (both
! in code units).
!------------------------------------------------------------------------
  use amr_commons
  implicit none
  real(dp)::tau, tproper
  integer::i
  if(.not. cosmo .or. tau .eq. 0d0) then ! this might happen quite often
     tproper = tau
     return
  endif
  i = 1
  do while( tau_frw(i) > tau .and. i < n_frw )
     i = i+1
  end do
  tproper = t_frw(i  )*(tau-tau_frw(i-1))/(tau_frw(i  )-tau_frw(i-1))+ &
          & t_frw(i-1)*(tau-tau_frw(i  ))/(tau_frw(i-1)-tau_frw(i  ))
END SUBROUTINE getProperTime
!------------------------------------------------------------------------
SUBROUTINE getAgeGyr(t_birth_proper, age)
! Calculate proper time passed, in Gyrs, since proper time t_birth_proper
! (given in code units) until the current time.
!------------------------------------------------------------------------
  use amr_commons
  use pm_commons
  use constants,only: Gyr2sec
  implicit none
  real(dp):: t_birth_proper, age
  real(dp),save:: scale_t_Gyr
  logical,save::scale_init=.false.
  real(dp):: scale_nH, scale_T2, scale_l, scale_d, scale_t, scale_v
  if( .not. scale_init) then
     ! The timescale has not been initialized
     call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
     scale_t_Gyr = (scale_t/aexp**2)/Gyr2sec
     scale_init=.true.
  endif
  age = (texp - t_birth_proper) * scale_t_Gyr
END SUBROUTINE getAgeGyr
!------------------------------------------------------------------------
SUBROUTINE getAgeSec(t_birth_proper, age)
! Calculate proper time passed, in sec, since proper time t_birth_proper
! (given in code units) until the current time.
!------------------------------------------------------------------------
  use amr_commons
  use pm_commons
  implicit none
  real(dp):: t_birth_proper, age
  real(dp),save:: scale_t_sec
  logical::scale_init=.false.
  real(dp):: scale_nH, scale_T2, scale_l, scale_d, scale_t, scale_v
  if( .not. scale_init) then
     ! The timescale has not been initialized
     call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
     scale_t_sec = (scale_t/aexp**2)
     scale_init=.true.
  endif
  age = (texp - t_birth_proper) * scale_t_sec
END SUBROUTINE getAgeSec
!------------------------------------------------------------------------
