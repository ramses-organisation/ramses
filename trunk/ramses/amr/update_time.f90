!=======================================================================
real(kind=8) function wallclock()
  implicit none
#ifdef WITHOUTMPI
  integer,      save :: tstart
  integer            :: tcur
  integer            :: count_rate
#else
  integer            :: info
  real(kind=8), save :: tstart
  real(kind=8)       :: tcur
#endif
  logical,      save :: first_call=.true.
  real(kind=8), save :: norm, offset=0.
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
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
  time(itimer) = 0.
end subroutine
end module
!=======================================================================
subroutine timer (label, cmd)
  use timer_m
  implicit none
  character(len=*) label, cmd
  real(kind=8) wallclock, current
  integer ierr
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
subroutine finalize_timer
  use amr_parameters
  use amr_commons
  use timer_m
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  real(kind=8) :: total, gtotal, avtime, rmstime
  real(kind=8), dimension(ncpu) :: vtime
  integer,      dimension(ncpu) :: all_ntimer
  logical,      dimension(ncpu) :: gprint_timer
  integer      :: imn, imx, mpi_err, icpu
  logical      :: o, print_timer
!-----------------------------------------------------------------------
  o = myid == 1
  total = 1e-9
  if (o .and. ncpu==1) write (*,'(/a,i7,a)') '     seconds         %    STEP (rank=',myid,')'
  do itimer = 1,ntimer
     total = total + time(itimer)
  end do
  if (ncpu==1) then
     do itimer = 1,ntimer
        if (o .and. time(itimer)/total > 0.001) write (*,'(f12.3,4x,f6.1,4x,a24)') &
          time(itimer), 100.*time(itimer)/total,labels(itimer)
     end do
     if (o) write (*,'(f12.3,4x,f6.1,4x,a)') total, 100., 'TOTAL'
  end if
#ifndef WITHOUTMPI
  if (ncpu > 1) then
     ! Check that timers are consistent across ranks
     call MPI_BARRIER(MPI_COMM_WORLD,mpi_err)
     call MPI_GATHER(ntimer,1,MPI_INTEGER,all_ntimer,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err)
     if (o) then
        if (maxval(all_ntimer) .ne. minval(all_ntimer)) then
           write (*,*)
           write (*,*) '--------------------------------------------------------------------'
           write (*,*) 'Error: Inconsistent number of timers on each rank. Min, max nr:', minval(all_ntimer), maxval(all_ntimer)
           write (*,*) 'Timing summary below can be misleading'
           write (*,*) 'Labels of timer on rank==1 :'
           write (*,*) '--------------------------------------------------------------------'
           do itimer=1,ntimer
              write(*,'(i3,1x,a)') itimer, labels(itimer)
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
        write (*,*)
        write (*,*) 'Labels of timer on rank==',myid
        write (*,*) '--------------------------------------------------------------------'
        do itimer=1,ntimer
           write(*,'(i3,1x,a)') itimer, labels(itimer)
        enddo
        write (*,*)
     endif

     call MPI_BARRIER(MPI_COMM_WORLD,mpi_err)
     call MPI_ALLREDUCE(total,gtotal,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,mpi_err)
     gtotal = gtotal / ncpu

     if (o) write (*,*) '--------------------------------------------------------------------'
     if (o) write (*,'(/a)') '     minimum       average       maximum' // &
                  '  standard dev        std/av       %   rmn   rmx  TIMER'
     do itimer = 1,ntimer
        call MPI_GATHER(real(time(itimer),kind=8),1,MPI_REAL8,vtime,1,MPI_REAL8,0,MPI_COMM_WORLD,mpi_err)
        if (o) then
           if (maxval(vtime)/gtotal > 0.001) then
              avtime  = sum(vtime) / ncpu ! average time used
              imn     = minloc(vtime,1)
              imx     = maxloc(vtime,1)
              rmstime = sqrt(sum((vtime - avtime)**2)/ncpu)
              write (*,'(5(f12.3,2x),f6.1,2x,2i4,4x,a24)') &
                 vtime(imn), avtime, vtime(imx), rmstime, rmstime/avtime, 100.*avtime/gtotal, imn, imx, labels(itimer)
           endif
        endif
     end do
     if (o) write (*,'(f12.3,4x,f6.1,4x,a)') total, 100., 'TOTAL'
  endif
#endif
end subroutine
!=======================================================================
subroutine reset_timer
   use timer_m
   implicit none
#ifndef WITHOUTMPI
   include 'mpif.h'
#endif
!-----------------------------------------------------------------------
   do itimer = 1,ntimer
      time(itimer)=0.0
   end do
end subroutine
!=======================================================================
subroutine update_time(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use cooling_module
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif  
  integer::ilevel

  real(dp)::dt,econs,mcons
  real(kind=8)::ttend
  real(kind=8),save::ttstart=0
  integer::i,itest,info

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
     if(const==0.0D0)then
        const=epot_tot+ekin_tot  ! initial total energy
        econs=0.0D0
     else
        econs=(ekin_tot+epot_tot-epot_tot_int-const) / &
             &(-(epot_tot-epot_tot_int-const)+ekin_tot)
     end if

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
        call clean_stop
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
  else
     aexp = 1.0
     hexp = 0.0
     texp = t
  end if

777 format(' Main step=',i7,' mcons=',1pe9.2,' econs=',1pe9.2, &
         & ' epot=',1pe9.2,' ekin=',1pe9.2)
778 format(' Main step=',i7,' mcons=',1pe9.2,' econs=',1pe9.2, &
         & ' epot=',1pe9.2,' ekin=',1pe9.2,' eint=',1pe9.2)
888 format(' Fine step=',i7,' t=',1pe12.5,' dt=',1pe10.3, &
         & ' a=',1pe10.3,' mem=',0pF4.1,'% ',0pF4.1,'%')
999 format(' Level ',I2,' has ',I10,' grids (',3(I8,','),')')
 
end subroutine update_time
  
subroutine clean_stop
  use amr_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::info

  call finalize_timer

#ifndef WITHOUTMPI
  call MPI_FINALIZE(info)
#endif
  stop
end subroutine clean_stop

subroutine writemem(usedmem)
  real::usedmem
  integer::getpagesize

#ifdef NOSYSTEM
!  call PXFSYSCONF(_SC_PAGESIZE,ipagesize,ierror)
  ipagesize=4096
#else
!  ipagesize = getpagesize()
  ipagesize=4096
#endif
  usedmem=dble(usedmem)*dble(ipagesize)

  if(usedmem>1024.**3.)then
     write(*,999)usedmem/1024.**3.
  else if (usedmem>1024.**2.) then
     write(*,998)usedmem/1024.**2
  else if (usedmem>1024.) then
     write(*,997)usedmem/1024.
  endif

997 format(' Used memory:',F6.1,' kb')
998 format(' Used memory:',F6.1,' Mb')
999 format(' Used memory:',F6.1,' Gb')

end subroutine writemem

subroutine getmem(outmem)
  use amr_commons,only:myid,IOGROUPSIZE,ncpu
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'  
#endif
  real(kind=4)::outmem
  character(len=300) :: dir, dir2,  cmd, file
  integer::read_status
  integer,parameter::tag=1134
  integer::dummy_io,info2
  integer::nmem,ind,j
  logical::file_exists
  
  file='/proc/self/stat'
#ifndef WITHOUTMPI
  if(IOGROUPSIZE>0) then
     if (mod(myid-1,IOGROUPSIZE)/=0) then
        call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
             & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
     end if
  endif
#endif

  inquire(file=file, exist=file_exists)
  if (file_exists) then
     open(unit=1,file=file,form='formatted')
     read(1,'(A300)',IOSTAT=read_status)dir
     close(1)
  else
     read_status=-1000
  endif

  ! Send the token
#ifndef WITHOUTMPI
  if(IOGROUPSIZE>0) then
     if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
        dummy_io=1
        call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
             & MPI_COMM_WORLD,info2)
     end if
  endif
#endif


  if (read_status < 0)then
     outmem=0.
     if (myid==1 .and. read_status .ne. -1000)write(*,*)'Problem in checking free memory'
  else
     ind=300
     j=0
     do while (j<23)
        ind=index(dir,' ')
        dir2=dir(ind+1:300)
        j=j+1
        dir=dir2
     end do
     ind=index(dir,' ')
     dir2=dir(1:ind)
     read(dir2,'(I12)')nmem
     outmem=real(nmem,kind=4)
  end if

end subroutine getmem

subroutine cmpmem(outmem)
  use amr_commons
  use hydro_commons
  implicit none

  real::outmem,outmem_int,outmem_dp,outmem_qdp
  outmem_int=0.0
  outmem_dp=0.0
  outmem_qdp=0.0

  outmem_dp =outmem_dp +ngridmax*ndim      ! xg
  outmem_int=outmem_int+ngridmax*twondim   ! nbor
  outmem_int=outmem_int+ngridmax           ! father
  outmem_int=outmem_int+ngridmax           ! next
  outmem_int=outmem_int+ngridmax           ! prev
  outmem_int=outmem_int+ngridmax*twotondim ! son 
  outmem_int=outmem_int+ngridmax*twotondim ! flag1
  outmem_int=outmem_int+ngridmax*twotondim ! flag2
  outmem_int=outmem_int+ngridmax*twotondim ! cpu_map1
  outmem_int=outmem_int+ngridmax*twotondim ! cpu_map2
  outmem_qdp=outmem_qdp+ngridmax*twotondim ! hilbert_key

  ! Add communicator variable here

  if(hydro)then
     
  outmem_dp =outmem_dp +ngridmax*twotondim*nvar ! uold
  outmem_dp =outmem_dp +ngridmax*twotondim*nvar ! unew

  if(pressure_fix)then

  outmem_dp =outmem_dp +ngridmax*twotondim ! uold
  outmem_dp =outmem_dp +ngridmax*twotondim ! uold

  endif

  endif

  write(*,*)'Estimated memory=',(outmem_dp*8.+outmem_int*4.+outmem_qdp*8.)/1024./1024.


end subroutine cmpmem
!------------------------------------------------------------------------
SUBROUTINE getProperTime(tau,tproper)
! Calculate proper time tproper corresponding to conformal time tau (both
! in code units).
!------------------------------------------------------------------------
  use amr_commons
  implicit none
  real(dp)::tau, tproper
  integer::i
  if(.not. cosmo .or. tau .eq. 0.d0) then ! this might happen quite often
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
  implicit none
  real(dp):: t_birth_proper, age
  real(dp), parameter:: yr = 3.15569d+07
  real(dp),save:: scale_t_Gyr
  logical,save::scale_init=.false.
  real(dp):: scale_nH, scale_T2, scale_l, scale_d, scale_t, scale_v
  if( .not. scale_init) then 
     ! The timescale has not been initialized
     call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
     scale_t_Gyr = (scale_t/aexp**2)/yr/1.e9
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






