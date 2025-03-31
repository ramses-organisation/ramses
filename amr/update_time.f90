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
  !---------------------------------------------------------------------------------
  ! This routine does several things.
  !
  ! It reports the progress of the simulation by printing following stats:
  !   - the AMR mesh structure
  !   - mass and energy conservation
  !   - time step number and size dt
  !   - current simulation time t and expension factor aexp
  !   - percentage of grid and particle memory in use
  !
  ! It checks whether the simulation reached the end, and calculates the total
  ! execution time.
  !
  ! It advances the global time variations: t, nstep, aexp, hexp, texp.
  ! If required, the magnetic field is adjusted to cosmic expansion.
  !
  ! The turbulence forcing is evolved.
  !--------------------------------------------------------------------------------

  real(dp)::dt,econs,mcons
#ifdef SOLVERmhd
  real(dp)::sqrt_aexp_prev
#endif
#if USE_TURB==1
  real(kind=dp) :: cur_turb_rms
#endif
  integer::i,itest

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
           if(.not.finish_run)write(*,*)'Run completed'
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
