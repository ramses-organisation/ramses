subroutine update_time(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use cooling_module
  use mpi_mod
  implicit none
  integer::ilevel

  real(dp)::dt,econs,mcons,ttend,momcons
  real(dp),save::ttstart=0
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
     ! Check conservation
     !--------------------------
     if(mass_tot_0==0.0D0)then
        mass_tot_0=mass_tot
        mcons=0.0D0
     else
        mcons=(mass_tot-mass_tot_0)/mass_tot_0
     end if

    if(e_tot_0==0.0D0)then
        e_tot_0=e_tot
        econs=0.0D0
     else
        econs=(e_tot-e_tot_0)/e_tot_0
     end if

    if(mom_tot_0==0.0D0)then
       mom_tot_0=mom_tot
        momcons=0.0D0
     else
        momcons=(mom_tot-mom_tot_0)/mom_tot_0
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
           if(scheme.eq.'induction')then
#ifdef SOLVERmhd
!              write(*,778)nstep_coarse,econs,epot_tot,ekin_tot,emag_tot
#endif
           else if(cooling.or.pressure_fix)then
!              write(*,778)nstep_coarse,econs,epot_tot,ekin_tot,eint_tot
           else
!              write(*,779)nstep_coarse,mcons,econs,ekin_tot,lor_max
               write(*,779)nstep_coarse,mcons,momcons,econs,lor_max
          end if
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
  end if

779 format(' Main step=',i6,' mcons=',1pe10.3,' momcons=',1pe10.3, &
         & ' econs=',1pe10.3,'  lor_max=',1pe10.3)
888 format(' Fine step=',i6,' t=',1pe12.5,' dt=',1pe10.3, &
         & ' a=',1pe10.3,' mem=',0pF4.1,'% ',0pF4.1,'%')
999 format(' Level ',I2,' has ',I10,' grids (',3(I8,','),')')

end subroutine update_time
