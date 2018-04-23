subroutine adaptive_loop
  use amr_commons
  use hydro_commons
  use pm_commons
  use poisson_commons
  use cooling_module
  use mpi_mod
  implicit none
  integer::ilevel,idim,ivar,info,i
  real(dp)::tt1,tt2


#ifndef WITHOUTMPI
  tt1=MPI_WTIME()
#endif

  call init_amr                      ! Initialize AMR variables
  call init_time                     ! Initialize time variables
  if(hydro)call init_hydro           ! Initialize hydro variables
  if(poisson)call init_poisson       ! Initialize poisson variables
  if(nrestart==0)call init_refine    ! Build initial AMR grid
  if(cooling)call set_table(dble(aexp))  ! Initialize cooling look up table
  if(pic)call init_part              ! Initialize particle variables
  if(pic)call init_tree              ! Initialize particle tree
  if(nrestart==0)call init_refine_2  ! Build initial AMR grid again

#ifndef WITHOUTMPI
  tt2=MPI_WTIME()
  if(myid==1)write(*,*)'Time elapsed since startup:',tt2-tt1
#endif

  i=0
  if(myid==1)then
     write(*,*)'Initial mesh structure'
     do ilevel=1,nlevelmax
        if(numbtot(1,ilevel)>0)write(*,999)ilevel,numbtot(1:4,ilevel)
     end do
  end if

  nstep_coarse_old=nstep_coarse

  if(myid==1)write(*,*)'Starting time integration'

  do ! Main time loop
     i=i+1
     if(verbose)write(*,*)'Entering amr_step_coarse'

     mass_tot=0.0D0  ! Reset total mass
     e_tot=0.0D0  ! Reset total energy
     mom_tot=0.0D0  ! Reset total momentum
#ifdef SOLVERmhd
     emag_tot=0.0D0  ! Reset total magnetic energy
#endif

     ! Make new refinements
     if(levelmin.lt.nlevelmax)then
        call refine_coarse
        do ilevel=1,levelmin
           call build_comm(ilevel)
           call make_virtual_fine_int(cpu_map(1),ilevel)
           if(hydro)then
#ifdef SOLVERmhd
              do ivar=1,nvar+3
#else
              do ivar=1,nvar
#endif
                 call make_virtual_fine_dp(uold(1,ivar),ilevel)
              end do
              if(simple_boundary)call make_boundary_hydro(ilevel)
              if(poisson)then
                 do idim=1,ndim
                    call make_virtual_fine_dp(f(1,idim),ilevel)
                 end do
              end if
           end if
           if(ilevel<levelmin)call refine_fine(ilevel)
        end do
     endif

     ! Call base level
     call amr_step(levelmin,1)

     if(levelmin.lt.nlevelmax)then
        ! Hydro book-keeping
        if(hydro)then
           do ilevel=levelmin-1,1,-1
              call upload_fine(ilevel)
#ifdef SOLVERmhd
              do ivar=1,nvar+3
#else
              do ivar=1,nvar
#endif
                 call make_virtual_fine_dp(uold(1,ivar),ilevel)
              end do
              if(simple_boundary)call make_boundary_hydro(ilevel)
              if(poisson)then
                 do idim=1,ndim
                    call make_virtual_fine_dp(f(1,idim),ilevel)
                 end do
              end if
           end do
        end if

        ! Build refinement map
        do ilevel=levelmin-1,1,-1
           call flag_fine(ilevel,2)
        end do
        call flag_coarse
     endif

     ! New coarse time-step
     nstep_coarse=nstep_coarse+1


  end do
close (26)
999 format(' Level ',I2,' has ',I10,' grids (',3(I8,','),')')

end subroutine adaptive_loop
