subroutine adaptive_loop
  use amr_commons
  use hydro_commons
  use pm_commons
  use poisson_commons
  use cooling_module
  use acc_commons
#ifdef RT
  use rt_hydro_commons
#endif
  implicit none
#ifndef WITHOUTMPI
#ifndef INCLUDEOK
#define INCLUDEOK
  include 'mpif.h'
#endif
#endif
  integer::ilevel,idim,ivar,info
  real(kind=8)::tt1,tt2
  real(kind=8)::tt1_amr_step,tt2_amr_step
  real(kind=4)::real_mem,real_mem_tot
#ifdef MPROFILING
#ifndef INCLUDEOK
#define INCLUDEOK
  include 'mpif.h'
#endif
#endif
#undef INCLUDEOK

#ifndef WITHOUTMPI
  tt1=MPI_WTIME(info)
#endif

  call init_amr                      ! Initialize AMR variables
  call init_time                     ! Initialize time variables
  if(hydro)call init_hydro           ! Initialize hydro variables
#ifdef RT
  if(rt.or.neq_chem) &
       call rt_init_hydro            ! Initialize radiation variables
#endif
  if(poisson)call init_poisson       ! Initialize poisson variables
#ifdef ATON
  if(aton)call init_radiation        ! Initialize radiation variables
#endif
  if(nrestart==0)call init_refine    ! Build initial AMR grid

#ifdef grackle
  if(cosmo)then
     ! Compute cooling table at current aexp
  endif
#else  
  if(cooling.and..not.neq_chem) &
       call set_table(dble(aexp))    ! Initialize cooling look up table
#endif
  if(pic)call init_part              ! Initialize particle variables
  if(pic)call init_tree              ! Initialize particle tree
  if(nrestart==0)call init_refine_2  ! Build initial AMR grid again

#if defined(_OPENACC) || defined(USE_ACC_VERSION)
  call init_acc                      ! Initialize OPENACC variables
#endif

#ifdef MPROFILING
  call init_mprof
#endif

  ! ACC communication buffers (see make_virtual)
  allocate(nemisspe_acc(1:ncpu,1:nlevelmax),nrecpe_acc(1:ncpu,1:nlevelmax))
  allocate(sum_nemisspe(1:ncpu,1:nlevelmax),sum_nrecpe(1:ncpu,1:nlevelmax))

#ifndef AMR_OPT_CPY
  !$acc data create(unew,uold,divu_acc,enew_acc,ucount0,gcount0,igrid_acc,divu,enew) &
  !$acc create(phi, phi_old, rho, f, flag2, father, xg, son, nbor, cpu_map) &
  !$acc create(pp_shfv,ind_cell_shfv,ind_grid_shfv) & !synchro_hydro_fine subroutine
  !$acc create(nemisspe_acc,nrecpe_acc,sum_nemisspe, sum_nrecpe) &
  !$acc copyin(state_diagram,ggg_acc,hhh_acc,oneontwotondim)
#else
  !$acc data create(unew,uold,f,divu_acc,enew_acc,ucount0,gcount0,igrid_acc,son,xg,nbor) &
  !$acc copyin(state_diagram,ggg_acc,hhh_acc,oneontwotondim)
#endif  

#if defined(_OPENACC) || defined(USE_ACC_VERSION)
  if(pressure_fix)then
  !$acc kernels
     sum_nemisspe = 0
     sum_nrecpe   = 0 
     divu_acc = 0.0d0
     enew_acc = 0.0d0
  !$acc end kernels
  endif
#endif
  
#ifndef WITHOUTMPI
  tt2=MPI_WTIME(info)
  if(myid==1)write(*,*)'Time elapsed since startup:',tt2-tt1
#endif

  if(myid==1)then
     write(*,*)'Initial mesh structure'
     do ilevel=1,nlevelmax
        if(numbtot(1,ilevel)>0)write(*,999)ilevel,numbtot(1:4,ilevel)
     end do
  end if

  nstep_coarse_old=nstep_coarse

  if(myid==1)write(*,*)'Starting time integration' 

  do ! Main time loop

#ifndef WITHOUTMPI
     tt1=MPI_WTIME(info)
#endif

#ifdef MPROFILING
     tt1=MPI_WTIME(info)
#endif

     if(verbose)write(*,*)'Entering amr_step_coarse'

     epot_tot=0.0D0  ! Reset total potential energy
     ekin_tot=0.0D0  ! Reset total kinetic energy
     mass_tot=0.0D0  ! Reset total mass
     eint_tot=0.0D0  ! Reset total internal energy
#ifdef SOLVERmhd
     emag_tot=0.0D0  ! Reset total magnetic energy
#endif

     ! Make new refinements
     if(levelmin.lt.nlevelmax .and..not.static)then
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
#ifdef SOLVERmhd
              end do
#else
              end do
#endif
              if(simple_boundary)call make_boundary_hydro(ilevel)
           endif
#ifdef RT
           if(rt)then
              do ivar=1,nrtvar
                 call make_virtual_fine_dp(rtuold(1,ivar),ilevel)
              end do
              if(simple_boundary)call rt_make_boundary_hydro(ilevel)
           endif
#endif
           if(poisson)then
              call make_virtual_fine_dp(phi(1),ilevel)
              do idim=1,ndim
                 call make_virtual_fine_dp(f(1,idim),ilevel)
              end do
           end if
           if(ilevel<levelmin)call refine_fine(ilevel)
        end do
     endif

#ifdef MPROFILING
  tt1_amr_step=MPI_WTIME(info)
#endif

     ! Call base level
     call amr_step(levelmin,1)
     !$acc exit data delete(active,reception,emission,boundary)

#ifdef MPROFILING
  tt2_amr_step=MPI_WTIME(info)
  acc_t_amr_step = acc_t_amr_step + (tt2_amr_step-tt1_amr_step)
#endif

     if(levelmin.lt.nlevelmax .and..not. static)then
        do ilevel=levelmin-1,1,-1
           ! Hydro book-keeping
           if(hydro)then
              call upload_fine(ilevel)
#ifdef SOLVERmhd
              do ivar=1,nvar+3
#else
              do ivar=1,nvar
#endif
                 call make_virtual_fine_dp(uold(1,ivar),ilevel)
#ifdef SOLVERmhd
              end do
#else
              end do
#endif
              if(simple_boundary)call make_boundary_hydro(ilevel)
           end if
#ifdef RT
           ! Radiation book-keeping
           if(rt)then
              call rt_upload_fine(ilevel)
              do ivar=1,nrtvar
                 call make_virtual_fine_dp(rtuold(1,ivar),ilevel)
              end do
              if(simple_boundary)call rt_make_boundary_hydro(ilevel)
           end if
#endif
           ! Gravity book-keeping
           if(poisson)then
              call make_virtual_fine_dp(phi(1),ilevel)
              do idim=1,ndim
                 call make_virtual_fine_dp(f(1,idim),ilevel)
              end do
           end if
        end do
        
        ! Build refinement map
        do ilevel=levelmin-1,1,-1
           call flag_fine(ilevel,2)
        end do
        call flag_coarse
     endif

     ! New coarse time-step
     nstep_coarse=nstep_coarse+1

#ifdef MPROFILING
        tt2=MPI_WTIME(info)
        acc_t_total_time = acc_t_total_time + (tt2-tt1)
        if(myid==1)then
           print *,""
           print *,"!**********************************************************!"
           print *,"!Manual Profiling {MPROF Flag (Makefile)} (time in seconds)!"
           print *,"!**********************************************************!"
           print *,""
           print *,"!********************* Communication **********************!"
           print '(A40,2X,F7.3)',"make_virtual_fine_dp",acc_t_make_virtual_fine_dp
           print '(A40,2X,F7.3)',"make_virtual_reverse_dp",acc_t_make_virtual_reverse_dp
           print '(A40,2X,F7.3)',"make_virtual_mg_dp",acc_t_make_virtual_mg_dp
           print '(A40,2X,F7.3)',"make_reverse_mg_dp",acc_t_make_reverse_mg_dp
           print *,""
           print *,"!********************* Hyrdo Kernel ***********************!"
           print '(A40,2X,F7.3)',"godunov_fine",acc_t_godunov_fine
           print '(A40,2X,F7.3)',"set_unew",acc_t_set_unew
           print '(A40,2X,F7.3)',"set_uold",acc_t_set_uold
           print '(A40,2X,F7.3)',"add_gravity_source_terms",acc_t_add_gravity_source_terms
           print '(A40,2X,F7.3)',"add_pdv_source_terms",acc_t_add_pdv_source_terms
           print '(A40,2X,F7.3)',"synchro_hydro_fine",acc_t_synchro_hydro_fine
           print *,""
           print *,"!********************* Gravity Kernel *********************!"
           print '(A40,2X,F7.3)',"save_phi_old",acc_t_save_phi_old
           print *,"commons"
           print '(A40,2X,F7.3)',"multigrid_fine",acc_t_multigrid_fine
           print '(A40,2X,F7.3)',"make_initial_phi",acc_t_make_initial_phi
           print '(A40,2X,F7.3)',"make_multipole_phi",acc_t_make_multipole_phi
           print '(A40,2X,F7.3)',"make_fine_mask",acc_t_make_fine_mask
           print '(A40,2X,F7.3)',"make_boundary_mask",acc_t_make_boundary_mask
           print '(A40,2X,F7.3)',"make_fine_bc_rhs",acc_t_make_fine_bc_rhs
           print *,"fine"
           print '(A40,2X,F7.3)',"restrict_mask_fine_reverse",acc_t_restrict_mask_fine_reverse
           print '(A40,2X,F7.3)',"cmp_residual_mg_fine",acc_t_cmp_residual_mg_fine
           print '(A40,2X,F7.3)',"cmp_residual_norm2_fine",acc_t_cmp_residual_norm2_fine
           print '(A40,2X,F7.3)',"gauss_seidel_mg_fine",acc_t_gauss_seidel_mg_fine
           print '(A40,2X,F7.3)',"restrict_residual_fine_reverse",acc_t_restrict_residual_fine_reverse
           print '(A40,2X,F7.3)',"interpolate_and_correct_fine",acc_t_interpolate_and_correct_fine
           print '(A40,2X,F7.3)',"set_scan_flag_fine",acc_t_set_scan_flag_fine
           print *,"coarse"
           print '(A40,2X,F7.3)',"restrict_mask_coarse_reverse",acc_t_restrict_mask_coarse_reverse
           print '(A40,2X,F7.3)',"cmp_residual_mg_coarse",acc_t_cmp_residual_mg_coarse
           print '(A40,2X,F7.3)',"gauss_seidel_mg_coarse",acc_t_gauss_seidel_mg_coarse
           print '(A40,2X,F7.3)',"restrict_residual_coarse_reverse",acc_t_restrict_residual_coarse_reverse
           print '(A40,2X,F7.3)',"interpolate_and_correct_coarse",acc_t_interpolate_and_correct_coarse
           print '(A40,2X,F7.3)',"set_scan_flag_coarse",acc_t_set_scan_flag_coarse
           print *,"force fine"
           print '(A40,2X,F7.3)',"make_boundary_phi",acc_t_make_boundary_phi
           print '(A40,2X,F7.3)',"force_fine",acc_t_force_fine
           print *,""
           print *,"!**********************************************************!"
           print '(A40,2X,F7.3)',"amr_step",acc_t_amr_step
           print '(A40,2X,F7.3)',"Total Time",acc_t_total_time
           print *,"!**********************************************************!"
           print *,""
        end if
#endif

#ifndef WITHOUTMPI
     tt2=MPI_WTIME(info)
     if(mod(nstep_coarse,ncontrol)==0)then
        call getmem(real_mem)
        call MPI_ALLREDUCE(real_mem,real_mem_tot,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,info)
        if(myid==1)then
           write(*,*)'Time elapsed since last coarse step:',tt2-tt1
           call writemem(real_mem_tot)
        endif
     endif
#endif

  end do
  
  !$acc end data
  deallocate(nemisspe_acc,nrecpe_acc)
  deallocate(sum_nemisspe,sum_nrecpe)
  
#ifdef WITHOUTMPI  
#ifdef MPROFILING
  call MPI_FINALIZE(info)
#endif
#endif

999 format(' Level ',I2,' has ',I10,' grids (',3(I8,','),')')

end subroutine adaptive_loop
