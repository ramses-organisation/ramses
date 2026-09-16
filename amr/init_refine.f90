!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_refine
  use amr_commons
  use pm_commons
#ifdef CRPHYS
  use cr_parameters, only: cr_advect
#endif
  implicit none
  !-------------------------------------------
  ! This routine builds the initial AMR grid
  !-------------------------------------------
  integer::ilevel

  if(myid==1)write(*,*)'Building initial AMR grid'
  init=.true.

  ! Base refinement
  do ilevel=1,levelmin
     call flag
     call refine
  end do

  ! Further refinements if necessary
  do ilevel=levelmin+1,nlevelmax
     if(initfile(levelmin).ne.' '.and.initfile(ilevel).eq.' ')exit
     if(hydro)call init_flow
#ifdef RT
     if(rt)call rt_init_flow
#endif
#ifdef CRPHYS
     if(cr_advect)call cr_init_flow
#endif
     if(ivar_refine==0)call init_refmap
     call flag
     call refine
     if(nremap>0)call load_balance
     if(numbtot(1,ilevel)==0)exit
  end do

  ! Final pass to initialize the flow
  init=.false.
  if(hydro)call init_flow
#ifdef RT
  if(rt)call rt_init_flow
#endif
#ifdef CRPHYS
  if(cr_advect)call cr_init_flow
#endif

end subroutine init_refine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_refine_2
  !--------------------------------------------------------------
  ! This routine builds additional refinements to the
  ! the initial AMR grid for filetype ne 'grafic'
  ! For DICE ICs, we ensure that all the particles are
  ! transfered down to level 1 before initialising the grid
  !--------------------------------------------------------------
  use amr_commons
  use hydro_commons
#ifdef RT
  use rt_hydro_commons
#endif
#ifdef CRPHYS
  use cr_parameters, only: cr_advect,ncrvar
  use cr_hydro_commons, only: cruold
#endif
  use pm_commons
  use poisson_commons
  use dice_commons
  implicit none
  integer::ilevel,i,ivar
  real(dp)::eps_star2

  if(filetype.ne.'grafic') then

     if(myid==1.and.dice_init) then
        write(*,*) "Initial conditions with AMR data structure"
        write(*,'(A50)')"__________________________________________________"
     end if
     do i=levelmin,nlevelmax+1
        do ilevel=levelmin-1,1,-1
           if(pic .and. dice_init)call merge_tree_fine(ilevel)
        enddo
        call refine_coarse
        do ilevel=1,nlevelmax
           call build_comm(ilevel)
           call make_virtual_fine_int(cpu_map(1),ilevel)
           call refine_fine(ilevel)
           if(pic.and.dice_init)call make_tree_fine(ilevel)
           if(hydro)call init_flow_fine(ilevel)
           if(pic.and.dice_init)then
              call kill_tree_fine(ilevel)
              call virtual_tree_fine(ilevel)
           endif
#ifdef RT
           if(rt)call rt_init_flow_fine(ilevel)
#endif
#ifdef CRPHYS
           if(cr_advect)call cr_init_flow_fine(ilevel)
#endif
        end do

        do ilevel=nlevelmax-1,levelmin,-1
           if(pic.and.dice_init)call merge_tree_fine(ilevel)
        enddo
        if(nremap>0)call load_balance

        do ilevel=levelmin,nlevelmax
           if(pic)call make_tree_fine(ilevel)
           if(poisson)call rho_fine(ilevel,2)
           if(hydro.and.dice_init)call init_flow_fine(ilevel)
           if(pic)then
              call kill_tree_fine(ilevel)
              call virtual_tree_fine(ilevel)
           endif
        end do

        do ilevel=nlevelmax,levelmin,-1
           if(pic)call merge_tree_fine(ilevel)
           if(hydro)then
              call upload_fine(ilevel)
              do ivar=1,nvar_all
                 call make_virtual_fine_dp(uold(1,ivar),ilevel)
              end do
              if(simple_boundary)call make_boundary_hydro(ilevel)
           endif
#ifdef RT
           if(rt)then
              call rt_upload_fine(ilevel)
              do ivar=1,nrtvar
                 call make_virtual_fine_dp(rtuold(1,ivar),ilevel)
              end do
              if(simple_boundary)call rt_make_boundary_hydro(ilevel)
           end if
#endif
#ifdef CRPHYS
           if(cr_advect)then
              call cr_upload_fine(ilevel)
              do ivar=1,ncrvar
                 call make_virtual_fine_dp(cruold(1,ivar),ilevel)
              end do
              if(simple_boundary)call cr_make_boundary_hydro(ilevel)
           end if
#endif
        end do

        do ilevel=nlevelmax,1,-1
           call flag_fine(ilevel,2)
        end do
        call flag_coarse

     end do
#if NDIM==3
     if(dice_init) then
        do ilevel=levelmin-1,1,-1
           if(pic)call merge_tree_fine(ilevel)
        enddo
        call kill_gas_part(1)
        do ilevel=1,nlevelmax
           if(pic)then
              call make_tree_fine(ilevel)
              call kill_tree_fine(ilevel)
              call virtual_tree_fine(ilevel)
           endif
        end do
        do ilevel=nlevelmax,levelmin,-1
           call merge_tree_fine(ilevel)
        end do
        deallocate(up)
        if(sf_virial)then
           eps_star2=eps_star
           eps_star=0d0
           do ilevel=nlevelmax,levelmin,-1
              call star_formation(ilevel)
           enddo
           eps_star=eps_star2
        endif
        dice_init=.false.
     end if
#endif
  endif ! if .not. 'grafic'


#ifdef RT
  if(neq_chem .and. rt_is_init_xion) then
     if(myid==1) write(*,*) 'Initializing ionization states from cell temperatures'
     do ilevel=nlevelmax,1,-1
        call rt_init_xion(ilevel)
        call upload_fine(ilevel)
        do ivar=1,nvar_all
           call make_virtual_fine_dp(uold(1,ivar),ilevel)
        end do
        if(simple_boundary)call make_boundary_hydro(ilevel)
     end do
  endif
#endif

end subroutine init_refine_2
!################################################################
!################################################################
!################################################################
!################################################################
