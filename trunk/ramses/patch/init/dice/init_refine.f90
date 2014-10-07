!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_refine
  use amr_commons
  use pm_commons
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

end subroutine init_refine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_refine_2
  !--------------------------------------------------------------
  ! This routine builds additional refinements to the
  ! the initial AMR grid for filetype ne 'grafic'
  !--------------------------------------------------------------
  use amr_commons
  use hydro_commons
#ifdef RT
  use rt_hydro_commons
#endif
  use pm_commons
  use poisson_commons
  implicit none
  integer::ilevel,i,ivar

  if(filetype.eq.'grafic')return

  do i=levelmin,nlevelmax+1

     !!! Patch DICE
     do ilevel=levelmin-1,1,-1
        if(pic)call merge_tree_fine(ilevel)
     enddo
     !!! ----------
     call refine_coarse
     if(myid==1) write(*,*) "Level ",i,"->Loading GAS"
     do ilevel=1,nlevelmax
        call build_comm(ilevel)
        call make_virtual_fine_int(cpu_map(1),ilevel)
        call refine_fine(ilevel)
        !!! Patch DICE
        if(pic)call make_tree_fine(ilevel)
        !!! ----------
        if(hydro)call init_flow_fine(ilevel)
        !!! Patch DICE
        if(pic)then
           call kill_tree_fine(ilevel)
           call virtual_tree_fine(ilevel)
        endif
        !!! ----------
#ifdef RT
        if(rt)call rt_init_flow_fine(ilevel)
#endif
     end do
      
     !!! Patch DICE
     do ilevel=nlevelmax-1,levelmin,-1
        if(pic)call merge_tree_fine(ilevel)
     enddo
     !!! ----------
     if(nremap>0)call load_balance

     do ilevel=levelmin,nlevelmax
        if(pic)call make_tree_fine(ilevel)
        if(poisson)call rho_fine(ilevel,2)
        if(hydro)call init_flow_fine(ilevel)
        if(pic)then
           call kill_tree_fine(ilevel)
           call virtual_tree_fine(ilevel)
        endif
     end do

     do ilevel=nlevelmax,levelmin,-1
        if(pic)call merge_tree_fine(ilevel)
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
     end do

     do ilevel=nlevelmax,1,-1
        call flag_fine(ilevel,2)
     end do
     call flag_coarse

  end do
  call remove_gas_particles

#ifdef RT
  if(rt_is_init_xion .and. rt_nregion .eq. 0) then
     if(myid==1) write(*,*) 'Initializing ionization states from T profile'
     do ilevel=nlevelmax,1,-1
        call rt_init_xion(ilevel)
        call upload_fine(ilevel)
     end do
  endif
#endif  

end subroutine init_refine_2
!################################################################
!################################################################
!################################################################
!################################################################
subroutine remove_gas_particles
  use amr_commons
  use pm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::ngas_tot_all
  integer,dimension(1:ncpu)::ngas_icpu_all
#endif
  !----------------------------------------------------------------------
  ! This subroutine removes the gas particles initially present in the gadget1 DICE output
  ! Valentin Perret
  !----------------------------------------------------------------------
  ! local constants
  integer::ip,icpu,igrid,jgrid,npart1,npart2,ipart,jpart,next_part
  integer::ngas,ngas_loc,ngas_tot,info,igas,ilevel,ivar
  integer,dimension(1:ncpu)::ngas_icpu
  real(dp)::vol_min,nISM,nCOM,d0
  integer,dimension(:),allocatable::ind_part,ind_grid
  logical,dimension(:),allocatable::ok_free
  integer,dimension(:),allocatable::indgas

  if(.not. hydro)return
  if(ndim.ne.3)return

  if(verbose)write(*,*)'Entering remove_gas_particles'
  !------------------------------------------------------
  ! Gather gas particles
  !------------------------------------------------------
  ngas_loc=0
  ! Loop over levels
  do icpu=1,ncpu
  ! Loop over cpus
     igrid=headl(icpu,levelmin)
     ! Loop over grids
     do jgrid=1,numbl(icpu,levelmin)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0
        ! Count old enough GMC particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if(idp(ipart).eq.1)then
                 npart2=npart2+1
              endif
              ipart=next_part  ! Go to next particle
            end do
        endif
        ngas_loc=ngas_loc+npart2   ! Add gase to the total
        igrid=next(igrid)   ! Go to next grid
     end do
  end do
  ! End loop over levels
  ngas_icpu=0
  ngas_icpu(myid)=ngas_loc
#ifndef WITHOUTMPI
  ! Give an array of number of gas on each cpu available to all cpus
  call MPI_ALLREDUCE(ngas_icpu,ngas_icpu_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  ngas_icpu=ngas_icpu_all
#endif

  ngas_tot=sum(ngas_icpu(1:ncpu))

  if (ngas_tot .eq. 0) return
  
  if(myid==1)then
     write(*,*)'-----------------------------------------------'
     write(*,*)'Number of gas particles to delete=',ngas_tot
     write(*,*)'-----------------------------------------------'
  endif

  ! Allocate arrays for particles index and parent grid
  if(ngas_loc>0)then
     allocate(ind_part(1:ngas_loc),ind_grid(1:ngas_loc),ok_free(1:ngas_loc))
  endif

  !------------------------------------------------------
  ! Flag the index of gas particles to be removed
  !------------------------------------------------------
  if(myid==1)then
     igas=0
  else
     igas=sum(ngas_icpu(1:myid-1))
  endif
  ! Loop over levels
  ip=0
  do icpu=1,ncpu
     igrid=headl(icpu,levelmin)
     ! Loop over grids
     do jgrid=1,numbl(icpu,levelmin)
        npart1=numbp(igrid)  ! Number of particles in the grid
        ! Count old enough star particles that have not exploded
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              if(idp(ipart).eq.1)then
                 igas=igas+1
                 ip=ip+1
                 ind_grid(ip)=igrid
                 ind_part(ip)=ipart
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif
        igrid=next(igrid)   ! Go to next grid
     end do
  end do 
  ! End loop over levels

  ! Remove gas particle
  if(ngas_loc>0)then
     ok_free=.true.
     call remove_list(ind_part,ind_grid,ok_free,ngas_loc)
     call add_free_cond(ind_part,ok_free,ngas_loc)
     deallocate(ind_part,ind_grid,ok_free)
  endif
  do ipart=1,npart
    idp(ipart) = idp(ipart)-1
  enddo

end subroutine
