#if NDIM==3
subroutine clump_finder(create_output,keep_alive)
  use amr_commons
  use poisson_commons, ONLY:rho
  use clfind_commons
  use hydro_commons
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  logical::create_output,keep_alive

  !----------------------------------------------------------------------------
  ! Description of clump_finder:
  ! The clumpfinder detect first all cells having a density above
  ! a given threshold. These cells are linked to their densest neighbors,
  ! defining a peak patch around each local density maximum.
  ! If a so called peak patch is considered irrelevant, it is merged to its
  ! neighboring peak patch with the highest saddle density.
  ! Parameters are read in a namelist and they are:
  ! - density_threshold: defines the cell population to consider
  ! - relevance_threshold: merge peaks that are considered as ``noise''
  ! - saddle_threshold: for cosmo runs, merge peaks into halos (HOP-style)
  ! - mass_threshold: output only clumps (or halos) above this mass
  ! Andreas Bleuler & Davide Martizzi & Romain Teyssier
  !----------------------------------------------------------------------------

  integer::istep,nskip,ilevel,icpu,nmove,nzero
  integer::i,levelmin_part
  integer(i8b)::ntest_all,nmove_tot,nzero_tot
  integer(i8b),dimension(1:ncpu)::ntest_cpu,ntest_cpu_all
  integer,dimension(1:ncpu)::npeaks_per_cpu_tot
  logical::verbose_all=.false.

#ifndef WITHOUTMPI
  integer(i8b)::nmove_all,nzero_all
#endif

  if (create_output) then
    if(nstep_coarse==nstep_coarse_old.and.nstep_coarse>0)return
    if(nstep_coarse==0.and.nrestart>0)return
  endif

  if(verbose.and.myid==1)write(*,*)' Entering clump_finder'

  ! When called from the create_sink, particles are all residing at level 1,
  ! otherwise at levelmin.

  if (create_output)then
     levelmin_part = levelmin
  else
     levelmin_part = 1
  end if

  !---------------------------------------------------------------
  ! Compute rho from gas density or dark matter particles
  !---------------------------------------------------------------
  if(ivar_clump==0 .or. ivar_clump==-1)then
     do ilevel=levelmin_part,nlevelmax
        if(pic)call make_tree_fine(ilevel)
        if(poisson)call rho_only(ilevel)
        if(pic)then
           call kill_tree_fine(ilevel)
           call virtual_tree_fine(ilevel)
        endif
     end do
     do ilevel=nlevelmax,levelmin_part,-1
        if(pic)call merge_tree_fine(ilevel)
     end do
  endif

  !------------------------------------------------------------------------
  ! count the number of cells with density above the threshold
  ! flag the cells, share info across processors
  !------------------------------------------------------------------------
  ntest=0
  do ilevel=levelmin,nlevelmax
     if(ivar_clump==0 .or. ivar_clump==-1)then ! action 1: count and flag
        call count_test_particle(rho(1),ilevel,0,1)
     else
        if(hydro)then      ! action 1: count and flag
           call count_test_particle(uold(1,ivar_clump),ilevel,0,1)
        endif
     end if
  end do
  ntest_cpu=0; ntest_cpu_all=0
  ntest_cpu(myid)=ntest
#ifndef WITHOUTMPI
#ifndef LONGINT
  call MPI_ALLREDUCE(ntest_cpu,ntest_cpu_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#else
  call MPI_ALLREDUCE(ntest_cpu,ntest_cpu_all,ncpu,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
#endif
  ntest_cpu(1)=ntest_cpu_all(1)
#endif
  do icpu=2,ncpu
     ntest_cpu(icpu)=ntest_cpu(icpu-1)+ntest_cpu_all(icpu)
  end do
  ntest_all=ntest_cpu(ncpu)
  if(myid==1)then
     if(ntest_all.gt.0.and.clinfo)then
        write(*,'(" Total number of cells above threshold=",I12)')ntest_all
     endif
  end if

  !------------------------------------------------------------------------
  ! Allocate arrays and create list of cells above the threshold
  !------------------------------------------------------------------------
  if (ntest>0) then
     allocate(denp(ntest),levp(ntest),imaxp(ntest),icellp(ntest))
     denp=0d0; levp=0; imaxp=0; icellp=0
  endif
  itest=0
  nskip=ntest_cpu(myid)-ntest
  do ilevel=levelmin,nlevelmax
     if(ivar_clump==0 .or. ivar_clump==-1)then
        call count_test_particle(rho(1),ilevel,nskip,2)
     else
        if(hydro)then
           call count_test_particle(uold(1,ivar_clump),ilevel,nskip,2)
        endif
     endif
  end do
  do ilevel=nlevelmax,levelmin,-1
     call make_virtual_fine_int(flag2(1),ilevel)
  end do

  !-----------------------------------------------------------------------
  ! Sort cells above threshold according to their density
  !-----------------------------------------------------------------------
  if (ntest>0) then
     allocate(testp_sort(ntest))
     do i=1,ntest
        denp(i)=-denp(i)
        testp_sort(i)=i
     end do
     call quick_sort_dp(denp(1),testp_sort(1),ntest)
     deallocate(denp)
  endif

  !-----------------------------------------------------------------------
  ! Count number of density peaks and share info across processors
  !-----------------------------------------------------------------------
  npeaks=0
  if(ntest>0)then
     if(ivar_clump==0 .or. ivar_clump==-1)then  ! case 1: count peaks
        call count_peaks(rho(1),npeaks)
     else
        if(hydro)then       ! case 1: count peaks
           call count_peaks(uold(1,ivar_clump),npeaks)
        endif
     endif
  endif
  allocate(npeaks_per_cpu(1:ncpu))
  allocate(ipeak_start(1:ncpu))
  npeaks_per_cpu=0
  npeaks_per_cpu(myid)=npeaks
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(npeaks_per_cpu,npeaks_per_cpu_tot,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  npeaks_per_cpu_tot=npeaks_per_cpu
#endif
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(npeaks,npeaks_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  npeaks_tot=npeaks
#endif
  if (myid==1.and.npeaks_tot>0) &
       & write(*,'(" Total number of density peaks found=",I10)')npeaks_tot

  !----------------------------------------------------------------------
  ! Determine peak-ids positions for each cpu
  !----------------------------------------------------------------------
  ipeak_start=0
  npeaks_per_cpu=npeaks_per_cpu_tot
  do icpu=2,ncpu
     ipeak_start(icpu)=ipeak_start(icpu-1)+npeaks_per_cpu(icpu-1)
  end do
  nskip=ipeak_start(myid)

  !----------------------------------------------------------------------
  ! Flag peaks with global peak id using flag2 array
  ! Compute peak density using max_dens array
  !----------------------------------------------------------------------

  ! Compute the size of the peak-based arrays
  npeaks_max=MAX(4*maxval(npeaks_per_cpu_tot),1000)
  allocate(max_dens(npeaks_max))
  allocate(peak_cell(npeaks_max))
  allocate(peak_cell_level(npeaks_max))
  max_dens=0d0; peak_cell=0; peak_cell_level=0;
  flag2=0
  if(ntest>0)then
     if(ivar_clump==0 .or. ivar_clump==-1)then
        call flag_peaks(rho(1),nskip)
     else
        if(hydro)then
           call flag_peaks(uold(1,ivar_clump),nskip)
        endif
     endif
  endif
  do ilevel=nlevelmax,levelmin,-1
     call make_virtual_fine_int(flag2(1),ilevel)
  end do

  !---------------------------------------------------------------------
  ! Determine peak-patches around each peak
  ! Main step:
  ! - order cells in descending density
  ! - get peak id from densest neighbor
  ! - nmove is number of peak id's passed along
  ! - done when nmove_tot=0 (for single core, only one sweep is necessary)
  !---------------------------------------------------------------------
  if (myid==1.and.ntest_all>0)write(*,*)'Finding peak patches'
  nmove_tot=1
  istep=0
  do while (nmove_tot.gt.0)
     nmove=0
     nzero=0
     if(ntest>0)then
        call propagate_flag(nmove,nzero)
     endif
     do ilevel=nlevelmax,levelmin,-1
        call make_virtual_fine_int(flag2(1),ilevel)
     end do
     istep=istep+1
     nmove_tot=nmove
     nzero_tot=nzero
#ifndef WITHOUTMPI
#ifndef LONGINT
     call MPI_ALLREDUCE(nmove_tot,nmove_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(nzero_tot,nzero_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#else
     call MPI_ALLREDUCE(nmove_tot,nmove_all,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(nzero_tot,nzero_all,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
#endif
     nmove_tot=nmove_all
     nzero_tot=nzero_all
#endif
     if(ntest_all>0.and.verbose)write(*,*)"istep=",istep,"nmove=",nmove_tot
  end do

  !------------------------------------
  ! Allocate peak-patch property arrays
  !------------------------------------
  call allocate_peak_patch_arrays
  call build_peak_communicator

  if(npeaks_tot > 0)then
     !------------------------------------------
     ! Compute the saddle point density matrix
     !------------------------------------------
     if(ivar_clump==0 .or. ivar_clump==-1)then
        call saddlepoint_search(rho(1))
     else
        if(hydro)then
           call saddlepoint_search(uold(1,ivar_clump))
        endif
     endif
     call build_peak_communicator

     !------------------------------------------
     ! Merge irrelevant peaks
     !------------------------------------------
     if(myid==1.and.clinfo)write(*,*)"Now merging irrelevant peaks."
     call merge_clumps('relevance')
     do ilevel=nlevelmax,levelmin,-1
        call make_virtual_fine_int(flag2(1),ilevel)
     end do

     !------------------------------------------
     ! Compute clumps properties
     !------------------------------------------
     if(myid==1.and.clinfo)write(*,*)"Computing relevant clump properties."
     if(ivar_clump==0 .or. ivar_clump==-1)then
        call compute_clump_properties(rho(1))
     else
        if(hydro)then
           call compute_clump_properties(uold(1,ivar_clump))
        endif
     endif

     !------------------------------------------
     ! Merge clumps into haloes
     !------------------------------------------
     if(saddle_threshold>0)then
        if(myid==1.and.clinfo)write(*,*)"Now merging peaks into halos."
        call merge_clumps('saddleden')
     endif

     !------------------------------------------
     ! Output clumps properties to file
     !------------------------------------------
     if(verbose)then
        write(*,*)"Output status of peak memory."
     endif

#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(verbose,verbose_all,1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,info)
#else
     verbose_all=verbose
#endif

     if(verbose_all)call analyze_peak_memory
     if(clinfo.and.saddle_threshold.LE.0)call write_clump_properties(.false.)
     if(create_output.and..not.unbind)then
        ! if unbind, output will be written in unbinding() routine
        if(myid==1)write(*,*)"Outputing clump properties to disc."
        call write_clump_properties(.true.)
        if(ivar_clump==0 .or. ivar_clump==-1)then
           if(pic)call output_part_clump_id()
        endif
        ! output the clump field
        if (output_clump_field)then
           if(myid==1)write(*,*)"Outputing clump field to disc"
           call write_clump_field
        end if
     endif

  end if

  !-----------------------------------------------------------------
  ! Call particle unbinding (and mergertree stuff)
  ! Call it even for npeaks_tot = 0, mergertrees need to know that
  ! there are no progenitors to work with
  !------------------------------------------------------------------
  if(unbind.and.create_output.and.pic) call unbinding()

  if (.not. keep_alive)then
     ! Deallocate test particle and peak arrays
     deallocate(npeaks_per_cpu)
     deallocate(ipeak_start)
     if (ntest>0)then
        deallocate(icellp)
        deallocate(levp)
        deallocate(testp_sort)
        deallocate(imaxp)
     endif
     call deallocate_all
  end if

end subroutine clump_finder
!################################################################
!################################################################
!################################################################
!################################################################
subroutine count_test_particle(xx,ilevel,nskip,action)
  use amr_commons
  use clfind_commons
  use mpi_mod
  implicit none
  integer::ilevel,nskip,action
  real(dp),dimension(1:ncoarse+ngridmax*twotondim)::xx

  !----------------------------------------------------------------------
  ! Description: This routine loops over all cells above and checks wether
  ! their density lies above the threshold. If so:
  ! case 1: count the new test particles and flag the cell
  ! case 2: create the test particle
  ! xx is on input the array containing the density field
  !----------------------------------------------------------------------

  integer ::ncache,ngrid
  integer ::igrid,ind,i,iskip
  integer ,dimension(1:nvector)::ind_grid,ind_cell
  logical ,dimension(1:nvector)::ok


  if(numbtot(1,ilevel)==0) return

  if(verbose .and. myid==1)then
     write(*,*)' Entering count test particle for level=',&
          & ilevel,' and action=',action
  endif

  ! Loop over grids
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     ! loop over cells
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do

        !checks
        do i=1,ngrid
           ok(i)=son(ind_cell(i))==0 !check if leaf cell
           ok(i)=ok(i).and.xx(ind_cell(i))>density_threshold !check density
        end do

        select case (action)
        case (1) !count and flag
           ! Compute test particle map
           do i=1,ngrid
              flag2(ind_cell(i))=0
              if(ok(i))then
                 flag2(ind_cell(i))=1
                 ntest=ntest+1
              endif
           end do
        case(2) !create 'testparticles'
           do i=1,ngrid
              if (ok(i))then
                 itest=itest+1                   ! Local 'test particle' index
                 levp(itest)=ilevel              ! Level
                 flag2(ind_cell(i))=itest+nskip  ! Initialize flag2 to GLOBAL test particle index
                 icellp(itest)=ind_cell(i)       ! Local cell index
                 denp(itest)=xx(ind_cell(i))     ! Save density values here!
              end if
           end do
        end select
     end do
  end do

end subroutine count_test_particle
!################################################################
!################################################################
!################################################################
!################################################################
subroutine count_peaks(xx,n)
  use amr_commons
  use clfind_commons
  use mpi_mod
  implicit none
  integer::n
  real(dp),dimension(1:ncoarse+ngridmax*twotondim)::xx
  !----------------------------------------------------------------------
  ! Count the peaks (cell with no denser neighbor)
  ! Store the index of the densest neighbor (can be the cell itself)
  ! for later usage.
  !----------------------------------------------------------------------
  integer::ilevel,next_level,ipart,jpart,ip
  integer,dimension(1:nvector)::ind_part,ind_cell,ind_max

  ! Group chunks of nvector cells (of the same level) and send them
  ! to the routine that constructs the neighboring cells
  ip=0
  do ipart=1,ntest
     ip=ip+1
     ilevel=levp(testp_sort(ipart))  ! level
     next_level=0                    ! level of next particle
     if(ipart<ntest)next_level=levp(testp_sort(ipart+1))
     ind_cell(ip)=icellp(testp_sort(ipart))
     ind_part(ip)=testp_sort(ipart)
     if(ip==nvector .or. next_level /= ilevel)then
        call neighborsearch(xx(1),ind_cell,ind_max,ip,n,ilevel,1)
        do jpart=1,ip
           imaxp(ind_part(jpart))=ind_max(jpart)
        end do
        ip=0
     endif
  end do
  if (ip>0)then
     call neighborsearch(xx(1),ind_cell,ind_max,ip,n,ilevel,1)
     do jpart=1,ip
        imaxp(ind_part(jpart))=ind_max(jpart)
     end do
  endif

end subroutine count_peaks
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine flag_peaks(xx,ipeak)
  use amr_commons
  use clfind_commons
  use mpi_mod
  implicit none
  integer::ipeak
  real(dp),dimension(1:ncoarse+ngridmax*twotondim)::xx
  !----------------------------------------------------------------------
  ! Flag (flag2 array) all cells that host a peak with the global peak id
  !----------------------------------------------------------------------
  integer::ipart,jpart
  do ipart=1,ntest
     jpart=testp_sort(ipart)
     if(imaxp(jpart).EQ.-1)then
        ipeak=ipeak+1
        flag2(icellp(jpart))=ipeak
        max_dens(ipeak-ipeak_start(myid))=xx(icellp(jpart))
        peak_cell(ipeak-ipeak_start(myid))=icellp(jpart)
        peak_cell_level(ipeak-ipeak_start(myid))=levp(jpart)
     endif
  end do
end subroutine flag_peaks
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine propagate_flag(nmove,nzero)
  use amr_commons
  use clfind_commons
  use mpi_mod
  implicit none
  integer::nmove,nzero
  !----------------------------------------------------------------------
  ! All cells above the threshold copy the flag2 value from their densest
  ! neighbors. Because of MPI boundaries, flag2=0 can be passed around.
  ! Once all flag2=0 values have disappiered (globally), we're done.
  !----------------------------------------------------------------------
  integer::ipart,jpart
  do ipart=1,ntest
     jpart=testp_sort(ipart)
     if(imaxp(jpart).NE.-1)then
        if(flag2(icellp(jpart)).ne.flag2(imaxp(jpart)))nmove=nmove+1
        flag2(icellp(jpart))=flag2(imaxp(jpart))
        if(flag2(icellp(jpart)).eq.0)nzero=nzero+1
     endif
  end do
end subroutine propagate_flag
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine saddlepoint_search(xx)
  use amr_commons
  use clfind_commons
  use mpi_mod
  implicit none
  real(dp),dimension(1:ncoarse+ngridmax*twotondim)::xx
  !---------------------------------------------------------------------------
  ! subroutine which creates a npeaks**2 sized array of saddlepoint densities
  ! by looping over all testparticles and passing them to neighborcheck with
  ! case 4, which means that saddlecheck will be called for each neighboring
  ! leaf cell. There it is checked, whether the two cells (original cell and
  ! neighboring cell) are connected by a new densest saddle.
  !---------------------------------------------------------------------------
  integer::ipart,ip,ilevel,next_level
  integer::dummyint
  integer,dimension(1:nvector)::ind_cell,ind_max

  ip=0
  do ipart=1,ntest
     ip=ip+1
     ilevel=levp(testp_sort(ipart)) ! level
     next_level=0 !level of next particle
     if(ipart<ntest)next_level=levp(testp_sort(ipart+1))
     ind_cell(ip)=icellp(testp_sort(ipart))
     if(ip==nvector .or. next_level /= ilevel)then
        call neighborsearch(xx(1),ind_cell,ind_max,ip,dummyint,ilevel,4)
        ip=0
     endif
  end do
  if (ip>0)call neighborsearch(xx(1),ind_cell,ind_max,ip,dummyint,ilevel,4)

end subroutine saddlepoint_search
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine neighborsearch(xx,ind_cell,ind_max,np,count,ilevel,action)
  use amr_commons
  implicit none
  integer::np,count,ilevel,action
  integer,dimension(1:nvector)::ind_max,ind_cell
  real(dp),dimension(1:ncoarse+ngridmax*twotondim)::xx

  !------------------------------------------------------------
  ! This routine constructs all neighboring leaf cells at levels
  ! ilevel-1, ilevel, ilevel+1.
  ! Depending on the action case value, fuctions performing
  ! further checks for the neighbor cells are called.
  ! xx is on input the array containing the density field
  !------------------------------------------------------------

  integer::j,ind,nx_loc,i1,j1,k1,i2,j2,k2,i3,j3,k3,ix,iy,iz
  integer::i1min,i1max,j1min,j1max,k1min,k1max
  integer::i2min,i2max,j2min,j2max,k2min,k2max
  integer::i3min,i3max,j3min,j3max,k3min,k3max
  real(dp)::dx,dx_loc,scale,vol_loc
  integer ,dimension(1:nvector)::clump_nr,indv,ind_grid,grid,ind_cell_coarse

  real(dp),dimension(1:twotondim,1:3)::xc
  integer ,dimension(1:99)::cell_index,cell_levl,test_levl
  real(dp),dimension(1:99,1:ndim)::xtest,xrel
  logical ,dimension(1:99)::ok
  real(dp),dimension(1:nvector)::density_max
  real(dp),dimension(1:3)::skip_loc
  logical ,dimension(1:nvector)::okpeak
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:threetondim)::nbors_father_cells_pass
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  integer::ntestpos,ntp,idim,ipos

#if NDIM==3
  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**3

  ! Integer constants
  i1min=0; i1max=0; i2min=0; i2max=0; i3min=0; i3max=0
  j1min=0; j1max=0; j2min=0; j2max=0; j3min=0; j3max=0
  k1min=0; k1max=0; k2min=0; k2max=0; k3min=0; k3max=0
  if(ndim>0)then
     i1max=1; i2max=2; i3max=3
  end if
  if(ndim>1)then
     j1max=1; j2max=2; j3max=3
  end if
  if(ndim>2)then
     k1max=1; k2max=2; k3max=3
  end if

  ! Cells center position relative to grid center position
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     xc(ind,1)=(dble(ix)-0.5D0)*dx
     xc(ind,2)=(dble(iy)-0.5D0)*dx
     xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  ! some preliminary action...
  do j=1,np
     indv(j)=(ind_cell(j)-ncoarse-1)/ngridmax+1 ! cell position in grid
     ind_grid(j)=ind_cell(j)-ncoarse-(indv(j)-1)*ngridmax ! grid index
     density_max(j)=xx(ind_cell(j))*1.0001d0 ! get cell density (1.0001 probably not necessary)
     ind_max(j)=ind_cell(j) !save cell index
     if (action.ge.4)clump_nr(j)=flag2(ind_cell(j)) ! save clump number
  end do

  ntestpos=3**ndim
  if(ilevel>levelmin)ntestpos=ntestpos+2**ndim
  if(ilevel<nlevelmax)ntestpos=ntestpos+4**ndim

  !================================
  ! generate neighbors level ilevel-1
  !================================
  ntp=0
  if(ilevel>levelmin)then
     ! Generate 2x2x2  neighboring cells at level ilevel-1
     do k1=k1min,k1max
        do j1=j1min,j1max
           do i1=i1min,i1max
              ntp=ntp+1
              xrel(ntp,1)=(2*i1-1)*dx_loc
              xrel(ntp,2)=(2*j1-1)*dx_loc
              xrel(ntp,3)=(2*k1-1)*dx_loc
              test_levl(ntp)=ilevel-1
           end do
        end do
     end do
  endif

  !================================
  ! generate neighbors at level ilevel
  !================================
  ! Generate 3x3x3 neighboring cells at level ilevel
  do k2=k2min,k2max
     do j2=j2min,j2max
        do i2=i2min,i2max
           ntp=ntp+1
           xrel(ntp,1)=(i2-1)*dx_loc
           xrel(ntp,2)=(j2-1)*dx_loc
           xrel(ntp,3)=(k2-1)*dx_loc
           test_levl(ntp)=ilevel
        end do
     end do
  end do

  !===================================
  ! generate neighbors at level ilevel+1
  !====================================
  if(ilevel<nlevelmax)then
     ! Generate 4x4x4 neighboring cells at level ilevel+1
     do k3=k3min,k3max
        do j3=j3min,j3max
           do i3=i3min,i3max
              ntp=ntp+1
              xrel(ntp,1)=(i3-1.5d0)*dx_loc/2
              xrel(ntp,2)=(j3-1.5d0)*dx_loc/2
              xrel(ntp,3)=(k3-1.5d0)*dx_loc/2
              test_levl(ntp)=ilevel+1
           end do
        end do
     end do
  endif



  ! Gather 27 neighboring father cells (should be present anytime !)
  do j=1,np
     ind_cell_coarse(j)=father(ind_grid(j))
  end do
  call get3cubefather(ind_cell_coarse,nbors_father_cells,nbors_father_grids,np,ilevel)


  ! initialze logical array
  okpeak=.true.

  do j=1,np
     ok=.false.
     do idim=1,ndim
        xtest(1:ntestpos,idim)=(xg(ind_grid(j),idim)+xc(indv(j),idim)-skip_loc(idim))*scale+xrel(1:ntestpos,idim)
        if(ilevel>levelmin)xtest(1:twotondim,idim)=xtest(1:twotondim,idim)+xc(indv(j),idim)*scale
     end do
     grid(1)=ind_grid(j)
     nbors_father_cells_pass=nbors_father_cells(j,1:threetondim)
     call get_cell_index_fast(cell_index,cell_levl,xtest,ind_grid(j),nbors_father_cells_pass,ntestpos,ilevel)

     do ipos=1,ntestpos
        if(son(cell_index(ipos))==0.and.cell_levl(ipos)==test_levl(ipos))ok(ipos)=.true.
     end do

     ! check those neighbors
     if (action==1)call peakcheck(xx(1),cell_index,okpeak(j),ok,density_max(j),ind_max(j),ntestpos)
     if (action==4)call saddlecheck(xx(1),ind_cell(j),cell_index,clump_nr(j),ok,ntestpos)

  end do

  ! Count peaks (only one action case left as case 2 and 3 are dealt with
  ! outside neighborsearch
  if (action==1) then
     do j=1,np
        if(okpeak(j))then
           count=count+1
           ind_max(j)=-1
        endif
     end do
  end if

#endif

end subroutine neighborsearch
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine peakcheck(xx,cell_index,okpeak,ok,density_max,ind_max,np)
  use amr_commons
  implicit none
  !----------------------------------------------------------------------
  ! routine to check wether neighbor is denser or not
  !----------------------------------------------------------------------
  logical,dimension(1:99)::ok
  logical::okpeak
  integer,dimension(1:99)::cell_index
  integer::ind_max
  real(dp)::density_max
  real(dp),dimension(1:ncoarse+ngridmax*twotondim)::xx
  integer::np,j


  do j=1,np
     ! only consider leaf-cells at correct level
     if(ok(j))then
        ! if cell is denser than densest neighbor
        if(xx(cell_index(j))>density_max)then
           okpeak=.false.                 ! cell is no peak
           density_max=xx(cell_index(j))  ! change densest neighbor dens
           ind_max=cell_index(j)          ! change densest neighbor index
        endif
     end if
  end do

end subroutine peakcheck
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine saddlecheck(xx,ind_cell,cell_index,clump_nr,ok,np)
  use amr_commons
  use clfind_commons, ONLY: sparse_saddle_dens
  use sparse_matrix
  implicit none
  !----------------------------------------------------------------------
  ! routine to check wether neighbor is connected through new densest saddle
  !----------------------------------------------------------------------
  logical,dimension(1:99)::ok
  integer,dimension(1:99)::cell_index,neigh_cl
  real(dp),dimension(1:99)::av_dens
  real(dp),dimension(1:ncoarse+ngridmax*twotondim)::xx
  integer::np,j,ipeak,jpeak,clump_nr,ind_cell

  do j=1,np
     neigh_cl(j)=flag2(cell_index(j))!index of the clump the neighboring cell is in
  end do
  do j=1,np
     ok(j)=ok(j).and. clump_nr/=0 ! temporary fix...
     ok(j)=ok(j).and. neigh_cl(j)/=0 !neighboring cell is in a clump
     ok(j)=ok(j).and. neigh_cl(j)/=clump_nr !neighboring cell is in another clump
     av_dens(j)=(xx(cell_index(j))+xx(ind_cell))/2 !average density of cell and neighbor cell
  end do
  do j=1,np
     if(ok(j))then ! if all criteria met, replace saddle density array value
        call get_local_peak_id(clump_nr,ipeak)
        call get_local_peak_id(neigh_cl(j),jpeak)
        if (get_value(ipeak,jpeak,sparse_saddle_dens) < av_dens(j))then
           call set_value(ipeak,jpeak,av_dens(j),sparse_saddle_dens)
        end if
        if (get_value(jpeak,ipeak,sparse_saddle_dens) < av_dens(j))then
           call set_value(jpeak,ipeak,av_dens(j),sparse_saddle_dens)
        end if
     end if
  end do

end subroutine saddlecheck
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine get_cell_index(cell_index,cell_levl,xpart,ilevel,n)
  use amr_commons
  implicit none

  integer::n,ilevel
  integer,dimension(1:nvector)::cell_index,cell_levl
  real(dp),dimension(1:nvector,1:3)::xpart

  !----------------------------------------------------------------------------
  ! This routine returns the index and level of the cell, (at maximum level
  ! ilevel), in which the input the position specified by xpart lies
  !----------------------------------------------------------------------------

  real(dp)::xx,yy,zz
  integer::i,j,ii,jj,kk,ind,iskip,igrid,ind_cell,igrid0

  if ((nx.eq.1).and.(ny.eq.1).and.(nz.eq.1)) then
  else if ((nx.eq.3).and.(ny.eq.3).and.(nz.eq.3)) then
  else
     write(*,*)"nx=ny=nz != 1,3 is not supported."
     call clean_stop
  end if

  ind_cell=0
  igrid0=son(1+icoarse_min+jcoarse_min*nx+kcoarse_min*nx*ny)
  do i=1,n
     xx = xpart(i,1)/boxlen + (nx-1)/2d0
     yy = xpart(i,2)/boxlen + (ny-1)/2d0
     zz = xpart(i,3)/boxlen + (nz-1)/2d0

     if(xx<0.)xx=xx+dble(nx)
     if(xx>dble(nx))xx=xx-dble(nx)
     if(yy<0.)yy=yy+dble(ny)
     if(yy>dble(ny))yy=yy-dble(ny)
     if(zz<0.)zz=zz+dble(nz)
     if(zz>dble(nz))zz=zz-dble(nz)

     igrid=igrid0
     do j=1,ilevel
        ii=1; jj=1; kk=1
        if(xx<xg(igrid,1))ii=0
        if(yy<xg(igrid,2))jj=0
        if(zz<xg(igrid,3))kk=0
        ind=1+ii+2*jj+4*kk
        iskip=ncoarse+(ind-1)*ngridmax
        ind_cell=iskip+igrid
        igrid=son(ind_cell)
        if(igrid==0.or.j==ilevel)exit
     end do
     cell_index(i)=ind_cell
     cell_levl(i)=j
  end do
end subroutine get_cell_index
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine read_clumpfind_params()
  use clfind_commons
  use mpi_mod
  implicit none

  namelist/clumpfind_params/ivar_clump,&
       & relevance_threshold,density_threshold,&
       & saddle_threshold,mass_threshold,clinfo,&
       & n_clfind,rho_clfind,age_cut_clfind,output_clump_field
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v

  ! Read namelist file
  rewind(1)
  read(1,NML=clumpfind_params,END=121)
  goto 122
121 if(myid==1)write(*,*)'You did not set up namelist &CLUMPFIND_PARAMS in parameter file.'

  if (.not. sink)then
     if(myid==1)write(*,*)'That block should a least contain a density '
     if(myid==1)write(*,*)'threshold n_clfind [parts/cc] or rho_clfind [g/cc]!'
     if(myid==1)write(*,*)'aborting...'
     call clean_stop
  end if

122 rewind(1)

  if (density_threshold>0.)then
     if (rho_clfind>0. .or. n_clfind >0.)then
        if(myid==1)write(*,*)'you provided the density threshold in code units.'
        if(myid==1)write(*,*)'Ignoring the input in physical units...'
     end if
  else
     if (cosmo) then
        if(myid==1)write(*,*)'For cosmological simulations you have to specify the'
        if(myid==1)write(*,*)'clumpfinder density threshold in code units.'
        call clean_stop
     end if
     call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
     if (rho_clfind>0. .and. n_clfind >0.)then   ! too much information...
        if(myid==1)write(*,*)'you set up the clumpfinder threshold in both, H/cc and g/cc, decide!'
        if(myid==1)write(*,*)'aborting...'
        call clean_stop
     else if (rho_clfind<0. .and. n_clfind <0.)then  !not enough information
        if (sink)then
           density_threshold=d_sink/10
           if(myid==1)write(*,*)'You did not specify a threshold for the clump finder. '
           if(myid==1)write(*,*)'Setting it to sink threshold / 10. !'
        else
           if(myid==1)write(*,*)'The &CLUMPFIND_PARAMS block should a least contain '
           if(myid==1)write(*,*)'density_threshold [code units], n_clfind [parts/cc]'
           if(myid==1)write(*,*)'or rho_clfind [g/cc]!'
           if(myid==1)write(*,*)'aborting...'
           call clean_stop
        end if
     else if (n_clfind>0.)then
        density_threshold=n_clfind/scale_nH
     else if(rho_clfind>0.)then
        density_threshold=rho_clfind/scale_d
     end if
  end if

end subroutine read_clumpfind_params
!################################################################
!################################################################
!################################################################
!################################################################
subroutine get_cell_index_fast(indp,cell_lev,xpart,ind_grid,nbors_father_cells,np,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
  integer::np,ilevel,ind_grid
  integer,dimension(1:99)::indp,cell_lev
  real(dp),dimension(1:99,1:ndim)::xpart
  integer ,dimension(1:threetondim)::nbors_father_cells

  !-----------------------------------------------------------------------
  ! This subroutine finds the leaf cell in which a particle sits
  !-----------------------------------------------------------------------

  integer::j,idim,nx_loc,ind,ix,iy,iz
  real(dp)::dx,dx_loc,scale,one_over_dx,one_over_scale
  ! Grid based arrays
  real(dp),dimension(1:ndim),save::x0
  real(dp),dimension(1:99,1:ndim),save::x
  integer ,dimension(1:99,1:ndim),save::id,igd,icd,icd_fine
  integer ,dimension(1:99),save::igrid,icell,kg,icell_fine
  real(dp),dimension(1:3),save::skip_loc
  real(dp),dimension(1:twotondim,1:3),save::xc
  logical,dimension(1:99),save::ok

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=0
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  one_over_dx=1/dx
  one_over_scale=1/scale

  ! Cells center position relative to grid center position
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     xc(ind,1)=(dble(ix)-0.5D0)*dx
     xc(ind,2)=(dble(iy)-0.5D0)*dx
     xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  ! Lower left corner of 3x3x3 grid-cube
  do idim=1,ndim
     x0(idim)=xg(ind_grid,idim)-3.0D0*dx
  end do

  ! Rescale position at level ilevel
  do idim=1,ndim
     do j=1,np
        x(j,idim)=xpart(j,idim)*one_over_scale+skip_loc(idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)-x0(idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)*one_over_dx
     end do
  end do

  ! NGP at level ilevel
  do idim=1,ndim
     do j=1,np
        id(j,idim)=int(x(j,idim))
     end do
  end do

   ! Compute parent grids
  do idim=1,ndim
     do j=1,np
        igd(j,idim)=id(j,idim)/2
     end do
  end do
#if NDIM==1
  do j=1,np
     kg(j)=1+igd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     kg(j)=1+igd(j,1)+3*igd(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
     kg(j)=1+igd(j,1)+3*igd(j,2)+9*igd(j,3)
  end do
#endif

  do j=1,np
     igrid(j)=son(nbors_father_cells(kg(j)))
  end do

  ! Check if particle has escaped to ilevel-1
  ok(1:np)=.true.
  do j=1,np
     if (igrid(j)==0)then
        ok(j)=.false.
        indp(j)=nbors_father_cells(kg(j))
        cell_lev(j)=ilevel-1
     end if
  end do

  ! Compute parent cell position
  do idim=1,ndim
     do j=1,np
        icd(j,idim)=id(j,idim)-2*igd(j,idim)
     end do
  end do

  call geticell99(icell,icd,np)

  ! Compute parent cell adress
  do j=1,np
     if(ok(j))then
        indp(j)=ncoarse+(icell(j)-1)*ngridmax+igrid(j)
     end if
  end do

  ! Check if particles have leaked into level ilevel+1
  do j=1,np
     if(ok(j))then
        if (son(indp(j))>0)then
           ok(j)=.false.
           cell_lev(j)=ilevel+1
           do idim=1,ndim
              icd_fine(1,idim)=int(2*(x(j,idim)-int(x(j,idim))))
           end do
           call geticell99(icell_fine,icd_fine,1)
           indp(j)=ncoarse+(icell_fine(1)-1)*ngridmax+son(indp(j))
        end if
     endif
  end do

  ! Cell center positions for particles which sit in the level ilevel
  do j=1,np
     if (ok(j))then
        cell_lev(j)=ilevel
     end if
  end do

end subroutine get_cell_index_fast
!################################################################
!################################################################
!################################################################
!################################################################
subroutine geticell99(icell,icd,np)
  use amr_parameters, only:ndim
  integer::np
  integer,dimension(1:99,1:ndim)::icd
  integer,dimension(1:99)::icell
  ! mini subroutine that gets the cell index (1 to 8)
  ! for certain coordinates (0,1 along each direction)
  ! put into a subroutine to make the code look less ugly
  integer::j

#if NDIM==1
  do j=1,np
     icell(j)=1+icd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     icell(j)=1+icd(j,1)+2*icd(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
     icell(j)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
  end do
#endif

end subroutine geticell99
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine rho_only(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  use mpi_mod
  implicit none
  integer::ilevel
  !------------------------------------------------------------------
  ! This routine computes the density field at level ilevel using
  ! the CIC scheme. Particles that are not entirely in
  ! level ilevel contribute also to the level density field
  ! (boundary particles) using buffer grids.
  !------------------------------------------------------------------
  integer::iskip,icpu,ind,i,nx_loc,ibound
  real(dp)::dx,scale,dx_loc

  if(.not. poisson)return
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  !-------------------------------------------------------
  ! Initialize rho to analytical and baryon density field
  !-------------------------------------------------------
  if(ilevel==levelmin)then
     do i=nlevelmax,ilevel,-1
        ! Compute mass multipole
        if(hydro)call multipole_fine(i)
        ! Perform TSC using pseudo-particle
#ifdef TSC
        if (ndim==3)then
           call tsc_from_multipole(i)
        else
           write(*,*)'TSC not supported for ndim neq 3'
           call clean_stop
        end if
#else
        ! Perform CIC using pseudo-particle
        call cic_from_multipole(i)
#endif
        ! Update boundaries
        call make_virtual_reverse_dp(rho(1),i)
        call make_virtual_fine_dp   (rho(1),i)
     end do
  end if

  !-------------------------------------------------------
  ! Initialize rho to zero in virtual boundaries
  !-------------------------------------------------------
  do icpu=1,ncpu
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,reception(icpu,ilevel)%ngrid
#ifdef LIGHT_MPI_COMM
           rho(reception(icpu,ilevel)%pcomm%igrid(i)+iskip)=0.0D0
#else
           rho(reception(icpu,ilevel)%igrid(i)+iskip)=0.0D0
#endif
        end do
     end do
  end do

  !---------------------------------------------------------
  ! Compute particle contribution to density field
  !---------------------------------------------------------
  ! Compute density due to current level particles
  if(pic)then
     call rho_only_level(ilevel)
  end if
  ! Update boudaries
  call make_virtual_reverse_dp(rho(1),ilevel)
  call make_virtual_fine_dp   (rho(1),ilevel)

  !----------------------------------------------------
  ! Reset rho in physical boundaries
  !----------------------------------------------------
  do ibound=1,nboundary
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,boundary(ibound,ilevel)%ngrid
           rho(boundary(ibound,ilevel)%igrid(i)+iskip)=0
        end do
     end do
  end do

111 format('   Entering rho_only for level ',I2)

end subroutine rho_only
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine rho_only_level(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  use clfind_commons
  implicit none
  integer::ilevel
  !------------------------------------------------------------------
  ! This routine computes the density field at level ilevel using
  ! the CIC scheme from particles that are not entirely in
  ! level ilevel (boundary particles).
  ! Arrays flag1 and flag2 are used as temporary work space.
  !------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,idim,icpu,next_part
  integer::i,ig,ip,npart1,npart2
  real(dp)::dx

  integer,dimension(1:nvector),save::ind_grid,ind_cell
  integer,dimension(1:nvector),save::ind_part,ind_grid_part
  real(dp),dimension(1:nvector,1:ndim),save::x0

  ! Mesh spacing in that level
  dx=0.5D0**ilevel

  ! Loop over cpus
  do icpu=1,ncpu
     ! Loop over grids
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0

        ! Count elligible particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              ! Select stars younger than age_cut_clfind
              if(age_cut_clfind>0d0 .and. star .and. use_proper_time) then
                 if((is_star(typep(ipart))).and.(t-tp(ipart).lt.age_cut_clfind).and.(tp(ipart).ne.0d0)) then
                    npart2=npart2+1
                 endif
              ! All particles
              else
                 npart2=npart2+1
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif

        ! Gather elligible particles
        if(npart2>0)then
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)

           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              ! Select stars younger than age_cut_clfind
              if(age_cut_clfind>0d0 .and. star .and. use_proper_time) then
                 if((is_star(typep(ipart))).and.(t-tp(ipart).lt.age_cut_clfind).and.(tp(ipart).ne.0d0)) then
                    if(ig==0)then
                       ig=1
                       ind_grid(ig)=igrid
                    end if
                    ip=ip+1
                    ind_part(ip)=ipart
                    ind_grid_part(ip)=ig
                 endif
              ! All particles
              else
                 if(ig==0)then
                    ig=1
                    ind_grid(ig)=igrid
                 end if
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig
              endif
              if(ip==nvector)then
                 ! Lower left corner of 3x3x3 grid-cube
                 do idim=1,ndim
                    do i=1,ig
                       x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
                    end do
                 end do
                 do i=1,ig
                    ind_cell(i)=father(ind_grid(i))
                 end do
#ifdef TSC
                 call tsc_only(ind_cell,ind_part,ind_grid_part,x0,ig,ip,ilevel)
#else
                 call cic_only(ind_cell,ind_part,ind_grid_part,x0,ig,ip,ilevel)
#endif

                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles

        end if

        igrid=next(igrid)   ! Go to next grid
     end do
     ! End loop over grids

     if(ip>0)then
        ! Lower left corner of 3x3x3 grid-cube
        do idim=1,ndim
           do i=1,ig
              x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
           end do
        end do
        do i=1,ig
           ind_cell(i)=father(ind_grid(i))
        end do
#ifdef TSC
        call tsc_only(ind_cell,ind_part,ind_grid_part,x0,ig,ip,ilevel)
#else
        call cic_only(ind_cell,ind_part,ind_grid_part,x0,ig,ip,ilevel)
#endif
     end if

  end do
  ! End loop over cpus

end subroutine rho_only_level
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine cic_only(ind_cell,ind_part,ind_grid_part,x0,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use poisson_commons
  use clfind_commons
  implicit none
  integer::ng,np,ilevel
  integer ,dimension(1:nvector)::ind_cell,ind_grid_part,ind_part
  real(dp),dimension(1:nvector,1:ndim)::x0
  !------------------------------------------------------------------
  ! This routine computes the density field at level ilevel using
  ! the CIC scheme. Only cells that are in level ilevel
  ! are updated by the input particle list.
  !------------------------------------------------------------------
  logical::error
  integer::j,ind,idim,nx_loc
  real(dp)::dx,dx_loc,scale,vol_loc
  ! Grid-based arrays
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle-based arrays
  logical ,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector),save::mmm
  real(dp),dimension(1:nvector),save::vol2
  real(dp),dimension(1:nvector,1:ndim),save::x,dd,dg
  integer ,dimension(1:nvector,1:ndim),save::ig,id,igg,igd,icg,icd
  real(dp),dimension(1:nvector,1:twotondim),save::vol
  integer ,dimension(1:nvector,1:twotondim),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim

  ! Gather neighboring father cells (should be present anytime !)
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ng,ilevel)

  ! Rescale particle position at level ilevel
  do idim=1,ndim
     do j=1,np
        x(j,idim)=xp(ind_part(j),idim)/scale+skip_loc(idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)-x0(ind_grid_part(j),idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)/dx
     end do
  end do

  ! Gather particle mass
  do j=1,np
     if(ivar_clump==-1)then
        if(is_star(typep(ind_part(j))))then
           mmm(j)=mp(ind_part(j))
        else
           mmm(j)=0d0
        end if
     else if(ivar_clump==0)then
        if(is_dm(typep(ind_part(j))))then
           mmm(j)=mp(ind_part(j))
        else
           mmm(j)=0d0
        end if
     end if
  end do

  ! Check for illegal moves
  error=.false.
  do idim=1,ndim
     do j=1,np
        if(x(j,idim)<0.5D0.or.x(j,idim)>5.5D0)error=.true.
     end do
  end do
  if(error)then
     write(*,*)'problem in cic_only'
     do idim=1,ndim
        do j=1,np
           if(x(j,idim)<0.5D0.or.x(j,idim)>5.5D0)then
              write(*,*)x(j,1:ndim)
           endif
        end do
     end do
     stop
  end if

  ! CIC at level ilevel (dd: right cloud boundary; dg: left cloud boundary)
  do idim=1,ndim
     do j=1,np
        dd(j,idim)=x(j,idim)+0.5D0
        id(j,idim)=int(dd(j,idim))
        dd(j,idim)=dd(j,idim)-id(j,idim)
        dg(j,idim)=1.0D0-dd(j,idim)
        ig(j,idim)=id(j,idim)-1
     end do
  end do

  ! Compute cloud volumes
#if NDIM==1
  do j=1,np
     vol(j,1)=dg(j,1)
     vol(j,2)=dd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     vol(j,1)=dg(j,1)*dg(j,2)
     vol(j,2)=dd(j,1)*dg(j,2)
     vol(j,3)=dg(j,1)*dd(j,2)
     vol(j,4)=dd(j,1)*dd(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
     vol(j,1)=dg(j,1)*dg(j,2)*dg(j,3)
     vol(j,2)=dd(j,1)*dg(j,2)*dg(j,3)
     vol(j,3)=dg(j,1)*dd(j,2)*dg(j,3)
     vol(j,4)=dd(j,1)*dd(j,2)*dg(j,3)
     vol(j,5)=dg(j,1)*dg(j,2)*dd(j,3)
     vol(j,6)=dd(j,1)*dg(j,2)*dd(j,3)
     vol(j,7)=dg(j,1)*dd(j,2)*dd(j,3)
     vol(j,8)=dd(j,1)*dd(j,2)*dd(j,3)
  end do
#endif

  ! Compute parent grids
  do idim=1,ndim
     do j=1,np
        igg(j,idim)=ig(j,idim)/2
        igd(j,idim)=id(j,idim)/2
     end do
  end do
#if NDIM==1
  do j=1,np
     kg(j,1)=1+igg(j,1)
     kg(j,2)=1+igd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     kg(j,1)=1+igg(j,1)+3*igg(j,2)
     kg(j,2)=1+igd(j,1)+3*igg(j,2)
     kg(j,3)=1+igg(j,1)+3*igd(j,2)
     kg(j,4)=1+igd(j,1)+3*igd(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
     kg(j,1)=1+igg(j,1)+3*igg(j,2)+9*igg(j,3)
     kg(j,2)=1+igd(j,1)+3*igg(j,2)+9*igg(j,3)
     kg(j,3)=1+igg(j,1)+3*igd(j,2)+9*igg(j,3)
     kg(j,4)=1+igd(j,1)+3*igd(j,2)+9*igg(j,3)
     kg(j,5)=1+igg(j,1)+3*igg(j,2)+9*igd(j,3)
     kg(j,6)=1+igd(j,1)+3*igg(j,2)+9*igd(j,3)
     kg(j,7)=1+igg(j,1)+3*igd(j,2)+9*igd(j,3)
     kg(j,8)=1+igd(j,1)+3*igd(j,2)+9*igd(j,3)
  end do
#endif
  do ind=1,twotondim
     do j=1,np
        igrid(j,ind)=son(nbors_father_cells(ind_grid_part(j),kg(j,ind)))
     end do
  end do

  ! Compute parent cell position
  do idim=1,ndim
     do j=1,np
        icg(j,idim)=ig(j,idim)-2*igg(j,idim)
        icd(j,idim)=id(j,idim)-2*igd(j,idim)
     end do
  end do
#if NDIM==1
  do j=1,np
     icell(j,1)=1+icg(j,1)
     icell(j,2)=1+icd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     icell(j,1)=1+icg(j,1)+2*icg(j,2)
     icell(j,2)=1+icd(j,1)+2*icg(j,2)
     icell(j,3)=1+icg(j,1)+2*icd(j,2)
     icell(j,4)=1+icd(j,1)+2*icd(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
     icell(j,1)=1+icg(j,1)+2*icg(j,2)+4*icg(j,3)
     icell(j,2)=1+icd(j,1)+2*icg(j,2)+4*icg(j,3)
     icell(j,3)=1+icg(j,1)+2*icd(j,2)+4*icg(j,3)
     icell(j,4)=1+icd(j,1)+2*icd(j,2)+4*icg(j,3)
     icell(j,5)=1+icg(j,1)+2*icg(j,2)+4*icd(j,3)
     icell(j,6)=1+icd(j,1)+2*icg(j,2)+4*icd(j,3)
     icell(j,7)=1+icg(j,1)+2*icd(j,2)+4*icd(j,3)
     icell(j,8)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
  end do
#endif

  ! Compute parent cell adress
  do ind=1,twotondim
     do j=1,np
        indp(j,ind)=ncoarse+(icell(j,ind)-1)*ngridmax+igrid(j,ind)
     end do
  end do

  ! Update mass density field
  do ind=1,twotondim
     do j=1,np
        ok(j)=(igrid(j,ind)>0) .and. is_not_tracer(typep(ind_part(j)))
     end do
     do j=1,np
        vol2(j)=mmm(j)*vol(j,ind)/vol_loc
     end do
     do j=1,np
        if(ok(j))then
           rho(indp(j,ind))=rho(indp(j,ind))+vol2(j)
        end if
     end do
  end do

end subroutine cic_only
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine tsc_only(ind_cell,ind_part,ind_grid_part,x0,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use poisson_commons
  implicit none
  integer::ng,np,ilevel
  integer ,dimension(1:nvector)::ind_cell,ind_grid_part,ind_part
  real(dp),dimension(1:nvector,1:ndim)::x0
  !------------------------------------------------------------------
  ! This routine computes the density field at level ilevel using
  ! the TSC scheme. Only cells that are in level ilevel
  ! are updated by the input particle list.
  !------------------------------------------------------------------
  integer::j,ind,idim,nx_loc
  real(dp)::dx,dx_loc,scale,vol_loc
  ! Grid-based arrays
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle-based arrays
  logical ,dimension(1:nvector),save::ok,abandoned
  real(dp),dimension(1:nvector),save::mmm
  real(dp),dimension(1:nvector),save::vol2
  real(dp),dimension(1:nvector,1:ndim),save::x,cl,cr,cc,wl,wr,wc
  integer ,dimension(1:nvector,1:ndim),save::igl,igr,igc,icl,icr,icc
  real(dp),dimension(1:nvector,1:threetondim),save::vol
  integer ,dimension(1:nvector,1:threetondim),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc

  if (ndim .ne. 3)then
     write(*,*)'TSC not supported for ndim neq 3'
     call clean_stop
  end if

#if NDIM==3

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim

  ! Gather neighboring father cells (should be present at anytime!)
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ng,ilevel)

  ! Rescale particle position at level ilevel
  do idim=1,ndim
     do j=1,np
        x(j,idim)=xp(ind_part(j),idim)/scale+skip_loc(idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)-x0(ind_grid_part(j),idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)/dx
     end do
  end do

  ! Gather particle mass
  do j=1,np
     mmm(j)=mp(ind_part(j))
  end do

  ! Check for illegal moves
  abandoned(1:np)=.false.
  do idim=1,ndim
     do j=1,np
        if(x(j,idim)<1.0D0.or.x(j,idim)>5.0D0) abandoned(j)=.true.
     end do
  end do

  ! TSC at level ilevel; a particle contributes
  !     to three cells in each dimension
  ! cl: position of leftmost cell centre
  ! cc: position of central cell centre
  ! cr: position of rightmost cell centre
  ! wl: weighting function for leftmost cell
  ! wc: weighting function for central cell
  ! wr: weighting function for rightmost cell
  do idim=1,ndim
     do j=1,np
        if(.not.abandoned(j)) then
           cl(j,idim)=dble(int(x(j,idim)))-0.5D0
           cc(j,idim)=dble(int(x(j,idim)))+0.5D0
           cr(j,idim)=dble(int(x(j,idim)))+1.5D0
           wl(j,idim)=0.50D0*(1.5D0-abs(x(j,idim)-cl(j,idim)))**2
           wc(j,idim)=0.75D0-          (x(j,idim)-cc(j,idim)) **2
           wr(j,idim)=0.50D0*(1.5D0-abs(x(j,idim)-cr(j,idim)))**2
        end if
     end do
  end do

  ! Compute cloud volumes
  do j=1,np
     if(.not.abandoned(j)) then
        vol(j,1 )=wl(j,1)*wl(j,2)*wl(j,3)
        vol(j,2 )=wc(j,1)*wl(j,2)*wl(j,3)
        vol(j,3 )=wr(j,1)*wl(j,2)*wl(j,3)
        vol(j,4 )=wl(j,1)*wc(j,2)*wl(j,3)
        vol(j,5 )=wc(j,1)*wc(j,2)*wl(j,3)
        vol(j,6 )=wr(j,1)*wc(j,2)*wl(j,3)
        vol(j,7 )=wl(j,1)*wr(j,2)*wl(j,3)
        vol(j,8 )=wc(j,1)*wr(j,2)*wl(j,3)
        vol(j,9 )=wr(j,1)*wr(j,2)*wl(j,3)
        vol(j,10)=wl(j,1)*wl(j,2)*wc(j,3)
        vol(j,11)=wc(j,1)*wl(j,2)*wc(j,3)
        vol(j,12)=wr(j,1)*wl(j,2)*wc(j,3)
        vol(j,13)=wl(j,1)*wc(j,2)*wc(j,3)
        vol(j,14)=wc(j,1)*wc(j,2)*wc(j,3)
        vol(j,15)=wr(j,1)*wc(j,2)*wc(j,3)
        vol(j,16)=wl(j,1)*wr(j,2)*wc(j,3)
        vol(j,17)=wc(j,1)*wr(j,2)*wc(j,3)
        vol(j,18)=wr(j,1)*wr(j,2)*wc(j,3)
        vol(j,19)=wl(j,1)*wl(j,2)*wr(j,3)
        vol(j,20)=wc(j,1)*wl(j,2)*wr(j,3)
        vol(j,21)=wr(j,1)*wl(j,2)*wr(j,3)
        vol(j,22)=wl(j,1)*wc(j,2)*wr(j,3)
        vol(j,23)=wc(j,1)*wc(j,2)*wr(j,3)
        vol(j,24)=wr(j,1)*wc(j,2)*wr(j,3)
        vol(j,25)=wl(j,1)*wr(j,2)*wr(j,3)
        vol(j,26)=wc(j,1)*wr(j,2)*wr(j,3)
        vol(j,27)=wr(j,1)*wr(j,2)*wr(j,3)
     end if
  end do

  ! Compute parent grids
  do idim=1,ndim
     do j=1,np
        if(.not.abandoned(j)) then
           igl(j,idim)=(int(cl(j,idim)))/2
           igc(j,idim)=(int(cc(j,idim)))/2
           igr(j,idim)=(int(cr(j,idim)))/2
        end if
     end do
  end do
  do j=1,np
     if(.not.abandoned(j)) then
        kg(j,1 )=1+igl(j,1)+3*igl(j,2)+9*igl(j,3)
        kg(j,2 )=1+igc(j,1)+3*igl(j,2)+9*igl(j,3)
        kg(j,3 )=1+igr(j,1)+3*igl(j,2)+9*igl(j,3)
        kg(j,4 )=1+igl(j,1)+3*igc(j,2)+9*igl(j,3)
        kg(j,5 )=1+igc(j,1)+3*igc(j,2)+9*igl(j,3)
        kg(j,6 )=1+igr(j,1)+3*igc(j,2)+9*igl(j,3)
        kg(j,7 )=1+igl(j,1)+3*igr(j,2)+9*igl(j,3)
        kg(j,8 )=1+igc(j,1)+3*igr(j,2)+9*igl(j,3)
        kg(j,9 )=1+igr(j,1)+3*igr(j,2)+9*igl(j,3)
        kg(j,10)=1+igl(j,1)+3*igl(j,2)+9*igc(j,3)
        kg(j,11)=1+igc(j,1)+3*igl(j,2)+9*igc(j,3)
        kg(j,12)=1+igr(j,1)+3*igl(j,2)+9*igc(j,3)
        kg(j,13)=1+igl(j,1)+3*igc(j,2)+9*igc(j,3)
        kg(j,14)=1+igc(j,1)+3*igc(j,2)+9*igc(j,3)
        kg(j,15)=1+igr(j,1)+3*igc(j,2)+9*igc(j,3)
        kg(j,16)=1+igl(j,1)+3*igr(j,2)+9*igc(j,3)
        kg(j,17)=1+igc(j,1)+3*igr(j,2)+9*igc(j,3)
        kg(j,18)=1+igr(j,1)+3*igr(j,2)+9*igc(j,3)
        kg(j,19)=1+igl(j,1)+3*igl(j,2)+9*igr(j,3)
        kg(j,20)=1+igc(j,1)+3*igl(j,2)+9*igr(j,3)
        kg(j,21)=1+igr(j,1)+3*igl(j,2)+9*igr(j,3)
        kg(j,22)=1+igl(j,1)+3*igc(j,2)+9*igr(j,3)
        kg(j,23)=1+igc(j,1)+3*igc(j,2)+9*igr(j,3)
        kg(j,24)=1+igr(j,1)+3*igc(j,2)+9*igr(j,3)
        kg(j,25)=1+igl(j,1)+3*igr(j,2)+9*igr(j,3)
        kg(j,26)=1+igc(j,1)+3*igr(j,2)+9*igr(j,3)
        kg(j,27)=1+igr(j,1)+3*igr(j,2)+9*igr(j,3)
     end if
  end do
  do ind=1,threetondim
     do j=1,np
        igrid(j,ind)=son(nbors_father_cells(ind_grid_part(j),kg(j,ind)))
     end do
  end do

  ! Compute parent cell position
  do idim=1,ndim
     do j=1,np
        if(.not.abandoned(j)) then
           icl(j,idim)=int(cl(j,idim))-2*igl(j,idim)
           icc(j,idim)=int(cc(j,idim))-2*igc(j,idim)
           icr(j,idim)=int(cr(j,idim))-2*igr(j,idim)
        end if
     end do
  end do
  do j=1,np
     if(.not.abandoned(j)) then
        icell(j,1 )=1+icl(j,1)+2*icl(j,2)+4*icl(j,3)
        icell(j,2 )=1+icc(j,1)+2*icl(j,2)+4*icl(j,3)
        icell(j,3 )=1+icr(j,1)+2*icl(j,2)+4*icl(j,3)
        icell(j,4 )=1+icl(j,1)+2*icc(j,2)+4*icl(j,3)
        icell(j,5 )=1+icc(j,1)+2*icc(j,2)+4*icl(j,3)
        icell(j,6 )=1+icr(j,1)+2*icc(j,2)+4*icl(j,3)
        icell(j,7 )=1+icl(j,1)+2*icr(j,2)+4*icl(j,3)
        icell(j,8 )=1+icc(j,1)+2*icr(j,2)+4*icl(j,3)
        icell(j,9 )=1+icr(j,1)+2*icr(j,2)+4*icl(j,3)
        icell(j,10)=1+icl(j,1)+2*icl(j,2)+4*icc(j,3)
        icell(j,11)=1+icc(j,1)+2*icl(j,2)+4*icc(j,3)
        icell(j,12)=1+icr(j,1)+2*icl(j,2)+4*icc(j,3)
        icell(j,13)=1+icl(j,1)+2*icc(j,2)+4*icc(j,3)
        icell(j,14)=1+icc(j,1)+2*icc(j,2)+4*icc(j,3)
        icell(j,15)=1+icr(j,1)+2*icc(j,2)+4*icc(j,3)
        icell(j,16)=1+icl(j,1)+2*icr(j,2)+4*icc(j,3)
        icell(j,17)=1+icc(j,1)+2*icr(j,2)+4*icc(j,3)
        icell(j,18)=1+icr(j,1)+2*icr(j,2)+4*icc(j,3)
        icell(j,19)=1+icl(j,1)+2*icl(j,2)+4*icr(j,3)
        icell(j,20)=1+icc(j,1)+2*icl(j,2)+4*icr(j,3)
        icell(j,21)=1+icr(j,1)+2*icl(j,2)+4*icr(j,3)
        icell(j,22)=1+icl(j,1)+2*icc(j,2)+4*icr(j,3)
        icell(j,23)=1+icc(j,1)+2*icc(j,2)+4*icr(j,3)
        icell(j,24)=1+icr(j,1)+2*icc(j,2)+4*icr(j,3)
        icell(j,25)=1+icl(j,1)+2*icr(j,2)+4*icr(j,3)
        icell(j,26)=1+icc(j,1)+2*icr(j,2)+4*icr(j,3)
        icell(j,27)=1+icr(j,1)+2*icr(j,2)+4*icr(j,3)
     end if
  end do

  ! Compute parent cell adress
  do ind=1,threetondim
     do j=1,np
        if(.not.abandoned(j)) then
           indp(j,ind)=ncoarse+(icell(j,ind)-1)*ngridmax+igrid(j,ind)
        end if
     end do
  end do

  ! Update mass density and number density fields
  do ind=1,threetondim

     do j=1,np
        if(.not.abandoned(j)) then
           ok(j)=(igrid(j,ind)>0) .and. is_not_tracer(typep(ind_part(j)))
        end if
     end do

     do j=1,np
        if(.not.abandoned(j)) then
           vol2(j)=mmm(j)*vol(j,ind)/vol_loc
        end if
     end do

     do j=1,np
        if(ok(j).and.(.not.abandoned(j))) then
           rho(indp(j,ind))=rho(indp(j,ind))+vol2(j)
        end if
     end do

  end do
#endif
end subroutine tsc_only
!###########################################################
!###########################################################
!###########################################################
!###########################################################
!######################################
!######################################
!######################################
subroutine output_part_clump_id()
  !---------------------------------------------------------------------------
  ! This subroutine loops over all test cells and assigns all particles in a
  ! testcell the peak ID the testcell has.
  !---------------------------------------------------------------------------
  use amr_commons
  use clfind_commons    ! unbinding stuff is all in here
  use pm_commons ! using mp
  use amr_parameters
  implicit none
  integer,dimension(:),allocatable::clump_ids
  character(len=80) :: fileloc
  character(len=5)  :: nchar,nchar2

  ! for looping over test cells and getting particle list
  integer   :: itestcell, ipart,this_part, global_peak_id, prtcls_in_grid

  ! getting particles per peak
  integer   :: ind, grid

  !getting in which cell of a grid a particle is
  integer   :: part_cell_ind, i, j, k

  if(verbose) write(*,*) "Entered get_clumpparticles"

  !-----------------------------------------------------------

  allocate(clmpidp(1:npartmax))
  allocate(clump_ids(1:npart))

  clmpidp=0

  do itestcell=1, ntest !loop over all test cells
     global_peak_id=flag2(icellp(itestcell))

     if (global_peak_id /= 0) then

        ind=(icellp(itestcell)-ncoarse-1)/ngridmax+1  ! get cell position
        grid=icellp(itestcell)-ncoarse-(ind-1)*ngridmax ! get grid index
        prtcls_in_grid = numbp(grid)          ! get number of particles in grid
        this_part=headp(grid)             ! get index of first particle

        ! loop over particles in grid
        do ipart = 1, prtcls_in_grid

            !check cell index of particle so you loop only once over each
            i=0
            j=0
            k=0
            if(xg(grid,1)-xp(this_part,1)/boxlen+(nx-1)/2.0 .le. 0) i=1
            if(xg(grid,2)-xp(this_part,2)/boxlen+(ny-1)/2.0 .le. 0) j=1
            if(xg(grid,3)-xp(this_part,3)/boxlen+(nz-1)/2.0 .le. 0) k=1

            part_cell_ind=i+2*j+4*k+1

            !If index is correct, assign clump id to particle
            if (part_cell_ind==ind) clmpidp(this_part)=global_peak_id
            !go to next particle in this grid
            this_part = nextp(this_part)

        end do
     end if   !global peak /=0
  end do   !loop over test cells

  call title(ifout, nchar)
  call title(myid, nchar2)
  fileloc=TRIM('output_'//TRIM(nchar)//'/id_clump.out'//TRIM(nchar2))
  open(unit=666,file=fileloc,form='unformatted')
  ipart=0
  do i=1,npartmax
    if(levelp(i)>0)then
      ipart=ipart+1
      clump_ids(ipart)=clmpidp(i)
    end if
  end do
  write(666) clump_ids
  close(666)

  deallocate(clmpidp)
  deallocate(clump_ids)

end subroutine output_part_clump_id
!########################################
!########################################
!########################################

#endif
