!! Patch written by Mladen Ivkovic (mladen.ivkovic@uzh.ch)
!! The main routine (unbinding) is called in the clump_finder routine
!! from the clump_finder.f90 after clump properties are written to file.
!! Changes in the clump_finder routines are marked by a "added for patch"
!! comment so they're easier to find and track.

!! New subroutines for this patch are:
!!
!! subroutine unbinding()
!! subroutine get_clumpparticles()
!! subroutine get_clump_properties_pb()
!! subroutine get_cmp
!! subroutine get_closest_border
!! subroutine unbinding_neighborsearch
!! subroutine bordercheck
!! subroutine particle_unbinding()
!!      contains function unbound
!!      contains function potential
!! subroutine compute_phi
!! subroutine allocate_unbinding_arrays()
!! subroutine deallocate_unbinding_arrays()
!! subroutine unbinding_write_formatted_output()
!! subroutine unbinding_formatted_particleoutput()



!! New namelist parameters for this pach:
!! (Can be set in the CLUMPFIND_PARAMS block)
!!
!! NAME                        DEFAULT VALUE        FUNCTION
!! unbind=                     .true.               Turn particle unbinding on 
!!                                                  or off
!!
!! nmassbins=                  50                   Number of bins for the mass 
!!                                                  binning of the cumulative
!!                                                  mass profile. Any integer >1.
!!
!! logbins=                    .true.               use logarithmic binning 
!!                                                  distances for cumulative mass
!!                                                  profiles (and gravitational 
!!                                                  potential of clumps).
!!                                                  If false, the code  will use 
!!                                                  linear binning distances.
!!
!! saddle_pot=                 .true.               Take neighbouring structures 
!!                                                  into account; Cut potentiall
!!                                                  off at closest saddle.
!!
!! unbinding_formatted_output= .false.              Create formatted output for 
!!                                                  particles, cumulative mass
!!                                                  profiles, binning distances, 
!!                                                  particle based clump
!!                                                  properties, gravitational 
!!                                                  potential of substructure
!!                                                  clumps 
!!
!! iter_properties=            .true.               whether to unbind multiple times 
!!                                                  with updated clump properties
!!                                                  determined by earlier unbindings
!!
!! conv_limit =                0.01                 convergence limit. If the 
!!                                                  v_clump_old/v_clump_new < conv_limit,
!!                                                  stop iterating for this clump. 
!!                                                  (only used when iter_properties=.true.)
!!
!! repeat_max =                100                  maximal number of loops per level
!!                                                  for iterative unbinding
!!                                                  (in case a clump doesn't converge)
!!                                                  (shouldn't happen)
!!                                                  (only used when iter_properties=.true.)









!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$                               $$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$   CLUMP FINDER STARTS HERE    $$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$                               $$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#if NDIM==3
subroutine clump_finder(create_output,keep_alive)
  use amr_commons
  use poisson_commons, ONLY:rho
  use clfind_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
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

#ifndef WITHOUTMPI
  integer(i8b)::nmove_all,nzero_all
#endif

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
  if(ivar_clump==0)then
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
     if(ivar_clump==0)then ! action 1: count and flag
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
     denp=0.d0; levp=0; imaxp=0; icellp=0
  endif
  itest=0
  nskip=ntest_cpu(myid)-ntest
  do ilevel=levelmin,nlevelmax
     if(ivar_clump==0)then
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
     if(ivar_clump==0)then  ! case 1: count peaks
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
  max_dens=0.d0; peak_cell=0; peak_cell_level=0;
  flag2=0
  if(ntest>0)then
     if(ivar_clump==0)then
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
     if(ivar_clump==0)then
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
     if(ivar_clump==0)then
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
     if(verbose)call analyze_peak_memory
     if(clinfo.and.saddle_threshold.LE.0)call write_clump_properties(.false.)
     if(create_output)then
        if(myid==1)write(*,*)"Outputing clump properties to disc."
        call write_clump_properties(.true.)
     endif

  end if




    !------------------------------------------
    ! Added for patch:
    ! Call particle unbinding
    !------------------------------------------
    
    if(unbind) call unbinding()





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
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
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
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
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
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
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
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
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
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
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
     density_max(j)=xx(ind_cell(j))*1.0001 ! get cell density (1.0001 probably not necessary)
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
              xrel(ntp,1)=(i3-1.5)*dx_loc/2.0
              xrel(ntp,2)=(j3-1.5)*dx_loc/2.0
              xrel(ntp,3)=(k3-1.5)*dx_loc/2.0
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
     call get_cell_index_fast(cell_index,cell_levl,xtest,ind_grid(j),nbors_father_cells(j,1:threetondim),ntestpos,ilevel)

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
     av_dens(j)=(xx(cell_index(j))+xx(ind_cell))*0.5 !average density of cell and neighbor cell
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
     xx = xpart(i,1)/boxlen + (nx-1)/2.0
     yy = xpart(i,2)/boxlen + (ny-1)/2.0
     zz = xpart(i,3)/boxlen + (nz-1)/2.0

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
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

! added for patch: unbinding parameters
  namelist/clumpfind_params/ivar_clump,& 
       & relevance_threshold,density_threshold,&
       & saddle_threshold,mass_threshold,clinfo,&
       & n_clfind,rho_clfind,&
       !unbinding parameters
       & unbind,nmassbins,logbins,unbinding_formatted_output, &
       & saddle_pot,iter_properties,conv_limit, repeat_max
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
           density_threshold=d_sink/10.
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
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  one_over_dx=1./dx
  one_over_scale=1./scale

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
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
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
           rho(reception(icpu,ilevel)%igrid(i)+iskip)=0.0D0
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
           rho(boundary(ibound,ilevel)%igrid(i)+iskip)=0.0
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
              if(age_cut_clfind>0.d0 .and. star) then
                 if((t-tp(ipart).lt.age_cut_clfind).and.(tp(ipart).ne.0.d0)) then
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
              if(age_cut_clfind>0.d0 .and. star) then
                 if((t-tp(ipart).lt.age_cut_clfind).and.(tp(ipart).ne.0.d0)) then
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
     mmm(j)=mp(ind_part(j))
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
        ok(j)=igrid(j,ind)>0
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
           ok(j)=igrid(j,ind)>0
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
!################                                ###########
!################    PARTICLE UNBINDING PATCH    ###########
!################    ROUTINES START HERE         ###########
!################                                ###########
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine unbinding()

  use amr_commons  !MPI stuff
  use pm_commons
  use clfind_commons !unbinding stuff

  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  !------------------------------------------------------------
  ! This subroutine assigns all particles that are in cells
  ! identified to be in a clump by the clump finder the peak
  ! ID of that clump and checks whether they are energetically
  ! bound to that structure. If not, they are passed on to the
  ! clump's parents.
  !------------------------------------------------------------

  integer           :: ipeak, info, ilevel, ipart, i
  integer           :: loop_counter=0
  integer, dimension(1:npart) :: clump_ids
  character(LEN=80)     :: fileloc
  character(LEN=5)      :: nchar,nchar2
  logical           :: loop_again_global, is_final_round, is_first



  !Logging/Announcing stuff
  if(myid==1) write(*,*) "Started unbinding."


  !------------------------
  ! Initial set-up
  !------------------------

  !update boundary relevance
  call build_peak_communicator
  call boundary_peak_dp(relevance)
  !call boundary_peak_int(new_peak) !already done in clump_finder, doesn't 
                    !need an update
  call boundary_peak_int(lev_peak)

  
  ! set up constants/
  GravConst=1d0   ! Gravitational constant
  if(cosmo)GravConst=3d0/8d0/3.1415926*omega_m*aexp

  periodical=(nx==1) !true if periodic
  nunbound=0     !count unbound particles
  candidates=0   !count unbinding candidates: Particles of child clumps
           !that will be tested. If a particle is passed on to the
           !parent clump, it will be counted twice.

  ! allocate necessary arrays
  call allocate_unbinding_arrays()



  !-------------------
  ! Gather particles 
  !-------------------

  !Get particles in substructrue, create linked lists
  call get_clumpparticles()

  !write output if required. This gives the particle clump assignment as given
  !by PHEW.
  if (unbinding_formatted_output) call unbinding_formatted_particleoutput(.true.) 



  !------------------
  ! Unbinding loop
  !------------------

  ! go level by level
  do ilevel=0, mergelevel_max
    !Start iteration loop 

    is_final_round=.true.
    if (iter_properties) is_final_round=.false.
    
    loop_again=.true.
    loop_counter=0

    !reset values
    to_iter = .true.
    hasatleastoneptcl=1


    do while(loop_again)
      loop_again = .false.  !set boolean whether to continue to false as default;
                  !will get reset if conditions are met

      loop_counter=loop_counter+1
      niterunbound=0

      !get particle based clump properties :
      !Center of Mass, bulk velocity, particle furthest away from CoM
      is_first = (loop_counter==1)
      call get_clump_properties_pb(ilevel,is_first)

      if (loop_counter==repeat_max) loop_again=.false.

#ifndef WITHOUTMPI
      !sync with other processors wheter you need to repeat loop
      call MPI_ALLREDUCE(loop_again, loop_again_global, 1, MPI_LOGICAL, MPI_LOR,MPI_COMM_WORLD, info)

      loop_again=loop_again_global
#endif

      ! get cumulative mass profiles
      call get_cmp(ilevel)

      ! get closest border to the center of mass
      if (saddle_pot) call get_closest_border(ilevel)

      !UNBINDING LOOP
      if (.not.loop_again) is_final_round=.true.
      do ipeak=1, hfree-1
        if (cmp_distances(ipeak,nmassbins)>0.0 .and. lev_peak(ipeak)==ilevel) then 
          ! use last cmp_distances as criterion because 
          ! that value is communicated across processing units
          call particle_unbinding(ipeak,is_final_round)
        end if
      end do

      if (loop_again) then
        !communicate whether peaks have remaining contributing particles
        call build_peak_communicator()
        call virtual_peak_int(hasatleastoneptcl,'max')
        call boundary_peak_int(hasatleastoneptcl)
        do ipeak=1,hfree-1
          !if peak has no contributing particles anymore
          if (hasatleastoneptcl(ipeak)==0) to_iter(ipeak)=.false. !dont loop anymore over this peak
        end do
      end if


      if (clinfo .and. iter_properties) then 
        !---------------------
        ! Talk to me.
        !---------------------

#ifndef WITHOUTMPI
        call MPI_ALLREDUCE(niterunbound, niterunbound_tot, 1, MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD, info)
#else
        niterunbound_tot=niterunbound
#endif
        if (myid==1) then
          write(*,'(A10,I10,A30,I5,A7,I5)') " Unbound", niterunbound_tot, "particles at level", ilevel, "loop", loop_counter
        end if
      end if


      if (.not. loop_again .and. clinfo .and. myid==1 .and. iter_properties .and. loop_counter < repeat_max) then
        write(*, '(A7,I5,A35,I5,A12)') "Level ", ilevel, "clump properties converged after ", loop_counter, "iterations."
      end if

      if (loop_counter==repeat_max) write(*,'(A7,I5,A20,I5,A35)') "Level ", ilevel, "not converged after ", repeat_max, "iterations. Moving on to next step."

    end do ! loop again for ilevel

  end do !loop over levels
  






  if (clinfo) then
#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(nunbound, nunbound_tot, 1, MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD, info)
    call MPI_ALLREDUCE(candidates, candidates_tot, 1, MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD, info)
#else
    nunbound_tot=nunbound
    candidates_tot=candidates
#endif
    if (myid==1) then
      write(*,'(A6,I10,A30,I10,A12)') " Found", nunbound_tot, "unbound particles out of ", candidates_tot, " candidates"
    end if
  end if
   


  !-----------------
  ! Write output
  !-----------------

  if(unbinding_formatted_output) call unbinding_write_formatted_output()
  
  call title(ifout-1, nchar)
  call title(myid, nchar2)
  fileloc=TRIM('output_'//TRIM(nchar)//'/unbinding.out'//TRIM(nchar2))

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


  !--------------------
  ! Deallocate arrays
  !--------------------
  call deallocate_unbinding_arrays()

  !--------------------
  ! Say good bye.
  !--------------------
  if(verbose.or.myid==1) write(*,*) "Finished unbinding."

end subroutine unbinding
!######################################
!######################################
!######################################
subroutine get_clumpparticles()

  !---------------------------------------------------------------------------
  ! This subroutine loops over all test cells and assigns all particles in a
  ! testcell the peak ID the testcell has. If the peak is not a namegiver
  ! (= not its own parent), the particles are added to a linked list for 
  ! unbinding later.
  !---------------------------------------------------------------------------

  use amr_commons
  use clfind_commons    !unbinding stuff is all in here
  use pm_commons !using mp
  use amr_parameters
  use hydro_commons !using mass_sph
  implicit none
  
#ifndef WITHOUTMPI
  integer :: info
  include 'mpif.h'
#endif
  ! for looping over test cells and getting particle list
  integer   :: itestcell, ipart,this_part, global_peak_id, local_peak_id, prtcls_in_grid 
  
  ! getting particles per peak
  integer   :: ind, grid

  !getting particle mass
  real(dp)  :: particle_mass, particle_mass_tot

  !getting in which cell of a grid a particle is
  integer   :: part_cell_ind,i,j,k

  !appending linked lists
  integer   :: ipeak, new_peak_local_id, ilevel


  if(verbose) write(*,*) "Entered get_clumpparticles"

  !get particle mass (copied from subroutine write_clump_properties)
  if(ivar_clump==0)then
    particle_mass=MINVAL(mp, MASK=(mp.GT.0.))
#ifndef WITHOUTMPI  
    call MPI_ALLREDUCE(particle_mass,particle_mass_tot,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,info)
    particle_mass=particle_mass_tot  
#endif
  else
    if(hydro)then
      particle_mass=mass_sph
    endif
  endif

  !if(myid==1) write(*,*) "--- Particle mass:", particle_mass


  !-----------------------------------------------------------
  ! Get particles from testcells into linked lists for clumps
  !-----------------------------------------------------------

  do itestcell=1, ntest !loop over all test cells
    global_peak_id=flag2(icellp(itestcell))

    if (global_peak_id /= 0) then

      ! get local peak id
      call get_local_peak_id(global_peak_id, local_peak_id)

      ! If peak ID is also a halo, there is no need for particle unbinding,
      ! so there is no need for a linked list of particles.
      ! This check is also important for the relevance conditions:
      ! There are clumps that satisfy the halo mass condition, but do not
      ! satisfy the clump relevance condition. If not separated, the 
      ! halos who are not relevant clumps will not be recognised.

      ! Check if peak is halo:      
      if (new_peak(local_peak_id)==global_peak_id) then
      
        !------------------
        !Found a halo cell
        !------------------

        ! if peak relevant:
        if(halo_mass(local_peak_id) > mass_threshold*particle_mass) then
        
          ind=(icellp(itestcell)-ncoarse-1)/ngridmax+1  ! get cell position
          grid=icellp(itestcell)-ncoarse-(ind-1)*ngridmax ! get grid index
          prtcls_in_grid = numbp(grid)          ! get number of particles in grid
          this_part=headp(grid)             ! get index of first particle
        
          ! If it is a halo: only assign particles the ID
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
        end if

      else != if not halo ID

        !---------------------
        !Found a subhalo cell
        !---------------------

        if (relevance(local_peak_id) > relevance_threshold) then
          ! If clump isn't also halo: assign ID to particles,
          ! create linked particle list
          
          ind=(icellp(itestcell)-ncoarse-1)/ngridmax+1  ! get cell position
          grid=icellp(itestcell)-ncoarse-(ind-1)*ngridmax ! get grid index
          prtcls_in_grid = numbp(grid)          ! get number of particles in grid
          this_part=headp(grid)             ! get index of first particle

           
          !loop over particles in grid
            do ipart=1, prtcls_in_grid
              !check cell index of particle so you loop only once over each
              i=0
              j=0
              k=0
              if(xg(grid,1)-xp(this_part,1)/boxlen+(nx-1)/2.0 .le. 0) i=1
              if(xg(grid,2)-xp(this_part,2)/boxlen+(ny-1)/2.0 .le. 0) j=1
              if(xg(grid,3)-xp(this_part,3)/boxlen+(nz-1)/2.0 .le. 0) k=1
  
              part_cell_ind=i+2*j+4*k+1
      
              ! If index is correct, assign clump id to particle
              if (part_cell_ind==ind) then
                !assign peak ID
                clmpidp(this_part)=global_peak_id
                !add particle to linked list of clumpparticles 
                !check if already particles are assigned
                if (nclmppart(local_peak_id)>0) then
                  !append to the last particle of the list
                  clmppart_next(clmppart_last(local_peak_id))=this_part
                else
                  !assign particle as first particle
                  !for this peak of linked list 
                  clmppart_first(local_peak_id)=this_part
                end if
                !update last particle for this clump
                nclmppart(local_peak_id)=nclmppart(local_peak_id)+1
                clmppart_last(local_peak_id)=this_part
              end if    
              !go to next particle in this grid
              this_part=nextp(this_part)
            end do
        end if !if clump is relevant
      end if     !if clump or halo
    end if       !global peak /=0
  end do         !loop over test cells




  !------------------------------------------------------
  !Append substructure particles to parents' linked list
  !------------------------------------------------------
  
  ! must be done level by level!
  do ilevel=0,mergelevel_max
    do ipeak=1, hfree-1
      ! append substructure linked lists to parent linked lists
      if(lev_peak(ipeak)==ilevel) then
        if (nclmppart(ipeak)>0) then
          ! get local id of parent
          call get_local_peak_id(new_peak(ipeak),new_peak_local_id)
          ! append particle LL to parent's LL if parent isn't a halo-namegiver
          if(new_peak(ipeak)/=new_peak(new_peak_local_id)) then !if peak isnt namegiver
            ! It might happen that the parent peak doesn't have a 
            ! particle linked list yet (on this processor).
            if (nclmppart(new_peak_local_id)>0) then !particle ll exists
              clmppart_next(clmppart_last(new_peak_local_id))=clmppart_first(ipeak)
            else
              clmppart_first(new_peak_local_id)=clmppart_first(ipeak)
            end if

            clmppart_last(new_peak_local_id)=clmppart_last(ipeak)
            nclmppart(new_peak_local_id)=nclmppart(new_peak_local_id)+nclmppart(ipeak)
          end if
        end if
      end if
    end do
  end do

end subroutine get_clumpparticles
!########################################
!########################################
!########################################
subroutine get_clump_properties_pb(ilevel,first)
  use amr_commons
  use pm_commons
  use clfind_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  logical, intent(in) :: first  !if it is the first time calculating
  integer, intent(in) :: ilevel !clump level currently on
   
  !--------------------------------------------------------------------------
  ! This subroutine computes the particle-based properties of the clumps:
  ! namely the center of mass and the clump's velocity.
  ! If it's called for the first time, it will compute the properties for
  ! all peak IDs. If not, it will go level by level.
  !--------------------------------------------------------------------------

  !particle furthest away
  real(dp) :: distance, biggest

  ! iterators
  integer :: ipeak, i, ipart
  integer :: thispart
  real(dp) :: vsq
  real(dp),dimension(1:3) :: period
  real(dp),dimension(1:npeaks_max) :: clmp_vel_pb_old
  logical :: check

  if (verbose) write(*,*)"Entered get_clump_properties (particle based)"


  !------------------------------------------------------------
  ! If iterative: Store old values, reset virtual peak values
  !------------------------------------------------------------

  do ipeak=1, hfree-1
    check = .not. first
    check = check .and. lev_peak(ipeak) == ilevel
    check = check .and. to_iter(ipeak)
    if (check) then
      clmp_vel_pb_old(ipeak)=clmp_vel_pb(ipeak,1)**2+clmp_vel_pb(ipeak,2)**2+clmp_vel_pb(ipeak,3)**2
      do i=1,3
        oldcom(ipeak,i)=clmp_com_pb(ipeak,i)
        oldvel(ipeak,i)=clmp_vel_pb(ipeak,i)
      end do
      oldcmpd(ipeak)=cmp_distances(ipeak,nmassbins)
      oldm(ipeak)=clmp_mass_pb(ipeak)
    end if

    if (iter_properties .and. ipeak>npeaks) then  
      !for communication: set virtual peak values=0
      !so they won't contribute in the communication sum
      !reset values
      do i=1,3
        clmp_vel_pb(ipeak,i)=0.0
        clmp_com_pb(ipeak,i)=0.0
      end do
      clmp_mass_pb(ipeak)=0.0
    end if
  end do  


  !------------------------------------------------------
  ! GET CENTER OF MASS, CENTER OF MOMENTUM FRAME VELOCITY
  !------------------------------------------------------

  do ipeak=1, hfree-1 !loop over all peaks
    if (to_iter(ipeak).and.lev_peak(ipeak)==ilevel) then !if peak has particles and needs to be iterated over

      !reset values
      do i=1,3
        clmp_vel_pb(ipeak,i)=0.0
        clmp_com_pb(ipeak,i)=0.0
      end do
      cmp_distances(ipeak,nmassbins)=0.0
      clmp_mass_pb(ipeak)=0.0


      if (hasatleastoneptcl(ipeak)>0 .and. nclmppart(ipeak)>0) then 
        !if there is work to do on this processing unit for this peak

        thispart=clmppart_first(ipeak)
        
        do ipart=1, nclmppart(ipeak)       !while there is a particle linked list
          if (contributes(thispart)) then    !if the particle should be considered
            if (periodical) then !determine periodic correction
              period=0.d0
              do i=1, 3
                if (xp(thispart,i)-peak_pos(ipeak,i) > 0.5*boxlen)    period(i)=(-1.0)*boxlen
                if (xp(thispart,i)-peak_pos(ipeak,i) < (-0.5*boxlen)) period(i)=boxlen
              end do
            end if

            clmp_mass_pb(ipeak)=clmp_mass_pb(ipeak)+mp(thispart)
            do i=1,3
              clmp_com_pb(ipeak,i)=clmp_com_pb(ipeak,i)+(xp(thispart,i)+period(i))*mp(thispart) !get center of mass sum       
              clmp_vel_pb(ipeak,i)=clmp_vel_pb(ipeak,i)+vp(thispart,i)*mp(thispart) !get velocity sum
            end do
          else
            contributes(thispart)=.true. !reset value
          end if
          thispart=clmppart_next(thispart) !go to next particle in linked list
        end do ! loop over particles
      end if     ! there is work for this peak on this processor
    end if       ! peak needs to be looked at
  end do         ! loop over peaks


  !----------------------------------------------------------------------
  ! communicate center of mass, clump mass and velocity across processors
  !----------------------------------------------------------------------
  call build_peak_communicator
  do i=1,3
    call virtual_peak_dp(clmp_com_pb(1,i),'sum')  !collect
    call boundary_peak_dp(clmp_com_pb(1,i))     !scatter
    call virtual_peak_dp(clmp_vel_pb(1,i),'sum')  !collect
    call boundary_peak_dp(clmp_vel_pb(1,i))     !scatter
  end do
  call virtual_peak_dp(clmp_mass_pb,'sum')      !collect
  call boundary_peak_dp(clmp_mass_pb)         !scatter


  do ipeak=1, hfree-1

    check = to_iter(ipeak)
    check = check .and. lev_peak(ipeak)==ilevel
    check = check .and. clmp_mass_pb(ipeak)>0

    if (check) then
      !calculate actual CoM and center of momentum frame velocity
      do i=1,3
        clmp_com_pb(ipeak,i)=clmp_com_pb(ipeak,i)/clmp_mass_pb(ipeak)
        clmp_vel_pb(ipeak,i)=clmp_vel_pb(ipeak,i)/clmp_mass_pb(ipeak)
      end do

      !------------------------------------
      !FIND PARTICLE FURTHEST AWAY FROM CoM
      !------------------------------------
      ! The maximal distance of a particle to the CoM is saved in the last
      ! cmp_distances array for every peak.
      if(nclmppart(ipeak)>0) then
        biggest=0.0
        thispart=clmppart_first(ipeak)
        do ipart=1, nclmppart(ipeak) !while there is a particle linked list

          period=0.d0
          
          if (periodical) then
            do i=1, 3
              if (xp(thispart,i)-peak_pos(ipeak,i)>0.5*boxlen) period(i)=(-1.0)*boxlen
              if (xp(thispart,i)-peak_pos(ipeak,i)<(-0.5*boxlen)) period(i)=boxlen
            end do
          end if

          distance=(xp(thispart,1)+period(1)-clmp_com_pb(ipeak,1))**2 + &
            (xp(thispart,2)+period(2)-clmp_com_pb(ipeak,2))**2 + &
            (xp(thispart,3)+period(3)-clmp_com_pb(ipeak,3))**2

          if(distance>biggest) biggest=distance ! save if it is biggest so far

          thispart=clmppart_next(thispart)

        end do
        if (biggest>0.0) cmp_distances(ipeak,nmassbins)=sqrt(biggest) !write if you have a result
      end if !to iterate
    end if
  end do ! over all peaks

  !-------------------------------------------------
  !communicate distance of particle furthest away
  !-------------------------------------------------
  call build_peak_communicator
  call virtual_peak_dp(cmp_distances(1,nmassbins), 'max')
  call boundary_peak_dp(cmp_distances(1,nmassbins))



  !-------------------------------------------------
  ! If iterative clump properties determination:
  ! Check whether bulk velocity converged
  !-------------------------------------------------

  if (iter_properties) then !if clump properties will be determined iteratively
    do ipeak=1, hfree-1

      check = to_iter(ipeak)
      check = check .and. lev_peak(ipeak)==ilevel
      check = check .and. cmp_distances(ipeak,nmassbins)>0.0

      if (check) then
        vsq=clmp_vel_pb(ipeak,1)**2+clmp_vel_pb(ipeak,2)**2+clmp_vel_pb(ipeak,3)**2

        if ( abs( sqrt(clmp_vel_pb_old(ipeak)/vsq) - 1.0) < conv_limit ) then
          to_iter(ipeak) = .false. ! consider bulk velocity as converged
          !write(*,'(A8,I3,A15,I8,A6,E15.6E2,A5,E15.6E2,A7,E15.6E2,A9,E15.6E2)') &
          !& "#####ID", myid, "clump CONVERGED", ipeak+ipeak_start(myid), "old:", &
          !& clmp_vel_pb_old(ipeak), "new:", vsq, "ratio",  abs( sqrt(clmp_vel_pb_old(ipeak)/vsq) - 1.0),&
          !& "v_bulk=", sqrt(vsq)
        else
          loop_again=.true. !repeat
        end if
      end if
    end do
  end if


end subroutine get_clump_properties_pb 
!###################################
!###################################
!###################################
subroutine get_cmp(ilevel)

  use amr_commons
  use pm_commons
  use clfind_commons
  implicit none
  integer, intent(in) :: ilevel 

  !----------------------------
  !Get cumulative mass profiles
  !----------------------------

  integer  :: ipeak, i, ipart, levelmax
  real(dp) :: r_null, distance
  integer  :: thispart
  real(dp),dimension(1:3) :: period
  logical  :: check
   
#ifndef WITHOUTMPI
  integer  :: levelmax_glob, info
  include 'mpif.h'
#endif

  if(verbose) write(*,*) "Entered get cumulative mass profiles"

  if (logbins) then
    !get minimal distance:
    levelmax=0
    do i=1,nlevelmax
       if(numbtot(1,i)>0) levelmax=levelmax+1
    end do

#ifndef WITHOUTMPI
    ! get system-wide levelmax
    call MPI_ALLREDUCE(levelmax,levelmax_glob, 1, MPI_INTEGER, MPI_MAX,MPI_COMM_WORLD, info)

    levelmax=levelmax_glob
#endif

    rmin=boxlen/2**levelmax
  end if




  do ipeak=1, hfree-1

    !peak must have need to be reiterated
    check=to_iter(ipeak)
    check=check.and.lev_peak(ipeak)==ilevel        !peak must have correct level
    check=check.and.nclmppart(ipeak)>0           !peak must have particles on this processor      

    !reset values
    if (check .or. ipeak > npeaks) then
      do i = 1, nmassbins
        cmp(ipeak,i) = 0.0
      end do
    end if


    if (check) then
      !------------------------------------------
      !Compute cumulative mass binning distances
      !------------------------------------------
      !The distances are not communicated later, but computed on each
      !processor independently, because each processor has all information it needs 
      !with cmp_distances(ipeak,nmassbins) and CoM

      if (logbins) then
        do i=1, nmassbins-1
          cmp_distances(ipeak,i)=rmin*(cmp_distances(ipeak,nmassbins)/rmin)**(real(i)/real(nmassbins))
        end do
      else !linear binnings
        r_null=cmp_distances(ipeak,nmassbins)/real(nmassbins)
        do i=0, nmassbins-1
          cmp_distances(ipeak,i)=r_null*i
        end do
      end if
      !The last bin must end with precicely with the maximal 
      !Distance of the particle. That is
      !needed because precision errors. The maximal distance is 
      !computed via particle data and the outermost bin is equal
      !to the distance of the outermost particle to the CoM.
      !Precision errors cause the code to crash here.

      !---------------------------------------------
      ! bin particles in cumulative mass profiles:
      ! get mass of each bin
      ! calculate particle distance to CoM
      !---------------------------------------------
      thispart=clmppart_first(ipeak)
      do ipart=1, nclmppart(ipeak)!while there is a particle linked list
        period=0.d0
        if (periodical) then
          do i=1, 3
            if (xp(thispart,i)-peak_pos(ipeak,i)>0.5*boxlen)  period(i)=(-1.0)*boxlen
            if (xp(thispart,i)-peak_pos(ipeak,i)<(-0.5*boxlen)) period(i)=boxlen
          end do
        end if
        distance=(xp(thispart,1)+period(1)-clmp_com_pb(ipeak,1))**2 + &
          (xp(thispart,2)+period(2)-clmp_com_pb(ipeak,2))**2 + &
          (xp(thispart,3)+period(3)-clmp_com_pb(ipeak,3))**2
        distance=sqrt(distance)


        i=1
        do 
          if (distance<=cmp_distances(ipeak,i)) then
            cmp(ipeak,i) = cmp(ipeak,i) + mp(thispart)
            exit
          else
            i=i+1
          end if
        end do

        thispart=clmppart_next(thispart)
      end do

      ! sum up masses to get profile instead of mass in shell
      do i=0,nmassbins-1
        cmp(ipeak,i+1)=cmp(ipeak,i+1)+cmp(ipeak,i) 
      end do

    end if  ! check
  end do    ! loop over peaks

  !--------------------------------------  
  !communicate cummulative mass profiles
  !--------------------------------------  
  call build_peak_communicator()
  do i=1,nmassbins
    call virtual_peak_dp(cmp(1,i), 'sum')
    call boundary_peak_dp(cmp(1,i)) 
  end do

end subroutine get_cmp
!########################################
!########################################
!########################################
subroutine get_closest_border(ilevel)
  use amr_commons
  use clfind_commons
  implicit none
  integer, intent(in) :: ilevel
  !---------------------------------------------------------------------------
  ! Find closest border to centre of mass. Modified subroutine saddlepoint_search
  !---------------------------------------------------------------------------
  integer             ::  ipart,ipeak,ip,jlevel,next_level
  integer             ::  local_peak_id,global_peak_id
  integer,dimension(1:nvector)  ::  ind_cell
  logical,dimension(1:npeaks_max) ::  check
  
  character(len=80) :: fileloc
  character(len=5)  :: nchar, nchar2

  if(verbose)write(*,*) "Entered get_closest_border"

  check=.false.
  do ipeak=1, hfree-1

    check(ipeak)=cmp_distances(ipeak,nmassbins)>0.0 !peak must have particles somewhere
    check(ipeak)=check(ipeak).and.to_iter(ipeak)
    check(ipeak)=check(ipeak).and.lev_peak(ipeak)==ilevel !peak must have correct level

    if(check(ipeak))  closest_border(ipeak) = 3.d0*boxlen**2 !reset value
  end do



  !-------------------------
  ! Loop over all testcells
  !-------------------------
  ip=0
  do ipart=1,ntest
    jlevel=levp(ipart) ! level
    next_level=0 !level of next particle
    if(ipart<ntest)next_level=levp(ipart+1)

    
    global_peak_id=flag2(icellp(ipart))
    if (global_peak_id/=0) then 

      call get_local_peak_id(global_peak_id,local_peak_id)

      if(check(local_peak_id)) then !if testcell is of interest:
        ip=ip+1
        ind_cell(ip)=icellp(ipart)
        if(ip==nvector .or. next_level /= jlevel)then
          call unbinding_neighborsearch(ind_cell,ip,jlevel)
          ip=0
        endif
      endif
    end if
  end do
  if (ip>0)call unbinding_neighborsearch(ind_cell,ip,jlevel)

  !------------------------
  ! Communicate results
  !------------------------

  call build_peak_communicator()
  call virtual_peak_dp(closest_border,'min')
  call boundary_peak_dp(closest_border)



  !------------------------
  ! Write output
  !------------------------

  if (unbinding_formatted_output) then
    call title(ifout-1, nchar)
    call title(myid, nchar2)

    fileloc=TRIM('output_'//TRIM(nchar)//'/unb_form_out_closestborders.txt'//nchar2)
    open(unit=666,file=fileloc,form='formatted')
    if (ilevel==mergelevel_max) write(666,'(4A18)') "peak id", "x", "y", "z"
    do ipart=1,npeaks
      if(cmp_distances(ipart,nmassbins)>0.and.ilevel==mergelevel_max) then
        write(666,'(I18,3E18.9E2)') ipart+ipeak_start(myid), closest_border(ipart)
      end if
    end do
    close(666)
  end if

end subroutine get_closest_border
!#####################################################
!#####################################################
!#####################################################
!#####################################################
subroutine unbinding_neighborsearch(ind_cell,np,jlevel)
  use amr_commons
  implicit none
  integer,dimension(1:nvector),intent(in) ::ind_cell  !array of indices of cells that I want to check
  integer,intent(in)            ::np      !number of actual cells in ind_cell
  integer,intent(in)            ::jlevel    !cell level

  !------------------------------------------------------------
  ! Modified subroutine neighborsearch
  ! This routine constructs all neighboring leaf cells at levels 
  ! jlevel-1, jlevel, jlevel+1.
  ! Then performs the check if the neighbors are a border
  ! in order to find the closest border to the center of mass
  !------------------------------------------------------------

  integer::j,ind,nx_loc,i1,j1,k1,i2,j2,k2,i3,j3,k3,ix,iy,iz
  integer::i1min,i1max,j1min,j1max,k1min,k1max
  integer::i2min,i2max,j2min,j2max,k2min,k2max
  integer::i3min,i3max,j3min,j3max,k3min,k3max
  real(dp)::dx,dx_loc,scale
  integer ,dimension(1:nvector)::clump_nr,indv,ind_grid,grid,ind_cell_coarse

  real(dp),dimension(1:twotondim,1:3)::xc
  integer ,dimension(1:99)::neigh_cell_index,cell_levl,test_levl
  real(dp),dimension(1:99,1:ndim)::xtest,xrel
  logical ,dimension(1:99)::ok
  real(dp),dimension(1:3)::skip_loc
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids 
  integer::ntestpos,ntp,idim,ipos

  real(dp),dimension(1:3)::this_cellpos



  ! Mesh spacing in that level
  dx=0.5D0**jlevel 
  nx_loc=(icoarse_max-icoarse_min+1)
  !skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  ! Integer constants
  i1min=0; i1max=1; i2min=0; i2max=2; i3min=0; i3max=3
  j1min=0; j1max=1; j2min=0; j2max=2; j3min=0; j3max=3
  k1min=0; k1max=1; k2min=0; k2max=2; k3min=0; k3max=3

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
    indv(j)   = (ind_cell(j)-ncoarse-1)/ngridmax+1     ! cell position in grid
    ind_grid(j) = ind_cell(j)-ncoarse-(indv(j)-1)*ngridmax ! grid index
    clump_nr(j) = flag2(ind_cell(j))             ! save clump number
  end do 




  ntestpos=3**ndim
  if(jlevel>levelmin)  ntestpos=ntestpos+2**ndim
  if(jlevel<nlevelmax) ntestpos=ntestpos+4**ndim

  !================================
  ! generate neighbors level jlevel-1
  !================================
  ntp=0
  if(jlevel>levelmin)then
    ! Generate 2x2x2  neighboring cells at level jlevel-1   
    do k1=k1min,k1max
      do j1=j1min,j1max
        do i1=i1min,i1max     
          ntp=ntp+1
          xrel(ntp,1)=(2*i1-1)*dx_loc
          xrel(ntp,2)=(2*j1-1)*dx_loc
          xrel(ntp,3)=(2*k1-1)*dx_loc
          test_levl(ntp)=jlevel-1
        end do
      end do
    end do
  endif

  !================================
  ! generate neighbors at level jlevel
  !================================
  ! Generate 3x3x3 neighboring cells at level jlevel
  do k2=k2min,k2max
    do j2=j2min,j2max
      do i2=i2min,i2max
        ntp=ntp+1
        xrel(ntp,1)=(i2-1)*dx_loc
        xrel(ntp,2)=(j2-1)*dx_loc
        xrel(ntp,3)=(k2-1)*dx_loc
        test_levl(ntp)=jlevel
      end do
    end do
  end do
  
  !===================================
  ! generate neighbors at level jlevel+1
  !====================================
  if(jlevel<nlevelmax)then
    ! Generate 4x4x4 neighboring cells at level jlevel+1
    do k3=k3min,k3max
      do j3=j3min,j3max
        do i3=i3min,i3max
          ntp=ntp+1
          xrel(ntp,1)=(i3-1.5)*dx_loc/2.0
          xrel(ntp,2)=(j3-1.5)*dx_loc/2.0
          xrel(ntp,3)=(k3-1.5)*dx_loc/2.0
          test_levl(ntp)=jlevel+1
        end do
      end do
    end do
  endif



  ! Gather 27 neighboring father cells (should be present anytime !)
  do j=1,np
    ind_cell_coarse(j)=father(ind_grid(j))
  end do
  call get3cubefather(ind_cell_coarse,nbors_father_cells,nbors_father_grids,np,jlevel)


  do j=1,np
    ok=.false.
    do idim=1,ndim
      !get real coordinates of neighbours
      xtest(1:ntestpos,idim)=(xg(ind_grid(j),idim)+xc(indv(j),idim)-skip_loc(idim))*scale+xrel(1:ntestpos,idim)
      if(jlevel>levelmin)xtest(1:twotondim,idim)=xtest(1:twotondim,idim)+xc(indv(j),idim)*scale
    end do
    grid(1)=ind_grid(j)
    call get_cell_index_fast(neigh_cell_index,cell_levl,xtest,ind_grid(j),nbors_father_cells(j,1:threetondim),ntestpos,jlevel)
   
    do ipos=1,ntestpos
      !make sure neighbour is a leaf cell
      if(son(neigh_cell_index(ipos))==0.and.cell_levl(ipos)==test_levl(ipos)) then
        ok(ipos)=.true. 
      end if
    end do
   
    !get position of the cell whose neighbours will be tested
    do idim=1,ndim
      this_cellpos(idim)=(xg(ind_grid(j),idim)+xc(indv(j),idim)-skip_loc(idim))*scale
    end do

    !check neighbors
!if (myid==1) write(*,*) "ID 1 calling bordercheck for peak", clump_nr(j)
    call bordercheck(this_cellpos,clump_nr(j),xtest,neigh_cell_index,ok,ntestpos)
    ! bordercheck (this_cellpos=position of cell to test;
    ! clump_nr(j)=peak ID of cell to test;
    ! xtest=positions of neighbour cells; 
    ! neigh_cell_index=index of neighbour cells;
    ! ok = if neighbour cell is leaf cell;
    ! ntestpos = how many neighbour cells there are
  end do

end subroutine unbinding_neighborsearch
!########################################
!########################################
!########################################
subroutine bordercheck(this_cellpos,clump_nr,xx,neigh_cell_index,ok,np)
  !----------------------------------------------------------------------
  ! routine to check wether neighbor belongs to another clump and is closer to the center
  ! of mass than all others before
  ! modified subroutine saddlecheck
  !----------------------------------------------------------------------
  use amr_commons
  use clfind_commons
  implicit none
  real(dp), dimension(1:np,1:ndim), intent(in)  :: xx     ! positions of neighbour cells
  real(dp), dimension(1:ndim),    intent(in)  :: this_cellpos   ! position of test cell whose neighbours 
                                  ! are to be tested
  integer,  dimension(1:99),      intent(in)  :: neigh_cell_index ! cell index of neighbours
  integer,              intent(in)  :: clump_nr   ! global peak ID of cell whose neighbours
                                  ! will be tested
  logical,  dimension(1:99),      intent(inout) :: ok     ! wether cell should be checkedre
  integer,              intent(in)  :: np     ! number of neighbours to be looped over


  real(dp), dimension(1:99,1:ndim)  :: pos    ! position of border for each neighbour
  integer,  dimension(1:99)     :: neigh_cl !clump number of neighbour,local peak id of neighbour
  real(dp)              :: newsum
  integer               :: i,j,ipeak
  real(dp), dimension(1:3)      :: period


  do j=1,np
    neigh_cl(j)=flag2(neigh_cell_index(j))!index of the clump the neighboring cell is in 

    ok(j)=ok(j).and. clump_nr/=0 ! temporary fix...
    ok(j)=ok(j).and. neigh_cl(j)/=0 !neighboring cell is in a clump. If neighbour not in clump, clump is still considered isolated.
    ok(j)=ok(j).and. neigh_cl(j)/=clump_nr !neighboring cell is in another clump 
  end do


  call get_local_peak_id(clump_nr,ipeak)

  do j=1,np
    if(ok(j))then ! if all criteria met, you've found a neighbour cell that belongs to a different clump 

      period=0.d0
      if (periodical) then
        do i=1, ndim
          if (xx(j,i)-clmp_com_pb(ipeak,i) > 0.5*boxlen)    period(i)=(-1.0)*boxlen
          if (xx(j,i)-clmp_com_pb(ipeak,i) < (-0.5*boxlen)) period(i)=boxlen
        end do
      end if

      do i=1, ndim
        !the cells will be nighbours, so no need to compute two different periodic corrections
        pos(j,i)=(xx(j,i)+period(i)+this_cellpos(i)+period(i))*0.5 
      end do

      newsum=0
      do i=1, ndim
        newsum=newsum+(pos(j,i)-clmp_com_pb(ipeak,i))**2
      end do

      if (newsum<closest_border(ipeak))  closest_border(ipeak)=newsum
      end if
  end do


end subroutine bordercheck
!########################################
!########################################
!########################################
subroutine particle_unbinding(ipeak,final_round)
  use amr_commons, only: dp
  use clfind_commons

  implicit none
  integer, intent(in) :: ipeak  !peak to loop over
  logical, intent(in) :: final_round !if it is the final round => whether to write
  !--------------------------------------------------------------
  !This subroutine loops over all particles in the linked list of
  !peak ipeak and checks if they are bound.
  !--------------------------------------------------------------

  integer :: thispart, ipeak_test, ipart
  real(dp):: phi_border !the potential at the border of the peak patch closest 
              !to the center of mass
  real(dp):: dist_border  !distance to the border


  !compute the potential for this peak on the points of the mass bin distances
  call compute_phi(ipeak)

  !compute potential at the closest border from the center of mass
  phi_border=0.d0
  if(saddle_pot) then
    dist_border=sqrt(closest_border(ipeak))
    if(dist_border<=cmp_distances(ipeak,nmassbins)) phi_border=potential(ipeak,dist_border)
  end if


  thispart=clmppart_first(ipeak)


  if (final_round) then
    do ipart=1, nclmppart(ipeak)    ! loop over particle LL
      call get_local_peak_id(clmpidp(thispart),ipeak_test)
      if (ipeak_test==ipeak) then   ! if this particle needs to be checked for unbinding
                      ! particle may be assigned to child clump
        candidates=candidates+1
        if(unbound(ipeak,thispart,phi_border)) then  
          !unbound(ipeak,thispart,phi_border) is a logical function. See below.

          nunbound=nunbound+1            !counter
          clmpidp(thispart)=new_peak(ipeak)    !update clump id
        end if
      end if
      thispart=clmppart_next(thispart)
    end do
  else 
    if (to_iter(ipeak)) then
      hasatleastoneptcl(ipeak)=0       ! set to false
      do ipart=1, nclmppart(ipeak)     ! loop over particle LL
        if(unbound(ipeak,thispart,phi_border)) then  
          !unbound(ipeak,thispart,phi_border) is a logical function. See below.

          niterunbound=niterunbound+1
          contributes(thispart) = .false.   ! particle doesn't contribute to
                            ! clump properties
        else
          hasatleastoneptcl(ipeak)=1 !there are contributing particles for this peak
        end if
        thispart=clmppart_next(thispart)
      end do

    end if 
  end if

!----------------------
!----------------------
!----------------------

contains

logical function unbound(ipeak, part_ind, phi_border)
  !-----------------------------------------------------------
  ! This function checks if the particle of clump ipeak
  ! with particle index part_ind is bound to the clump or not.
  ! It returns TRUE if the particle is not energetically bound.
  !-----------------------------------------------------------

  use pm_commons
  use clfind_commons
  implicit none

  integer, intent(in) :: ipeak, part_ind
  real(dp), intent(in) :: phi_border
  real(dp) :: distance, kinetic_energy, minusphi
  real(dp),dimension(1:3) :: period
  integer :: i

  period=0.d0
  if (periodical) then
    do i=1, ndim
      if (xp(part_ind,i)-clmp_com_pb(ipeak,i) > 0.5*boxlen   ) period(i)=(-1.0)*boxlen
      if (xp(part_ind,i)-clmp_com_pb(ipeak,i) < (-0.5*boxlen)) period(i)=boxlen
    end do
  end if

  distance=(xp(part_ind,1)+period(1)-clmp_com_pb(ipeak,1))**2 + &
      (xp(part_ind,2)+period(2)-clmp_com_pb(ipeak,2))**2 + &
      (xp(part_ind,3)+period(3)-clmp_com_pb(ipeak,3))**2
  distance=sqrt(distance)


  kinetic_energy=0.5*((vp(part_ind,1)-clmp_vel_pb(ipeak,1))**2 + &
      (vp(part_ind,2)-clmp_vel_pb(ipeak,2))**2 + &
      (vp(part_ind,3)-clmp_vel_pb(ipeak,3))**2)
  
  minusphi=potential(ipeak,distance)

  unbound=(kinetic_energy>=minusphi-phi_border)


end function unbound
!----------------------
!----------------------
!----------------------
real(dp) function potential(ipeak,distance)
  !------------------------------------------------------------------
  !This function interpolates the potential of a particle for given distance
  !It returns (-1)*phi
  !------------------------------------------------------------------

  integer, intent(in) :: ipeak
  real(dp),intent(in) :: distance !is computed in function 'unbound', then passed

  integer :: ibin, thisbin
  real(dp) :: a,b

  ibin=1
  !thisbin: the first cmp_distance which is greater than particle distance
  thisbin=1
  do 
    if (distance<=cmp_distances(ipeak,ibin)) then
      thisbin=ibin
      exit
    else
      ibin=ibin+1
    end if
  end do

  a=(phi_unb(thisbin)-phi_unb(thisbin-1))/(cmp_distances(ipeak,thisbin)-cmp_distances(ipeak,thisbin-1))
  b=phi_unb(thisbin-1)-a*cmp_distances(ipeak,thisbin-1)
  potential=(-1)*a*distance-b

end function potential 
end subroutine particle_unbinding
!###############################################
!###############################################
!###############################################
subroutine compute_phi(ipeak)
  !-----------------------------------------------------------
  ! This subroutine computes the potential on each massbin
  ! It writes potential[r=ibin] into the array phi_unb[ibin]
  !-----------------------------------------------------------
  use clfind_commons
  use amr_commons!, only: dp
  integer, intent(in) :: ipeak
  real(dp) :: delta,add
  integer  :: i

  !Writing unformatted output
  character(len=5)  :: bins
  character(len=5)  :: peak
  !character(len=10) :: peak !in case there are too many (>1E6) peaks
  character(len=80) :: fileloc
  character(len=5)  :: nchar

  !compute part of integral/sum for each bin
  phi_unb(nmassbins)=0.0
  do i=2,nmassbins
    delta=cmp_distances(ipeak,i)-cmp_distances(ipeak,i-1)
    phi_unb(i-1)=-0.5*GravConst*(cmp(ipeak,i)/cmp_distances(ipeak,i)**2+cmp(ipeak,i-1)/cmp_distances(ipeak,i-1)**2)*delta
  end do
  delta=cmp_distances(ipeak,1)-cmp_distances(ipeak,0)
  phi_unb(0)=-0.5*GravConst*(cmp(ipeak,1)/cmp_distances(ipeak,1)**2)*delta

  !sum bins up
  !does not need to be done for i=nmassbins!
  add=-cmp(ipeak,nmassbins)/cmp_distances(ipeak,nmassbins)*GravConst !-G*M_tot/r_max
  do i=nmassbins-1,0,-1
    phi_unb(i)=phi_unb(i)+phi_unb(i+1) !stopps at phi_unb(1)
    phi_unb(i+1)=phi_unb(i+1)+add
  end do
  phi_unb(0)=phi_unb(0)+add !bypass division by 0, needed for interpolation.


  if (unbinding_formatted_output) then
  !WARNING:
  !Crashes if there are more than 10^6 peaks because it creates a file with ****** in the name,
  !and multiple processes try to write into same file.
  !This happens because the subroutine title() only goes up to 5 digits.
  !Instead, switch to commented out part. (Also the other declaration for char peaks)
    if (ipeak<=npeaks) then
      !generate filename integers
      call title(ifout-1, nchar)
      call title(nmassbins,bins)
      call title(ipeak+ipeak_start(myid),peak)
      !write(peak,'(i10)') ipeak+ipeak_start ! in case of too many peaks

      fileloc=TRIM('output_'//TRIM(nchar)//'/unb_form_out_phi_'//TRIM(peak)//'-'//TRIM(bins)//'.txt')
      open(unit=667,file=fileloc,form='formatted')
      write(667,'(2A20)') "distance", "potential"
      do i=0,nmassbins
        write(667,'(2E20.8E2)') cmp_distances(ipeak,i), phi_unb(i)
      end do
      close(667)
      
      fileloc=TRIM('output_'//TRIM(nchar)//'/unb_form_out_cmp_'//TRIM(peak)//'-'//TRIM(bins)//'.txt')
      open(unit=666,file=fileloc,form='formatted')
      write(666,'(2A20)') "distance", "cumulative mass"
      do i=0,nmassbins
        write(666,'(2E20.8E2)') cmp_distances(ipeak,i), cmp(ipeak,i)
      end do
      close(666)

      fileloc=TRIM('output_'//TRIM(nchar)//'/unb_form_out_clumpproperties_'//TRIM(peak)//'.txt')
      open(unit=666, file=fileloc, form='formatted')
      write(666, '(A10,7A18)') "clmp id", "x", "y ", "z", "vx", "vy", "vz", "max_dist"
      if (clmp_mass_pb(ipeak)>0) then
        write(666,'(I10,7E18.9E2)') ipeak+ipeak_start(myid), clmp_com_pb(ipeak,1), clmp_com_pb(ipeak,2), clmp_com_pb(ipeak,3), clmp_vel_pb(ipeak,1), clmp_vel_pb(ipeak,2),clmp_vel_pb(ipeak,3), cmp_distances(ipeak, nmassbins)
      end if
      close(666)
    end if
  end if

end subroutine compute_phi 
!###############################################
!###############################################
!###############################################
subroutine allocate_unbinding_arrays()
  use clfind_commons
  use pm_commons, only:npartmax
  implicit none

  !----------------------------------------------
  ! This subroutine allocates the necessary 
  ! arrays and gives them initial values.
  !----------------------------------------------

  !-------------------
  ! Clump properties
  !-------------------
  allocate(clmp_com_pb(1:npeaks_max,1:3))
  clmp_com_pb=0.0
  allocate(clmp_vel_pb(1:npeaks_max,1:3))
  clmp_vel_pb=0.0
  allocate(clmp_mass_pb(1:npeaks_max))
  clmp_mass_pb=0.0
  allocate(cmp_distances(1:npeaks_max,0:nmassbins))
  cmp_distances=0.0
  allocate(cmp(1:npeaks_max,0:nmassbins))
  cmp=0.d0
  ! careful with this! The first index of the second subscript
  ! of the cumulative mass aray (index 0) is there for reference
  ! for the enclosed mass interpolation.

  allocate(phi_unb(0:nmassbins)) ! array where to store the potential
  phi_unb=0.d0

  if (saddle_pot) then
    allocate(closest_border(1:npeaks_max)) !point of the closest border to CoM
    closest_border=3.d0*boxlen**2
  end if

  allocate(to_iter(1:npeaks_max)) ! peak needs to be checked or not
  to_iter=.true.

  if (iter_properties) then
    allocate(oldcom(1:npeaks_max,1:3))
    allocate(oldvel(1:npeaks_max,1:3))
    allocate(oldcmpd(1:npeaks_max))
    allocate(oldm(1:npeaks_max))
  end if

  allocate(hasatleastoneptcl(1:npeaks_max))
  hasatleastoneptcl=1 !initiate to yes


  !----------------------
  ! Particle linked list
  !----------------------
  allocate(clmpidp(1:npartmax))
  clmpidp=0

  allocate(clmppart_first(1:npeaks_max)) !linked lists containing particles 
  clmppart_first=0

  allocate(clmppart_last(1:npeaks_max)) !linked lists containing particles 
  clmppart_last=0

  allocate(clmppart_next(1:npartmax)) !linked lists containing particles 
  clmppart_next=0

  allocate(nclmppart(1:npeaks_max)) !linked lists containing particles 
  nclmppart=0

  allocate(contributes(1:npartmax)) !particle contributes to clump properties or not
  contributes=.true.



end subroutine allocate_unbinding_arrays
!########################################
!########################################
!########################################
subroutine deallocate_unbinding_arrays()
  use clfind_commons
  implicit none

  deallocate(clmp_com_pb)
  deallocate(clmp_vel_pb)
  deallocate(clmp_mass_pb)
  deallocate(cmp_distances)
  deallocate(cmp)

  deallocate(phi_unb)

  if(saddle_pot) deallocate(closest_border)

  deallocate(to_iter)

  if (iter_properties) deallocate(oldcom,oldvel,oldcmpd,oldm)

  deallocate(hasatleastoneptcl)


  deallocate(clmpidp)
  deallocate(clmppart_last)
  deallocate(clmppart_first)
  deallocate(clmppart_next)
  deallocate(nclmppart)
  deallocate(contributes)

end subroutine deallocate_unbinding_arrays
!############################################
!############################################
!############################################
subroutine unbinding_write_formatted_output()

  !------------------------------------------------------------------
  ! This subroutine outputs all the interesting particle attributes.
  !------------------------------------------------------------------
  
  use amr_commons
  use pm_commons
  use clfind_commons
  implicit none
  !select which output(s)
  logical :: particles, CoM, mass_prof, dist


  !iterators
  integer:: ipeak, i

  !filename
  character(len=80) :: fileloc
  character(len=5) :: nchar, nchar2


  ! set which output you want here:
  particles=.true.  ! all particles: coordinates, velocity, clump ID
  CoM=.true.      ! for all clumps which are not halo-namegivers:
            ! clump id, PHEW peak position, particlebased Center
            ! of Mass, maximal distance of particles to CoM
  dist=.false.     ! mass bin distances
  mass_prof=.false.  ! mass profiles; binning width is bin distances


  !generate filename integers
  call title(ifout-1, nchar)
  call title(myid, nchar2)
  
  if(particles) call unbinding_formatted_particleoutput(.false.)


  if (CoM) then
    fileloc=TRIM('output_'//TRIM(nchar)//'/unb_form_out_COM.txt'//nchar2)
  
    open(unit=666, file=fileloc, form='formatted')
    write(666, '(A10,7A18)') "clmp id", "x", "y ", "z", "vx", "vy", "vz", "max_dist"
    do ipeak=1, npeaks
      if (clmp_mass_pb(ipeak)>0) then
        write(666,'(I10,7E18.9E2)') ipeak+ipeak_start(myid), clmp_com_pb(ipeak,1), clmp_com_pb(ipeak,2), clmp_com_pb(ipeak,3), clmp_vel_pb(ipeak,1), clmp_vel_pb(ipeak,2),clmp_vel_pb(ipeak,3), cmp_distances(ipeak, nmassbins)
      end if
    end do
    close(666)
  end if


  !CUMULATIVE MASS PROFILES
  if (mass_prof) then
    if (logbins) then
      fileloc=TRIM('output_'//TRIM(nchar)//'/unb_form_out_CMP-log.txt'//nchar2)
    else
      fileloc=TRIM('output_'//TRIM(nchar)//'/unb_form_out_CMP-lin.txt'//nchar2)
    end if
    
    open(unit=666, file=fileloc, form='formatted')
    write(666,'(A10)',advance='no') ipeak+ipeak_start(myid)
    do i=1, nmassbins
      write(666,'(I18)',advance='no') i
    end do
    write(666,*)
    do ipeak=1, npeaks
      if (cmp_distances(ipeak,nmassbins)>0.0) then
        write(666,'(I10)',advance='no') ipeak+ipeak_start(myid)
        do i=1, nmassbins
          write(666,'(E18.9E2)',advance='no') cmp(ipeak,i)
        end do
        write(666,*)
      end if
    end do
    close(666)
  end if



  !MASS PROFILE BIN DISTANCES
  if (dist) then
    fileloc=TRIM('output_'//TRIM(nchar)//'/unb_form_out_distances.txt'//nchar2)

    open(unit=666, file=fileloc, form='formatted')
    write(666,'(A10)',advance='no') ipeak+ipeak_start(myid)
    do i=1, nmassbins
      write(666,'(I18)',advance='no') i
    end do
    write(666,*)
    do ipeak=1, npeaks
      if (cmp_distances(ipeak,nmassbins)>0.0) then
        write(666,'(I10)',advance='no') ipeak+ipeak_start(myid)
        do i=1, nmassbins
          write(666,'(E18.9E2)',advance='no') cmp_distances(ipeak,i)
        end do
        write(666,*)
      end if
    end do
    close(666)
  end if


end subroutine unbinding_write_formatted_output
!###############################################
!###############################################
!###############################################
subroutine unbinding_formatted_particleoutput(before)

  !--------------------------------------------------------------------
  ! This subroutine writes all the interesting particle attributes to
  ! file. 
  ! If before = .true. (called before the unbinding starts), it will 
  ! create a new directory "before" in the output directory and
  ! write the particle attributes as found by PHEW to file.
  !--------------------------------------------------------------------

  use amr_commons
  use pm_commons
  use clfind_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer :: info
#endif

  logical,intent(in) :: before
  
  !filename
  character(len=80) :: fileloc
  character(len=5)  :: nchar, nchar2

  !local vars
  integer       :: i
  character(len=80) :: cmnd

  if (before) then

    if (myid==1) then ! create before dir
      call title(ifout-1,nchar)
      cmnd='mkdir -p output_'//TRIM(nchar)//'/before'
      call system(TRIM(cmnd))
    end if

#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif

  end if


  !generate filename
  call title(ifout-1, nchar)
  call title(myid, nchar2)

  if (before) then
    fileloc=TRIM('output_'//TRIM(nchar)//'/before/unb_form_out_particleoutput.txt'//nchar2)
  else
    fileloc=TRIM('output_'//TRIM(nchar)//'/unb_form_out_particleoutput.txt'//nchar2)
  end if

 

  open(unit=666, file=fileloc, form='formatted')
  write(666, '(9A18)') "x", "y", "z", "vx", "vy", "vz", "clmp id", "mass", "pid"
  do i=1, npartmax
    if(levelp(i)>0) then
      write(666, '(6E18.9E2,I18,E18.9E2,I18)') xp(i,1), xp(i,2), xp(i,3), vp(i,1), vp(i,2), vp(i,3), clmpidp(i),mp(i),idp(i)
    end if 
  end do

  close(666)
end subroutine unbinding_formatted_particleoutput
!#############################################
!#############################################
!#############################################

! endif: NDIM == 3
#endif


