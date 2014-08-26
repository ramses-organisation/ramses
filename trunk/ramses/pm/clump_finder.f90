subroutine clump_finder(create_output,keep_alive)
  use amr_commons
  use poisson_commons, ONLY:phi,rho
  use clfind_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
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

  integer::istep,nskip,ilevel,info,icpu,nmove,nmove_all,nzero,nzero_all
  integer::i,j,ntest_all,peak_nr
  integer,dimension(1:ncpu)::ntest_cpu,ntest_cpu_all
  integer,dimension(1:ncpu)::npeaks_per_cpu_tot
  logical::all_bound

  if(verbose.and.myid==1)write(*,*)' Entering clump_finder'
  first_pass=.true.

  !---------------------------------------------------------------
  ! Compute rho from gas density or dark matter particles
  !---------------------------------------------------------------
  if(ivar_clump==0)then
     do ilevel=levelmin,nlevelmax
        if(pic)call make_tree_fine(ilevel)
        if(poisson)call rho_fine(ilevel,2)
        if(pic)then
           call kill_tree_fine(ilevel)
           call virtual_tree_fine(ilevel)
        endif
     end do
     do ilevel=nlevelmax,levelmin,-1
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
  call MPI_ALLREDUCE(ntest_cpu,ntest_cpu_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  ntest_cpu(1)=ntest_cpu_all(1)
#endif
  do icpu=2,ncpu
     ntest_cpu(icpu)=ntest_cpu(icpu-1)+ntest_cpu_all(icpu)
  end do
  ntest_all=ntest_cpu(ncpu)
  if(myid==1)then
     if(ntest_all.gt.0.and.clinfo)then
        write(*,'(" Total number of cells above threshold=",I10)')ntest_all
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
  npeaks=0; nmove=0; nzero=0
  if(ntest>0)then
     if(ivar_clump==0)then  ! case 1: count peaks
        call scan_for_peaks(rho(1),npeaks,nzero,1)
     else
        if(hydro)then       ! case 1: count peaks
           call scan_for_peaks(uold(1,ivar_clump),npeaks,nzero,1)
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
  peak_nr=ipeak_start(myid)

  !----------------------------------------------------------------------
  ! Flag peaks with global peak id using flag2 array
  ! Compute peak density using max_dens array
  !----------------------------------------------------------------------
  nmove=0
  nskip=peak_nr
  ! Compute the size of the peak-based arrays
  npeaks_max=MAX(4*maxval(npeaks_per_cpu_tot),1000)
  allocate(max_dens(npeaks_max))
  max_dens=0.
  flag2=0
  if(ntest>0)then
     if(ivar_clump==0)then  ! case 2: flag peaks
        call scan_for_peaks(rho(1),nskip,nzero,2)
     else
        if(hydro)then       ! case 2: flag peaks
           call scan_for_peaks(uold(1,ivar_clump),nskip,nzero,2)
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
  ! - done when nmove=0 (for single core, only one sweep is necessary)
  !---------------------------------------------------------------------
  if (myid==1.and.ntest_all>0)write(*,*)'Finding peak patches'
  nmove=1
  istep=0
  do while (nmove.gt.0)
     nmove=0
     nzero=0
     nskip=peak_nr
     if(ntest>0)then
        if(ivar_clump==0)then
           call scan_for_peaks(rho(1),nmove,nzero,3)
        else
           if(hydro)then
              call scan_for_peaks(uold(1,ivar_clump),nmove,nzero,3)
           endif
        endif
     endif
     do ilevel=nlevelmax,levelmin,-1
        call make_virtual_fine_int(flag2(1),ilevel)
     end do
     istep=istep+1
#ifndef WITHOUTMPI 
     call MPI_ALLREDUCE(nmove,nmove_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
     nmove=nmove_all
     call MPI_ALLREDUCE(nzero,nzero_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
     nzero=nzero_all
#endif   
     if(myid==1.and.ntest_all>0.and.clinfo)write(*,*)"istep=",istep,"nmove=",nmove
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
     if(myid==1.and.clinfo)then
        write(*,*)"Output status of peak memory."
     endif
     if(clinfo)call analyze_peak_memory
     if(clinfo.and.saddle_threshold.LE.0)call write_clump_properties(.false.)
     if(create_output)then
        if(myid==1)write(*,*)"Outputing clump properties to disc."
        call write_clump_properties(.true.)
     endif
     
  end if

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
subroutine scan_for_peaks(xx,n,nzero,action)
  use amr_commons
  use clfind_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h' 
#endif
  integer::n,nzero,action
  real(dp),dimension(1:ncoarse+ngridmax*twotondim)::xx
  !----------------------------------------------------------------------
  ! vectorization of the neighborsearch for the action cases
  ! 1: count the peaks (no denser neighbor)
  ! 2: count and flag peaks with global peak index number
  ! 3: get global clump index from densest neighbor
  !----------------------------------------------------------------------
  integer::ilevel,next_level,ipart,jpart,ip
  integer,dimension(1:nvector)::ind_part,ind_cell,ind_max
!  logical,save::first_pass=.true. !potential problem here when clumpfinder is called multiple times?
  !Deal with when running clump_finder for sinks...

  ! shortcut neighborsearch when it has already been done before
  ! by saving each cells densest neighbor in imaxp 
  if(.not. first_pass)then
     select case (action)
     case (1)   ! Count peaks  
        do ipart=1,ntest
           jpart=testp_sort(ipart)
           if(imaxp(jpart).EQ.-1)n=n+1
        end do
     case (2)   ! Initialize flag2 to peak global index
        do ipart=1,ntest
           jpart=testp_sort(ipart)
           if(imaxp(jpart).EQ.-1)then
              n=n+1
              flag2(icellp(jpart))=n
              max_dens(n-ipeak_start(myid))=xx(icellp(jpart))
           endif
        end do
     case (3) ! Propagate flag2
        do ipart=1,ntest
           jpart=testp_sort(ipart)
           if(imaxp(jpart).NE.-1)then
              if(flag2(icellp(jpart)).ne.flag2(imaxp(jpart)))n=n+1
              flag2(icellp(jpart))=flag2(imaxp(jpart))
              if(flag2(icellp(jpart)).eq.0)nzero=nzero+1
           endif
        end do
     end select
     return
  endif

  ip=0
  do ipart=1,ntest
     ip=ip+1
     ilevel=levp(testp_sort(ipart)) ! level
     next_level=0 !level of next particle
     if(ipart<ntest)next_level=levp(testp_sort(ipart+1))
     ind_cell(ip)=icellp(testp_sort(ipart))
     ind_part(ip)=testp_sort(ipart)
     if(ip==nvector .or. next_level /= ilevel)then
        call neighborsearch(xx(1),ind_cell,ind_max,ip,n,nzero,ilevel,action)
        do jpart=1,ip
           imaxp(ind_part(jpart))=ind_max(jpart)
        end do
        ip=0
     endif
  end do
  if (ip>0)then
     call neighborsearch(xx(1),ind_cell,ind_max,ip,n,nzero,ilevel,action)
     do jpart=1,ip
        imaxp(ind_part(jpart))=ind_max(jpart)
     end do
  endif

  first_pass=.false.

end subroutine scan_for_peaks
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
  integer::i,j,info,dummyint,dummyzero
  integer,dimension(1:nvector)::ind_cell,ind_max

  ip=0
  do ipart=1,ntest
     ip=ip+1
     ilevel=levp(testp_sort(ipart)) ! level
     next_level=0 !level of next particle
     if(ipart<ntest)next_level=levp(testp_sort(ipart+1))
     ind_cell(ip)=icellp(testp_sort(ipart))
     if(ip==nvector .or. next_level /= ilevel)then
        call neighborsearch(xx(1),ind_cell,ind_max,ip,dummyint,dummyzero,ilevel,4)
        ip=0
     endif
  end do
  if (ip>0)call neighborsearch(xx(1),ind_cell,ind_max,ip,dummyint,dummyzero,ilevel,4)

end subroutine saddlepoint_search
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine neighborsearch(xx,ind_cell,ind_max,np,count,count_zero,ilevel,action)
  use amr_commons
  use clfind_commons,ONLY:max_dens,ipeak_start
  implicit none
  integer::np,count,count_zero,ilevel,action
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
  logical ,dimension(1:nvector)::okpeak,okdummy
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
     if (action<4)call peakcheck(xx(1),cell_index,okpeak(j),ok,density_max(j),ind_max(j),ntestpos)
     if (action==4)call saddlecheck(xx(1),ind_cell(j),cell_index,clump_nr(j),ok,ntestpos)
     
  end do
     



  !===================================
  ! choose action for different cases
  !====================================
  select case (action)
  case (1)   ! Count peaks  
     do j=1,np
        if(okpeak(j))then
           count=count+1
           ind_max(j)=-1
        endif
     end do  
  case (2)   ! Initialize flag2 to peak global index
     do j=1,np
        if(okpeak(j))then 
           count=count+1
           ind_max(j)=-1
           flag2(ind_cell(j))=count
           max_dens(count-ipeak_start(myid))=xx(ind_cell(j))
        end if
     end do
  case (3) ! Propagate flag2
     do j=1,np
        if(flag2(ind_cell(j)).ne.flag2(ind_max(j)))count=count+1
        flag2(ind_cell(j))=flag2(ind_max(j))
        if(flag2(ind_cell(j)).eq.0)count_zero=count_zero+1
        if(okpeak(j))ind_max(j)=-1
     end do
  end select

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
     if(ok(j))then !so if there is a denser neighbor
        if(xx(cell_index(j))>density_max)then           
           okpeak=.false. !no peak
           density_max=xx(cell_index(j)) !change densest neighbor dens
           ind_max=cell_index(j) !change densest neighbor index
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

  namelist/clumpfind_params/ivar_clump,& 
       & relevance_threshold,density_threshold,&
       & saddle_threshold,mass_threshold,clinfo,&
       & n_clfind,rho_clfind,merge_unbound
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m  
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  
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
  integer::ng,np,ilevel,ind_grid
  integer,dimension(1:99)::indp,cell_lev
  real(dp),dimension(1:99,1:ndim)::xpart
  integer ,dimension(1:threetondim)::nbors_father_cells


  !-----------------------------------------------------------------------
  ! Very similar as get_cell_index_for_particle but optimized for the usage
  ! inside neighborsearch.
  !-----------------------------------------------------------------------

  integer::i,j,idim,nx_loc,ind,ix,iy,iz
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
!           icd_fine(1,1:ndim)=int(2*(xpart(j,1:ndim)-xg(son(indp(j)),1:ndim)+0.5*dx)/dx)
           icd_fine(1,1:ndim)=int(2*(xpart(j,1:ndim)*one_over_scale+skip_loc(1:ndim)-xg(son(indp(j)),1:ndim)+0.5*dx)/dx)
           call geticell99(icell_fine,icd_fine,1)
           indp(j)=ncoarse+(icell_fine(1)-1)*ngridmax+son(indp(j))
        end if
     endif
  end do

  !cell center positions for particles which sit in the level ilevel
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
  ! same as geticell but for input vector size of 99 instead of nvector
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
