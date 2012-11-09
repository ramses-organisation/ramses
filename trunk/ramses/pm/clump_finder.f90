subroutine clump_finder(create_output)
  use amr_commons
  use pm_commons
  use hydro_commons
  use clfind_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  logical::create_output
  !----------------------------------------------------------------------------
  ! Description of clump_finder:
  ! The clumpfinder assigns a test particle to each cell having a density above 
  ! a given threshold. These particles are moved to the densest neighbors until 
  ! particles sit in a local density maximum. The particles (now containing the
  ! peak_nr they belong to) are moved back to their original position and all
  ! the relevant properties are computed. If a so called peak patch is 
  ! considered irrelevant, it is merged to the neighbor which it is connected 
  ! to through the saddle point with the highest density.
  ! Andreas Bleuler & Romain Teyssier 10/2010 - ?
  ! Davide Martizzi & Romain Teyssier 10/2012 - ?
  !----------------------------------------------------------------------------
  ! local constants
  integer::ipart,itest,istep,nskip,ilevel,info,icpu,igrid,nmove,nmove_all
  character(LEN=5)::nchar
  character(LEN=80)::filename
  integer::jgrid

  !new variables for clump/sink comb
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  integer::j,jj,i
  real(kind=8),dimension(1:nvector,1:3)::pos
  integer,dimension(1:nvector)::cell_index,cell_levl,cc

  integer,allocatable,dimension(:)::tempr
  real(kind=8),allocatable,dimension(:)::tempi

  integer::ntest,ntest_all
  integer,dimension(1:ncpu)::ntest_cpu,ntest_cpu_all

  integer::peak_nr
  integer,dimension(1:ncpu)::npeaks_per_cpu,npeaks_per_cpu_tot

  if(verbose)write(*,*)' Entering clump_finder'

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  !-------------------------------------------------------------------------------
  ! Count test particles
  !-------------------------------------------------------------------------------
  ntest=0
  do ilevel=levelmin,nlevelmax
     call count_test_particle(ilevel,ntest)
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
     if(ntest_all.gt.0)then
        write(*,'(" Total number of test particles=",I6)')ntest_all
     endif
  end if

  !-------------------------------------------------------------------------------
  ! Allocate test particle arrays
  !-------------------------------------------------------------------------------
  allocate(denp(ntest),levp(ntest),iglobalp(ntest),icellp(ntest))
  denp=0.d0; levp=0; iglobalp=0; icellp=0
  !-------------------------------------------------------------------------------
  ! Compute test particle properties
  !-------------------------------------------------------------------------------
  itest=0
  nskip=ntest_cpu(myid)-ntest
  do ilevel=levelmin,nlevelmax
     call create_test_particle(ilevel,itest,nskip) 
  end do
  do ilevel=nlevelmax,levelmin,-1
     call make_virtual_fine_int(flag2(1),ilevel)
  end do
  !-------------------------------------------------------------------------------
  ! Sort particles according to density
  !-------------------------------------------------------------------------------
  allocate(testp_sort(ntest)) 
  do i=1,ntest
     denp(i)=-denp(i)
     testp_sort(i)=i
  end do
  if(ntest>0)call quick_sort(denp(1),testp_sort(1),ntest) 
  deallocate(denp)

  !-------------------------------------------------------------------------------               
  ! Count number of density peaks
  !-------------------------------------------------------------------------------
  npeaks=0; nmove=0
  if(ntest>0) then 
     call scan_for_peaks(ntest,nmove,npeaks,1)
  end if
  npeaks_per_cpu=0
  npeaks_per_cpu(myid)=npeaks
  write(*,*)'n_peaks on processor number',myid,'= ',npeaks

  !----------------------------------------------------------------------------                       
  ! Share number of peaks per cpu and create a list  
  !----------------------------------------------------------------------------
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
  if (myid==1.and.npeaks_tot>0)write(*,'(" Total number of density peaks found=",I6)')npeaks_tot

  !----------------------------------------------------------------------------
  ! Determine peak-ids positions for each cpu
  !----------------------------------------------------------------------------
  peak_nr=0
  do icpu=1,myid-1
     peak_nr=peak_nr+npeaks_per_cpu_tot(icpu)
  end do

  !----------------------------------------------------------------------------
  ! Change the value of flag2 at the peak positions
  !----------------------------------------------------------------------------
  nmove=0
  nskip=peak_nr
  flag2=0
  if(ntest>0)call scan_for_peaks(ntest,nmove,nskip,2)
  do ilevel=nlevelmax,levelmin,-1
     call make_virtual_fine_int(flag2(1),ilevel)
  end do

  !-------------------------------------------------------------------------------
  ! Compute position of the peaks in global peak array peak_pos_tot
  !-------------------------------------------------------------------------------
  nskip=peak_nr
  if(nstep>0)call assign_part_to_peak(ntest,nskip)

  !-------------------------------------------------------------------------------               
  ! Identify peak patches using density ordering
  !-------------------------------------------------------------------------------
  nmove=1
  istep=0
  do while (nmove.gt.0)
     nmove=0
     nskip=peak_nr
     if(ntest>0)call scan_for_peaks(ntest,nmove,nskip,3)
     do ilevel=nlevelmax,levelmin,-1
        call make_virtual_fine_int(flag2(1),ilevel)
     end do
     istep=istep+1
#ifndef WITHOUTMPI 
     call MPI_ALLREDUCE(nmove,nmove_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
     nmove=nmove_all
#endif   
     if(myid==1)write(*,*)"istep=",istep,"nmove=",nmove   
  end do

  !-------------------------------------------------------------------------------
  ! Allocate peak-patch property arrays
  !-------------------------------------------------------------------------------
  call allocate_peak_patch_arrays()

  !-------------------------------------------------------------------------------
  ! Compute peak-patch mass etc. and output these properties before merging 
  !-------------------------------------------------------------------------------
  call compute_clump_properties(ntest,ntest_all) 
  if (verbose)call write_clump_properties(.false.)

  !-------------------------------------------------------------------------------
  ! Find the saddle point densities and merge irrelevant clumps
  !-------------------------------------------------------------------------------
  if (npeaks_tot > 0)then
     !communicate across boundaries
     do ilevel=nlevelmax,levelmin,-1
        call make_virtual_fine_int(flag2(1),ilevel)
        call make_virtual_fine_dp(phi(1),ilevel)
     end do
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     call saddlepoint_search(ntest) 
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !do i=1,npeaks_tot
     !   if (myid==1)write(*,'(50(F4.3,X))'),saddle_dens_tot(i,1:npeaks_tot)
     !end do
     call merge_clumps(ntest) 
  end if

  !-------------------------------------------------------------------------------
  ! output to file  clump_properties and a complete map of all the cell-centers 
  ! together with the peak the cell belongs to
  !-------------------------------------------------------------------------------
  if (npeaks_tot > 0)then
     call clump_phi
     call compute_clump_properties_round2(ntest,ntest_all)
!     if ((sink .eqv. .false.) .or. (mod(nstep_coarse,ncontrol)==0))
     call write_clump_properties(.false.)
     if(create_output)then
        call write_peak_map(ntest)
        call write_clump_properties(.true.)
     end if
  end if

  !------------------------------------------------------------------------------
  ! if the clumpfinder is used to produce sinks, flag all the cells which contain
  ! a relevant density peak whose peak patch doesn't yet contain a sink.
  !------------------------------------------------------------------------------
  if(sink)then
     allocate(occupied(1:npeaks_tot),occupied_all(1:npeaks_tot))
     occupied=0; occupied_all=0;
     ! loop over sinks and mark all clumps containing a sink
     pos=0.0
     if(myid==1 .and. verbose)write(*,*)'looping over ',nsink,' sinks and marking their clumps'
     do j=1,nsink
        pos(1,1:3)=xsink(j,1:3)
        call cmp_cpumap(pos,cc,1)
        if (cc(1) .eq. myid)then
           call get_cell_index(cell_index,cell_levl,pos,nlevelmax,1)
           if (flag2(cell_index(1))>0)then
              occupied(flag2(cell_index(1)))=1
              if(verbose)write(*,*)'CPU # ',myid,'blocked clump # ',flag2(cell_index(1)),' for sink production because of sink # ',j
           end if
        end if
     end do
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(occupied,occupied_all,npeaks_tot,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
     occupied_all=occupied
#endif
     !------------------------------------------------------------------------------
     ! determine whether a peak patch is eligible to form a new sink.
     ! if a new sink has to be created, flag2 is set to 1 at the peak position
     !------------------------------------------------------------------------------     
     pos=0.0
     flag2=0
     call heapsort_index(max_dens_tot,sort_index,npeaks_tot)
     do j=npeaks_tot,1,-1
        jj=sort_index(j)
        if (verbose .and. myid==1)write(*,*)'clump number: ',jj
        if (relevance_tot(jj) > 1.0d-1 .and. occupied_all(jj)==0 .and. minmatch_tot(jj)==1)then           
           if (verbose .and. myid==1)write(*,*)'relevance occupied and minmatch ok'
           if (e_bind_tot4(jj)/(e_thermal_tot4(jj)+e_kin_int_tot4(jj)) > 1.)then
              if (verbose .and. myid==1)write(*,*)'bound'
              if(max_dens_tot(jj)>(n_sink/scale_nH))then
                 if (verbose .and. myid==1)write(*,*)'peak_density ok'
                 pos(1,1:3)=peak_pos_tot(jj,1:3)
                 call cmp_cpumap(pos,cc,1)
                 if (cc(1) .eq. myid)then
                    call get_cell_index(cell_index,cell_levl,pos,nlevelmax,1)
                    flag2(cell_index(1))=1
                 end if
              end if
           end if
        end if
     end do
     deallocate(occupied,occupied_all)
  endif
  
  ! Deallocate test particle and peak arrays
  deallocate(icellp)
  deallocate(levp)
  deallocate(testp_sort)
  deallocate(iglobalp)
  call deallocate_all
  
end subroutine clump_finder
!################################################################
!################################################################
!################################################################
!################################################################
subroutine count_test_particle(ilevel,ntot)
  use amr_commons
  use pm_commons
  use hydro_commons
  use cooling_module, ONLY: XH=>X, rhoc, mH 
  use random
  use clfind_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  integer::ntot
  !----------------------------------------------------------------------
  ! Description: This routine creates a particle in each cell which lies 
  ! above the density threshold.
  ! Yann Rasera  10/2002-01/2003
  !----------------------------------------------------------------------
  ! local constants
  real(dp)::d0
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:twotondim,1:3)::xc
  ! other variables
  integer ::ncache,nnew,ngrid,icpu,index_star
  integer ::igrid,ix,iy,iz,ind,i,iskip,nx_loc
  integer ::ntot_all,info
  logical ::ok_free
  real(dp),dimension(1:3)::skip_loc
  real(dp)::d,x,y,z,dx,dx_loc,scale,vol_loc,dx_min,vol_min
  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  integer ,dimension(1:nvector),save::ind_grid_new,ind_cell_new,ind_part
  logical ,dimension(1:nvector),save::ok,ok_new=.true.
  integer ,dimension(1:ncpu)::ntot_star_cpu,ntot_star_all

  if(numbtot(1,ilevel)==0) return
  if(.not. hydro)return

  if(verbose)write(*,*)' Entering count test particle'

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Clump density threshold from H/cc to code units
  d0 = density_threshold/scale_nH

  !------------------------------------------------
  ! Compute number of new test particles in the level ilevel
  !------------------------------------------------
  ! Loop over grids
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     ! Test particle formation ---> logical array ok(i)
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        ! Flag leaf cells
        do i=1,ngrid
           ok(i)=son(ind_cell(i))==0
        end do
        ! Density criterion
        do i=1,ngrid
           d=uold(ind_cell(i),1)
           if(d<=d0)ok(i)=.false. 
        end do
        ! Compute test particle map
        do i=1,ngrid
           flag2(ind_cell(i))=0
           if(ok(i))then
              flag2(ind_cell(i))=1 
              ntot=ntot+1
           endif
        end do
     end do
  end do

end subroutine count_test_particle
!################################################################
!################################################################
!################################################################
!################################################################
subroutine create_test_particle(ilevel,ntot,nskip)
  use amr_commons
  use pm_commons
  use hydro_commons
  use random
  use clfind_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  integer::ntot,nskip
  !----------------------------------------------------------------------
  ! Compute test particle properties
  !----------------------------------------------------------------------
  ! local constants
  real(dp)::d0
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:twotondim,1:3)::xc
  ! other variables
  integer ::ncache,nnew,ngrid,icpu,index_star
  integer ::igrid,ix,iy,iz,ind,i,iskip,nx_loc
  integer ::ntot_all,info
  logical ::ok_free
  real(dp),dimension(1:3)::skip_loc
  real(dp)::d,x,y,z,dx,dx_loc,scale,vol_loc,dx_min,vol_min
  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  integer ,dimension(1:nvector),save::ind_grid_new,ind_cell_new,ind_part
  logical ,dimension(1:nvector),save::ok,ok_new=.true.
  integer ,dimension(1:ncpu)::ntot_star_cpu,ntot_star_all

  if(numbtot(1,ilevel)==0) return
  if(.not. hydro)return

  if(verbose)write(*,*)' Entering test particle creation'

  ! Loop over grids
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     ! Loop over cells
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        ! Flag cells with test particle
        do i=1,ngrid
           ok(i)=flag2(ind_cell(i))>0
        end do
        ! Calculate new test particle positions
        do i=1,ngrid
           if (ok(i))then
              ntot=ntot+1                    ! Local test particle index
              levp(ntot)=ilevel              ! Level
              iglobalp(ntot)=ntot+nskip      ! Global test particle index
              flag2(ind_cell(i))=ntot+nskip  ! Initialize flag2 to GLOBAL test particle index
              icellp(ntot)=ind_cell(i)       ! Local cell index
              denp(ntot)=uold(ind_cell(i),1) ! Save density values here!
           end if
        end do
        ! End loop over new test particles
     end do
     ! End loop over cells
  end do
  ! End loop over grids

end subroutine create_test_particle
!################################################################
!################################################################
!################################################################
!################################################################
subroutine scan_for_peaks(npartt,nmove,counter,action)
  use amr_commons
  use pm_commons
  use clfind_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h' 
#endif
  integer::npartt,ilevel,nmove,action,counter
  !----------------------------------------------------------------------
  ! Assign particle/cell on ilevel to the densest peak
  !----------------------------------------------------------------------
  integer::nv
  integer::igrid,jgrid,ipart,jpart,next_part,ig,ip,npart1,npartmin
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part,indv

  do ipart=1,npartt
     nv=1
     ilevel=levp(testp_sort(ipart)) ! level
     indv(nv)=(icellp(testp_sort(ipart))-ncoarse-1)/ngridmax+1 ! cell position
     ind_grid(nv)=icellp(testp_sort(ipart))-ncoarse-(indv(nv)-1)*ngridmax ! grid index
     ind_part(nv)=testp_sort(ipart)
     ind_grid_part(nv)=1
     ig=1
     ip=1
     call flag_peak(indv,ind_grid,ind_part,ind_grid_part,ig,ip,nmove,ilevel,counter,action)
  end do

  if(verbose)write(*,*)'   Exiting scan_for_peaks',nmove

end subroutine scan_for_peaks
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine flag_peak(indv,ind_grid,ind_part,ind_grid_part,ng,np,nm,ilevel,counter,action)
  use amr_commons
  use pm_commons
  use poisson_commons
  use clfind_commons, ONLY: icellp,npeaks
  use hydro_commons, ONLY: uold
  implicit none
  integer::ng,np,nm,ilevel,counter,action
  integer,dimension(1:nvector)::ind_grid,indv
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !------------------------------------------------------------
  ! This routine moves the particles in the arrays of length 
  ! nvector one step to the densest neighbor. It returns the
  ! number of particles which have effectively moved.
  !------------------------------------------------------------
  integer::nv=1
  logical::error
  logical,dimension(1:nvector)::okpeak
  integer::i,j,ind,idim,nx_loc,i1,j1,k1,i2,j2,k2,i3,j3,k3,ix,iy,iz
  real(dp)::dx,dx_loc,scale,vol_loc
  integer::i1min,i1max,j1min,j1max,k1min,k1max
  integer::i2min,i2max,j2min,j2max,k2min,k2max
  integer::i3min,i3max,j3min,j3max,k3min,k3max
  real(dp),dimension(1:twotondim,1:3)::xc
  ! Grid-based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
 ! Particle-based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x,xtest,xmax
  real(dp),dimension(1:nvector),save::density_max
  integer ,dimension(1:nvector,1:ndim),save::ig,id
  integer ,dimension(1:nvector),save::cell_index,cell_levl,ind_max
  real(dp),dimension(1:nvector,1:ndim,1:twotondim),save::xpart
  real(dp),dimension(1:3)::skip_loc

  okpeak=.true.

  ! Meshspacing in that level
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
  
  !==============================
  ! Get particle density and cell
  !==============================
  do j=1,np
     xtest(j,1)=(xg(ind_grid(j),1)+xc(indv(j),1)-skip_loc(1))*scale
     xtest(j,2)=(xg(ind_grid(j),2)+xc(indv(j),2)-skip_loc(2))*scale
     xtest(j,3)=(xg(ind_grid(j),3)+xc(indv(j),3)-skip_loc(3))*scale
  end do
  call get_cell_index(cell_index,cell_levl,xtest,ilevel,np)
  do j=1,np
     density_max(j)=uold(cell_index(j),1)*1.0001
     ind_max(j)=cell_index(j)
  end do

  !====================================================
  ! Check for potential new positions at level ilevel-1
  !====================================================
  if(ilevel>levelmin)then
     ! Generate 2x2x2 neighboring cells at level ilevel-1
     do k1=k1min,k1max
        do j1=j1min,j1max
           do i1=i1min,i1max
              do j=1,np
                 xtest(j,1)=(xg(ind_grid(j),1)+2*xc(indv(j),1)-skip_loc(1))*scale+(2*i1-1)*dx_loc
                 xtest(j,2)=(xg(ind_grid(j),2)+2*xc(indv(j),2)-skip_loc(2))*scale+(2*j1-1)*dx_loc
                 xtest(j,3)=(xg(ind_grid(j),3)+2*xc(indv(j),3)-skip_loc(3))*scale+(2*k1-1)*dx_loc
              end do
              call get_cell_index(cell_index,cell_levl,xtest,ilevel,np)
              do j=1,np
                 if(son(cell_index(j))==0.and.cell_levl(j)==(ilevel-1))then
                    if(uold(cell_index(j),1)>density_max(j))then
                       okpeak(j)=.false.
                       density_max(j)=uold(cell_index(j),1)
                       ind_max(j)=cell_index(j)
                    endif
                 endif
              end do
           end do
        end do
     end do
  endif

  !====================================================
  ! Check for potential new positions at level ilevel
  !====================================================
  ! Generate 3x3x3 neighboring cells at level ilevel
  do k2=k2min,k2max
     do j2=j2min,j2max
        do i2=i2min,i2max
           do j=1,np
              xtest(j,1)=(xg(ind_grid(j),1)+xc(indv(j),1)-skip_loc(1))*scale+(i2-1)*dx_loc
              xtest(j,2)=(xg(ind_grid(j),2)+xc(indv(j),2)-skip_loc(2))*scale+(j2-1)*dx_loc
              xtest(j,3)=(xg(ind_grid(j),3)+xc(indv(j),3)-skip_loc(3))*scale+(k2-1)*dx_loc
           end do
           call get_cell_index(cell_index,cell_levl,xtest,ilevel,np)
           do j=1,np
              if(son(cell_index(j))==0.and.cell_levl(j)==ilevel)then
                 if(uold(cell_index(j),1)>density_max(j))then
                    okpeak(j)=.false.
                    density_max(j)=uold(cell_index(j),1)
                    ind_max(j)=cell_index(j)
                 endif
              endif
           end do
        end do
     end do
  end do

  !====================================================
  ! Check for potential new positions at level ilevel+1
  !====================================================
  if(ilevel<nlevelmax)then
     ! Generate 4x4x4 neighboring cells at level ilevel+1
     do k3=k3min,k3max
        do j3=j3min,j3max
           do i3=i3min,i3max
              do j=1,np
                 xtest(j,1)=(xg(ind_grid(j),1)+xc(indv(j),1)-skip_loc(1))*scale+(i3-1.5)*dx_loc/2.0
                 xtest(j,2)=(xg(ind_grid(j),3)+xc(indv(j),2)-skip_loc(2))*scale+(j3-1.5)*dx_loc/2.0
                 xtest(j,3)=(xg(ind_grid(j),3)+xc(indv(j),3)-skip_loc(3))*scale+(k3-1.5)*dx_loc/2.0
              end do
              call get_cell_index(cell_index,cell_levl,xtest,ilevel+1,np)
              do j=1,np
                 if(son(cell_index(j))==0.and.cell_levl(j)==(ilevel+1))then
                    if(uold(cell_index(j),1)>density_max(j))then
                       okpeak(j)=.false.
                       density_max(j)=uold(cell_index(j),1)
                       ind_max(j)=cell_index(j)
                    endif
                 endif
              end do
           end do
        end do
     end do
  endif

  select case (action)
  case (1)
     !====================================================
     !  Count peaks
     !====================================================
     do j=1,np
        if(okpeak(j))counter=counter+1
     end do
  case (2)
     !====================================================
     !  Initialize flag2 to peak global index
     !====================================================
     do j=1,np
        if(okpeak(j))then 
           counter=counter+1
           flag2(icellp(ind_part(j)))=counter
        end if
     end do     
  case (3)
     !====================================================
     !  Propagate flag2
     !====================================================
     do j=1,np
        if(flag2(icellp(ind_part(j))).ne.flag2(ind_max(j)))nm=nm+1
        flag2(icellp(ind_part(j)))=flag2(ind_max(j))
     end do
  end select

end subroutine flag_peak
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine assign_part_to_peak(ntest,peak_nr)
  use amr_commons
  use hydro_commons
  use pm_commons
  use clfind_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer ntest,peak_nr
  !----------------------------------------------------------------------------
  ! This subroutine loops over all particles and marks every cell containing a
  ! particle (flag 2). Every cell containig a particle is a peak.
  ! The number of peaks on each cpu is counted. Using MPI communication, a GLOBAL
  ! peak index is given to each peak. 
  ! In a second loop over all particles, the peak index a particle belongs to,
  ! is written into the mass variable of each particle.
  !----------------------------------------------------------------------------

  integer::igrid,jgrid,ipart,jpart,next_part,ig,ip,ilevel,npart1,info
  integer::jj,icpu,peak_loc
  integer::n_cls,nv

  integer::minf2,maxf2,cum
  integer,dimension(1:nvector)::ind_grid,ind_cell,init_ind_cell,init_cell_lev,cell_lev
  integer,dimension(1:nvector)::ind_part
  integer,dimension(:),allocatable::flip
  real(dp),dimension(1:nvector,1:ndim)::pos,init_pos
  real(kind=8),allocatable,dimension(:,:)::peak_pos
  integer,dimension(1:ncpu)::npeaks_per_cpu,npeaks_per_cpu_tot

  integer ::ix,iy,iz,ind,i,iskip,nx_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:3)::skip_loc
  real(dp)::d,x,y,z,dx,dx_loc,scale,vol_loc,dx_min,vol_min

  ! Set local constants
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_min=(0.5D0**nlevelmax)*scale
  vol_min=dx_min**ndim
  ! Cells center position relative to grid center position
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     xc(ind,1)=(dble(ix)-0.5D0)
     xc(ind,2)=(dble(iy)-0.5D0)
     xc(ind,3)=(dble(iz)-0.5D0)
  end do

  !----------------------------------------------------------------------------
  ! allocate arrays where the postiton of the peaks is stored
  !----------------------------------------------------------------------------
  allocate(peak_pos_tot(1:npeaks_tot,1:ndim))
  allocate(peak_pos(1:npeaks_tot,1:ndim))
  peak_pos=0.
  
  !---------------------------------------------------------------------------- 
  ! determine peak-ids positions for each cpu
  !----------------------------------------------------------------------------
  peak_nr=peak_nr+1
  peak_loc=1

  do i=1,ntest
     if(flag2(icellp(testp_sort(i)))==peak_nr)then
        ilevel=levp(testp_sort(i)) ! cell level
        dx=0.5D0**ilevel ! mesh spacing at that level
        ind=(icellp(testp_sort(i))-ncoarse-1)/ngridmax+1 ! cell position
        igrid=icellp(testp_sort(i))-ncoarse-(ind-1)*ngridmax ! grid index 
        peak_pos(peak_nr,1)=(xg(igrid,1)+xc(ind,1)*dx-skip_loc(1))*scale
        peak_pos(peak_nr,2)=(xg(igrid,2)+xc(ind,2)*dx-skip_loc(2))*scale
        peak_pos(peak_nr,3)=(xg(igrid,3)+xc(ind,3)*dx-skip_loc(3))*scale
        ! jump to next peak
        peak_loc=peak_loc+1
        peak_nr=peak_nr+1
     end if
  end do
  ! reset peak skip
  peak_nr=peak_nr-peak_loc

  !----------------------------------------------------------------------------
  ! create global list of peak positions
  !----------------------------------------------------------------------------
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(peak_pos,peak_pos_tot,3*npeaks_tot,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  peak_pos_tot=peak_pos
#endif

  deallocate(peak_pos) ! from here on only peak_pos_tot is used

end subroutine assign_part_to_peak
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine get_cell_index(cell_index,cell_levl,xpart,ilevel,np)
  use amr_commons
  implicit none

  integer::np,ilevel
  integer,dimension(1:nvector)::cell_index,cell_levl
  real(dp),dimension(1:nvector,1:3)::xpart

  !----------------------------------------------------------------------------
  ! This routine returns the index of the cell, at maximum level
  ! ilevel, in which the input particle sits
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
  do i=1,np
     xx = xpart(i,1)/boxlen + (nx-1)/2.0
     yy = xpart(i,2)/boxlen + (ny-1)/2.0
     zz = xpart(i,3)/boxlen + (nz-1)/2.0
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
  use amr_commons
  use hydro_commons

  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  
  !--------------------------------------------------                           
  ! Namelist definitions                                                        
  !--------------------------------------------------                           
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  namelist/clumpfind_params/relevance_threshold,density_threshold,mass_threshold

  ! Read namelist file 
  rewind(1)
  read(1,NML=clumpfind_params,END=101)
  goto 102
101 write(*,*)' You need to set up namelist &CLUMPFIND_PARAMS in parameter file'
  call clean_stop
102 rewind(1)
  
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  !convert mass_threshold from solar masses to grams
  !...and to user units
  mass_threshold=mass_threshold*1.98892d33 /(scale_l**3. * scale_d)
end subroutine read_clumpfind_params
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
