! Ifront_writer patch (see !IF comments): 
!        Added subroutine write_ifront to find ifront position in a
!        stromgren sphere/shadow experiment and write it to std output.
! RT patch: Plenty of change here. To see it, do a diff.
!*************************************************************************
subroutine adaptive_loop
  use amr_commons
  use hydro_commons
  use pm_commons
  use poisson_commons
  use cooling_module
#ifdef RT
  use rt_hydro_commons
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer(kind=8)::n_step
  integer::ilevel,idim,ivar,info,tot_pt
  real(kind=8)::tt1,tt2,muspt,muspt_this_step
  real(kind=4)::real_mem,real_mem_tot

#ifndef WITHOUTMPI
  tt1=MPI_WTIME()
#endif

  call init_amr                      ! Initialize AMR variables
  call init_time                     ! Initialize time variables
  if(hydro)call init_hydro           ! Initialize hydro variables
#ifdef RT
  if(rt.or.neq_chem) &
       & call rt_init_hydro          ! Initialize radiation variables
#endif
  if(poisson)call init_poisson       ! Initialize poisson variables
#ifdef ATON
  if(aton)call init_radiation        ! Initialize radiation variables
#endif
  if(nrestart==0)call init_refine    ! Build initial AMR grid

#ifdef grackle
  if(use_grackle==0)then
     if(cooling.and..not.neq_chem) &
        call set_table(dble(aexp))    ! Initialize cooling look up table
  endif
#else  
  if(cooling.and..not.neq_chem) &
       call set_table(dble(aexp))    ! Initialize cooling look up table
#endif
  if(pic)call init_part              ! Initialize particle variables
  if(pic)call init_tree              ! Initialize particle tree
  if(nrestart==0)call init_refine_2  ! Build initial AMR grid again

#ifndef WITHOUTMPI
  muspt=0.
  tot_pt=-1
  tt2=MPI_WTIME()
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
                               call timer('coarse levels','start')

#ifndef WITHOUTMPI
     tt1=MPI_WTIME()
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
     if(levelmin.lt.nlevelmax.and.(.not.static.or.(nstep_coarse_old.eq.nstep_coarse.and.restart_remap)))then
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

     !call write_ifront_stromgren !-------------------------------------!IF !sln use this
     !call write_ifront_shadow !----------------------------------------!IF
 
    ! Call base level
     call amr_step(levelmin,1)
     call timer('coarse levels','start')

     if(levelmin.lt.nlevelmax.and.(.not.static.or.(nstep_coarse_old.eq.nstep_coarse.and.restart_remap)))then
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

#ifndef WITHOUTMPI
     tt2=MPI_WTIME()
     if(mod(nstep_coarse,ncontrol)==0)then
        call getmem(real_mem)
        call MPI_ALLREDUCE(real_mem,real_mem_tot,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,info)
        if(myid==1)then
           if (tot_pt==0) muspt=0. ! dont count first timestep
           n_step = int(numbtot(1,levelmin),kind=8)*twotondim
           do ilevel=levelmin+1,nlevelmax
             n_step = n_step + int(numbtot(1,ilevel),kind=8)*product(nsubcycle(levelmin:ilevel-1))*(twotondim-1)
           enddo
           muspt_this_step = (tt2-tt1)*1e6/n_step*ncpu
           muspt = muspt + muspt_this_step
           tot_pt = tot_pt + 1
           write(*,'(a,f8.2,a,f12.2,a,f12.2,a)')' Time elapsed since last coarse step:',tt2-tt1 &
          ,' s',muspt_this_step,' mus/pt'  &
          ,muspt / max(tot_pt,1), ' mus/pt (av)'
           call writemem(real_mem_tot)
        endif
     endif
#endif

  end do

999 format(' Level ',I2,' has ',I10,' grids (',3(I8,','),')')

end subroutine adaptive_loop

!************************************************************************
SUBROUTINE write_ifront_stromgren

! Find and write to standard output the position of an ionization front
! from the x=y=z=0 origin of the box.
!-------------------------------------------------------------------------
  use amr_commons
  use hydro_commons
  use rt_parameters
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::igrid,jgrid,ind,info
  integer::i,j,ilevel,ncache,ix,iy,iz
  integer::nx_loc,nbins_loc,j0
  real(dp)::dx,scale
  real(dp),dimension(1:twotondim,1:3)::xc
  integer,dimension(:),allocatable::ind_grid,ind_cell
  real(kind=8)::xion,rr
  real(kind=8),dimension(:,:),allocatable::bins,bins_all  !  r,r_avg,x_avg
  integer,dimension(:),allocatable::bin_count,bin_count_all
  real(dp)::ifront_rad
!-------------------------------------------------------------------------

#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
  rr=0.0D0; xion=0.0D0
  nbins_loc=100
  allocate(bin_count(nbins_loc), bin_count_all(nbins_loc))  ; bin_count=0
  allocate(bins(3,nbins_loc)   , bins_all(3,nbins_loc))     ; bins=0
  do i=1,nbins_loc !initialize the bins
     bins(1,i) = (i-1) * boxlen/(nbins_loc) ! Lower radial values for bins
  end do

  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  do ilevel=1,nlevelmax
     ncache=numbl(myid,ilevel)
     if(ncache > 0)then
        dx=0.5D0**ilevel
        allocate(ind_grid(1:ncache),ind_cell(1:ncache))
        ! Gather all grids
        igrid=headl(myid,ilevel)
        do jgrid=1,ncache
           ind_grid(jgrid)=igrid
           igrid=next(igrid)
        end do
        ! Gather variables
        do ind=1,twotondim
           iz=(ind-1)/4
           iy=(ind-1-4*iz)/2
           ix=(ind-1-2*iy-4*iz)
           xc(ind,1)=(dble(ix)-0.5D0)*dx
           xc(ind,2)=(dble(iy)-0.5D0)*dx
           xc(ind,3)=(dble(iz)-0.5D0)*dx
           do i=1,ncache
              ind_cell(i)=ncoarse+(ind-1)*ngridmax+ind_grid(i)
           end do
           do i=1,ncache
              if(son(ind_cell(i))==0)then
                 ! Now need to bin the cell's ionization state according 
                 ! to it's radius
                 rr = sqrt(   &
                      ((xg(ind_grid(i),1)+xc(ind,1)-dble(icoarse_min))   &
                       *scale)**2 &
                      +((xg(ind_grid(i),2)+xc(ind,2)-dble(jcoarse_min))  &
                       *scale)**2 &
                      +((xg(ind_grid(i),3)+xc(ind,3)-dble(kcoarse_min))  &
                       *scale)**2 )
                 xion = uold(ind_cell(i),iIons)/uold(ind_cell(i),1)
                 do j=1,nbins_loc
                    if(bins(1,j) .ge. rr) then
                       bins(2,j) = bins(2,j) + rr
                       bins(3,j) = bins(3,j) + xion
                       bin_count(j)=bin_count(j)+1
                       exit
                    endif
                 end do

              end if
           end do
        end do
        deallocate(ind_grid, ind_cell)
     end if
  end do

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(bins,      bins_all,      nbins_loc*3, &
                     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(bin_count, bin_count_all, nbins_loc, &
                     MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, info)
#endif

  if(myid==1) then
     bins(2:3,:)=bins_all(2:3,:); bin_count=bin_count_all
     bins(2,:)=bins(2,:)/bin_count ! Average radius per radius bin
     bins(3,:)=bins(3,:)/bin_count ! Average xion per radius bin
     ifront_rad=0.
     j0=0
     do j=1,nbins_loc
        if(bin_count(j) .gt. 0) then
           if(bins(3,j) .le. 0.5) then
              exit ! Found the front. j points at first bin w x<=0.5
           endif
           j0=j ! Adjacent nonempty bin to the left for interpolation
        endif
     end do
     ! Interpolate
     if(j0 .eq. 0) then
        ifront_rad = 0.
     else
        ifront_rad = bins(2,j0) + (0.5d0-bins(3,j0))                     &
                         / (bins(3,j)-bins(3,j0)) * (bins(2,j)-bins(2,j0))
     endif
     
     write(*, 111) t, ifront_rad

  end if

  ! Deallocate local arrays
  deallocate(bin_count,bins)
  deallocate(bin_count_all,bins_all)

 
#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif

111 format('IFRONT     ', f21.6, f21.6)

END SUBROUTINE write_ifront_stromgren

!*************************************************************************
SUBROUTINE write_ifront_shadow
! Find and write to standard output the position of an ionization front
! along the x-axis at the center of the yz-plane, i.e. at the center of 
! the box.
!-------------------------------------------------------------------------
  use amr_commons
  use hydro_commons
  use rt_parameters
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::igrid,jgrid,ind,info
  integer::i,j,k,ilevel,ncache,ix,iy,iz
  integer::nx_loc,nbins_loc,j0
  real(dp)::dx,scale
  real(dp),dimension(1:twotondim,1:3)::xc
  integer,dimension(:),allocatable::ind_grid,ind_cell
  real(kind=8)::xHII, rx, ry, rz, dr_max
  real(kind=8),dimension(:,:),allocatable::bins,bins_all!x,x_avg,xHII_avg
  integer,dimension(:),allocatable::bin_count,bin_count_all
  real(dp)::ifront_x
!-------------------------------------------------------------------------
#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
  rx=0.0D0; xHII=0.0D0
  nbins_loc=200
  allocate(bin_count(nbins_loc), bin_count_all(nbins_loc))  ; bin_count=0
  allocate(bins(3,nbins_loc)   , bins_all(3,nbins_loc))     ; bins=0
  do i=1,nbins_loc !initialize the bins
     bins(1,i) = (i-1) * boxlen/(nbins_loc) ! Lower x values for bins
  end do

  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  ! width of a cell
  !dr_max = sqrt(2.)*scale/2**(levelmin)                   ! Use this in Iliev3
  dr_max = sqrt(2.)*scale/2**(nlevelmax) ! Careful here  ! Use this in Iliev7
  do ilevel=1,nlevelmax
     ncache=numbl(myid,ilevel)
     if(ncache > 0)then
        dx=0.5D0**ilevel
        allocate(ind_grid(1:ncache),ind_cell(1:ncache))
        ! Gather all grids
        igrid=headl(myid,ilevel)
        do jgrid=1,ncache
           ind_grid(jgrid)=igrid
           igrid=next(igrid)
        end do
        ! Gather variables
        do ind=1,twotondim
           iz=(ind-1)/4
           iy=(ind-1-4*iz)/2
           ix=(ind-1-2*iy-4*iz)
           xc(ind,1)=(dble(ix)-0.5D0)*dx
           xc(ind,2)=(dble(iy)-0.5D0)*dx
           xc(ind,3)=(dble(iz)-0.5D0)*dx
           do i=1,ncache
              ind_cell(i)=ncoarse+(ind-1)*ngridmax+ind_grid(i)
           end do
           do i=1,ncache
              if(son(ind_cell(i))==0)then ! only bin leaf cells
                 ! bin the cell's ionization state according to it's x-pos
                 rx= (xg(ind_grid(i),1)+xc(ind,1)-dble(icoarse_min))*scale
                 ry= (xg(ind_grid(i),2)+xc(ind,2)-dble(jcoarse_min))*scale
                 rz= (xg(ind_grid(i),3)+xc(ind,3)-dble(kcoarse_min))*scale
                 if(sqrt((ry-boxlen/2)**2+(rz-boxlen/2)**2) .gt. dr_max) &
                      cycle ! Only consider cells along x in box middle
                 xHII = uold(ind_cell(i),iIons)/uold(ind_cell(i),1)
                 do j=1,nbins_loc
                    if(bins(1,j) .ge. rx) then
                       bins(2,j) = bins(2,j) + rx
                       bins(3,j) = bins(3,j) + xHII
                       bin_count(j)=bin_count(j)+1
                       exit
                    endif
                 end do
              end if ! End leaf cells
           end do ! End ncache loop
        end do ! End twotondim loop
        deallocate(ind_grid, ind_cell)
     end if ! End ncache>0?
  end do ! End ilevel loop

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(bins,      bins_all,      nbins_loc*3,                  &
                     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(bin_count, bin_count_all, nbins_loc,                    &
                     MPI_INTEGER,          MPI_SUM, MPI_COMM_WORLD, info)
#endif

  if(myid==1) then
     bins(2:3,:)=bins_all(2:3,:); bin_count=bin_count_all
     bins(2,:)=bins(2,:)/bin_count ! Average x per x bin
     bins(3,:)=bins(3,:)/bin_count ! Average xHII per x bin
     ifront_x=0.
     j0=0
     do j=1,nbins_loc
        if(bin_count(j) .gt. 0) then
           if(bins(3,j) .le. 0.5) then
              exit ! Found the front. j points at first bin w xHII<=0.5
           endif
           j0=j ! Adjacent nonempty bin to the left for interpolation
        endif
     end do
     ! Interpolate
     if(j0 .eq. 0) then
        ifront_x = 0.
     else
        ifront_x = bins(2,j0) + (0.5d0-bins(3,j0))                       &
                         / (bins(3,j)-bins(3,j0)) * (bins(2,j)-bins(2,j0))
     endif
    
     write(*, 111) t, ifront_x
  end if

  ! Deallocate local arrays
  deallocate(bin_count,bins)
  deallocate(bin_count_all,bins_all)

 
#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif

111 format('IFRONT     ', f21.6, f21.6)

END SUBROUTINE write_ifront_shadow
