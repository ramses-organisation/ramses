!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine force_fine(ilevel,icount)
  use amr_commons
  use pm_commons
  use poisson_commons
  use ffv
  use ffv_subroutines
  implicit none
#ifndef WITHOUTMPI
#ifndef INCLUDEOK
#define INCLUDEOK
  include 'mpif.h'
#endif
#endif
  integer::ilevel,icount
  !----------------------------------------------------------
  ! This routine computes the gravitational acceleration,
  ! the maximum density rho_max, and the potential energy
  !----------------------------------------------------------
  integer::igrid,ngrid,ncache,i,ind,iskip,ix,iy,iz
  integer::info,ibound,nx_loc,idim
  real(dp)::dx,dx_loc,scale,fact,fourpi
  real(kind=8)::rho_loc,rho_all,epot_loc,epot_all
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:3)::skip_loc

  integer ,dimension(1:nvector),save::ind_grid,ind_cell,ind_cell_father
  real(dp),dimension(1:nvector,1:ndim),save::xx,ff

  integer::nbin,ibin
  integer,dimension(:,:,:),allocatable::ix_aux,iy_aux,iz_aux,iix_aux,iiy_aux,iiz_aux
  real(dp)::a,b
#ifdef MPROFILING
#ifndef INCLUDEOK
#define INCLUDEOK
  include 'mpif.h'
#endif
  integer::info
  real(kind=8)::tt1,tt2
  tt1=MPI_WTIME(info)
#endif
#undef INCLUDEOK
 
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Mesh size at level ilevel in coarse cell units
  dx=0.5D0**ilevel
  
  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  a=0.50D0*4.0D0/3.0D0/dx
  b=0.25D0*1.0D0/3.0D0/dx
  !   |dim
  !   | |node
  !   | | |cell
  !   v v v
  ggg(1,1,1:8)=(/1,0,1,0,1,0,1,0/); hhh(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,2,1:8)=(/0,2,0,2,0,2,0,2/); hhh(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  ggg(1,3,1:8)=(/1,1,1,1,1,1,1,1/); hhh(1,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(1,4,1:8)=(/2,2,2,2,2,2,2,2/); hhh(1,4,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(2,1,1:8)=(/3,3,0,0,3,3,0,0/); hhh(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,2,1:8)=(/0,0,4,4,0,0,4,4/); hhh(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  ggg(2,3,1:8)=(/3,3,3,3,3,3,3,3/); hhh(2,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(2,4,1:8)=(/4,4,4,4,4,4,4,4/); hhh(2,4,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(3,1,1:8)=(/5,5,5,5,0,0,0,0/); hhh(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,2,1:8)=(/0,0,0,0,6,6,6,6/); hhh(3,2,1:8)=(/5,6,7,8,1,2,3,4/)
  ggg(3,3,1:8)=(/5,5,5,5,5,5,5,5/); hhh(3,3,1:8)=(/1,2,3,4,5,6,7,8/)
  ggg(3,4,1:8)=(/6,6,6,6,6,6,6,6/); hhh(3,4,1:8)=(/1,2,3,4,5,6,7,8/)

  ! CIC method constants
  aa = 1.0D0/4.0D0**ndim
  bb = 3.0D0*aa
  cc = 9.0D0*aa
  dd = 27.D0*aa
  bbbb(:)  =(/aa ,bb ,bb ,cc ,bb ,cc ,cc ,dd/)

  !sampling positions in the 3x3x3 father cell cube
  ccc(:,1)=(/1 ,2 ,4 ,5 ,10,11,13,14/)
  ccc(:,2)=(/3 ,2 ,6 ,5 ,12,11,15,14/)
  ccc(:,3)=(/7 ,8 ,4 ,5 ,16,17,13,14/)
  ccc(:,4)=(/9 ,8 ,6 ,5 ,18,17,15,14/)
  ccc(:,5)=(/19,20,22,23,10,11,13,14/)
  ccc(:,6)=(/21,20,24,23,12,11,15,14/)
  ccc(:,7)=(/25,26,22,23,16,17,13,14/)
  ccc(:,8)=(/27,26,24,23,18,17,15,14/)

  iii_aux=(/1,2,1,2,1,2,1,2/)
  jjj_aux=(/3,3,4,4,3,3,4,4/)
  kkk_aux=(/5,5,5,5,6,6,6,6/)
  lll=0; mmm=0
  lll(1:3,1,1)=(/2,1,1/)
  mmm(1:3,1,1)=(/2,1,2/)
  lll(1:3,2,1)=(/1,1,2/)
  mmm(1:3,2,1)=(/1,2,1/)
  lll(1:9,1,2)=(/4,3,3,2,1,1,2,1,1/)
  mmm(1:9,1,2)=(/4,3,4,2,1,2,4,3,4/)
  lll(1:9,2,2)=(/3,3,4,1,1,2,1,1,2/)
  mmm(1:9,2,2)=(/3,4,3,1,2,1,3,4,3/)
  lll(1:9,3,2)=(/2,1,1,2,1,1,4,3,3/)
  mmm(1:9,3,2)=(/2,1,2,4,3,4,2,1,2/)
  lll(1:9,4,2)=(/1,1,2,1,1,2,3,3,4/)
  mmm(1:9,4,2)=(/1,2,1,3,4,3,1,2,1/)
  lll(1:27,1,3)=(/8,7,7,6,5,5,6,5,5,4,3,3,2,1,1,2,1,1,4,3,3,2,1,1,2,1,1/)
  mmm(1:27,1,3)=(/8,7,8,6,5,6,8,7,8,4,3,4,2,1,2,4,3,4,8,7,8,6,5,6,8,7,8/)
  lll(1:27,2,3)=(/7,7,8,5,5,6,5,5,6,3,3,4,1,1,2,1,1,2,3,3,4,1,1,2,1,1,2/)
  mmm(1:27,2,3)=(/7,8,7,5,6,5,7,8,7,3,4,3,1,2,1,3,4,3,7,8,7,5,6,5,7,8,7/)
  lll(1:27,3,3)=(/6,5,5,6,5,5,8,7,7,2,1,1,2,1,1,4,3,3,2,1,1,2,1,1,4,3,3/)
  mmm(1:27,3,3)=(/6,5,6,8,7,8,6,5,6,2,1,2,4,3,4,2,1,2,6,5,6,8,7,8,6,5,6/)
  lll(1:27,4,3)=(/5,5,6,5,5,6,7,7,8,1,1,2,1,1,2,3,3,4,1,1,2,1,1,2,3,3,4/)
  mmm(1:27,4,3)=(/5,6,5,7,8,7,5,6,5,1,2,1,3,4,3,1,2,1,5,6,5,7,8,7,5,6,5/)
  lll(1:27,5,3)=(/4,3,3,2,1,1,2,1,1,4,3,3,2,1,1,2,1,1,8,7,7,6,5,5,6,5,5/)
  mmm(1:27,5,3)=(/4,3,4,2,1,2,4,3,4,8,7,8,6,5,6,8,7,8,4,3,4,2,1,2,4,3,4/)
  lll(1:27,6,3)=(/3,3,4,1,1,2,1,1,2,3,3,4,1,1,2,1,1,2,7,7,8,5,5,6,5,5,6/)
  mmm(1:27,6,3)=(/3,4,3,1,2,1,3,4,3,7,8,7,5,6,5,7,8,7,3,4,3,1,2,1,3,4,3/)
  lll(1:27,7,3)=(/2,1,1,2,1,1,4,3,3,2,1,1,2,1,1,4,3,3,6,5,5,6,5,5,8,7,7/)
  mmm(1:27,7,3)=(/2,1,2,4,3,4,2,1,2,6,5,6,8,7,8,6,5,6,2,1,2,4,3,4,2,1,2/)
  lll(1:27,8,3)=(/1,1,2,1,1,2,3,3,4,1,1,2,1,1,2,3,3,4,5,5,6,5,5,6,7,7,8/)
  mmm(1:27,8,3)=(/1,2,1,3,4,3,1,2,1,5,6,5,7,8,7,5,6,5,1,2,1,3,4,3,1,2,1/)

  ncache=active(ilevel)%ngrid
  nbin = int(ncache/nvector)+1

  allocate(ngrid_aux(1:nbin),ind_grid_aux(1:nbin,1:nvector))
  allocate(ind_cell_aux(1:nbin,1:nvector))
  allocate(ind_left(1:nbin,1:nvector,1:ndim),ind_right(1:nbin,1:nvector,1:ndim))
  allocate(igridn(1:nbin,1:nvector,0:twondim))
  allocate(phi1(1:nbin,1:nvector),phi2(1:nbin,1:nvector), & 
  &        phi3(1:nbin,1:nvector),phi4(1:nbin,1:nvector))
  allocate(phi_left(1:nbin,1:nvector,1:twotondim,1:ndim), &
  &        phi_right(1:nbin,1:nvector,1:twotondim,1:ndim))
  allocate(nbors_father_grids(1:nbin,1:ndim,1:nvector,1:twotondim), &
  &        nbors_father_cells(1:nbin,1:ndim,1:nvector,1:threetondim))
  allocate(ix_aux(1:nbin,1:ndim,1:nvector),iy_aux(1:nbin,1:ndim,1:nvector), &
  & iz_aux(1:nbin,1:ndim,1:nvector),iix_aux(1:nbin,1:ndim,1:nvector), & 
  & iiy_aux(1:nbin,1:ndim,1:nvector),iiz_aux(1:nbin,1:ndim,1:nvector))
  allocate(pos(1:nbin,1:ndim,1:nvector), & 
  & ind_grid_father(1:nbin,1:ndim,1:nvector),ind_grid_ok(1:nbin,1:ndim,1:nvector))
  allocate(nbors_father_ok(1:nbin,1:ndim,1:nvector,1:threetondim))
  allocate(nbors_grids_ok(1:nbin,1:ndim,1:nvector,1:twotondim))
  allocate(ind_grid1(1:nbin,1:ndim,1:nvector),ind_grid2(1:nbin,1:ndim,1:nvector), &
  &        ind_grid3(1:nbin,1:ndim,1:nvector))
  allocate(nbors_grids(1:nbin,1:ndim,1:nvector,1:twotondim))

  phi_left = 0d0
  phi_right = 0d0
  
  !------------------------------------------------------------------------------------
  ! If in the future this subroutine is implemented by the !$acc routine directive,
  ! then this data region must be transformed to !$acc declare from !$acc data
  ! and be moved in the ffv module after the declaration part of the vars!
  !------------------------------------------------------------------------------------

  !$acc data create(ngrid_aux,ind_grid_aux,ind_cell_aux,ind_left,ind_right) &
  !$acc create(igridn,phi1,phi2,phi3,phi4,phi_left,phi_right,nbors_grids_ok) &
  !$acc create(nbors_father_grids,nbors_father_cells,nbors_father_ok) &
  !$acc create(pos,ind_grid_father,ind_grid_ok) &
  !$acc create(ind_grid1,ind_grid2,ind_grid3,nbors_grids) &
  !$acc create(ix_aux,iy_aux,iz_aux,iix_aux,iiy_aux,iiz_aux) &
  !$acc copyin(iii_aux,jjj_aux,kkk_aux,lll,mmm,dtnew,dtold,ggg,hhh,ccc,bbbb) &
  !$acc present(active,rho,f,phi,phi_old,son,nbor)
  
  !-------------------------------------
  ! Compute analytical gravity force
  !-------------------------------------
  if(gravity_type>0)then 

     ! Loop over myid grids by vector sweeps
     ncache=active(ilevel)%ngrid
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do
     
        ! Loop over cells
        do ind=1,twotondim

           ! Gather cell indices
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do
           ! Gather cell centre positions
           do idim=1,ndim
              do i=1,ngrid
                 xx(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
              end do
           end do
           ! Rescale position from code units to user units
           do idim=1,ndim
              do i=1,ngrid
                 xx(i,idim)=(xx(i,idim)-skip_loc(idim))*scale
              end do
           end do
           
           ! Call analytical gravity routine
           call gravana(xx,ff,dx_loc,ngrid)
           
           ! Scatter variables
           do idim=1,ndim
              do i=1,ngrid
                 f(ind_cell(i),idim)=ff(i,idim)
              end do
           end do

        end do
        ! End loop over cells

     end do
     ! End loop over grids

     ! Update boundaries
     do idim=1,ndim
        call make_virtual_fine_dp(f(1,idim),ilevel)
     end do
     if(simple_boundary)call make_boundary_force(ilevel)

  !------------------------------
  ! Compute gradient of potential
  !------------------------------
  else

     ! Update physical boundaries
     call make_boundary_phi(ilevel)

     ! Loop over myid grids by vector sweeps
     !$acc parallel loop gang worker independent
     do igrid=1,ncache,nvector
      ibin = int(igrid/nvector)+1
        ngrid_aux(ibin)=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid_aux(ibin)
           ind_grid_aux(ibin,i)=active(ilevel)%igrid(igrid+i-1)
        end do
        ! Compute gradient of potential
        call gradient_phi(ilevel,icount,nbin,ibin, & 
                        & ix_aux,iy_aux,iz_aux,iix_aux,iiy_aux,iiz_aux)
     end do
     ! End loop over grids
    
     if (sink)then
        call f_gas_sink(ilevel)
     end if
     
     ! Update boundaries
     do idim=1,ndim
        call make_virtual_fine_dp_acc(f(1,idim),ilevel)
     end do

     if(simple_boundary)call make_boundary_force(ilevel)

  endif

  !----------------------------------------------
  ! Compute gravity potential and maximum density
  !----------------------------------------------
  rho_loc =0.0; rho_all =0.0
  epot_loc=0.0; epot_all=0.0
  fourpi=4.0D0*ACOS(-1.0D0)
  if(cosmo)fourpi=1.5D0*omega_m*aexp
  fact=-dx_loc**ndim/fourpi/2.0D0

  !$acc data copy(epot_loc,rho_loc)

  ! Loop over myid grids by vector sweeps
  !$acc parallel loop gang vector reduction(+:epot_loc) reduction(max:rho_loc)
  do igrid=1,ncache,nvector    
     ! Loop over cells
     do ind=1,twotondim
        ibin = int(igrid/nvector)+1
        ngrid_aux(ibin)=MIN(nvector,ncache-igrid+1)
        ! Gather cell indices
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid_aux(ibin)
           ind_grid_aux(ibin,i)=active(ilevel)%igrid(igrid+i-1)
           ind_cell_aux(ibin,i)=iskip+ind_grid_aux(ibin,i)
        end do
        ! Loop over dimensions
        do idim=1,ndim
           do i=1,ngrid_aux(ibin)
              if(son(ind_cell_aux(ibin,i))==0)then
                 epot_loc=epot_loc+fact*f(ind_cell_aux(ibin,i),idim)**2
              end if
           end do
        end do
        ! End loop over dimensions
        do i=1,ngrid_aux(ibin)
           rho_loc=MAX(rho_loc,dble(abs(rho(ind_cell_aux(ibin,i)))))
        end do
     end do
     ! End loop over cells
  end do
  ! End loop over grids

  !$acc end data

#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(epot_loc,epot_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(rho_loc ,rho_all ,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
     epot_loc=epot_all
     rho_loc =rho_all
#endif
     epot_tot=epot_tot+epot_loc
     rho_max(ilevel)=rho_loc

     !$acc end data

     deallocate(ind_grid_aux,ngrid_aux)
     deallocate(ind_cell_aux)
     deallocate(ind_left,ind_right)
     deallocate(igridn)
     deallocate(phi1,phi2,phi3,phi4)
     deallocate(phi_left,phi_right)
     deallocate(nbors_father_grids,nbors_father_cells)
     deallocate(ix_aux,iy_aux,iz_aux,iix_aux,iiy_aux,iiz_aux)
     deallocate(pos,ind_grid_father,ind_grid_ok)
     deallocate(nbors_father_ok)
     deallocate(nbors_grids_ok)
     deallocate(ind_grid1,ind_grid2,ind_grid3)
     deallocate(nbors_grids)
     
#ifdef MPROFILING
  tt2=MPI_WTIME(info)
  acc_t_force_fine = acc_t_force_fine + (tt2-tt1)
#endif
  
111 format('   Entering force_fine for level ',I2)

end subroutine force_fine
!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine gradient_phi(ilevel,icount,nbin,ibin,ix_aux,iy_aux,iz_aux,iix_aux,iiy_aux,iiz_aux)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  use ffv
  use ffv_subroutines
  implicit none
  integer::ilevel,icount,nbin,ibin
  integer,dimension(1:nbin,1:ndim,1:nvector)::ix_aux,iy_aux,iz_aux,iix_aux,iiy_aux,iiz_aux
  !-------------------------------------------------
  ! This routine compute the 3-force for all cells
  ! in grids ind_grid(:) at level ilevel, using a
  ! 5 nodes kernel (5 points FDA).
  !-------------------------------------------------
  integer::i,idim,ind,iskip,nx_loc
  integer::id1,id2,id3,id4
  integer::ig1,ig2,ig3,ig4
  integer::ih1,ih2,ih3,ih4
  real(dp)::dx,a,b,scale,dx_loc

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel

  a=0.50D0*4.0D0/3.0D0/dx
  b=0.25D0*1.0D0/3.0D0/dx

  ! Rescaling factor
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  ! Gather neighboring grids
  do i=1,ngrid_aux(ibin)
     igridn(ibin,i,0)=ind_grid_aux(ibin,i)
  end do
  do idim=1,ndim
     do i=1,ngrid_aux(ibin)
        ind_left (ibin,i,idim)=nbor(ind_grid_aux(ibin,i),2*idim-1)
        ind_right(ibin,i,idim)=nbor(ind_grid_aux(ibin,i),2*idim  )
        igridn(ibin,i,2*idim-1)=son(ind_left (ibin,i,idim))
        igridn(ibin,i,2*idim  )=son(ind_right(ibin,i,idim))
     end do
  end do

  ! Interpolate potential from upper level
  if (ilevel>levelmin)then
     do idim=1,ndim
        call interpol_phi_ffv(ind_left ,phi_left ,ilevel,icount,ibin,nbin,idim,ix_aux,iy_aux,iz_aux,iix_aux,iiy_aux,iiz_aux)
        call interpol_phi_ffv(ind_right,phi_right,ilevel,icount,ibin,nbin,idim,ix_aux,iy_aux,iz_aux,iix_aux,iiy_aux,iiz_aux)
     end do
  end if

  ! Loop over cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid_aux(ibin)
        ind_cell_aux(ibin,i)=iskip+ind_grid_aux(ibin,i)
     end do

     ! Loop over dimensions
     do idim=1,ndim

        ! Loop over nodes
        id1=hhh(idim,1,ind); ig1=ggg(idim,1,ind); ih1=ncoarse+(id1-1)*ngridmax
        id2=hhh(idim,2,ind); ig2=ggg(idim,2,ind); ih2=ncoarse+(id2-1)*ngridmax
        id3=hhh(idim,3,ind); ig3=ggg(idim,3,ind); ih3=ncoarse+(id3-1)*ngridmax
        id4=hhh(idim,4,ind); ig4=ggg(idim,4,ind); ih4=ncoarse+(id4-1)*ngridmax

        ! Gather potential
        do i=1,ngrid_aux(ibin)

           if(igridn(ibin,i,ig1)>0)then
              phi1(ibin,i)=phi(igridn(ibin,i,ig1)+ih1)
           else
              phi1(ibin,i)=phi_left(ibin,i,id1,idim)
           end if

           if(igridn(ibin,i,ig2)>0)then
              phi2(ibin,i)=phi(igridn(ibin,i,ig2)+ih2)
           else
              phi2(ibin,i)=phi_right(ibin,i,id2,idim)
           end if

           if(igridn(ibin,i,ig3)>0)then
              phi3(ibin,i)=phi(igridn(ibin,i,ig3)+ih3)
           else
              phi3(ibin,i)=phi_left(ibin,i,id3,idim)
           end if

           if(igridn(ibin,i,ig4)>0)then
              phi4(ibin,i)=phi(igridn(ibin,i,ig4)+ih4)
           else
              phi4(ibin,i)=phi_right(ibin,i,id4,idim)
           end if

           f(ind_cell_aux(ibin,i),idim)=a*(phi1(ibin,i)-phi2(ibin,i)) - b*(phi3(ibin,i)-phi4(ibin,i))
        end do
     end do
  end do

end subroutine gradient_phi
