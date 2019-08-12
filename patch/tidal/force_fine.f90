!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine force_fine(ilevel,icount)
  use amr_commons
  use pm_commons
  use poisson_commons
  use mpi_mod
  use constants, only : twopi
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  integer::ilevel,icount
  !----------------------------------------------------------
  ! This routine computes the gravitational acceleration,
  ! the maximum density rho_max, and the potential energy
  !----------------------------------------------------------
  integer::igrid,ngrid,ncache,i,ind,iskip,ix,iy,iz
  integer::nx_loc,idim
  real(dp)::dx,dx_loc,scale,fact,fourpi
  real(kind=8)::rho_loc,rho_all,epot_loc,epot_all
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:3)::skip_loc

  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp),dimension(1:nvector,1:ndim),save::xx,ff

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
     ncache=active(ilevel)%ngrid
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do
        ! Compute gradient of potential
        call gradient_phi(ind_grid,ngrid,ilevel,icount)
     end do
     ! End loop over grids

#if NDIM==3
     if (sink)then
        call f_gas_sink(ilevel)
     end if
#endif
     ! Update boundaries
     do idim=1,ndim
        call make_virtual_fine_dp(f(1,idim),ilevel)
     end do
     if(simple_boundary)call make_boundary_force(ilevel)

  endif

  !----------------------------------------------
  ! Compute gravity potential and maximum density
  !----------------------------------------------
  rho_loc =0; rho_all =0
  epot_loc=0; epot_all=0
  fourpi=2*twopi
  if(cosmo)fourpi=1.5D0*omega_m*aexp
  fact=-dx_loc**ndim/fourpi/2.0D0

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
        ! Loop over dimensions
        do idim=1,ndim
           do i=1,ngrid
              if(son(ind_cell(i))==0)then
                 epot_loc=epot_loc+fact*f(ind_cell(i),idim)**2
              end if
           end do
        end do
        ! End loop over dimensions
        do i=1,ngrid
           rho_loc=MAX(rho_loc,dble(abs(rho(ind_cell(i)))))
        end do
     end do
     ! End loop over cells
  end do
  ! End loop over grids

#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(epot_loc,epot_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(rho_loc ,rho_all ,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
     epot_loc=epot_all
     rho_loc =rho_all
#endif
     epot_tot=epot_tot+epot_loc
     rho_max(ilevel)=rho_loc

111 format('   Entering force_fine for level ',I2)

end subroutine force_fine
!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine gradient_phi(ind_grid,ngrid,ilevel,icount)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ngrid,ilevel,icount
  integer,dimension(1:nvector)::ind_grid
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
  integer,dimension(1:3,1:4,1:8)::ggg,hhh

  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:ndim),save::ind_left,ind_right
  integer ,dimension(1:nvector,0:twondim),save::igridn
  real(dp),dimension(1:nvector),save::phi1,phi2,phi3,phi4
  real(dp),dimension(1:nvector,1:twotondim,1:ndim),save::phi_left,phi_right

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel

  ! Rescaling factor
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

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

  ! Gather neighboring grids
  do i=1,ngrid
     igridn(i,0)=ind_grid(i)
  end do
  do idim=1,ndim
     do i=1,ngrid
        ind_left (i,idim)=nbor(ind_grid(i),2*idim-1)
        ind_right(i,idim)=nbor(ind_grid(i),2*idim  )
        igridn(i,2*idim-1)=son(ind_left (i,idim))
        igridn(i,2*idim  )=son(ind_right(i,idim))
     end do
  end do

  ! Interpolate potential from upper level
  if (ilevel>levelmin)then
     do idim=1,ndim
        call interpol_phi(ind_left (1,idim),phi_left (1,1,idim),ngrid,ilevel,icount)
        call interpol_phi(ind_right(1,idim),phi_right(1,1,idim),ngrid,ilevel,icount)
     end do
  end if
  ! Loop over cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
     end do

     ! Loop over dimensions
     do idim=1,ndim

        ! Loop over nodes
        id1=hhh(idim,1,ind); ig1=ggg(idim,1,ind); ih1=ncoarse+(id1-1)*ngridmax
        id2=hhh(idim,2,ind); ig2=ggg(idim,2,ind); ih2=ncoarse+(id2-1)*ngridmax
        id3=hhh(idim,3,ind); ig3=ggg(idim,3,ind); ih3=ncoarse+(id3-1)*ngridmax
        id4=hhh(idim,4,ind); ig4=ggg(idim,4,ind); ih4=ncoarse+(id4-1)*ngridmax

        ! Gather potential
        do i=1,ngrid
           if(igridn(i,ig1)>0)then
              phi1(i)=phi(igridn(i,ig1)+ih1)
           else
              phi1(i)=phi_left(i,id1,idim)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig2)>0)then
              phi2(i)=phi(igridn(i,ig2)+ih2)
           else
              phi2(i)=phi_right(i,id2,idim)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig3)>0)then
              phi3(i)=phi(igridn(i,ig3)+ih3)
           else
              phi3(i)=phi_left(i,id3,idim)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig4)>0)then
              phi4(i)=phi(igridn(i,ig4)+ih4)
           else
              phi4(i)=phi_right(i,id4,idim)
           end if
        end do
        do i=1,ngrid
           f(ind_cell(i),idim)=a*(phi1(i)-phi2(i)) &
                &             -b*(phi3(i)-phi4(i))
        end do
     end do
  end do

end subroutine gradient_phi

!EDIT TINE
subroutine calc_tidal_field(ilevel,icount)
  use amr_commons
  use poisson_commons
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  integer::ilevel,icount
  !-------------------------------------------------
  ! this routine calculates the eigenvalues of the tidal field tensor
  !-------------------------------------------------
  integer::igrid,ngrid,ncache,i,ind,iskip,ncell,idim
  integer ,dimension(1:nvector),save::ind_grid
  real(dp),allocatable,dimension(:,:,:)::tidal_field    ! tidal tensor
  real(dp)::abs_err,A1=0, A2=0, A3=0
  real(dp),allocatable,dimension(:,:)::eigenv,a
  integer ,dimension(1:nvector),save::ind_cell

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ncell=ncoarse+twotondim*ngridmax
  allocate(tidal_field(1:ncell,1:ndim,1:ndim))
  allocate(tidal_eigval(1:ncell,1:ndim))
  allocate(eigenv(1:ndim,1:ndim))
  allocate(a(1:ndim,1:ndim))
  tidal_field=0; tidal_eigval=0; eigenv=0; a=0

  ! Loop over myid grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     ! Calculate the components of the local tidal tensor
     do idim=1,ndim
        call gradient_f(ind_grid,ngrid,ilevel,icount,idim, tidal_field)
     end do
     !call gradient_f(ind_grid,ngrid,ilevel,icount,2, tidal_field)
     !call gradient_f(ind_grid,ngrid,ilevel,icount,3, tidal_field)
     ! Calculate the eigenvalues
     ! Loop over cells
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i) = iskip+ind_grid(i)
           a=tidal_field(ind_cell(i),:,:)
           !abs_err=1d-8*Icl(j)**2+1d-40 # what to pick?
           !call jacobi(a,eigenv,abs_err) !in flagformationsites
           do idim=1,ndim
               !tidal_eigval(int_cell(i),idim)=a(idim,idim)
               tidal_eigval(ind_cell(i),idim)=DBLE(idim)
           end do
           !tidal_eigval(ind_cell(i),1)=1.d+0
           !tidal_eigval(ind_cell(i),2)=2.d+0
           !tidal_eigval(ind_cell(i),3)=3.d+0
        end do
     end do
  end do
  ! End loop over grids
  deallocate(tidal_field)

111 format('   Entering calc_tidal_field for level ',I2)

end subroutine calc_tidal_field

!EDIT TINE
subroutine gradient_f(ind_grid,ngrid,ilevel,icount,direction, tidal_field)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ngrid,ilevel,icount,direction
  integer,dimension(1:nvector)::ind_grid
  real(dp),dimension(1:nvector,1:ndim,1:ndim)::tidal_field
  !-------------------------------------------------
  ! This routine compute the components of the tidal tensor for all cells
  ! in grids ind_grid(:) at level ilevel, using a
  ! 5 nodes kernel (5 points FDA).
  !-------------------------------------------------
  integer::i,idim,ind,iskip,nx_loc
  integer::id1,id2,id3,id4
  integer::ig1,ig2,ig3,ig4
  integer::ih1,ih2,ih3,ih4
  real(dp)::dx,a,b,scale,dx_loc
  integer,dimension(1:3,1:4,1:8)::ggg,hhh

  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:ndim),save::ind_left,ind_right
  integer ,dimension(1:nvector,0:twondim),save::igridn
  real(dp),dimension(1:nvector),save::fi1,fi2,fi3,fi4
  real(dp),dimension(1:nvector,1:twotondim,1:ndim),save::fi_left,fi_right

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel

  ! Rescaling factor
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

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

  ! Gather neighboring grids
  do i=1,ngrid
     igridn(i,0)=ind_grid(i)
  end do
  do idim=1,ndim
     do i=1,ngrid
        ind_left (i,idim)=nbor(ind_grid(i),2*idim-1)
        ind_right(i,idim)=nbor(ind_grid(i),2*idim  )
        igridn(i,2*idim-1)=son(ind_left (i,idim))
        igridn(i,2*idim  )=son(ind_right(i,idim))
     end do
  end do

  ! Interpolate acceleration from upper level
  if (ilevel>levelmin)then
     do idim=1,ndim
        call interpol_f(ind_left (1,idim),fi_left (1,1,idim),ngrid,ilevel,icount,direction)
        call interpol_f(ind_right(1,idim),fi_right(1,1,idim),ngrid,ilevel,icount,direction)
     end do
  end if
  ! Loop over cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
     end do

     ! Loop over dimensions
     do idim=1,ndim

        ! Loop over nodes
        id1=hhh(idim,1,ind); ig1=ggg(idim,1,ind); ih1=ncoarse+(id1-1)*ngridmax
        id2=hhh(idim,2,ind); ig2=ggg(idim,2,ind); ih2=ncoarse+(id2-1)*ngridmax
        id3=hhh(idim,3,ind); ig3=ggg(idim,3,ind); ih3=ncoarse+(id3-1)*ngridmax
        id4=hhh(idim,4,ind); ig4=ggg(idim,4,ind); ih4=ncoarse+(id4-1)*ngridmax

        ! Gather acceleration
        do i=1,ngrid
           if(igridn(i,ig1)>0)then
              fi1(i)=f(igridn(i,ig1)+ih1,direction)
           else
              fi1(i)=fi_left(i,id1,idim)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig2)>0)then
              fi2(i)=f(igridn(i,ig2)+ih2,direction)
           else
              fi2(i)=fi_right(i,id2,idim)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig3)>0)then
              fi3(i)=f(igridn(i,ig3)+ih3,direction)
           else
              fi3(i)=fi_left(i,id3,idim)
           end if
        end do
        do i=1,ngrid
           if(igridn(i,ig4)>0)then
              fi4(i)=f(igridn(i,ig4)+ih4,direction)
           else
              fi4(i)=fi_right(i,id4,idim)
           end if
        end do
        do i=1,ngrid
           !TODO: I think the - sign is already included?
           tidal_field(ind_cell(i),idim,direction)=a*(fi1(i)-fi2(i)) &
                &             -b*(fi3(i)-fi4(i))
        end do
     end do
  end do

end subroutine gradient_f

subroutine interpol_f(ind_cell,fi_int,ncell,ilevel,icount,direction)
  use amr_commons
  use poisson_commons, only:f
  implicit none
  integer::ncell,ilevel,icount,direction
  integer ,dimension(1:nvector)::ind_cell
  real(dp),dimension(1:nvector,1:twotondim)::fi_int

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Routine for interpolation at level-boundaries. Interpolation is used for
  ! - computing tidal field (gradient_f) at fine level for cells close to boundary
  ! Interpolation is performed in space (CIC)
  ! DOES NOT TAKE INTO ACCOUNT ADAPTIVE TIMESTEPPING
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer::i,ind,indice,ind_average,ind_father
  real(dp)::dx,tfrac
  real(dp)::aa,bb,cc,dd,coeff,add
  integer,dimension(1:8,1:8)::ccc
  real(dp),dimension(1:8)::bbbb

  ! CIC method constants
  aa = 1d0/4d0**ndim
  bb = 3*aa
  cc = 9*aa
  dd = 27*aa
  bbbb(:)  =(/aa ,bb ,bb ,cc ,bb ,cc ,cc ,dd/)

  ! Sampling positions in the 3x3x3 father cell cube
  ccc(:,1)=(/1 ,2 ,4 ,5 ,10,11,13,14/)
  ccc(:,2)=(/3 ,2 ,6 ,5 ,12,11,15,14/)
  ccc(:,3)=(/7 ,8 ,4 ,5 ,16,17,13,14/)
  ccc(:,4)=(/9 ,8 ,6 ,5 ,18,17,15,14/)
  ccc(:,5)=(/19,20,22,23,10,11,13,14/)
  ccc(:,6)=(/21,20,24,23,12,11,15,14/)
  ccc(:,7)=(/25,26,22,23,16,17,13,14/)
  ccc(:,8)=(/27,26,24,23,18,17,15,14/)

  if (icount .ne. 1 .and. icount .ne. 2)then
     write(*,*)'icount has bad value'
     call clean_stop
  endif

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ncell,ilevel)

  ! Third order f interpolation
  do ind=1,twotondim
     do i=1,ncell
        fi_int(i,ind)=0
     end do
     do ind_average=1,twotondim
        ind_father=ccc(ind_average,ind)
        coeff=bbbb(ind_average)
        do i=1,ncell
           indice=nbors_father_cells(i,ind_father)
           if (indice==0) then
              add=coeff*(f(ind_cell(i),direction))!+(phi(ind_cell(i))-phi_old(ind_cell(i)))*tfrac)
           else
              add=coeff*(f(indice,direction))!+(phi(indice)-phi_old(indice))*tfrac)
           endif
           fi_int(i,ind)=fi_int(i,ind)+add
        end do
     end do
  end do

 end subroutine interpol_f
