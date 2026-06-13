!###############################################################
!###############################################################
subroutine crmom_step(ilevel)
  ! Two-moment cosmic-ray HYPERBOLIC transport driver: subcycled
  ! free-streaming P1/M1 advection, ported from ramses_cral feat/CR_tests
  ! cr/cr_godunov_fine.f90:875-950 (crmom_step). SOURCE TERMS and cooling
  ! are deferred to Phase 2: cral's add_cr_source_terms (line 922) and
  ! cr_cooling_fine (line 926) are NOT called here.
  !
  ! Either do one transport step on ilevel, with CR field updates in the
  ! coarser-level neighbours, or, if cr_nsubcycle>1, do many substeps in
  ! ilevel only, using Dirichlet boundary conditions for the level
  ! boundaries.
  !
  ! Central transformation vs cral: CR state lives in the separate arrays
  ! cruold/crunew(1:ncell,1:ncrvars) with iCRu=1 (cral embedded it in
  ! uold/unew at nvar+4). Virtual-boundary exchange loops therefore run
  ! over 1:ncrvars on crunew/cruold, never over nvar+3+1..nvar+3+ncrvars.
  use amr_commons
  use cr_parameters
  use cr_hydro_commons
  use constants
  use mpi_mod
  implicit none
  integer,intent(in)::ilevel
  real(dp)::dt_hydro,t_left,dt_cr,t_save
  integer::i_substep,ivar
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  dt_hydro=dtnew(ilevel)                     ! Store hydro timestep length
  t_left=dt_hydro
  ! Shift the time backwards one hydro-dt, in case of CR subcycling. Not
  ! strictly necessary as CR advection has no explicit time dependence.
  t_save=t ; t=t-t_left

  ! Initialise crunew=cruold for the first subcycle (zeroes virtual cells).
  ! In cral this was the external cr_set_unew call in amr_step right after
  ! set_unew; sno has no such call, so it is done here once before the loop.
  call cr_set_unew(ilevel)

  call get_crmom_courant(dt_cr,ilevel)
  i_substep=0
  do while(t_left>0)                         !                CR sub-cycle
     i_substep=i_substep+1
     ! Temporarily change timestep length to the CR step:
     dtnew(ilevel)=MIN(t_left,dt_cr)
     if(i_substep.gt.cr_nsubcycle)then
        print*,'Doing CR substeps but should not!! ',i_substep,cr_nsubcycle,dt_hydro,t_left
        call clean_stop
     endif
     t=t+dtnew(ilevel)                        ! Shift the time forwards one dt_cr

     if(i_substep>1)call cr_set_unew(ilevel)

     call cr_godunov_fine(ilevel)
     ! Pull in updates from finer level on other cpus -- only needed when
     ! updating coarser level (i.e. not subcycling).
     if(cr_nsubcycle==1)then
        do ivar=1,ncrvars
           call make_virtual_reverse_dp(crunew(1,ivar),ilevel)
        end do
     endif

     ! Set cruold equal to crunew, but only for the CR vars
     call cr_set_uold(ilevel)

     do ivar=1,ncrvars
        call make_virtual_fine_dp(cruold(1,ivar),ilevel)
     end do

     if(simple_boundary)call cr_make_boundary_hydro(ilevel)

     t_left=t_left-dtnew(ilevel)
  end do                                     !          End CR subcycle loop
  dtnew(ilevel)=dt_hydro                     ! Restore hydro timestep length
  t=t_save        ! Restore original time (otherwise tiny roundoff error)

  ! Restriction operator to update this level split cells (gas only;
  ! gas is frozen in Phase 1 so this is a no-op restriction on it).
  call upload_fine(ilevel)

  if(myid==1 .and. mod(nstep_coarse,ncontrol)==0)then
     call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
     write(*,901)ilevel,i_substep,cr_vmax(ilevel)*scale_v/1e5, &
          cr_vmax(ilevel)*scale_v/c_cgs,dt_cr
  endif

111 format('   Entering crmom_step for level ',I2)
901 format(' Performed level',I3,' CR-step with ',I5,' subcycles, cr_vmax(km/s)=', &
         1pe9.2,' cr_c_fraction=',1pe9.2,',  dt_cr=',1pe9.2)

end subroutine crmom_step

!************************************************************************
SUBROUTINE cr_godunov_fine(ilevel)

! This routine is a wrapper to the grid solver for cosmic-ray transport.
! Small grids (2x2x2) are gathered from level ilevel and sent to the CR
! solver. On entry, CR variables are gathered from array cruold. On exit,
! crunew has been updated. Ported from cral cr/cr_godunov_fine.f90:2-33.
!------------------------------------------------------------------------
   use amr_commons
   use cr_hydro_commons
   implicit none
   integer::ilevel
   integer::i,igrid,ncache,ngrid
   integer,dimension(1:nvector),save::ind_grid
!------------------------------------------------------------------------
   if(numbtot(1,ilevel)==0)return  ! # of grids at ilevel
   if(verbose)write(*,111)ilevel

   ncache=active(ilevel)%ngrid  ! total # of grids at level ilevel
   do igrid=1,ncache,nvector    ! take steps of nvector grids up to ncache
      ngrid=MIN(nvector,ncache-igrid+1) ! # of grids in each sweep
      do i=1,ngrid              ! collect grid indices for one sweep
         ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
      end do
      call cr_godfine1(ind_grid,ngrid,ilevel)
   end do

   if(verbose)write(*,112)ilevel

111 format('   Entering cr_godunov_fine for level ',i2)
112 format('   Exiting cr_godunov_fine for level ',i2)

END SUBROUTINE cr_godunov_fine

!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE cr_set_unew(ilevel)
   use amr_commons
   use cr_hydro_commons
   use cr_parameters
   implicit none
   integer::ilevel
   !--------------------------------------------------------------------------
   ! This routine sets array crunew to its initial value cruold before
   ! calling the CR advection scheme. crunew is set to zero in virtual
   ! boundaries. Ported from cral cr/cr_godunov_fine.f90:39-96, operating on
   ! the separate cruold/crunew arrays (iCRu=1) over 1:ncrvars.
   !--------------------------------------------------------------------------
   integer::i,ivar,ind,icpu,iskip
   real(dp)::fred,sqrt3

   if(numbtot(1,ilevel)==0)return
   if(verbose)write(*,111)ilevel

   if(isotropic_pressure)then
      sqrt3=sqrt(3d0)
   else
      sqrt3=1.0d0
   endif

   if(reduced_CR_flux_correction)then
      do ind=1,twotondim
         iskip=ncoarse+(ind-1)*ngridmax
         do i=1,active(ilevel)%ngrid
            fred=sqrt(sum(cruold(active(ilevel)%igrid(i)+iskip,iCRu+1:iCRu+ndim)**2)) &
                 /(cr_vmax(ilevel)*cruold(active(ilevel)%igrid(i)+iskip,iCRu))*sqrt3
            if(fred>1.0)then
               cruold(active(ilevel)%igrid(i)+iskip,iCRu+1:iCRu+ndim) = &
                    cruold(active(ilevel)%igrid(i)+iskip,iCRu+1:iCRu+ndim)/fred
            endif
         end do
      end do
   endif

   ! Set crunew to cruold for myid cells
   do ind=1,twotondim
      iskip=ncoarse+(ind-1)*ngridmax
      do ivar=1,ncrvars
         do i=1,active(ilevel)%ngrid
            crunew(active(ilevel)%igrid(i)+iskip,ivar)=cruold(active(ilevel)%igrid(i)+iskip,ivar)
         end do
      end do
   end do

   ! Set crunew to 0 for virtual boundary cells
   do icpu=1,ncpu
   do ind=1,twotondim
      iskip=ncoarse+(ind-1)*ngridmax
      do ivar=1,ncrvars
         do i=1,reception(icpu,ilevel)%ngrid
            crunew(reception(icpu,ilevel)%igrid(i)+iskip,ivar)=0.0
         end do
      end do
   end do
   end do

111 format('   Entering cr_set_unew for level ',i2)

END SUBROUTINE cr_set_unew

!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE cr_set_uold(ilevel)
   use amr_commons
   use cr_hydro_commons
   use cr_parameters
   implicit none
   integer::ilevel
   !--------------------------------------------------------------------------
   ! This routine sets array cruold to its new value crunew after the CR
   ! step. Ported from cral cr/cr_godunov_fine.f90:102-148, operating on the
   ! separate cruold/crunew arrays (iCRu=1) over 1:ncrvars.
   !
   ! Transport-robustness floor: the per-group CR energy is floored at
   ! smallcr after the update (mirroring rt_set_uold's max(...,smallNp)
   ! "prevent beam induced crash" fix in rt/rt_godunov_fine.f90). In cral
   ! this floor lives in add_cr_source_terms, which Phase 1 omits; without
   ! it the explicit transport of a degenerate state can drive E_cr slightly
   ! negative and trip the Ecr<0 guard in cmp_cr_flux_tensors. This is pure
   ! transport robustness -- no scattering/streaming/cooling source physics.
   !--------------------------------------------------------------------------
   integer::i,ivar,ind,iskip,iGrp,icrE
   real(dp)::fred,sqrt3

   if(numbtot(1,ilevel)==0)return
   if(verbose)write(*,111)ilevel

   if(isotropic_pressure)then
      sqrt3=sqrt(3d0)
   else
      sqrt3=1.0d0
   endif

   if(reduced_CR_flux_correction)then
      do ind=1,twotondim
         iskip=ncoarse+(ind-1)*ngridmax
         do i=1,active(ilevel)%ngrid
            fred=sqrt(sum(crunew(active(ilevel)%igrid(i)+iskip,iCRu+1:iCRu+ndim)**2)) &
                 /(cr_vmax(ilevel)*crunew(active(ilevel)%igrid(i)+iskip,iCRu))*sqrt3
            if(fred>1.0)then
               crunew(active(ilevel)%igrid(i)+iskip,iCRu+1:iCRu+ndim) = &
                    crunew(active(ilevel)%igrid(i)+iskip,iCRu+1:iCRu+ndim)/fred
            endif
         end do
      end do
   endif

   ! Set cruold to crunew for myid cells
   do ind=1,twotondim
      iskip=ncoarse+(ind-1)*ngridmax
      do ivar=1,ncrvars
         do i=1,active(ilevel)%ngrid
            cruold(active(ilevel)%igrid(i)+iskip,ivar)=crunew(active(ilevel)%igrid(i)+iskip,ivar)
         end do
      end do

      ! Floor each group's CR energy at smallcr (no negative CR densities),
      ! mirroring the rt_set_uold photon-density fix.
      do iGrp=1,ncr
         icrE=iCRu+(ndim+1)*(iGrp-1)
         do i=1,active(ilevel)%ngrid
            cruold(active(ilevel)%igrid(i)+iskip,icrE)= &
                 max(cruold(active(ilevel)%igrid(i)+iskip,icrE),smallcr)
         end do
      end do
   end do

111 format('   Entering cr_set_uold for level ',i2)

END SUBROUTINE cr_set_uold

!*************************************************************************
SUBROUTINE cr_godfine1(ind_grid, ncache, ilevel)

! This routine first gathers CR variables from neighboring grids to set
! initial conditions in a 6x6x6 grid. It interpolates from coarser level to
! missing grid variables. It then calls the solver that computes fluxes.
! These fluxes are zeroed at coarse-fine boundaries, since contribution from
! finer levels has already been taken into account. Conservative variables
! are updated and stored in array crunew(:), both at the current level and at
! the coarser level if necessary.
!
! Ported from cral cr/cr_godunov_fine.f90:151-495. Central transformation:
! the single embedded uloc(...,1:nvar+3+ncrvars) stencil is split into TWO
! stencils -- uin_gas(...,1:nvar+3) gathered from uold (rho/vel/B, needed by
! the M1 closure), and uin_cr(...,1:ncrvars) gathered from cruold (E,F with
! iCRu=1). The conservative update writes crunew(:,1:ncrvars); there are NO
! gas writes (gas momentum/energy back-reaction lives only in the omitted
! add_cr_source_terms). The buffer-neighbour gather uses the same
! get3cubefather/getnborfather/interpol machinery as gas: gas via
! interpol_hydro, CR via the local cr_interpol_hydro (mirroring RT).
!
! in ind_grid: Indexes of grids/octs to solve in
! in ncache:   Length of ind_grid (i.e. number of grids)
! in ilevel:   Level at which the grids are
!-------------------------------------------------------------------------
   use amr_commons
   use hydro_commons
   use cr_hydro_commons
   use cr_parameters
   use cr_flux_module
   use amr_constants, only:i1min,i1max,j1min,j1max,k1min,k1max, &
                        &  i2min,i2max,j2min,j2max,k2min,k2max, &
                        &  i3min,i3max,j3min,j3max,k3min,k3max
   implicit none
   integer::ilevel,ncache
   real(dp)::dt
   integer   ,dimension(1:nvector)::ind_grid
   integer   ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
   integer   ,dimension(1:nvector,0:twondim    ),save::ibuffer_father
   integer   ,dimension(1:nvector,0:twondim    ),save::ind1
   ! Split buffers for the coarser-level neighbour interpolation:
   real(dp),dimension(1:nvector,0:twondim  ,1:nvar+3),save::u1_gas
   real(dp),dimension(1:nvector,1:twotondim,1:nvar+3),save::u2_gas
   real(dp),dimension(1:nvector,0:twondim  ,1:ncrvars),save::u1_cr
   real(dp),dimension(1:nvector,1:twotondim,1:ncrvars),save::u2_cr

   ! Split 6x6x6 stencils:
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3),save::uin_gas
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ncrvars),save::uin_cr
   real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:ncrvars,1:ndim),save::flux
   logical,dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::ok
   integer,dimension(1:nvector),save::igrid_nbor,ind_cell,ind_buffer
   integer,dimension(1:nvector),save::ind_exist,ind_nexist

   integer::i,j,idim,ind_son,ind_father,iskip,nbuffer
   integer::i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,k3,nx_loc,nb_noneigh,nexist
   integer::ivar,iGrp
   real(dp)::dx,scale,oneontwotondim
   real(dp)::fred,sqrt3
!------------------------------------------------------------------------
   oneontwotondim = 1d0/dble(twotondim) ! 1/8 in 3D
   ! Mesh spacing
   nx_loc=icoarse_max-icoarse_min+1    ! =1
   scale=boxlen/dble(nx_loc)           ! length per coarse oct (=boxlen)
   dx=0.5D0**ilevel*scale              ! length per oct/grid at ilevel
   dt=dtnew(ilevel)
   if(isotropic_pressure)then
      sqrt3=sqrt(3d0)
   else
      sqrt3=1.0d0
   endif

   !------------------------------------------
   ! Gather 3^ndim neighboring father cells
   ! of grid (ilevel-1). ind_cell are indexes
   ! in uold/cruold.
   !------------------------------------------
   do i=1,ncache
      ind_cell(i)=father(ind_grid(i))
   end do
   call get3cubefather(ind_cell,nbors_father_cells,ncache,ilevel)

   !---------------------------
   ! Gather 6x6x6 cells stencil
   !---------------------------
   ! Loop over 3x3x3 neighboring father cells (ilevel-1)
   do k1=k1min,k1max
   do j1=j1min,j1max
   do i1=i1min,i1max

      ! Check if neighbor has a grid (or is just a cell)
      nbuffer=0
      nexist=0

      ind_father=1+i1+3*j1+9*k1
      do i=1,ncache
         igrid_nbor(i)=son(nbors_father_cells(i,ind_father))
         if(igrid_nbor(i)>0) then
            nexist=nexist+1
            ind_exist(nexist)=i
         else
            nbuffer=nbuffer+1
            ind_nexist(nbuffer)=i
            ind_buffer(nbuffer)=nbors_father_cells(i,ind_father)
         end if
      end do

      if(nbuffer>0)then
         ! For those nongrid cells we interpolate variables from parent cells.
         call getnborfather(ind_buffer,ibuffer_father,nbuffer,ilevel)
         do j=0,twondim
            ! Gas (rho/vel/B) from uold, interpolated with the MHD-aware
            ! interpol_hydro (preserves div B).
            do ivar=1,nvar+3
               do i=1,nbuffer
                  u1_gas(i,j,ivar)=uold(ibuffer_father(i,j),ivar)
               end do
            end do
            ! CR (E,F) from cruold, interpolated as cell-centered scalars.
            do ivar=1,ncrvars
               do i=1,nbuffer
                  u1_cr(i,j,ivar)=cruold(ibuffer_father(i,j),ivar)
               end do
            end do
            do i=1,nbuffer
               ind1(i,j)=son(ibuffer_father(i,j))
            end do
         end do
         call interpol_hydro(u1_gas,ind1,u2_gas,nbuffer)
         call cr_interpol_hydro(u1_cr,u2_cr,nbuffer)
         if(reduced_CR_flux_correction)then
            do j=1,twotondim
               do i=1,nbuffer
                  ! F<e*c/sqrt(3) for P1 closure
                  fred=sqrt(sum(u2_cr(i,j,iCRu+1:iCRu+ndim)**2)) &
                       /(cr_vmax(ilevel)*u2_cr(i,j,iCRu))*sqrt3
                  if(fred.gt.1.0)then
                     u2_cr(i,j,iCRu+1:iCRu+ndim)=u2_cr(i,j,iCRu+1:iCRu+ndim)/fred
                  endif
               end do
            end do
         endif
      endif

      ! Loop 2x2x2 cells within father cell and add them to the stencils
      do k2=k2min,k2max
      do j2=j2min,j2max
      do i2=i2min,i2max

         ind_son=1+i2+2*j2+4*k2
         iskip=ncoarse+(ind_son-1)*ngridmax
         do i=1,nexist
            ind_cell(i)=iskip+igrid_nbor(ind_exist(i))
         end do

         i3=1; j3=1; k3=1
         if(ndim>0)i3=1+2*(i1-1)+i2 ! From -1 to 4 over outer loop, but
         if(ndim>1)j3=1+2*(j1-1)+j2 ! only over 2x2x2 indexes in inner loop
         if(ndim>2)k3=1+2*(k1-1)+k2

         ! Gather gas variables (rho/vel/B) from uold into uin_gas
         do ivar=1,nvar+3
            do i=1,nexist
               uin_gas(ind_exist(i),i3,j3,k3,ivar)=uold(ind_cell(i),ivar)
            end do
            do i=1,nbuffer
               uin_gas(ind_nexist(i),i3,j3,k3,ivar)=u2_gas(i,ind_son,ivar)
            end do
         end do

         ! Gather CR variables (E,F) from cruold into uin_cr
         do ivar=1,ncrvars
            do i=1,nexist
               uin_cr(ind_exist(i),i3,j3,k3,ivar)=cruold(ind_cell(i),ivar)
            end do
            do i=1,nbuffer
               uin_cr(ind_nexist(i),i3,j3,k3,ivar)=u2_cr(i,ind_son,ivar)
            end do
         end do

         ! Gather refinement flag
         do i=1,nexist
            ok(ind_exist(i),i3,j3,k3)=son(ind_cell(i))>0
         end do
         do i=1,nbuffer
            ok(ind_nexist(i),i3,j3,k3)=.false.
         end do

      end do
      end do
      end do
      ! End loop over cells

   end do
   end do
   end do
   ! End loop over neighboring grids

   !----------------------------------------------------------------------
   ! Compute fluxes of each CR group, using the Eddington tensor
   !----------------------------------------------------------------------
   do iGrp=1,ncr
      call cmp_cr_faces(uin_gas,uin_cr,flux,dx,dt,iGrp,ncache,ilevel)
   end do

   !----------------------------------------------------------------------
   ! Reset flux along direction at refined interface, if not cr-subcycling
   !----------------------------------------------------------------------
   if(cr_nsubcycle==1)then
   do idim=1,ndim
      i0=0; j0=0; k0=0
      if(idim==1)i0=1
      if(idim==2)j0=1
      if(idim==3)k0=1
      do k3=k3min,k3max+k0
      do j3=j3min,j3max+j0
      do i3=i3min,i3max+i0
         do i=1,ncache
            if(ok(i,i3-i0,j3-j0,k3-k0) .or. ok(i,i3,j3,k3))then
               flux(i,i3,j3,k3,:,idim)=0.0d0
            end if
         end do
      end do
      end do
      end do
   end do
   endif

   !--------------------------------------
   ! Conservative update at level ilevel
   !--------------------------------------
   do idim=1,ndim
      i0=0; j0=0; k0=0
      if(idim==1)i0=1
      if(idim==2)j0=1
      if(idim==3)k0=1
      do k2=k2min,k2max      ! all from 0 to 1
      do j2=j2min,j2max      ! => update 2x2x2 cells
      do i2=i2min,i2max      ! i.e. update one grid at level ilevel
         ind_son=1+i2+2*j2+4*k2
         iskip=ncoarse+(ind_son-1)*ngridmax
         do i=1,ncache
            ind_cell(i)=iskip+ind_grid(i)
         end do
         i3=1+i2
         j3=1+j2   ! just because the flux indexes are (1:3), not (0:2)
         k3=1+k2
         ! Update conservative CR variables (new state vector)
         do i=1,ncache
            crunew(ind_cell(i),1:ncrvars)= &
               crunew(ind_cell(i),1:ncrvars)   &
               & +(flux(i,i3   ,j3   ,k3   ,1:ncrvars,idim)  &
               & - flux(i,i3+i0,j3+j0,k3+k0,1:ncrvars,idim))
         end do
      end do
      end do
      end do
   end do

   !--------------------------------------
   ! Conservative update at level ilevel-1
   !--------------------------------------
   if(cr_nsubcycle==1)then
      ! Loop over dimensions
      do idim=1,ndim
         i0=0; j0=0; k0=0
         if(idim==1)i0=1
         if(idim==2)j0=1
         if(idim==3)k0=1

         !----------------------
         ! Left flux at boundary
         !----------------------
         ! Check if grid sits near left boundary
         ! and gather neighbor father cells index
         nb_noneigh=0
         do i=1,ncache
            if(son(nbor(ind_grid(i),2*idim-1))==0)then
               nb_noneigh=nb_noneigh+1
               ind_buffer(nb_noneigh)=nbor(ind_grid(i),2*idim-1)
               ind_cell(nb_noneigh)=i
            end if
         end do
         ! Conservative update of new state variables
         ! Loop over boundary cells
         do k3=k3min,k3max-k0 ! 1 to 1 if dim=3, 1 to 2 otherwise
            do j3=j3min,j3max-j0 ! 1 to 1 if dim=2, 1 to 2 otherwise
               do i3=i3min,i3max-i0 ! 1 to 1 if dim=1, 1 to 2 otherwise
                  do i=1,nb_noneigh
                     crunew(ind_buffer(i),1:ncrvars)= &
                          crunew(ind_buffer(i),1:ncrvars)   &
                          & -flux(ind_cell(i),i3,j3,k3,:,idim)  &
                          & *oneontwotondim
                  end do
               end do
            end do
         end do

         !-----------------------
         ! Right flux at boundary
         !-----------------------
         ! Check if grid sits near right boundary
         ! and gather neighbor father cells index
         nb_noneigh=0
         do i=1,ncache
            if(son(nbor(ind_grid(i),2*idim))==0)then
               nb_noneigh=nb_noneigh+1
               ind_buffer(nb_noneigh)=nbor(ind_grid(i),2*idim)
               ind_cell(nb_noneigh)=i
            end if
         end do
         ! Conservative update of new state variables
         ! Loop over boundary cells
         do k3=k3min+k0,k3max
            do j3=j3min+j0,j3max
               do i3=i3min+i0,i3max
                  do i=1,nb_noneigh
                     crunew(ind_buffer(i),1:ncrvars)= &
                          crunew(ind_buffer(i),1:ncrvars)   &
                          & +flux(ind_cell(i),i3+i0,j3+j0,k3+k0,:,idim)  &
                          & *oneontwotondim
                  end do
               end do
            end do
         end do

      end do
      ! End loop over dimensions
   end if
   ! End if-clause for cr-subcycling

END SUBROUTINE cr_godfine1

!************************************************************************
SUBROUTINE cr_interpol_hydro(u1,u2,nn)
! Interpolate the CR buffer (E,F over 1:ncrvars) from a coarse father cell to
! its 2^ndim children, treating every CR variable as a cell-centered scalar
! with the chosen slope limiter. This mirrors rt/rt_interpol_hydro.f90 (the
! separate-array analog) and reproduces the CR block of cral's combined
! interpol_hydro (mhd/interpol_hydro.f90:686-713), which the sno split
! stencil cannot route through interpol_hydro (sized 1:nvar+3).
!------------------------------------------------------------------------
  use amr_commons
  use cr_parameters
  use cr_hydro_commons
  use hydro_parameters, only: interpol_type
  implicit none
  integer::nn
  real(dp),dimension(1:nvector,0:twondim  ,1:ncrvars)::u1
  real(dp),dimension(1:nvector,1:twotondim,1:ncrvars)::u2
  integer::i,j,ivar,idim,ind,ix,iy,iz

  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,0:twondim),save::a
  real(dp),dimension(1:nvector,1:ndim),save::w
!------------------------------------------------------------------------
  ! Set position of cell centers relative to grid (oct) center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)
  end do

  ! Loop over CR interpolation variables
  do ivar=1,ncrvars

     ! Load father variable
     do j=0,twondim
        do i=1,nn
           a(i,j)=u1(i,j,ivar)
        end do
     end do

     ! Reset gradient
     w(1:nn,1:ndim)=0.0D0

     ! Compute gradient with chosen limiter
     if(interpol_type==1)call compute_limiter_minmod(a,w,nn)
     if(interpol_type==2)call compute_limiter_central(a,w,nn)
     if(interpol_type==3)call compute_central(a,w,nn)

     ! Interpolate over children cells
     do ind=1,twotondim
        u2(1:nn,ind,ivar)=a(1:nn,0)  ! center value
        do idim=1,ndim
           do i=1,nn
              u2(i,ind,ivar)=u2(i,ind,ivar)+w(i,idim)*xc(ind,idim)
           end do
        end do
     end do

  end do
  ! End loop over variables

END SUBROUTINE cr_interpol_hydro
