!###########################################################
!###########################################################
!###########################################################
SUBROUTINE cr_courant_fine(ilevel)

   ! CR contribution to the level timestep, called after courant_fine has set
   ! dtnew(ilevel) from the gas: accumulate ecrs_tot and, when the CRs evolve,
   ! reduce dtnew by the gas Courant condition with CR pressure (cr_cmpdt),
   ! refresh cr_vmax/Dcr_code, apply the cr_varvmax cap, clamp to c, and
   ! reduce dtnew by the CR Courant condition.
   !-------------------------------------------------------------------------
   use amr_commons
   use hydro_commons
   use poisson_commons, only: f
   use cr_hydro_commons, only: cruold
   use cr_parameters, only: cr_advect,Ecr_idx,ecrs_tot, &
        & cr_vmax,cr_c_code,cr_nsubcycle,cr_varvmax, &
        & cr_varvmax_fudge,cr_varvmax_vdvs,cr_va_max, &
        & gamma_cr,Dcr_code,cr_streaming_diffusion,ncr_groups
   use mpi_mod
#if USE_TURB==1
   use turb_commons
#endif
   implicit none
#ifndef WITHOUTMPI
   integer::info
   real(kind=8)::cr_va_max_all
#endif
   integer::ilevel
   integer::i,ivar,idim,ind,ncache,igrid,iskip,igrp
   integer::nleaf,ngrid,nx_loc
   integer,dimension(1:nvector),save::ind_grid,ind_cell,ind_leaf
   real(dp)::dt_cr,dt_lev,dx,vol,scale
   real(kind=8)::ecrs_loc,ecrs_all,dt_loc,dt_all
   real(dp),dimension(1:nvector,1:nvar_all),save::uu
   real(dp),dimension(1:nvector,1:ndim),save::gg
   real(dp),dimension(1:nvector,1:ncr_groups),save::ecr
   !-------------------------------------------------------------------------
   if(numbtot(1,ilevel)==0)return

   ecrs_loc=0.0d0
   cr_va_max=0.0d0
   dt_loc=dtnew(ilevel)

   ! Mesh spacing at this level
   nx_loc=icoarse_max-icoarse_min+1
   scale=boxlen/dble(nx_loc)
   dx=0.5D0**ilevel*scale
   vol=dx**ndim

   ! Loop over active grids by vector sweeps
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
            ind_cell(i)=ind_grid(i)+iskip
         end do

         ! Gather leaf cells
         nleaf=0
         do i=1,ngrid
            if(son(ind_cell(i))==0)then
               nleaf=nleaf+1
               ind_leaf(nleaf)=ind_cell(i)
            end if
         end do

         ! Compute total CR energy
         do igrp=1,ncr_groups
            do i=1,nleaf
               ecrs_loc=ecrs_loc+cruold(ind_leaf(i),Ecr_idx(igrp))*vol
            end do
         end do

         if(.not.cr_advect .or. nleaf==0)cycle

         ! Gather hydro variables
         do ivar=1,nvar_all
            do i=1,nleaf
               uu(i,ivar)=uold(ind_leaf(i),ivar)
            end do
         end do

         ! Gather CR energies
         do igrp=1,ncr_groups
            do i=1,nleaf
               ecr(i,igrp)=cruold(ind_leaf(i),Ecr_idx(igrp))
            end do
         end do

         ! Gather gravitational acceleration
         gg=0.0d0
         if(poisson)then
            do idim=1,ndim
               do i=1,nleaf
                  gg(i,idim)=f(ind_leaf(i),idim)
               end do
            end do
         end if

#if USE_TURB==1
         ! Gather turbulence driving acceleration
         if (driven_turb) then
            do idim=1,ndim
               do i=1,nleaf
                  gg(i,idim)=gg(i,idim)+fturb(ind_leaf(i),idim)
               end do
            end do
         end if
#endif

         ! Compute CFL time-step with CR pressure
         call cr_cmpdt(uu,gg,ecr,dx,dt_lev,nleaf)
         dt_loc=min(dt_loc,dt_lev)

         if(cr_varvmax .and. cr_varvmax_vdvs .and. cr_streaming_diffusion) &
              & call cr_cmpdt_vamax(uu,nleaf)

      end do
      ! End loop over cells

   end do
   ! End loop over grids

   ! Compute global quantities
#ifndef WITHOUTMPI
   call MPI_ALLREDUCE(ecrs_loc,ecrs_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
        &MPI_COMM_WORLD,info)
   if(cr_advect)then
      call MPI_ALLREDUCE(dt_loc,dt_all,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
           &MPI_COMM_WORLD,info)
      call MPI_ALLREDUCE(cr_va_max,cr_va_max_all,1,MPI_DOUBLE_PRECISION,MPI_MAX,&
           &MPI_COMM_WORLD,info)
      cr_va_max=cr_va_max_all
   endif
#endif
#ifdef WITHOUTMPI
   ecrs_all=ecrs_loc
   dt_all=dt_loc
#endif

   ecrs_tot=ecrs_tot+ecrs_all
   if(.not.cr_advect)return
   dtnew(ilevel)=MIN(dtnew(ilevel),dt_all)

   call update_cr_vmax_and_Dcr_code(cr_vmax(ilevel))

   ! Adaptive reduced light speed (cr_varvmax): keep cr_vmax at least
   ! cr_varvmax_fudge times larger than the gas signal speed dx/3/dt so the
   ! CR wave always outruns the moving gas. Inert when cr_varvmax=.false.
   if(cr_varvmax)then
      cr_vmax(ilevel) = max(cr_vmax(ilevel), dx/3d0/dtnew(ilevel) * cr_varvmax_fudge)
      if(cr_varvmax_vdvs)then
         if(cr_streaming_diffusion)then
            do igrp=1,ncr_groups
               cr_vmax(ilevel) = max(cr_vmax(ilevel), gamma_cr(igrp)*cr_va_max * cr_varvmax_fudge)
            end do
         endif
         do igrp=1,ncr_groups
            cr_vmax(ilevel) = max(cr_vmax(ilevel), Dcr_code(igrp)/dx * (gamma_cr(igrp)-1d0) * cr_varvmax_fudge)
         end do
      endif
   endif

   ! Finally, make sure vmax <= c
   if(cr_vmax(ilevel) .gt. cr_c_code) cr_vmax(ilevel)=cr_c_code

   call get_crmom_courant(dt_cr,ilevel)
   dtnew(ilevel) = MIN(dtnew(ilevel), dt_cr * cr_nsubcycle*0.99999d0)

END SUBROUTINE cr_courant_fine
