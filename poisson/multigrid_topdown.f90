! ------------------------------------------------------------------------
! Multigrid Poisson solver for refined AMR levels
! ------------------------------------------------------------------------
! This file contains top-down MG-fine-level and MG-coarse-level related routines
! (OBSOLETE, UNUSED)
!
! Used variables:
!                       finest(AMR)level     coarse(MG)levels
!     -----------------------------------------------------------------
!     potential            phi            active_mg(myid,ilevel)%u(:,1)
!     physical RHS         rho            active_mg(myid,ilevel)%u(:,2)
!     residual             f(:,1)         active_mg(myid,ilevel)%u(:,3)
!     BC-modified RHS      f(:,2)                  N/A
!     mask                 f(:,3)         active_mg(myid,ilevel)%u(:,4)
!
! ------------------------------------------------------------------------


! ------------------------------------------------------------------------
! Mask restriction (top-down, OBSOLETE, UNUSED)
! ------------------------------------------------------------------------

subroutine restrict_mask_fine(ifinelevel,allmasked)
   use amr_commons
   use poisson_commons
   implicit none
   integer, intent(in)  :: ifinelevel
   logical, intent(out) :: allmasked

   integer :: ind_c_cell,ind_f_cell

   integer :: iskip_f_amr, iskip_c_amr, iskip_c_mg
   integer :: igrid_f_amr, igrid_c_amr, igrid_c_mg
   integer :: icell_f_amr, icell_c_amr, icell_c_mg

   real(dp) :: ngpmask
   real(dp) :: dtwotondim = (twotondim)

   integer :: icoarselevel
   icoarselevel=ifinelevel-1
   allmasked=.true.

   if(ifinelevel==1) return

   ! Loop over coarse cells of the coarse active comm for myid
   do ind_c_cell=1,twotondim
      iskip_c_amr=ncoarse+(ind_c_cell-1)*ngridmax
      iskip_c_mg =(ind_c_cell-1)*active_mg(myid,icoarselevel)%ngrid

      ! Loop over coarse grids
      do igrid_c_mg=1,active_mg(myid,icoarselevel)%ngrid
#ifdef LIGHT_MPI_COMM
         igrid_c_amr=active_mg(myid,icoarselevel)%pcomm%igrid(igrid_c_mg)
#else
         igrid_c_amr=active_mg(myid,icoarselevel)%igrid(igrid_c_mg)
#endif
         icell_c_amr=iskip_c_amr+igrid_c_amr
         icell_c_mg =iskip_c_mg +igrid_c_mg

         if(son(icell_c_amr)==0) then
            ! Cell is not refined
            ngpmask      = -1.0d0
         else
            igrid_f_amr=son(icell_c_amr)
            ngpmask = 0.0d0
            do ind_f_cell=1,twotondim
               iskip_f_amr=ncoarse+(ind_f_cell-1)*ngridmax
               icell_f_amr=iskip_f_amr+igrid_f_amr
               ngpmask=ngpmask+f(icell_f_amr,3)
            end do
            ngpmask=ngpmask/dtwotondim
         end if
         ! Store cell mask
#ifdef LIGHT_MPI_COMM
         active_mg(myid,icoarselevel)%pcomm%u(icell_c_mg,4)=ngpmask
#else
         active_mg(myid,icoarselevel)%u(icell_c_mg,4)=ngpmask
#endif
         allmasked=allmasked .and. (ngpmask<=0.0)
      end do
   end do

end subroutine restrict_mask_fine


! ------------------------------------------------------------------------
! Residual restriction (top-down, OBSOLETE, UNUSED)
! ------------------------------------------------------------------------

subroutine restrict_residual_fine(ifinelevel)
   ! Restrict fine (AMR) residual at level ifinelevel using injection
   ! into coarser residual at level ifinelevel-1
   ! Restricted residual is stored into the RHS at the coarser level
   use amr_commons
   use pm_commons
   use poisson_commons
   implicit none
   integer, intent(in) :: ifinelevel

   real(dp) :: val
   real(dp) :: dtwotondim = (twotondim)

   integer  :: icoarselevel
   integer  :: ngrid_c, ind_c, iskip_c_amr, iskip_c_mg
   integer  :: igrid_c_amr, icell_c_amr, icell_c_mg, igrid_c_mg
   integer  :: ind_f, igrid_f_amr, iskip_f_amr, icell_f_amr

   icoarselevel=ifinelevel-1

   ! Loop over coarse MG cells
   ngrid_c=active_mg(myid,icoarselevel)%ngrid
   do ind_c=1,twotondim
      iskip_c_amr = ncoarse + (ind_c-1)*ngridmax
      iskip_c_mg  = (ind_c-1)*ngrid_c

      do igrid_c_mg=1,ngrid_c
#ifdef LIGHT_MPI_COMM
         igrid_c_amr = active_mg(myid,icoarselevel)%pcomm%igrid(igrid_c_mg)
#else
         igrid_c_amr = active_mg(myid,icoarselevel)%igrid(igrid_c_mg)
#endif
         icell_c_amr = igrid_c_amr + iskip_c_amr
         icell_c_mg  = igrid_c_mg  + iskip_c_mg

         ! Get AMR child grid
         igrid_f_amr = son(icell_c_amr)
         if(igrid_f_amr==0) then
            ! Nullify residual (coarser RHS)
#ifdef LIGHT_MPI_COMM
            active_mg(myid,icoarselevel)%pcomm%u(icell_c_mg,2) = 0.0d0
#else
            active_mg(myid,icoarselevel)%u(icell_c_mg,2) = 0.0d0
#endif
            cycle
         end if

         val = 0.0d0
         ! Loop over child (fine MG) cells
         do ind_f=1,twotondim
            iskip_f_amr = ncoarse + (ind_f-1)*ngridmax
            icell_f_amr = igrid_f_amr + iskip_f_amr

            if (f(icell_f_amr,3)<=0.0) cycle
            val = val + f(icell_f_amr,1)
         end do
         ! Store restricted residual into RHS of coarse level
#ifdef LIGHT_MPI_COMM
         active_mg(myid,icoarselevel)%pcomm%u(icell_c_mg,2) = val/dtwotondim
#else
         active_mg(myid,icoarselevel)%u(icell_c_mg,2) = val/dtwotondim
#endif
      end do
   end do
end subroutine restrict_residual_fine


! ------------------------------------------------------------------------
! Mask restriction (top-down, OBSOLETE, UNUSED)
! ------------------------------------------------------------------------

subroutine restrict_mask_coarse(ifinelevel,allmasked)
   use amr_commons
   use poisson_commons
   implicit none
   integer, intent(in) :: ifinelevel
   logical, intent(out) :: allmasked

   integer :: ind_c_cell, ind_f_cell, cpu_amr

   integer :: iskip_c_amr, iskip_c_mg
   integer :: igrid_c_amr, igrid_c_mg
   integer :: icell_c_amr, icell_c_mg

   integer :: iskip_f_mg
   integer :: igrid_f_amr, igrid_f_mg
   integer :: icell_f_mg

   real(dp) :: ngpmask
   real(dp) :: dtwotondim = (twotondim)

   integer :: icoarselevel
   icoarselevel=ifinelevel-1
   allmasked=.true.

   ! Loop over coarse cells of the myid active comm
   do ind_c_cell=1,twotondim
      iskip_c_amr=ncoarse+(ind_c_cell-1)*ngridmax
      iskip_c_mg =(ind_c_cell-1)*active_mg(myid,icoarselevel)%ngrid

      ! Loop over coarse grids of myid
      do igrid_c_mg=1,active_mg(myid,icoarselevel)%ngrid
#ifdef LIGHT_MPI_COMM
         igrid_c_amr=active_mg(myid,icoarselevel)%pcomm%igrid(igrid_c_mg)
#else
         igrid_c_amr=active_mg(myid,icoarselevel)%igrid(igrid_c_mg)
#endif
         icell_c_amr=iskip_c_amr+igrid_c_amr
         icell_c_mg =iskip_c_mg +igrid_c_mg
         igrid_f_amr=son(icell_c_amr)
         cpu_amr=cpu_map(icell_c_amr)
         if(igrid_f_amr==0) then
            ! Cell is not refined
            ngpmask      = -1.0d0
         else
            ! Cell is refined
            ! Check if son grid is in MG hierarchy
            igrid_f_mg=lookup_mg(igrid_f_amr)
            if(igrid_f_mg<=0) then
               ! Child oct is not in multigrid hierarchy
               ngpmask=-1.0d0
            else
               ! Child oct is within MG hierarchy
               ! Loop over fine cells and gather ngpmask
               ngpmask=0.0d0
               do ind_f_cell=1,twotondim
                  ! Extract fine mask value in the corresponding MG comm
                  iskip_f_mg=(ind_f_cell-1)*active_mg(cpu_amr,ifinelevel)%ngrid
                  icell_f_mg=iskip_f_mg+igrid_f_mg
#ifdef LIGHT_MPI_COMM
                  ngpmask=ngpmask+active_mg(cpu_amr,ifinelevel)%pcomm%u(icell_f_mg,4)
#else
                  ngpmask=ngpmask+active_mg(cpu_amr,ifinelevel)%u(icell_f_mg,4)
#endif
               end do
               ngpmask=ngpmask/dtwotondim
            end if
         end if
         ! Store cell mask
#ifdef LIGHT_MPI_COMM
         active_mg(myid,icoarselevel)%pcomm%u(icell_c_mg,4)=ngpmask
#else
         active_mg(myid,icoarselevel)%u(icell_c_mg,4)=ngpmask
#endif
         allmasked=allmasked .and. (ngpmask<=0.0)
      end do
   end do

end subroutine restrict_mask_coarse


! ------------------------------------------------------------------------
! Residual restriction (top-down, OBSOLETE, UNUSED)
! ------------------------------------------------------------------------

subroutine restrict_residual_coarse(ifinelevel)
   ! Restrict coarser (MG) residual at level ifinelevel using NGP into coarser residual at level
   ! ifinelevel-1
   ! Restricted residual is stored into the RHS at the coarser level
   use amr_commons
   use pm_commons
   use poisson_commons
   implicit none
   integer, intent(in) :: ifinelevel

   real(dp) :: val, w
   integer  :: icoarselevel, cpu_amr
   integer  :: ngrid_c, ind_c, iskip_c_amr, iskip_c_mg, igrid_c_amr, icell_c_amr, icell_c_mg, igrid_c_mg
   integer  :: ind_f, igrid_f_amr, igrid_f_mg, icell_f_mg

   icoarselevel=ifinelevel-1

   ! Loop over coarse MG cells
   ngrid_c=active_mg(myid,icoarselevel)%ngrid
   do ind_c=1,twotondim
      iskip_c_amr = ncoarse + (ind_c-1)*ngridmax
      iskip_c_mg  = (ind_c-1)*ngrid_c

      do igrid_c_mg=1,ngrid_c
#ifdef LIGHT_MPI_COMM
         igrid_c_amr = active_mg(myid,icoarselevel)%pcomm%igrid(igrid_c_mg)
#else
         igrid_c_amr = active_mg(myid,icoarselevel)%igrid(igrid_c_mg)
#endif
         icell_c_amr = igrid_c_amr + iskip_c_amr
         cpu_amr     = cpu_map(icell_c_amr)
         icell_c_mg  = igrid_c_mg  + iskip_c_mg

         ! Get AMR child grid
         igrid_f_amr = son(icell_c_amr)
         if(igrid_f_amr==0) then
#ifdef LIGHT_MPI_COMM
            active_mg(myid,icoarselevel)%pcomm%u(icell_c_mg,2) = 0.0d0    ! Nullify residual (coarser RHS)
#else
            active_mg(myid,icoarselevel)%u(icell_c_mg,2) = 0.0d0    ! Nullify residual (coarser RHS)
#endif
            cycle
         end if

         ! Get child MG grid id
         igrid_f_mg = lookup_mg(igrid_f_amr)
         if(igrid_f_mg<=0) then
            ! Son grid is not in MG hierarchy
#ifdef LIGHT_MPI_COMM
            active_mg(myid,icoarselevel)%pcomm%u(icell_c_mg,2) = 0.0d0    ! Nullify residual (coarser RHS)
#else
            active_mg(myid,icoarselevel)%u(icell_c_mg,2) = 0.0d0    ! Nullify residual (coarser RHS)
#endif
            cycle
         end if

         ! Loop over child (fine MG) cells
         val = 0
         w = 0
         do ind_f=1,twotondim
            icell_f_mg = igrid_f_mg + (ind_f-1)*active_mg(cpu_amr,ifinelevel)%ngrid

#ifdef LIGHT_MPI_COMM
            if (active_mg(cpu_amr,ifinelevel)%pcomm%u(icell_f_mg,4)<=0.0) cycle
            val = val + active_mg(cpu_amr,ifinelevel)%pcomm%u(icell_f_mg,3)
#else
            if (active_mg(cpu_amr,ifinelevel)%u(icell_f_mg,4)<=0.0) cycle
            val = val + active_mg(cpu_amr,ifinelevel)%u(icell_f_mg,3)
#endif
            w = w + 1d0
         end do
         ! Store restricted residual into RHS of coarse level
#ifdef LIGHT_MPI_COMM
         if(w>0) then
            active_mg(myid,icoarselevel)%pcomm%u(icell_c_mg,2) = val/w
         else
            active_mg(myid,icoarselevel)%pcomm%u(icell_c_mg,2) = 0
         end if
#else
         if(w>0) then
            active_mg(myid,icoarselevel)%u(icell_c_mg,2) = val/w
         else
            active_mg(myid,icoarselevel)%u(icell_c_mg,2) = 0
         end if
#endif
      end do
   end do
end subroutine restrict_residual_coarse
