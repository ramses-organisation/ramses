! ------------------------------------------------------------------------
! Multigrid Poisson solver for refined AMR levels
! ------------------------------------------------------------------------
! This file contains all MG-coarse-level related routines
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
         igrid_c_amr=active_mg(myid,icoarselevel)%igrid(igrid_c_mg)
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
                  ngpmask=ngpmask+active_mg(cpu_amr,ifinelevel)%u(icell_f_mg,4)
               end do
               ngpmask=ngpmask/dtwotondim
            end if
         end if
         ! Store cell mask
         active_mg(myid,icoarselevel)%u(icell_c_mg,4)=ngpmask
         allmasked=allmasked .and. (ngpmask<=0.0)
      end do
   end do

end subroutine restrict_mask_coarse

! ------------------------------------------------------------------------
! Mask restriction (bottom-up)
! ------------------------------------------------------------------------

subroutine restrict_mask_coarse_reverse(ifinelevel)
   use amr_commons
   use poisson_commons
   use poisson_commons_acc
   use ind_accessing
   implicit none
   integer, intent(in) :: ifinelevel

   integer :: ind_c_cell, ind_f_cell, cpu_amr
   integer :: iskip_c_mg
   integer :: igrid_c_amr, igrid_c_mg
   integer :: icell_c_amr, icell_c_mg
   integer :: iskip_f_mg
   integer :: igrid_f_amr, igrid_f_mg
   integer :: icell_f_mg
   real(dp) :: dtwotondim = (twotondim)
   integer :: icoarselevel
   real(dp),dimension(1:twotondim,1:ngrid_a_mg(myid,ifinelevel)) :: ngpmask
   integer,dimension(1:twotondim,1:ngrid_a_mg(myid,ifinelevel))  :: active_ind_aux
#ifdef MPROFILING
  include 'mpif.h'
  integer::info
  real(kind=8)::tt1,tt2
  tt1=MPI_WTIME(info)
#endif

   icoarselevel=ifinelevel-1

   !$acc data present(father,cpu_map,lookup_mg,active_mg_u,ngrid_a_mg,active_mg_igrid) &
   !$acc create(ngpmask,active_ind_aux)

   ! Loop over fine cells of the myid active comm
   !$acc parallel loop collapse(2) gang vector independent
   ! Loop over fine grids of myid
   do igrid_f_mg=1,ngrid_a_mg(myid,ifinelevel)
      do ind_f_cell=1,twotondim
         iskip_f_mg =(ind_f_cell-1)*ngrid_a_mg(myid,ifinelevel)
         icell_f_mg=iskip_f_mg+igrid_f_mg
         active_ind = ind_acc_a_mg_igrid(ifinelevel,myid)+igrid_f_mg
         igrid_f_amr=active_mg_igrid(active_ind)
         ! Get coarse grid AMR index and CPU id
         icell_c_amr=father(igrid_f_amr)
         ind_c_cell=(icell_c_amr-ncoarse-1)/ngridmax+1
         igrid_c_amr=icell_c_amr-ncoarse-(ind_c_cell-1)*ngridmax
         cpu_amr=cpu_map(father(igrid_c_amr))
         ! Convert to MG index, get MG coarse cell id
         igrid_c_mg=lookup_mg(igrid_c_amr)
         iskip_c_mg=(ind_c_cell-1)*ngrid_a_mg(cpu_amr,icoarselevel)
         icell_c_mg=iskip_c_mg+igrid_c_mg
         ! Stack cell volume fraction in coarse cell
         active_ind=ind_acc_a_mg_u(ifinelevel,myid,4)+icell_f_mg
         ngpmask(ind_f_cell,igrid_f_mg)=(1d0+active_mg_u(active_ind))/2d0/dtwotondim
         active_ind_aux(ind_f_cell,igrid_f_mg) = &
         & ind_acc_a_mg_u(icoarselevel,cpu_amr,4)+icell_c_mg
      end do
   end do

   !$acc parallel
   do ind_f_cell=1,twotondim
      !$acc loop
      do igrid_f_mg=1,ngrid_a_mg(myid,ifinelevel)
         active_mg_u(active_ind_aux(ind_f_cell,igrid_f_mg))= & 
         & active_mg_u(active_ind_aux(ind_f_cell,igrid_f_mg))+ & 
         & ngpmask(ind_f_cell,igrid_f_mg)
      end do
   end do
   !$acc end parallel

   !$acc end data

#ifdef MPROFILING
  tt2=MPI_WTIME(info)
  acc_t_restrict_mask_coarse_reverse = acc_t_restrict_mask_coarse_reverse + (tt2-tt1)
#endif   

end subroutine restrict_mask_coarse_reverse

! ------------------------------------------------------------------------
! Residual computation
! ------------------------------------------------------------------------

subroutine cmp_residual_mg_coarse(ilevel)
! Computes the residual for pure MG levels, and stores it into active_mg(myid,ilevel)%u(:,3)
   use amr_commons
   use poisson_commons
   use poisson_commons_acc
   use ind_accessing
   implicit none
   integer, intent(in) :: ilevel

   real(dp) :: dx, oneoverdx2, phi_c, nb_sum
   integer  :: ngrid
   integer  :: ind, igrid_mg, idim, inbor
   integer  :: icell_mg, iskip_mg, igrid_nbor_mg, icell_nbor_mg
   integer  :: igrid_amr, iskip_amr, cpu_nbor_amr 
   integer  :: igshift, igrid_nbor_amr

   real(dp) :: dtwondim = (twondim)
#ifdef MPROFILING
  include 'mpif.h'
  integer::info
  real(kind=8)::tt1,tt2
  tt1=MPI_WTIME(info)
#endif

   ! Set constants
   dx  = 0.5d0**ilevel
   oneoverdx2 = 1.0d0/(dx*dx)
   ngrid=ngrid_a_mg(myid,ilevel)
   nb_sum=0.0d0 

   !$acc data present(son,nbor,iii,jjj,cpu_map,lookup_mg) &
   !$acc present(active_mg_u,active_mg_f,ngrid_a_mg,active_mg_igrid)

   ! Loop over cells myid & active grids myid
   !$acc parallel loop collapse(2) gang worker vector independent &
   !$acc private(active_ind,tmp_ind) 
   do igrid_mg=1,ngrid
      do ind=1,twotondim
      iskip_mg  = (ind-1)*ngrid
      iskip_amr = ncoarse+(ind-1)*ngridmax

      active_ind = ind_acc_a_mg_igrid(ilevel,myid)+igrid_mg
      igrid_amr = active_mg_igrid(active_ind)
      icell_mg = igrid_mg + iskip_mg

      active_ind = ind_acc_a_mg_u(ilevel,myid,1)+icell_mg
      phi_c = active_mg_u(active_ind)
      nb_sum=0.0d0  ! Sum of phi on neighbors

!         ! SCAN FLAG TEST
!         active_ind = ind_acc_a_mg_f(ilevel,myid)+icell_mg
!         if(.not. btest(active_mg_f(active_ind),0)) then ! NO SCAN
!            do inbor=1,2
!               do idim=1,ndim
!                  ! Get neighbor grid
!                  igshift = iii(idim,inbor,ind)
!                  if(igshift==0) then
!                     igrid_nbor_amr = igrid_amr
!                     cpu_nbor_amr   = myid
!                  else
!                     igrid_nbor_amr = son(nbor(igrid_amr,igshift))
!                    cpu_nbor_amr   = cpu_map(nbor(igrid_amr,igshift))
!                  end if
!                  igrid_nbor_mg = lookup_mg(igrid_nbor_amr) 
!                  ! Add up
!                  icell_nbor_mg = igrid_nbor_mg + &
!                      (jjj(idim,inbor,ind)-1)*ngrid_a_mg(cpu_nbor_amr,ilevel)
!
!                  active_ind = ind_acc_a_mg_u(ilevel,cpu_nbor_amr,1)+icell_nbor_mg
!                  nb_sum = nb_sum + active_mg_u(active_ind)
!               end do
!            end do
!         else ! PERFORM SCAN
            active_ind = ind_acc_a_mg_u(ilevel,myid,4)+icell_mg
            if(active_mg_u(active_ind)<=0.0) then
               active_ind = ind_acc_a_mg_u(ilevel,myid,3)+icell_mg
               active_mg_u(active_ind)=0.0
               cycle
            end if
            do idim=1,ndim
               do inbor=1,2
                  ! Get neighbor grid
                  igshift = iii(idim,inbor,ind)
                  if(igshift==0) then
                     igrid_nbor_amr = igrid_amr
                     cpu_nbor_amr   = myid
                  else
                     igrid_nbor_amr = son(nbor(igrid_amr,igshift))
                     cpu_nbor_amr   = cpu_map(nbor(igrid_amr,igshift))
                  end if

                  if(igrid_nbor_amr==0) then
                     ! No neighbor cell !
                     ! Virtual phi value on unrefnd neighbor cell : -phi_c/mask_c
                     ! (simulates mask=-1.0 for the nonexistent refined cell)
                     active_ind = ind_acc_a_mg_u(ilevel,myid,4)+icell_mg
                     nb_sum = nb_sum + (- phi_c/active_mg_u(active_ind))
                  else
                     ! Fetch neighbor cell
                     igrid_nbor_mg  = lookup_mg(igrid_nbor_amr)
                     if(igrid_nbor_mg<=0) then
                        active_ind = ind_acc_a_mg_u(ilevel,myid,4)+icell_mg
                        nb_sum=nb_sum + (-phi_c/active_mg_u(active_ind))
                        cycle
                     end if

                     icell_nbor_mg  = igrid_nbor_mg + &
                       (jjj(idim,inbor,ind)-1)*ngrid_a_mg(cpu_nbor_amr,ilevel)

                     active_ind = ind_acc_a_mg_u(ilevel,cpu_nbor_amr,4)+icell_nbor_mg
                     if(active_mg_u(active_ind)<=0.0) then
                        active_ind = ind_acc_a_mg_u(ilevel,cpu_nbor_amr,4)+icell_nbor_mg
                        tmp_ind = ind_acc_a_mg_u(ilevel,myid,4)+icell_mg
                        ! Neighbor cell is masked : compute its virtual phi with the mask
                        nb_sum = nb_sum+phi_c*(active_mg_u(active_ind)/active_mg_u(tmp_ind))
                     else
                        active_ind = ind_acc_a_mg_u(ilevel,cpu_nbor_amr,1)+icell_nbor_mg
                        ! Neighbor cell is active, use its true potential
                        nb_sum = nb_sum + active_mg_u(active_ind)
                     end if
                  end if
               end do
            end do
!         end if ! END SCAN TEST
         active_ind = ind_acc_a_mg_u(ilevel,myid,3)+icell_mg
         tmp_ind    = ind_acc_a_mg_u(ilevel,myid,2)+icell_mg
         ! Store ***MINUS THE RESIDUAL***
         active_mg_u(active_ind)=-oneoverdx2*(nb_sum-dtwondim*phi_c)+active_mg_u(tmp_ind)
      end do
   end do

   !$acc end data
   
#ifdef MPROFILING
  tt2=MPI_WTIME(info)
  acc_t_cmp_residual_mg_coarse = acc_t_cmp_residual_mg_coarse + (tt2-tt1)
#endif   

end subroutine cmp_residual_mg_coarse

! ##################################################################
! ##################################################################

subroutine cmp_uvar_norm2_coarse(ivar, ilevel, norm2)
   use amr_commons
   use poisson_commons
   use ind_accessing
   implicit none

   integer,  intent(in)  :: ilevel, ivar
   real(dp), intent(out) :: norm2

   real(dp) :: dx2
   integer  :: ngrid
   integer  :: ind, igrid_mg, icell_mg, iskip_mg

   ! Set constants
   dx2  = (0.5d0**ilevel)**ndim
   ngrid=active_mg(myid,ilevel)%ngrid

   norm2 = 0.0d0
   ! Loop over cells
   do ind=1,twotondim
      iskip_mg = (ind-1)*ngrid
      ! Loop over active grids
      do igrid_mg=1,ngrid
         icell_mg = iskip_mg + igrid_mg
         if(active_mg(myid,ilevel)%u(icell_mg,4)<=0.0 .and. ivar/=4) cycle
         norm2 = norm2 + active_mg(myid,ilevel)%u(icell_mg,ivar)**2
      end do
   end do
   norm2 = dx2*norm2
end subroutine cmp_uvar_norm2_coarse

! ##################################################################
! ##################################################################

subroutine cmp_fvar_norm2_coarse(ivar, ilevel, norm2)
   use amr_commons
   use poisson_commons
   implicit none

   integer,  intent(in)  :: ilevel, ivar
   real(dp), intent(out) :: norm2

   real(dp) :: dx2
   integer  :: ngrid
   integer  :: ind, igrid_mg, icell_mg, iskip_mg

   ! Set constants
   dx2  = (0.5d0**ilevel)**ndim
   ngrid=active_mg(myid,ilevel)%ngrid

   norm2 = 0.0d0
   ! Loop over cells
   do ind=1,twotondim
      iskip_mg = (ind-1)*ngrid
      ! Loop over active grids
      do igrid_mg=1,ngrid
         icell_mg = iskip_mg + igrid_mg
         if(active_mg(myid,ilevel)%u(icell_mg,4)<=0.0) cycle
         norm2 = norm2 + active_mg(myid,ilevel)%f(icell_mg,ivar)**2
      end do
   end do
   norm2 = dx2*norm2
end subroutine cmp_fvar_norm2_coarse

! ------------------------------------------------------------------------
! Gauss-Seidel smoothing
! ------------------------------------------------------------------------

subroutine gauss_seidel_mg_coarse(ilevel,safe,redstep)
   use amr_commons
   use pm_commons
   use poisson_commons
   use poisson_commons_acc
   use ind_accessing
   implicit none
   integer, intent(in) :: ilevel

   logical, intent(in) :: safe
   logical, intent(in) :: redstep
   real(dp) :: dx2, nb_sum, weight
   integer  :: ngrid
   integer  :: ind, ind0, igrid_mg, idim, inbor
   integer  :: igrid_amr, cpu_nbor_amr 
   integer  :: iskip_mg, igrid_nbor_mg, icell_mg, icell_nbor_mg
   integer  :: igshift, igrid_nbor_amr
   real(dp) :: dtwondim = (twondim)
#ifdef MPROFILING
  include 'mpif.h'
  integer::info
  real(kind=8)::tt1,tt2
  tt1=MPI_WTIME(info)
#endif

   ! Set constants
   dx2  = (0.5d0**ilevel)**2
   ngrid=ngrid_a_mg(myid,ilevel)

   !$acc data present(son,nbor,cpu_map,lookup_mg,ired,iblack,iii,jjj) &
   !$acc present(active_mg_u,active_mg_f,ngrid_a_mg,active_mg_igrid)

   ! Loop over cells, with red/black ordering
   !$acc parallel loop collapse(2) gang worker vector independent &
   !$acc private(active_ind,tmp_ind)
   ! Loop over active grids
   do igrid_mg=1,ngrid
      do ind0=1,twotondim/2      ! Only half of the cells for a red or black sweep
      if(redstep) then
         ind = ired  (ndim,ind0)
      else
         ind = iblack(ndim,ind0)
      end if

      iskip_mg  = (ind-1)*ngrid

      active_ind = ind_acc_a_mg_igrid(ilevel,myid)+igrid_mg
      igrid_amr = active_mg_igrid(active_ind)
      icell_mg  = iskip_mg  + igrid_mg

      nb_sum=0.0d0                       ! Sum of phi on neighbors
!      ! Read scan flag
!      active_ind = ind_acc_a_mg_f(ilevel,myid)+icell_mg
!      if(.not. btest(active_mg_f(active_ind),0)) then
!      ! Use max-speed "dumb" Gauss-Seidel for "inner" cells
!      ! Those cells are active, have all their neighbors active
!      ! and all neighbors are in the AMR+MG trees
!        do inbor=1,2
!           do idim=1,ndim
!              ! Get neighbor grid shift
!              igshift = iii(idim,inbor,ind)
!              ! Get neighbor grid
!              if(igshift==0) then
!                 igrid_nbor_amr = igrid_amr
!                 cpu_nbor_amr   = myid
!              else
!                 igrid_nbor_amr = son(nbor(igrid_amr,igshift))
!                 cpu_nbor_amr   = cpu_map(nbor(igrid_amr,igshift))
!              end if
!              ! Get neighbor cpu
!              igrid_nbor_mg  = lookup_mg(igrid_nbor_amr)
!              icell_nbor_mg  = igrid_nbor_mg + &
!                      (jjj(idim,inbor,ind)-1)*ngrid_a_mg(cpu_nbor_amr,ilevel)
!
!              active_ind = ind_acc_a_mg_u(ilevel,cpu_nbor_amr,1)+icell_nbor_mg
!              nb_sum = nb_sum + active_mg_u(active_ind)
!           end do
!        end do
!        weight=0.0d0
!      else

! ***** Use directly this branch of the if with more checks ***** !

         ! Use the finer "solve" Gauss-Seidel near boundaries,
         ! with all necessary checks

         active_ind = ind_acc_a_mg_u(ilevel,myid,4)+icell_mg
         if(active_mg_u(active_ind)<=0.0) cycle
         if(safe .and. active_mg_u(active_ind)<1.0) cycle

         weight=0.0d0                       ! Central weight for "Solve G-S"

         do inbor=1,2
            do idim=1,ndim
               ! Get neighbor grid shift
               igshift = iii(idim,inbor,ind)

               ! Get neighbor grid
               if(igshift==0) then
                  igrid_nbor_amr = igrid_amr
                  cpu_nbor_amr   = myid
               else
                  igrid_nbor_amr = son(nbor(igrid_amr,igshift))
                  cpu_nbor_amr   = cpu_map(nbor(igrid_amr,igshift))
               end if

               if(igrid_nbor_amr==0) then
                  ! No neighbor cell, set mask=-1 on nonexistent neighbor cell
                  active_ind = ind_acc_a_mg_u(ilevel,myid,4)+icell_mg
                  weight = weight + (- 1.0d0/active_mg_u(active_ind))
               else
                  ! Fetch neighbor cell
                  igrid_nbor_mg  = lookup_mg(igrid_nbor_amr)
                  if(igrid_nbor_mg<=0) then
                     active_ind = ind_acc_a_mg_u(ilevel,myid,4)+icell_mg
                     ! No MG neighbor
                     weight = weight + (- 1.0d0/active_mg_u(active_ind))
                  else
                     icell_nbor_mg  = igrid_nbor_mg  + (jjj(idim,inbor,ind)-1)*ngrid_a_mg(cpu_nbor_amr,ilevel)
                     active_ind = ind_acc_a_mg_u(ilevel,cpu_nbor_amr,4)+icell_nbor_mg
                     if(active_mg_u(active_ind)<=0.0) then
                        active_ind = ind_acc_a_mg_u(ilevel,cpu_nbor_amr,4)+icell_nbor_mg
                        tmp_ind    = ind_acc_a_mg_u(ilevel,myid,4)+icell_mg
                        ! Neighbor cell is masked
                        weight = weight + active_mg_u(active_ind)/active_mg_u(tmp_ind)
                      else
                        ! Neighbor cell is active, increment neighbor sum
                        active_ind = ind_acc_a_mg_u(ilevel,cpu_nbor_amr,1)+icell_nbor_mg
                        nb_sum = nb_sum + active_mg_u(active_ind)
                     end if
                  end if
               end if
            end do
         end do
!      end if
       ! Update the potential, solving for potential on icell_amr
       active_ind = ind_acc_a_mg_u(ilevel,myid,1)+icell_mg
       tmp_ind    = ind_acc_a_mg_u(ilevel,myid,2)+icell_mg
       active_mg_u(active_ind)=(nb_sum-dx2*active_mg_u(tmp_ind))/(dtwondim - weight)
      end do
   end do

   !$acc end data
   
#ifdef MPROFILING
  tt2=MPI_WTIME(info)
  acc_t_gauss_seidel_mg_coarse = acc_t_gauss_seidel_mg_coarse + (tt2-tt1)
#endif

end subroutine gauss_seidel_mg_coarse

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
         igrid_c_amr = active_mg(myid,icoarselevel)%igrid(igrid_c_mg)
         icell_c_amr = igrid_c_amr + iskip_c_amr
         cpu_amr     = cpu_map(icell_c_amr)
         icell_c_mg  = igrid_c_mg  + iskip_c_mg

         ! Get AMR child grid
         igrid_f_amr = son(icell_c_amr)
         if(igrid_f_amr==0) then
            active_mg(myid,icoarselevel)%u(icell_c_mg,2) = 0.0d0    ! Nullify residual (coarser RHS)
            cycle
         end if

         ! Get child MG grid id
         igrid_f_mg = lookup_mg(igrid_f_amr)
         if(igrid_f_mg<=0) then
            ! Son grid is not in MG hierarchy
            active_mg(myid,icoarselevel)%u(icell_c_mg,2) = 0.0d0    ! Nullify residual (coarser RHS)
            cycle
         end if

         ! Loop over child (fine MG) cells
         val = 0.0d0
         w = 0d0
         do ind_f=1,twotondim
            icell_f_mg = igrid_f_mg + (ind_f-1)*active_mg(cpu_amr,ifinelevel)%ngrid

            if (active_mg(cpu_amr,ifinelevel)%u(icell_f_mg,4)<=0.0) cycle
            val = val + active_mg(cpu_amr,ifinelevel)%u(icell_f_mg,3)
            w = w + 1d0
         end do
         ! Store restricted residual into RHS of coarse level
         if(w>0) then
            active_mg(myid,icoarselevel)%u(icell_c_mg,2) = val/w
         else
            active_mg(myid,icoarselevel)%u(icell_c_mg,2) = 0d0
         end if
      end do
   end do
end subroutine restrict_residual_coarse

! ------------------------------------------------------------------------
! Residual restriction (bottom-up)
! ------------------------------------------------------------------------

subroutine restrict_residual_coarse_reverse(ifinelevel)
   use amr_commons
   use poisson_commons
   use poisson_commons_acc
   use ind_accessing
   implicit none
   integer, intent(in) :: ifinelevel

   integer :: ind_c_cell, ind_f_cell, cpu_amr
   integer :: iskip_c_mg
   integer :: igrid_c_amr, igrid_c_mg
   integer :: icell_c_amr, icell_c_mg
   integer :: iskip_f_mg
   integer :: igrid_f_amr, igrid_f_mg
   integer :: icell_f_mg
   real(dp) :: res
   real(dp) :: dtwotondim = (twotondim)

   integer :: icoarselevel
#ifdef MPROFILING
  include 'mpif.h'
  integer::info
  real(kind=8)::tt1,tt2
  tt1=MPI_WTIME(info)
#endif
   icoarselevel=ifinelevel-1

   !$acc data present(father,cpu_map,lookup_mg,active_mg_u,ngrid_a_mg,active_mg_igrid) &
   !$acc create(active_ind)

   ! Loop over fine cells of the myid active comm
   !$acc parallel loop collapse(2) gang vector
   do ind_f_cell=1,twotondim
      ! Loop over fine grids of myid
      do igrid_f_mg=1,ngrid_a_mg(myid,ifinelevel)
         iskip_f_mg =(ind_f_cell-1)*ngrid_a_mg(myid,ifinelevel)
         icell_f_mg=iskip_f_mg+igrid_f_mg

         ! Is fine cell masked?
         active_ind = ind_acc_a_mg_u(ifinelevel,myid,4)+icell_f_mg
         if(active_mg_u(active_ind)<=0d0) cycle

         ! Get coarse grid AMR index and CPU id
         active_ind = ind_acc_a_mg_igrid(ifinelevel,myid)+igrid_f_mg
         igrid_f_amr=active_mg_igrid(active_ind)
         icell_c_amr=father(igrid_f_amr)
         ind_c_cell=(icell_c_amr-ncoarse-1)/ngridmax+1
         igrid_c_amr=icell_c_amr-ncoarse-(ind_c_cell-1)*ngridmax
         cpu_amr=cpu_map(father(igrid_c_amr))

         ! Convert to MG index, get MG coarse cell id
         igrid_c_mg=lookup_mg(igrid_c_amr)
         iskip_c_mg=(ind_c_cell-1)*ngrid_a_mg(cpu_amr,icoarselevel)
         icell_c_mg=iskip_c_mg+igrid_c_mg

         ! Is coarse cell masked?
         active_ind = ind_acc_a_mg_u(icoarselevel,cpu_amr,4)+icell_c_mg
         if(active_mg_u(active_ind)<=0d0) cycle

         ! Stack fine cell residual in coarse cell rhs
         active_ind = ind_acc_a_mg_u(ifinelevel,myid,3)+icell_f_mg
         res=active_mg_u(active_ind)/dtwotondim

         active_ind = ind_acc_a_mg_u(icoarselevel,cpu_amr,2)+icell_c_mg
         !$acc atomic update
         active_mg_u(active_ind)=active_mg_u(active_ind)+res
      end do
   end do

   !$acc end data
   
#ifdef MPROFILING
  tt2=MPI_WTIME(info)
  acc_t_restrict_residual_coarse_reverse = acc_t_restrict_residual_coarse_reverse + (tt2-tt1)
#endif

end subroutine restrict_residual_coarse_reverse

! ------------------------------------------------------------------------
! Interpolation and correction
! ------------------------------------------------------------------------

subroutine interpolate_and_correct_coarse(ifinelevel)
   use amr_commons
   use poisson_commons
   use poisson_commons_acc
   use ind_accessing
   use iacv
   use iacv_subroutines
   implicit none
   integer, intent(in) :: ifinelevel

   integer  :: i, ind_father, ind_average, ind_f, iskip_f_amr, ngrid_f, istart
   integer  :: icell_c_amr, igrid_c_amr, igrid_c_mg, icell_c_mg, iskip_f_mg, icell_f_mg
   integer  :: icoarselevel, ind_c, cpu_c_amr
   real(dp) :: coeff,corr_aux
   integer  ::nbin,ibin
   integer ,dimension(:,:),allocatable::ix,iy,iz,iix,iiy,iiz
#ifdef MPROFILING
  include 'mpif.h'
  integer::info
  real(kind=8)::tt1,tt2
  tt1=MPI_WTIME(info)
#endif

   icoarselevel=ifinelevel-1
   ngrid_f=ngrid_a_mg(myid,ifinelevel)
   nbin = int(ngrid_f/nvector)+1

   allocate(nbatch(1:nbin))
   allocate(igrid_f_amr(1:nbin,1:nvector))
   allocate(icell_amr(1:nbin,1:nvector),cpu_amr_aux(1:nbin,1:nvector))
   allocate(nbors_father_cells(1:nbin,1:nvector,1:threetondim))
   allocate(nbors_father_grids(1:nbin,1:nvector,1:twotondim))
   allocate(pos(1:nbin,1:nvector),ind_grid_father(1:nbin,1:nvector),ind_grid_ok(1:nbin,1:nvector))
   allocate(nbors_father_ok(1:nbin,1:nvector,1:threetondim),nbors_grids_ok(1:nbin,1:nvector,1:twotondim))
   allocate(ix(1:nbin,1:nvector),iy(1:nbin,1:nvector),iz(1:nbin,1:nvector),iix(1:nbin,1:nvector),iiy(1:nbin,1:nvector),iiz(1:nbin,1:nvector))
   allocate(ind_grid1(1:nbin,1:nvector),ind_grid2(1:nbin,1:nvector),ind_grid3(1:nbin,1:nvector))
   allocate(nbors_grids(1:nbin,1:nvector,1:twotondim))

   !$acc data present(father,cpu_map,lookup_mg,nbor,son) &
   !$acc present(active_mg_u,ngrid_a_mg,active_mg_igrid,bbb,ccc) &
   !$acc create(nbatch,igrid_f_amr,icell_amr,nbors_father_cells,nbors_father_grids) &
   !$acc create(pos,ind_grid_father,ind_grid_ok,nbors_father_ok,nbors_grids_ok) &
   !$acc create(ind_grid1,ind_grid2,ind_grid3,nbors_grids,cpu_amr_aux) &
   !$acc create(ix,iy,iz,iix,iiy,iiz) &
   !$acc present(iii_aux,jjj_aux,kkk_aux,lll,mmm)

   ! Loop over fine grids by vector sweeps
   !$acc parallel loop gang worker independent private(active_ind)
   do istart=1,ngrid_f,nvector
      ibin = int(istart/nvector)+1
      nbatch(ibin)=MIN(nvector,ngrid_f-istart+1)

      !$acc loop vector
      do i=1,nbatch(ibin)
      ! Gather nvector grids
         active_ind = ind_acc_a_mg_igrid(ifinelevel,myid)+(istart+i-1)
         igrid_f_amr(ibin,i)=active_mg_igrid(active_ind)
      ! Compute father (coarse) cell index
         icell_amr(ibin,i)=father(igrid_f_amr(ibin,i))
         cpu_amr_aux(ibin,i)=cpu_map(icell_amr(ibin,i))
      end do

      ! Gather 3x3x3 neighboring parent cells
      call get3cubefather_iacv(ifinelevel,nbin,ibin,ix,iy,iz,iix,iiy,iiz)
   end do

   corr_aux = 0.0d0

   !$acc parallel loop gang independent private(active_ind)
   do istart=1,ngrid_f,nvector
   ibin = int(istart/nvector)+1
      !$acc loop collapse(2) vector reduction(+:corr_aux)
      do ind_f=1,twotondim
      do i=1,nbatch(ibin)
      corr_aux = 0.0d0
      iskip_f_mg  = (ind_f-1)*ngrid_f
      iskip_f_amr = ncoarse+(ind_f-1)*ngridmax
      icell_amr(ibin,i) = iskip_f_amr + igrid_f_amr(ibin,i)
         ! Loop over relevant parent cells
         do ind_average=1,twotondim
            ind_father = ccc(ind_average,ind_f)
            coeff      = bbb(ind_average)
            icell_f_mg  = iskip_f_mg + istart+i-1
            active_ind = ind_acc_a_mg_u(ifinelevel,cpu_amr_aux(ibin,i),4)+icell_f_mg
            if(active_mg_u(active_ind)<=0.0) then
               corr_aux=0.0d0        ! Fine cell is masked : no correction
               cycle
            end if
            icell_c_amr = nbors_father_cells(ibin,i,ind_father)
            ind_c       = (icell_c_amr-ncoarse)/ngridmax + 1
            igrid_c_amr = icell_c_amr - ncoarse - (ind_c-1)*ngridmax
            igrid_c_mg  = lookup_mg(igrid_c_amr)
            cpu_c_amr   = cpu_map(father(igrid_c_amr))
            if(igrid_c_mg<=0) cycle

            icell_c_mg  = (ind_c-1)*ngrid_a_mg(cpu_c_amr,icoarselevel) + igrid_c_mg
            active_ind = ind_acc_a_mg_u(icoarselevel,cpu_c_amr,1)+icell_c_mg
            corr_aux=corr_aux+coeff*active_mg_u(active_ind)
         end do
         ! Correct potential        
         icell_f_mg  = iskip_f_mg + istart+i-1
         active_ind = ind_acc_a_mg_u(ifinelevel,cpu_amr_aux(ibin,i),1)+icell_f_mg
         active_mg_u(active_ind)=active_mg_u(active_ind)+corr_aux
      end do
      end do
      ! End loop over cells
   end do
   ! End loop over grids

   !$acc end data

   deallocate(nbatch)
   deallocate(igrid_f_amr,icell_amr,cpu_amr_aux)
   deallocate(nbors_father_cells,nbors_father_grids)
   deallocate(pos,ind_grid_father,ind_grid_ok)
   deallocate(nbors_father_ok,nbors_grids_ok)
   deallocate(ind_grid1,ind_grid2,ind_grid3)
   deallocate(nbors_grids)
   deallocate(ix,iy,iz,iix,iiy,iiz)
   
#ifdef MPROFILING
  tt2=MPI_WTIME(info)
  acc_t_interpolate_and_correct_coarse = acc_t_interpolate_and_correct_coarse + (tt2-tt1)
#endif

end subroutine interpolate_and_correct_coarse

! ------------------------------------------------------------------------
! Flag setting
! ------------------------------------------------------------------------

subroutine set_scan_flag_coarse(ilevel)
   use amr_commons
   use poisson_commons
   use poisson_commons_acc
   use ind_accessing
   implicit none
   integer, intent(in) :: ilevel

   integer :: ind, ngrid, scan_flag
   integer :: igrid_mg, inbor, idim, igshift
   integer :: igrid_amr, igrid_nbor_amr, cpu_nbor_amr, icell_nbor_amr

   integer :: iskip_mg, icell_mg, igrid_nbor_mg, icell_nbor_mg
#ifdef MPROFILING
  include 'mpif.h'
  integer::info
  real(kind=8)::tt1,tt2
  tt1=MPI_WTIME(info)
#endif

   ngrid = ngrid_a_mg(myid,ilevel)
   if(ngrid==0) return

   !$acc data present(lookup_mg,son,nbor,cpu_map,iii,jjj) &
   !$acc      present(active_mg_u,ngrid_a_mg,active_mg_igrid)

   ! Loop over cells and set coarse SCAN flag
   !$acc parallel loop independent
   do ind=1,twotondim
      !$acc loop vector
      do igrid_mg=1,ngrid
         iskip_mg  = (ind-1)*ngrid
         active_ind = ind_acc_a_mg_igrid(ilevel,myid)+igrid_mg
         igrid_amr = active_mg_igrid(active_ind)
         icell_mg  = iskip_mg  + igrid_mg

         active_ind = ind_acc_a_mg_u(ilevel,myid,4)+icell_mg
         if(active_mg_u(active_ind)==1d0) then
            scan_flag=0       ! Init flag to 'no scan needed'
            scan_flag_loop: do inbor=1,2
               do idim=1,ndim
                  igshift = iii(idim,inbor,ind)
                  if(igshift==0) then
                     igrid_nbor_amr = igrid_amr
                     cpu_nbor_amr   = myid
                  else
                     igrid_nbor_amr = son(nbor(igrid_amr,igshift))
                     cpu_nbor_amr   = cpu_map(nbor(igrid_amr,igshift))
                  end if

                  if(igrid_nbor_amr==0) then
                     scan_flag=1
                     exit scan_flag_loop
                  else
                     igrid_nbor_mg = lookup_mg(igrid_nbor_amr)
                     if(igrid_nbor_mg<=0) then
                        scan_flag=1
                        exit scan_flag_loop
                     else
                        icell_nbor_mg  = igrid_nbor_mg  + &
                                 (jjj(idim,inbor,ind)-1)*ngrid_a_mg(cpu_nbor_amr,ilevel)

                        active_ind = ind_acc_a_mg_u(ilevel,cpu_nbor_amr,4)+icell_nbor_mg
                        if(active_mg_u(active_ind)<=0.0) then
                           scan_flag=1
                           exit scan_flag_loop
                        end if
                     end if
                  end if
               end do
            end do scan_flag_loop
         else
            scan_flag=1
         end if
         active_ind = ind_acc_a_mg_f(ilevel,myid)+icell_mg
         active_mg_f(active_ind)=scan_flag
      end do
   end do

   !$acc end data
 
#ifdef MPROFILING
  tt2=MPI_WTIME(info)
  acc_t_set_scan_flag_coarse = acc_t_set_scan_flag_coarse + (tt2-tt1)
#endif  
   
end subroutine set_scan_flag_coarse
