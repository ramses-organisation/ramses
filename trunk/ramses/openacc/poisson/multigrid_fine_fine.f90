! ------------------------------------------------------------------------
! Multigrid Poisson solver for refined AMR levels
! ------------------------------------------------------------------------
! This file contains all MG-fine-level related routines
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
         igrid_c_amr=active_mg(myid,icoarselevel)%igrid(igrid_c_mg)
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
         active_mg(myid,icoarselevel)%u(icell_c_mg,4)=ngpmask
         allmasked=allmasked .and. (ngpmask<=0.0)
      end do
   end do

end subroutine restrict_mask_fine

! ------------------------------------------------------------------------
! Mask restriction (bottom-up)
! ------------------------------------------------------------------------

subroutine restrict_mask_fine_reverse(ifinelevel)
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
   integer :: iskip_f_amr
   integer :: igrid_f_amr, igrid_f_mg
   integer :: icell_f_amr
   real(dp) :: dtwotondim = (twotondim)
   integer :: icoarselevel
   real(dp),dimension(1:twotondim,1:active(ifinelevel)%ngrid) :: ngpmask
   integer,dimension(1:twotondim,1:active(ifinelevel)%ngrid)   :: active_ind_aux
#ifdef MPROFILING
  include 'mpif.h'
  integer::info
  real(kind=8)::tt1,tt2
  tt1=MPI_WTIME(info)
#endif

   icoarselevel=ifinelevel-1

   !$acc data present(active,father,cpu_map,f,lookup_mg) &
   !$acc present(active_mg_u,ngrid_a_mg) create(ngpmask,active_ind_aux)

   ! Loop over fine cells of the myid active comm
   !$acc parallel loop collapse(2) gang vector
   do ind_f_cell=1,twotondim
      ! Loop over fine grids of myid
      do igrid_f_mg=1,active(ifinelevel)%ngrid
         iskip_f_amr=ncoarse+(ind_f_cell-1)*ngridmax
         igrid_f_amr=active(ifinelevel)%igrid(igrid_f_mg)
         icell_f_amr=igrid_f_amr+iskip_f_amr

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
         ngpmask(ind_f_cell,igrid_f_mg)=(1d0+f(icell_f_amr,3))/2d0/dtwotondim
         active_ind_aux(ind_f_cell,igrid_f_mg) = & 
         & ind_acc_a_mg_u(icoarselevel,cpu_amr,4)+icell_c_mg
      end do
   end do

   ! Loop over fine cells of the myid active comm
   !$acc parallel
   do ind_f_cell=1,twotondim
      ! Loop over fine grids of myid
      !$acc loop
      do igrid_f_mg=1,active(ifinelevel)%ngrid
         active_mg_u(active_ind_aux(ind_f_cell,igrid_f_mg)) = & 
         & active_mg_u(active_ind_aux(ind_f_cell,igrid_f_mg))+ngpmask(ind_f_cell,igrid_f_mg)
      end do
   end do
   !$acc end parallel

  !$acc end data
  
#ifdef MPROFILING
  tt2=MPI_WTIME(info)
  acc_t_restrict_mask_fine_reverse = acc_t_restrict_mask_fine_reverse + (tt2-tt1)
#endif

end subroutine restrict_mask_fine_reverse

! ------------------------------------------------------------------------
! Residual computation
! ------------------------------------------------------------------------

subroutine cmp_residual_mg_fine(ilevel)
   ! Computes the residual the fine (AMR) level, and stores it into f(:,1)
   use amr_commons
   use poisson_commons
   use ind_accessing
   implicit none
   integer, intent(in) :: ilevel

   integer, dimension(1:3,1:2,1:8) :: iii, jjj
   real(dp) :: dx, oneoverdx2, phi_c, nb_sum
   integer  :: ngrid
   integer  :: ind, igrid_mg, idim, inbor
   integer  :: igrid_amr, icell_amr, iskip_amr
   integer  :: igshift, igrid_nbor_amr, icell_nbor_amr

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

   iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
   iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

   !$acc data present(f,phi,nbor,son,active,flag2) copyin(iii,jjj)

   ngrid=active(ilevel)%ngrid

   ! Loop over cells
   !$acc parallel loop collapse(2) gang worker vector independent
   ! Loop over active grids
   do igrid_mg=1,ngrid
      do ind=1,twotondim
         iskip_amr = ncoarse+(ind-1)*ngridmax
         igrid_amr = active(ilevel)%igrid(igrid_mg)
         icell_amr = iskip_amr + igrid_amr
         phi_c = phi(icell_amr)           ! Value of potential on center cell

         nb_sum=0.0d0                     ! Sum of phi on neighbors

!         ! SCAN FLAG TEST
!         if(flag2(icell_amr)/ngridmax==0) then
!            do inbor=1,2
!               do idim=1,ndim
!                  ! Get neighbor grid
!                  igshift = iii(idim,inbor,ind)
!                  if(igshift==0) then
!                     igrid_nbor_amr = igrid_amr
!                  else
!                     igrid_nbor_amr = son(nbor(igrid_amr,igshift))
!                  end if
!                  icell_nbor_amr = igrid_nbor_amr + &
!                      (ncoarse + (jjj(idim,inbor,ind)-1)*ngridmax)
!                  nb_sum = nb_sum + phi(icell_nbor_amr)
!               end do
!            end do
!         else ! PERFORM SCAN
            if(f(icell_amr,3)<=0.0) then
               f(icell_amr,1)=0.0d0
               cycle
            end if
            do idim=1,ndim
               do inbor=1,2
                  ! Get neighbor grid
                  igshift = iii(idim,inbor,ind)
                  if(igshift==0) then
                     igrid_nbor_amr = igrid_amr
                  else
                     igrid_nbor_amr = son(nbor(igrid_amr,igshift))
                  end if

                  if(igrid_nbor_amr==0) then
                     ! No neighbor cell !
                     ! Virtual phi value on unrefined neighbor cell :
                     ! -phi_c/mask_c
                     ! (simulates mask=-1.0 for the nonexistent refined cell)
                     nb_sum = nb_sum - phi_c/f(icell_amr,3)
                  else
                     ! Fetch neighbor cell
                     icell_nbor_amr = igrid_nbor_amr + &
                         (ncoarse + (jjj(idim,inbor,ind)-1)*ngridmax)
                     if(f(icell_nbor_amr,3)<=0.0) then
                        ! Neighbor cell is masked :
                        ! compute its virtual phi with the mask
                        nb_sum = nb_sum + &
                            phi_c*(f(icell_nbor_amr,3)/f(icell_amr,3))
                     else
                        ! Neighbor cell is active, use its true potential
                        nb_sum = nb_sum + phi(icell_nbor_amr)
                     end if
                  end if
               end do
            end do
!         end if ! END SCAN TEST

         ! Store ***MINUS THE RESIDUAL*** in f(:,1), using BC-modified RHS
         f(icell_amr,1) = -oneoverdx2*( nb_sum - dtwondim*phi_c )+f(icell_amr,2)
      end do
   end do
 
   !$acc end data
   
#ifdef MPROFILING
  tt2=MPI_WTIME(info)
  acc_t_cmp_residual_mg_fine = acc_t_cmp_residual_mg_fine + (tt2-tt1)
#endif   
   
end subroutine cmp_residual_mg_fine

! ##################################################################
! ##################################################################

subroutine cmp_residual_norm2_fine(ilevel, norm2)
   use amr_commons
   use poisson_commons
   use ind_accessing
   implicit none
   integer,  intent(in)  :: ilevel
   real(kind=8), intent(out) :: norm2

   real(kind=8) :: dx2
   integer  :: ngrid
   integer  :: ind, igrid_mg
   integer  :: igrid_amr, icell_amr, iskip_amr
#ifdef MPROFILING
  include 'mpif.h'
  integer::info
  real(kind=8)::tt1,tt2
  tt1=MPI_WTIME(info)
#endif

   ! Set constants
   dx2  = (0.5d0**ilevel)**2
   ngrid=active(ilevel)%ngrid
   norm2 = 0.0d0

  !$acc data present(active,f) copy(norm2)

   ! Race conditions for collapse
   ! Loop over cells
   !$acc parallel loop collapse(2) gang worker vector independent reduction(+:norm2)
   ! Loop over active grids
   do igrid_mg=1,ngrid
      do ind=1,twotondim
         iskip_amr = ncoarse+(ind-1)*ngridmax
         igrid_amr = active(ilevel)%igrid(igrid_mg)
         icell_amr = iskip_amr + igrid_amr
         if(f(icell_amr,3)<=0.0) then      ! Do not count masked cells
            cycle
         end if
         norm2 = norm2 + f(icell_amr,1)**2
      end do
   end do

   !$acc end data
   norm2 = dx2*norm2

#ifdef MPROFILING
  tt2=MPI_WTIME(info)
  acc_t_cmp_residual_norm2_fine = acc_t_cmp_residual_norm2_fine + (tt2-tt1)
#endif 

end subroutine cmp_residual_norm2_fine

! ------------------------------------------------------------------------
! Gauss-Seidel smoothing
! ------------------------------------------------------------------------

subroutine gauss_seidel_mg_fine(ilevel,redstep)
   use amr_commons
   use pm_commons
   use poisson_commons
   use poisson_commons_acc
   use ind_accessing
   implicit none
   integer, intent(in) :: ilevel
   logical, intent(in) :: redstep

   real(dp) :: dx2, nb_sum, weight
   integer  :: ngrid
   integer  :: ind, ind0, igrid_mg, idim, inbor
   integer  :: igrid_amr, icell_amr, iskip_amr
   integer  :: igshift, igrid_nbor_amr, icell_nbor_amr

   real(dp) :: dtwondim = (twondim)
#ifdef MPROFILING
  include 'mpif.h'
  integer::info
  real(kind=8)::tt1,tt2
  tt1=MPI_WTIME(info)
#endif

   ! Set constants
   dx2  = (0.5d0**ilevel)**2

   ngrid=active(ilevel)%ngrid
   weight=0.0d0
   nb_sum=0.0d0

   !$acc data present(phi,f,active,nbor,son,flag2,ired,iblack,iii,jjj)

   ! Loop over cells, with red/black ordering
   !$acc parallel loop collapse(2) gang worker vector independent
   ! Loop over active grids
   do igrid_mg=1,ngrid
      do ind0=1,twotondim/2      ! Only half of the cells for a red or black sweep
      if(redstep) then
         ind = ired  (ndim,ind0)
      else
         ind = iblack(ndim,ind0)
      end if

         iskip_amr = ncoarse+(ind-1)*ngridmax
         igrid_amr = active(ilevel)%igrid(igrid_mg)
         icell_amr = iskip_amr + igrid_amr

         nb_sum=0.0d0                       ! Sum of phi on neighbors

!         ! Read scan flag
!         if(flag2(icell_amr)/ngridmax==0) then
!            ! Use max-speed "dumb" Gauss-Seidel for "inner" cells
!            ! Those cells are active, have all their neighbors active
!           ! and all neighbors are in the AMR+MG trees
!            do inbor=1,2
!               do idim=1,ndim
!                  ! Get neighbor grid shift
!                  igshift = iii(idim,inbor,ind)
!                  ! Get neighbor grid
!                  if(igshift==0) then
!                     igrid_nbor_amr = igrid_amr
!                  else
!                     igrid_nbor_amr = son(nbor(igrid_amr,igshift))
!                  end if
!                  icell_nbor_amr = igrid_nbor_amr + &
!                      (ncoarse + (jjj(idim,inbor,ind)-1)*ngridmax)
!                  nb_sum = nb_sum + phi(icell_nbor_amr)
!               end do
!            end do
!            weight=0.0d0
!         else
            ! Use the finer "solve" Gauss-Seidel near boundaries,
            ! with all necessary checks
            if (f(icell_amr,3)<=0.0) cycle
            if (safe_mode(ilevel) .and. f(icell_amr,3)<1.0) cycle

            weight=0.0d0 ! Central weight for "Solve G-S"
            do inbor=1,2
               do idim=1,ndim
                  ! Get neighbor grid shift
                  igshift = iii(idim,inbor,ind)

                  ! Get neighbor grid
                  if(igshift==0) then
                     igrid_nbor_amr = igrid_amr
                  else
                     igrid_nbor_amr = son(nbor(igrid_amr,igshift))
                  end if

                  if(igrid_nbor_amr==0) then
                     ! No neighbor cell,
                     ! set mask=-1 on nonexistent neighbor cell
                     weight = weight + (- 1.0d0/f(icell_amr,3))
                  else
                     ! Fetch neighbor cell
                     icell_nbor_amr = igrid_nbor_amr + &
                         (ncoarse + (jjj(idim,inbor,ind)-1)*ngridmax)
                     if(f(icell_nbor_amr,3)<=0.0) then
                        ! Neighbor cell is masked
                        weight = weight + f(icell_nbor_amr,3)/f(icell_amr,3)
                     else
                        ! Neighbor cell is active, increment neighbor sum
                        nb_sum = nb_sum + phi(icell_nbor_amr)
                     end if
                  end if
               end do
            end do
!         end if
            ! Update the potential, solving for potential on icell_amr
            phi(icell_amr) = (nb_sum - dx2*f(icell_amr,2)) / (dtwondim - weight)
      end do
   end do

   !$acc end data
   
#ifdef MPROFILING
  tt2=MPI_WTIME(info)
  acc_t_gauss_seidel_mg_fine = acc_t_gauss_seidel_mg_fine + (tt2-tt1)
#endif  
   
end subroutine gauss_seidel_mg_fine

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
         igrid_c_amr = active_mg(myid,icoarselevel)%igrid(igrid_c_mg)
         icell_c_amr = igrid_c_amr + iskip_c_amr
         icell_c_mg  = igrid_c_mg  + iskip_c_mg

         ! Get AMR child grid
         igrid_f_amr = son(icell_c_amr)
         if(igrid_f_amr==0) then
            ! Nullify residual (coarser RHS)
            active_mg(myid,icoarselevel)%u(icell_c_mg,2) = 0.0d0
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
         active_mg(myid,icoarselevel)%u(icell_c_mg,2) = val/dtwotondim
      end do
   end do
end subroutine restrict_residual_fine


! ------------------------------------------------------------------------
! Residual restriction (bottom-up)
! ------------------------------------------------------------------------

subroutine restrict_residual_fine_reverse(ifinelevel)
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
   integer :: iskip_f_amr
   integer :: igrid_f_amr, igrid_f_mg
   integer :: icell_f_amr
   real(dp):: dtwotondim = (twotondim)
   integer :: icoarselevel
   integer ,dimension(1:twotondim,1:active(ifinelevel)%ngrid):: active_ind_aux
   real(dp),dimension(1:twotondim,1:active(ifinelevel)%ngrid):: res
#ifdef MPROFILING
  include 'mpif.h'
  integer::info
  real(kind=8)::tt1,tt2
  tt1=MPI_WTIME(info)
#endif

   icoarselevel=ifinelevel-1

   !$acc data present(active,father,cpu_map,lookup_mg,f) &
   !$acc present(active_mg_u,ngrid_a_mg) create(active_ind_aux,active_ind,res)

   ! Loop over fine cells of the myid active comm
   !$acc parallel loop collapse(2) gang vector
   do ind_f_cell=1,twotondim
      ! Loop over fine grids of myid
      do igrid_f_mg=1,active(ifinelevel)%ngrid
         iskip_f_amr=ncoarse+(ind_f_cell-1)*ngridmax
         igrid_f_amr=active(ifinelevel)%igrid(igrid_f_mg)
         icell_f_amr=igrid_f_amr+iskip_f_amr
         ! Is fine cell masked?
         if(f(icell_f_amr,3)<=0d0) cycle

         ! Get coarse grid AMR index and CPU id
         icell_c_amr=father(igrid_f_amr)
         ind_c_cell=(icell_c_amr-ncoarse-1)/ngridmax+1
         igrid_c_amr=icell_c_amr-ncoarse-(ind_c_cell-1)*ngridmax
         cpu_amr=cpu_map(father(igrid_c_amr))

         ! Convert to MG index, get MG coarse cell id
         igrid_c_mg=lookup_mg(igrid_c_amr)
         iskip_c_mg=(ind_c_cell-1)*ngrid_a_mg(cpu_amr,icoarselevel)
         icell_c_mg=iskip_c_mg+igrid_c_mg

         active_ind = ind_acc_a_mg_u(icoarselevel,cpu_amr,4)+icell_c_mg
         ! Is coarse cell masked?
         if(active_mg_u(active_ind)<=0d0) cycle

         ! Stack fine cell residual in coarse cell rhs
         res(ind_f_cell,igrid_f_mg)=f(icell_f_amr,1)/dtwotondim

         active_ind_aux(ind_f_cell,igrid_f_mg) = &
         & ind_acc_a_mg_u(icoarselevel,cpu_amr,2)+icell_c_mg
      end do
   end do

   !$acc parallel
   do ind_f_cell=1,twotondim
      !$acc loop
      do igrid_f_mg=1,active(ifinelevel)%ngrid
         active_mg_u(active_ind_aux(ind_f_cell,igrid_f_mg))= & 
         & active_mg_u(active_ind_aux(ind_f_cell,igrid_f_mg))+res(ind_f_cell,igrid_f_mg)
      end do
   end do
   !$acc end parallel

   !$acc end data
   
#ifdef MPROFILING
  tt2=MPI_WTIME(info)
  acc_t_restrict_residual_fine_reverse = acc_t_restrict_residual_fine_reverse + (tt2-tt1)
#endif   

end subroutine restrict_residual_fine_reverse

! ------------------------------------------------------------------------
! Interpolation and correction
! ------------------------------------------------------------------------

subroutine interpolate_and_correct_fine(ifinelevel)
   use amr_commons
   use poisson_commons
   use poisson_commons_acc
   use iacv
   use iacv_subroutines
   use ind_accessing
   implicit none
   integer, intent(in) :: ifinelevel

   integer  :: i, ind_father, ind_average, ind_f, iskip_f_amr
   integer  :: ngrid_f, istart
   integer  :: icell_c_amr, igrid_c_amr, igrid_c_mg, icell_c_mg
   integer  :: icoarselevel, ind_c, cpu_amr
   real(dp) :: coeff, corr_aux
   integer::nbin,ibin
   integer ,dimension(:,:),allocatable::ix,iy,iz,iix,iiy,iiz
#ifdef MPROFILING
  include 'mpif.h'
  integer::info
  real(kind=8)::tt1,tt2
  tt1=MPI_WTIME(info)
#endif

   icoarselevel=ifinelevel-1
   ngrid_f=active(ifinelevel)%ngrid
   nbin = int(ngrid_f/nvector)+1

   allocate(nbatch(1:nbin))
   allocate(igrid_f_amr(1:nbin,1:nvector))
   allocate(icell_amr(1:nbin,1:nvector))
   allocate(nbors_father_cells(1:nbin,1:nvector,1:threetondim))
   allocate(nbors_father_grids(1:nbin,1:nvector,1:twotondim))
   allocate(pos(1:nbin,1:nvector),ind_grid_father(1:nbin,1:nvector),ind_grid_ok(1:nbin,1:nvector))
   allocate(nbors_father_ok(1:nbin,1:nvector,1:threetondim),nbors_grids_ok(1:nbin,1:nvector,1:twotondim))
   allocate(ix(1:nbin,1:nvector),iy(1:nbin,1:nvector),iz(1:nbin,1:nvector),iix(1:nbin,1:nvector),iiy(1:nbin,1:nvector),iiz(1:nbin,1:nvector))
   allocate(ind_grid1(1:nbin,1:nvector),ind_grid2(1:nbin,1:nvector),ind_grid3(1:nbin,1:nvector))
   allocate(nbors_grids(1:nbin,1:nvector,1:twotondim))

   !$acc data present(phi,f,active,father,cpu_map,lookup_mg,bbb,ccc,nbor,son) &
   !$acc create(nbatch,igrid_f_amr,icell_amr,nbors_father_cells,nbors_father_grids) &
   !$acc create(pos,ind_grid_father,ind_grid_ok,nbors_father_ok,nbors_grids_ok) &
   !$acc create(ind_grid1,ind_grid2,ind_grid3,nbors_grids) &
   !$acc create(ix,iy,iz,iix,iiy,iiz) &
   !$acc present(active_mg_u,ngrid_a_mg) &
   !$acc present(iii_aux,jjj_aux,kkk_aux,lll,mmm)

   ! Loop over fine grids by vector sweeps
   !$acc parallel loop gang worker independent
   do istart=1,ngrid_f,nvector
      ibin = int(istart/nvector)+1
      nbatch(ibin)=MIN(nvector,ngrid_f-istart+1)

      !$acc loop vector
      do i=1,nbatch(ibin)
      ! Gather nvector grids
         igrid_f_amr(ibin,i)=active(ifinelevel)%igrid(istart+i-1)
      ! Compute father (coarse) cell index
         icell_amr(ibin,i)=father(igrid_f_amr(ibin,i))
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
      iskip_f_amr = ncoarse+(ind_f-1)*ngridmax
      icell_amr(ibin,i) = iskip_f_amr + igrid_f_amr(ibin,i)
         ! Loop over relevant parent cells
         do ind_average=1,twotondim
            ind_father = ccc(ind_average,ind_f)
            coeff      = bbb(ind_average)
            if(f(icell_amr(ibin,i),3)<=0.0) then
               corr_aux=0.0d0        ! Fine cell is masked : no correction
               cycle
            end if
            icell_c_amr = nbors_father_cells(ibin,i,ind_father)
            ind_c       = (icell_c_amr-ncoarse-1)/ngridmax + 1
            igrid_c_amr = icell_c_amr - ncoarse - (ind_c-1)*ngridmax
            cpu_amr     = cpu_map(father(igrid_c_amr))
            igrid_c_mg  = lookup_mg(igrid_c_amr)
            if(igrid_c_mg<=0) cycle

            icell_c_mg=(ind_c-1)*ngrid_a_mg(cpu_amr,icoarselevel)+igrid_c_mg
            active_ind = ind_acc_a_mg_u(icoarselevel,cpu_amr,1)+icell_c_mg
            corr_aux=corr_aux+coeff*active_mg_u(active_ind)
         end do
         ! Correct potential
         phi(icell_amr(ibin,i))=phi(icell_amr(ibin,i))+corr_aux
      end do
      end do
      ! End loop over cells
   end do
   ! End loop over grids

   !$acc end data

   deallocate(nbatch)
   deallocate(igrid_f_amr,icell_amr)
   deallocate(nbors_father_cells,nbors_father_grids)
   deallocate(pos,ind_grid_father,ind_grid_ok)
   deallocate(nbors_father_ok,nbors_grids_ok)
   deallocate(ind_grid1,ind_grid2,ind_grid3)
   deallocate(nbors_grids)
   deallocate(ix,iy,iz,iix,iiy,iiz)
   
#ifdef MPROFILING
  tt2=MPI_WTIME(info)
  acc_t_interpolate_and_correct_fine = acc_t_interpolate_and_correct_fine + (tt2-tt1)
#endif   

end subroutine interpolate_and_correct_fine


! ------------------------------------------------------------------------
! Flag setting
! ------------------------------------------------------------------------

subroutine set_scan_flag_fine(ilevel)
   use amr_commons
   use poisson_commons
   use poisson_commons_acc
   use ind_accessing
   implicit none
   integer, intent(in) :: ilevel

   integer :: ind, ngrid, scan_flag
   integer :: igrid_mg, inbor, idim, igshift
   integer :: igrid_amr, igrid_nbor_amr
   integer :: iskip_amr, icell_amr, icell_nbor_amr
#ifdef MPROFILING
  include 'mpif.h'
  integer::info
  real(kind=8)::tt1,tt2
  tt1=MPI_WTIME(info)
#endif

   ngrid = active(ilevel)%ngrid

   !$acc data present(active,son,nbor,f,flag2,iii,jjj)

   ! Loop over cells and set fine SCAN flag
   !$acc parallel loop collapse(2) independent
   do ind=1,twotondim
      do igrid_mg=1,ngrid
         iskip_amr = ncoarse+(ind-1)*ngridmax
         igrid_amr = active(ilevel)%igrid(igrid_mg)
         icell_amr = iskip_amr + igrid_amr

         if(f(icell_amr,3)==1.0) then
            scan_flag=0       ! Init flag to 'no scan needed'
            scan_flag_loop: do inbor=1,2
               do idim=1,ndim
                  igshift = iii(idim,inbor,ind)
                  if(igshift==0) then
                     igrid_nbor_amr = igrid_amr
                  else
                     igrid_nbor_amr = son(nbor(igrid_amr,igshift))
                  end if

                  if(igrid_nbor_amr==0) then
                     scan_flag=1
                     exit scan_flag_loop
                  else
                     icell_nbor_amr = igrid_nbor_amr + &
                           ncoarse+(jjj(idim,inbor,ind)-1)*ngridmax
                     if(f(icell_nbor_amr,3)<=0.0) then
                        scan_flag=1
                        exit scan_flag_loop
                     end if
                  end if
               end do
            end do scan_flag_loop
         else
            scan_flag=1
         end if
         ! Update flag2 with scan flag,
         ! BEWARE as lookup_mg backups are stored in flag2
         ! Safety init:
         if(flag2(icell_amr)>ngridmax .or. flag2(icell_amr)<0) flag2(icell_amr)=0 
         ! Do NOT overwrite flag2 !
         flag2(icell_amr)=flag2(icell_amr)+ngridmax*scan_flag
      end do
   end do
   
#if defined(_OPENACC)
   call update_globalvar_int_to_host(flag2,ilevel)
#endif
   !$acc end data
   
#ifdef MPROFILING
  tt2=MPI_WTIME(info)
  acc_t_set_scan_flag_fine = acc_t_set_scan_flag_fine + (tt2-tt1)
#endif

end subroutine set_scan_flag_fine
