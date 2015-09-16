! ------------------------------------------------------------------------
! Multigrid Poisson solver for refined AMR levels
! ------------------------------------------------------------------------
! This file contains all generic fine multigrid routines, such as
!   * multigrid iterations @ fine and coarse MG levels
!   * communicator building
!   * MPI routines
!   * helper functions
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
! Main multigrid routine, called by amr_step
! ------------------------------------------------------------------------

subroutine multigrid_fine(ilevel,icount)
   use amr_commons
   use poisson_commons
   use poisson_parameters
   use poisson_commons_acc
   use ind_accessing
   implicit none
#ifndef WITHOUTMPI
#ifndef INCLUDEOK
#define INCLUDEOK
  include 'mpif.h'
#endif
#endif

   integer, intent(in) :: ilevel,icount

   integer, parameter  :: MAXITER  = 10
   real(dp), parameter :: SAFE_FACTOR = 0.5

   integer  :: ifine, i, iter, info, icpu, j, step
   real(kind=8) :: res_norm2, i_res_norm2, i_res_norm2_tot, res_norm2_tot
   real(kind=8) :: err, last_err,maxval_aux
   integer,dimension(1:ilevel):: tmp_max

   logical :: allmasked, allmasked_tot
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

   if(gravity_type>0)return
   if(numbtot(1,ilevel)==0)return

   if(verbose) print '(A,I2)','Entering fine multigrid at level ',ilevel

   ! ---------------------------------------------------------------------
   ! Prepare first guess, mask and BCs at finest level
   ! ---------------------------------------------------------------------

   if(ilevel>levelmin)then
      call make_initial_phi(ilevel,icount)  ! Interpolate phi down
   else
      call make_multipole_phi(ilevel)       ! Fill with simple initial guess
   endif
   call make_virtual_fine_dp_acc(phi(1),ilevel) ! Update boundaries
   call make_boundary_phi(ilevel)           ! Update physical boundaries

   call make_fine_mask  (ilevel)            ! Fill the fine mask
   call make_virtual_fine_dp_acc(f(:,3),ilevel) ! Communicate mask
   call make_boundary_mask(ilevel)          ! Set mask to -1 in phys bounds

   call make_fine_bc_rhs(ilevel,icount)     ! Fill BC-modified RHS

   ! ---------------------------------------------------------------------
   ! Build communicators up
   ! ---------------------------------------------------------------------

   ! @ finer level
   call build_parent_comms_mg(active(ilevel),ilevel)
   ! @ coarser levels
   do ifine=(ilevel-1),2,-1
      call build_parent_comms_mg(active_mg(myid,ifine),ifine)
   end do

   ! ---------------------------------------------------------------------
   ! Replacement of derived data types 
   ! ---------------------------------------------------------------------

   allocate(ngrid_a_mg(1:ncpu,1:ilevel),ngrid_e_mg(1:ncpu,1:ilevel))

   do ifine=1,ilevel
      do icpu=1,ncpu
         ngrid_a_mg(icpu,ifine) = active_mg(icpu,ifine)%ngrid
         ngrid_e_mg(icpu,ifine) = emission_mg(icpu,ifine)%ngrid
      end do
   tmp_max(ifine) = maxval(ngrid_a_mg(:,ifine))
   end do

   allocate(active_mg_u      (1:sum(ngrid_a_mg(:,:)*twotondim*4)))
   allocate(active_mg_f      (1:sum(ngrid_a_mg(:,:)*twotondim  )))
   allocate(emission_mg_u    (1:sum(ngrid_e_mg(:,:)*twotondim*4)))
   allocate(active_mg_igrid  (1:sum(ngrid_a_mg(:,:)            )))
   allocate(emission_mg_igrid(1:sum(ngrid_e_mg(:,:)            )))

   ! copy active_mg%u to active_mg_u
   do lvl_ind=1,ilevel
   do component=1,4
   do cpu_ind=1,ncpu
   do j=1,twotondim
   do i=1,ngrid_a_mg(cpu_ind,lvl_ind)
      step=(j-1)*ngrid_a_mg(cpu_ind,lvl_ind)
      active_ind = ind_acc_a_mg_u(lvl_ind,cpu_ind,component)+step+i
      active_mg_u(active_ind) = active_mg(cpu_ind,lvl_ind)%u(i+step,component)
   end do
   end do
   end do
   end do
   end do

   ! copy active_mg%f to active_mg_f
   do lvl_ind=1,ilevel
   do cpu_ind=1,ncpu
   do j=1,twotondim
   do i=1,ngrid_a_mg(cpu_ind,lvl_ind)
      step=(j-1)*ngrid_a_mg(cpu_ind,lvl_ind)
      active_ind = ind_acc_a_mg_f(lvl_ind,cpu_ind)+step+i
      active_mg_f(active_ind) = active_mg(cpu_ind,lvl_ind)%f(i+step,1)
   end do
   end do
   end do
   end do

   ! copy emission_mg%u to emission_mg_u
   do lvl_ind=1,ilevel
   do component=1,4
   do cpu_ind=1,ncpu
   do j=1,twotondim
   do i=1,ngrid_e_mg(cpu_ind,lvl_ind)
      step=(j-1)*ngrid_e_mg(cpu_ind,lvl_ind)
      emission_ind = ind_acc_e_mg_u(lvl_ind,cpu_ind,component)+step+i
      emission_mg_u(emission_ind) = emission_mg(cpu_ind,lvl_ind)%u(i+step,component)
   end do
   end do
   end do
   end do
   end do

   ! copy active_mg%igrid to active_mg_igrid
   do lvl_ind=1,ilevel
   do cpu_ind=1,ncpu
   do i=1,ngrid_a_mg(cpu_ind,lvl_ind)
      active_ind = ind_acc_a_mg_igrid(lvl_ind,cpu_ind)+i
      active_mg_igrid(active_ind) = active_mg(cpu_ind,lvl_ind)%igrid(i)
   end do
   end do
   end do

   ! copy emission_mg%igrid to emission_mg_igrid
   do lvl_ind=1,ilevel
   do cpu_ind=1,ncpu
   do i=1,ngrid_e_mg(cpu_ind,lvl_ind)
      emission_ind = ind_acc_e_mg_igrid(lvl_ind,cpu_ind)+i
      emission_mg_igrid(emission_ind) = emission_mg(cpu_ind,lvl_ind)%igrid(i)
   end do
   end do
   end do

   ! ---------------------------------------------------------------------
   ! Prepare Gravity update (reduce data transfer GPU <--> CPU)
   ! ---------------------------------------------------------------------

   ired  (1,1:4)=(/1,0,0,0/)
   iblack(1,1:4)=(/2,0,0,0/)
   ired  (2,1:4)=(/1,4,0,0/)
   iblack(2,1:4)=(/2,3,0,0/)
   ired  (3,1:4)=(/1,4,6,7/)
   iblack(3,1:4)=(/2,3,5,8/)
   iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
   iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

   a = 1.0D0/4.0D0**ndim
   b = 3.0D0*a
   c = 9.0D0*a
   d = 27.D0*a
   bbb(:)  =(/a ,b ,b ,c ,b ,c ,c ,d/)
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

   ! ---------------------------------------------------------------------
   ! Prepare Device for the compuatations
   ! ---------------------------------------------------------------------

   !$acc data copyin(ngrid_a_mg,ngrid_e_mg,active_mg_u,active_mg_f,emission_mg_u) &
   !$acc copyin(active_mg_igrid,emission_mg_igrid,tmp_max) &
   !$acc copyin(lookup_mg,ired,iblack,iii,jjj,bbb,ccc,iii_aux,jjj_aux,kkk_aux,lll,mmm)
#if defined(_OPENACC)
   call update_globalvar_int_to_device(flag2,ilevel)
#endif

   ! ---------------------------------------------------------------------
   ! Restrict mask up, then set scan flag
   ! ---------------------------------------------------------------------
   ! @ finer level

   if(ilevel>1) then
      ! Restrict and communicate mask
      call restrict_mask_fine_reverse(ilevel)
      call make_reverse_mg_dp_acc(4,ilevel-1)
      call make_virtual_mg_dp_acc(4,ilevel-1)

      ! Convert volume fraction to mask value
      !$acc parallel loop collapse(3) gang worker vector independent
      do icpu=1,ncpu
         do j=1,twotondim
         do i=1,tmp_max(ilevel-1) !ngrid_a_mg(icpu,ilevel-1)
         if(ngrid_a_mg(icpu,ilevel-1)==0)   cycle
         if(i.GT.ngrid_a_mg(icpu,ilevel-1)) cycle
         step=(j-1)*ngrid_a_mg(icpu,ilevel-1)
         active_ind = ind_acc_a_mg_u(ilevel-1,icpu,4)+step+i
         active_mg_u(active_ind)=2d0*active_mg_u(active_ind)-1d0
         end do
         end do
      end do

      ! Check active mask state
      maxval_aux=0
      if(ngrid_a_mg(myid,ilevel-1)>0) then
         !$acc parallel loop collapse(2) gang vector independent reduction(max:maxval_aux) &
         !$acc copy(maxval_aux)
         do j=1,twotondim
         do i=1,ngrid_a_mg(myid,ilevel-1)
         step=(j-1)*ngrid_a_mg(myid,ilevel-1)
         active_ind = ind_acc_a_mg_u(ilevel-1,myid,4)+step+i
         if (active_mg_u(active_ind).GT.maxval_aux) maxval_aux=active_mg_u(active_ind)
         end do
         end do
         allmasked=maxval_aux<=0d0
      else
         allmasked=.true.
      end if

      ! Allreduce on mask state
#ifndef WITHOUTMPI
      call MPI_ALLREDUCE(allmasked, allmasked_tot, 1, MPI_LOGICAL, &
           & MPI_LAND, MPI_COMM_WORLD, info)
      allmasked=allmasked_tot
#endif
   else
      allmasked=.true.
   endif

   ! @ coarser levels
   ! Restrict mask and compute levelmin_mg in the process
   if (.not. allmasked) then
      levelmin_mg=1
      do ifine=(ilevel-1),2,-1

         ! Restrict and communicate mask
         call restrict_mask_coarse_reverse(ifine)
         call make_reverse_mg_dp_acc(4,ifine-1)
         call make_virtual_mg_dp_acc(4,ifine-1)

         ! Convert volume fraction to mask value
         !$acc parallel loop collapse(3) gang worker vector independent
         do icpu=1,ncpu
            do j=1,twotondim
            do i=1,tmp_max(ifine-1)
            if(ngrid_a_mg(icpu,ifine-1)==0)   cycle
            if(i.GT.ngrid_a_mg(icpu,ifine-1)) cycle
            step=(j-1)*ngrid_a_mg(icpu,ifine-1)
            active_ind = ind_acc_a_mg_u(ifine-1,icpu,4)+step+i
            active_mg_u(active_ind)=2d0*active_mg_u(active_ind)-1d0
            end do
            end do
         end do

         ! Check active mask state
         if(ngrid_a_mg(myid,ifine-1)>0) then
         maxval_aux=0
         !$acc parallel loop collapse(2) gang vector independent reduction(max:maxval_aux) &
         !$acc copy(maxval_aux)
         do j=1,twotondim
         do i=1,ngrid_a_mg(myid,ifine-1)
         step=(j-1)*ngrid_a_mg(myid,ifine-1)
         active_ind = ind_acc_a_mg_u(ifine-1,myid,4)+step+i
         if (active_mg_u(active_ind).GT.maxval_aux) maxval_aux=active_mg_u(active_ind)
         end do
         end do
            allmasked=maxval_aux<=0d0
         else
            allmasked=.true.
         end if

         ! Allreduce on mask state
#ifndef WITHOUTMPI
         call MPI_ALLREDUCE(allmasked,allmasked_tot,1,MPI_LOGICAL, &
                 & MPI_LAND,MPI_COMM_WORLD,info)
         allmasked=allmasked_tot
#endif

         if(allmasked) then ! Coarser level is fully masked: stop here
            levelmin_mg=ifine
            exit
         end if
      end do
   else
      levelmin_mg=ilevel
   end if
   if(nboundary>0)levelmin_mg=max(levelmin_mg,2)

   ! Update flag with scan flag
   call set_scan_flag_fine(ilevel)
   do ifine=levelmin_mg,ilevel-1
      call set_scan_flag_coarse(ifine)
   end do

   ! ---------------------------------------------------------------------
   ! Initiate solve at fine level
   ! ---------------------------------------------------------------------
   iter = 0
   err = 1.0d0
   main_iteration_loop: do
      iter=iter+1

      ! Pre-smoothing
      do i=1,ngs_fine
         call gauss_seidel_mg_fine(ilevel,.true. )  ! Red step
         call make_virtual_fine_dp_acc(phi(1),ilevel)   ! Communicate phi
         call gauss_seidel_mg_fine(ilevel,.false.)  ! Black step
         call make_virtual_fine_dp_acc(phi(1),ilevel)   ! Communicate phi
      end do

      ! Compute residual and restrict into upper level RHS
      call cmp_residual_mg_fine(ilevel)
      call make_virtual_fine_dp_acc(f(1,1),ilevel) ! communicate residual
     
      if(iter==1) then
         call cmp_residual_norm2_fine(ilevel,i_res_norm2)
#ifndef WITHOUTMPI
         call MPI_ALLREDUCE(i_res_norm2,i_res_norm2_tot,1, &
                 & MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
         i_res_norm2=i_res_norm2_tot
#endif
      end if

      ! First clear the rhs in coarser reception comms
      !$acc parallel loop collapse(3) gang worker vector independent
      do icpu=1,ncpu
         do j=1,twotondim
         do i=1,tmp_max(ilevel-1)
         if(ngrid_a_mg(icpu,ilevel-1)==0) cycle
         if(i.GT.ngrid_a_mg(icpu,ilevel-1)) cycle
         step=(j-1)*ngrid_a_mg(icpu,ilevel-1)
         active_ind = ind_acc_a_mg_u(ilevel-1,icpu,2)+step+i
         active_mg_u(active_ind)=0.0d0
         end do
         end do
      end do

      ! Restrict and do reverse-comm
      call restrict_residual_fine_reverse(ilevel)
      call make_reverse_mg_dp_acc(2,ilevel-1) ! communicate rhs

      if(ilevel>1) then
         ! Reset correction at upper level before solve
         !$acc parallel loop collapse(3) gang worker vector independent
         do icpu=1,ncpu
         do j=1,twotondim
         do i=1,tmp_max(ilevel-1)
         if(ngrid_a_mg(icpu,ilevel-1)==0) cycle
         if(i.GT.ngrid_a_mg(icpu,ilevel-1)) cycle
         step=(j-1)*ngrid_a_mg(icpu,ilevel-1)
         active_ind = ind_acc_a_mg_u(ilevel-1,icpu,1)+step+i
         active_mg_u(active_ind)=0.0d0
         end do
         end do
         end do

         ! Multigrid-solve the upper level
         call recursive_multigrid_coarse(ilevel-1, safe_mode(ilevel),tmp_max,ilevel)

         ! Interpolate coarse solution and correct fine solution
         call interpolate_and_correct_fine(ilevel)
         call make_virtual_fine_dp_acc(phi(1),ilevel)   ! Communicate phi       
      end if

      ! Post-smoothing
      do i=1,ngs_fine
         call gauss_seidel_mg_fine(ilevel,.true. )  ! Red step
         call make_virtual_fine_dp_acc(phi(1),ilevel)   ! Communicate phi
         call gauss_seidel_mg_fine(ilevel,.false.)  ! Black step
         call make_virtual_fine_dp_acc(phi(1),ilevel)   ! Communicate phi
      end do

      ! Update fine residual
      call cmp_residual_mg_fine(ilevel)
      call make_virtual_fine_dp_acc(f(1,1),ilevel) ! communicate residual

      call cmp_residual_norm2_fine(ilevel,res_norm2)
#ifndef WITHOUTMPI
      call MPI_ALLREDUCE(res_norm2,res_norm2_tot,1, &
              & MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
      res_norm2=res_norm2_tot
#endif

      last_err = err
      err = sqrt(res_norm2/i_res_norm2)

      ! Verbosity
      if(verbose) print '(A,I5,A,1pE10.3)','   ==> Step=', &
         & iter,' Error=',err

      ! Converged?
      if(err<epsilon .or. iter>=MAXITER) exit

      ! Not converged, check error and possibly enable safe mode for the level
      if(err > last_err*SAFE_FACTOR .and. (.not. safe_mode(ilevel))) then
         if(verbose)print *,'CAUTION: Switching to safe MG mode for level ',ilevel
         safe_mode(ilevel) = .true.
      end if

   end do main_iteration_loop

   !$acc end data

   if(myid==1) print '(A,I5,A,I5,A,1pE10.3)','   ==> Level=',ilevel, ' Step=', &
            iter,' Error=',err
   if(myid==1 .and. iter==MAXITER) print *,'WARN: Fine multigrid &
      &Poisson failed to converge...'

   ! ---------------------------------------------------------------------
   ! Cleanup MG levels after solve complete
   ! ---------------------------------------------------------------------
   do ifine=1,ilevel-1
      call cleanup_mg_level(ifine)
   end do

   ! clean up
   deallocate(ngrid_a_mg,ngrid_e_mg)
   deallocate(active_mg_u,active_mg_f,emission_mg_u)
   deallocate(active_mg_igrid,emission_mg_igrid)

#ifdef MPROFILING
  tt2=MPI_WTIME(info)
  acc_t_multigrid_fine = acc_t_multigrid_fine + (tt2-tt1)
#endif

end subroutine multigrid_fine

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Recursive multigrid routine for coarse MG levels
! ------------------------------------------------------------------------
recursive subroutine recursive_multigrid_coarse(ifinelevel, safe, tmp_max, base_level)
   use amr_commons
   use poisson_commons
   use poisson_commons_acc
   use ind_accessing
   implicit none
#ifndef WITHOUTMPI
   include "mpif.h"
#endif
   integer             :: base_level
   integer, intent(in) :: ifinelevel
   logical, intent(in) :: safe
   integer,dimension(1:base_level) :: tmp_max
   integer :: i, icpu, info, icycle, ncycle
   integer :: j, step

   if(ifinelevel<=levelmin_mg) then
      ! Solve 'directly' :
      do i=1,2*ngs_coarse
         call gauss_seidel_mg_coarse(ifinelevel,safe,.true. )  ! Red step
         call make_virtual_mg_dp_acc(1,ifinelevel)  ! Communicate solution
         call gauss_seidel_mg_coarse(ifinelevel,safe,.false.)  ! Black step
         call make_virtual_mg_dp_acc(1,ifinelevel)  ! Communicate solution
      end do
      return
   end if

   if(safe) then
      ncycle=ncycles_coarse_safe
   else
      ncycle=1
   endif

   !$acc data present(active_mg_u,ngrid_a_mg)

   do icycle=1,ncycle

      ! Pre-smoothing
      do i=1,ngs_coarse
         call gauss_seidel_mg_coarse(ifinelevel,safe,.true. )  ! Red step
         call make_virtual_mg_dp_acc(1,ifinelevel)  ! Communicate solution
         call gauss_seidel_mg_coarse(ifinelevel,safe,.false.)  ! Black step
         call make_virtual_mg_dp_acc(1,ifinelevel)  ! Communicate solution
      end do

      ! Compute residual and restrict into upper level RHS
      call cmp_residual_mg_coarse(ifinelevel)
      call make_virtual_mg_dp_acc(3,ifinelevel)  ! Communicate residual

      ! First clear the rhs in coarser reception comms
      !$acc parallel loop collapse(3) gang worker vector independent
      do icpu=1,ncpu
         do j=1,twotondim
         do i=1,tmp_max(ifinelevel-1)
         if(ngrid_a_mg(icpu,ifinelevel-1)==0) cycle
         if(i.GT.ngrid_a_mg(icpu,ifinelevel-1)) cycle
         step=(j-1)*ngrid_a_mg(icpu,ifinelevel-1)
         active_ind = ind_acc_a_mg_u(ifinelevel-1,icpu,2)+step+i
         active_mg_u(active_ind)=0.0d0
         end do
         end do
      end do
         
      ! Restrict and do reverse-comm
      call restrict_residual_coarse_reverse(ifinelevel)
      call make_reverse_mg_dp_acc(2,ifinelevel-1) ! communicate rhs

      ! Reset correction from upper level before solve
      !$acc parallel loop collapse(3) gang worker vector independent
      do icpu=1,ncpu
         do j=1,twotondim
         do i=1,tmp_max(ifinelevel-1)
         if(ngrid_a_mg(icpu,ifinelevel-1)==0) cycle
         if(i.GT.ngrid_a_mg(icpu,ifinelevel-1)) cycle
         step=(j-1)*ngrid_a_mg(icpu,ifinelevel-1)
         active_ind = ind_acc_a_mg_u(ifinelevel-1,icpu,1)+step+i
         active_mg_u(active_ind)=0.0d0
         end do
         end do
      end do

      ! Multigrid-solve the upper level
      call recursive_multigrid_coarse(ifinelevel-1, safe, tmp_max, base_level)

      ! Interpolate coarse solution and correct back into fine solution
      call interpolate_and_correct_coarse(ifinelevel)
      call make_virtual_mg_dp_acc(1,ifinelevel)  ! Communicate solution

      ! Post-smoothing
      do i=1,ngs_coarse
         call gauss_seidel_mg_coarse(ifinelevel,safe,.true. )  ! Red step
         call make_virtual_mg_dp_acc(1,ifinelevel)  ! Communicate solution
         call gauss_seidel_mg_coarse(ifinelevel,safe,.false.)  ! Black step
         call make_virtual_mg_dp_acc(1,ifinelevel)  ! Communicate solution
      end do

   end do

   !$acc end data

end subroutine recursive_multigrid_coarse


! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Multigrid communicator building
! ------------------------------------------------------------------------
subroutine build_parent_comms_mg(active_f_comm, ifinelevel)
   use amr_commons
   use poisson_commons
   implicit none

#ifndef WITHOUTMPI
   include "mpif.h"
   integer, dimension (MPI_STATUS_SIZE, ncpu) :: statuses
#endif

   integer, intent(in) :: ifinelevel
   type(communicator), intent(in) :: active_f_comm

   integer :: icoarselevel
   integer :: ngrids, cur_grid, cur_cpu, cur_cell, newgrids
   integer :: i, nbatch, ind, icpu, istart, info

   integer :: nact_tot, nreq_tot, nreq_tot2
   integer, dimension(1:ncpu) :: nreq, nreq2

   integer, dimension(1:nvector), save :: ind_cell_father
   integer, dimension(1:nvector,1:twotondim),   save :: nbors_father_grids
   integer, dimension(1:nvector,1:threetondim), save :: nbors_father_cells

   type(communicator), dimension(1:ncpu) :: comm_send, comm_receive
   type(communicator), dimension(1:ncpu) :: comm_send2, comm_receive2

   integer, dimension(1:ncpu) :: indx, recvbuf, recvbuf2
   integer, dimension(1:ncpu) :: reqsend, reqrecv
   integer :: countrecv, countsend
   integer :: tag = 777


   icoarselevel=ifinelevel-1

   nact_tot=0
   nreq_tot=0; nreq=0
   indx=0; recvbuf=0

   ! ---------------------------------------------------------------------
   ! STAGE 1 : Coarse grid MG activation for local grids (1st pass)
   ! ---------------------------------------------------------------------

   ! Loop over the AMR active communicator first
   ngrids = active_f_comm%ngrid
   do istart=1,ngrids,nvector
      nbatch=min(nvector,ngrids-istart+1)
      ! Gather grid indices and retrieve parent cells
      do i=1,nbatch
         ind_cell_father(i)=father( active_f_comm%igrid(istart+i-1) )
      end do

      ! Compute neighbouring father cells and grids
      call get3cubefather(ind_cell_father,nbors_father_cells, &
         & nbors_father_grids,nbatch,ifinelevel)

      ! Now process the twotondim father grids
      do ind=1,twotondim
         do i=1,nbatch
            cur_grid = nbors_father_grids(i,ind)
            if(lookup_mg(cur_grid)>0) cycle ! Grid already active

            cur_cpu=cpu_map(father(cur_grid))
            if(cur_cpu==0) cycle

            if(cur_cpu==myid) then
               ! Stack grid for local activation
               ! We own the grid: fill lookup_mg with its final value
               nact_tot=nact_tot+1
               flag2(nact_tot)=cur_grid
               lookup_mg(cur_grid)=nact_tot
            else
               ! Stack grid into 2nd part of flag2 for remote activation
               ! Here, lookup_mg(cur_grid) should be MINUS the AMR index
               ! of the grid in the corresponding CPU AMR structure.
               nreq_tot=nreq_tot+1
               nreq(cur_cpu)=nreq(cur_cpu)+1
               flag2(ngridmax+nreq_tot)=cur_grid ! "home" index in 2nd part
               lookup_mg(cur_grid)=abs(lookup_mg(cur_grid)) ! Flag visited (>0)
            end if
         end do
      end do
   end do


   ! ---------------------------------------------------------------------
   ! STAGE 2 : Coarse grid MG activation request
   ! ---------------------------------------------------------------------

#ifndef WITHOUTMPI
   ! Share number of requests and replies
   call MPI_ALLTOALL(nreq, 1, MPI_INTEGER, recvbuf, 1, MPI_INTEGER, &
      & MPI_COMM_WORLD, info)

   ! Allocate inbound comms
   do icpu=1,ncpu
      comm_receive(icpu)%ngrid=recvbuf(icpu)
      if(recvbuf(icpu)>0) allocate(comm_receive(icpu)%igrid(1:recvbuf(icpu)))
   end do

   ! Receive to-be-activated grids
   countrecv=0; reqrecv=0
   do icpu=1,ncpu
      if(icpu==myid) cycle
      ngrids=comm_receive(icpu)%ngrid
      if(ngrids>0) then
         countrecv=countrecv+1
         call MPI_IRECV(comm_receive(icpu)%igrid, ngrids, MPI_INTEGER, &
            & icpu-1, tag, MPI_COMM_WORLD, reqrecv(countrecv), info)
      end if
   end do

   ! Allocate and then fill outbound (activation request) communicators
   do icpu=1,ncpu
      comm_send(icpu)%ngrid=nreq(icpu)
      if(nreq(icpu)>0) allocate(comm_send(icpu)%igrid(1:nreq(icpu)))
   end do
   nreq=0
   do i=1,nreq_tot
      cur_grid=flag2(ngridmax+i) ! Local AMR index
      cur_cpu =cpu_map(father(cur_grid))
      nreq(cur_cpu)=nreq(cur_cpu)+1
      comm_send(cur_cpu)%igrid(nreq(cur_cpu))=lookup_mg(cur_grid) ! Remote
   end do

   ! Send to-be-activated grids
   countsend=0; reqsend=0
   do icpu=1,ncpu
      if(icpu==myid) cycle
      ngrids=comm_send(icpu)%ngrid
      if(ngrids>0) then
         countsend=countsend+1
         call MPI_ISEND(comm_send(icpu)%igrid, ngrids, MPI_INTEGER, &
            & icpu-1, tag, MPI_COMM_WORLD, reqsend(countsend), info)
      end if
   end do

   ! Wait for completion of receives
   call MPI_WAITALL(countrecv, reqrecv, statuses, info)

   ! Activate requested grids
   do icpu=1,ncpu
      if(icpu==myid) cycle
      ngrids=comm_receive(icpu)%ngrid
      if(ngrids>0) then
         do i=1,ngrids
            cur_grid=comm_receive(icpu)%igrid(i) ! Local AMR index
            if(lookup_mg(cur_grid)>0) cycle      ! Already active: cycle
            ! Activate grid
            nact_tot=nact_tot+1
            flag2(nact_tot)=cur_grid
            lookup_mg(cur_grid)=nact_tot
         end do
      end if
   end do

   ! Wait for completion of sends
   call MPI_WAITALL(countsend, reqsend, statuses, info)
#endif

   ! ---------------------------------------------------------------------
   ! STAGE 3 : Coarse grid MG active comm gathering
   ! ---------------------------------------------------------------------

   active_mg(myid,icoarselevel)%ngrid=nact_tot
   if(nact_tot>0) then
      allocate( active_mg(myid,icoarselevel)%igrid(1:nact_tot) )
      allocate( active_mg(myid,icoarselevel)%u(1:nact_tot*twotondim,1:4) )
      allocate( active_mg(myid,icoarselevel)%f(1:nact_tot*twotondim,1:1) )
      active_mg(myid,icoarselevel)%igrid=0
      active_mg(myid,icoarselevel)%u=0.0d0
      active_mg(myid,icoarselevel)%f=0
   end if
   do i=1,nact_tot
      active_mg(myid,icoarselevel)%igrid(i)=flag2(i)
   end do

   ! ---------------------------------------------------------------------
   ! STAGE 4 : Screen active grid neighbors for new reception grids
   ! ---------------------------------------------------------------------
   ngrids = active_mg(myid,icoarselevel)%ngrid
   nreq2 = 0
   nreq_tot2 = 0
   do istart=1,ngrids,nvector
      nbatch=min(nvector,ngrids-istart+1)
      ! Gather grid indices and retrieve parent cells
      do i=1,nbatch
         ind_cell_father(i)=father( active_mg(myid,icoarselevel)%igrid(istart+i-1) )
      end do

      ! Compute neighbouring father cells
      call get3cubefather(ind_cell_father,nbors_father_cells,nbors_father_grids,nbatch,icoarselevel)

      ! Now process the father grids
      do ind=1,threetondim
         do i=1,nbatch
            cur_cell = nbors_father_cells(i,ind)
            cur_cpu  = cpu_map(cur_cell)
            if(cur_cpu==0) cycle
            cur_grid = son(cur_cell)
            if(cur_cpu/=myid) then
               ! Neighbor cell is not managed by current CPU
               if (cur_grid==0) cycle              ! No grid there
               if (lookup_mg(cur_grid)>0) cycle    ! Already selected
               ! Add grid to request
               nreq_tot2=nreq_tot2+1
               nreq2(cur_cpu)=nreq2(cur_cpu)+1
               flag2(ngridmax+nreq_tot+nreq_tot2)=cur_grid
               lookup_mg(cur_grid)=abs(lookup_mg(cur_grid))  ! Mark visited
            end if
         end do
      end do
   end do

   ! ---------------------------------------------------------------------
   ! STAGE 5 : Share new reception grid requests, build emission comms
   ! ---------------------------------------------------------------------

#ifndef WITHOUTMPI
   ! Share number of requests and replies
   recvbuf2=0
   call MPI_ALLTOALL(nreq2, 1, MPI_INTEGER, recvbuf2, 1, MPI_INTEGER, &
      & MPI_COMM_WORLD, info)

   ! Allocate inbound comms
   do icpu=1,ncpu
      comm_receive2(icpu)%ngrid=recvbuf2(icpu)
      if(recvbuf2(icpu)>0) allocate(comm_receive2(icpu)%igrid(1:recvbuf2(icpu)))
   end do

   ! Receive potential reception grids
   countrecv=0; reqrecv=0
   do icpu=1,ncpu
      if(icpu==myid) cycle
      ngrids=comm_receive2(icpu)%ngrid
      if(ngrids>0) then
         countrecv=countrecv+1
         call MPI_IRECV(comm_receive2(icpu)%igrid, ngrids, MPI_INTEGER, &
            & icpu-1, tag, MPI_COMM_WORLD, reqrecv(countrecv), info)
      end if
   end do

   ! Allocate and then fill outbound (reception request) communicators
   do icpu=1,ncpu
      comm_send2(icpu)%ngrid=nreq2(icpu)
      if(nreq2(icpu)>0) allocate(comm_send2(icpu)%igrid(1:nreq2(icpu)))
   end do
   nreq2=0
   do i=1,nreq_tot2
      cur_grid=flag2(ngridmax+nreq_tot+i) ! Local AMR index
      cur_cpu =cpu_map(father(cur_grid))
      nreq2(cur_cpu)=nreq2(cur_cpu)+1
      comm_send2(cur_cpu)%igrid(nreq2(cur_cpu))=lookup_mg(cur_grid) ! Remote AMR index
      ! Restore negative lookup_mg
      lookup_mg(cur_grid)=-abs(lookup_mg(cur_grid))
   end do

   ! Send reception request grids
   countsend=0; reqsend=0
   do icpu=1,ncpu
      if(icpu==myid) cycle
      ngrids=comm_send2(icpu)%ngrid
      if(ngrids>0) then
         countsend=countsend+1
         call MPI_ISEND(comm_send2(icpu)%igrid,ngrids,MPI_INTEGER,icpu-1, &
            & tag, MPI_COMM_WORLD, reqsend(countsend), info)
      end if
   end do

   ! Wait for completion of receives
   call MPI_WAITALL(countrecv, reqrecv, statuses, info)

   ! Compute local MG indices of inbound grids, alloc & fill emission comms
   do icpu=1,ncpu
      if(icpu==myid) cycle
      newgrids=0
      do i=1,recvbuf2(icpu)
         ! MAP AMR -> MG INDICES IN PLACE
         comm_receive2(icpu)%igrid(i)=lookup_mg(comm_receive2(icpu)%igrid(i))
         if(comm_receive2(icpu)%igrid(i)>0) newgrids=newgrids+1
      end do
      ! Allocate emission communicators
      ngrids=recvbuf(icpu)+newgrids
      emission_mg(icpu,icoarselevel)%ngrid=ngrids
      if(ngrids>0) then
         allocate(emission_mg(icpu,icoarselevel)%igrid(1:ngrids))
         allocate(emission_mg(icpu,icoarselevel)%u(1:ngrids*twotondim,1:4) )
         allocate(emission_mg(icpu,icoarselevel)%f(1:ngrids*twotondim,1:1))
         emission_mg(icpu,icoarselevel)%igrid=0
         emission_mg(icpu,icoarselevel)%u=0.0d0
         emission_mg(icpu,icoarselevel)%f=0
      end if
      ! First part: activation request emission grids
      do i=1,recvbuf(icpu)
         emission_mg(icpu,icoarselevel)%igrid(i)=lookup_mg(comm_receive(icpu)%igrid(i))
      end do
      ! Second part: new emission grids
      cur_grid=recvbuf(icpu)
      do i=1,recvbuf2(icpu)
         if(comm_receive2(icpu)%igrid(i)>0) then
            cur_grid=cur_grid+1
            emission_mg(icpu,icoarselevel)%igrid(cur_grid)=comm_receive2(icpu)%igrid(i)
         end if
      end do
   end do

   ! Wait for completion of sends
   call MPI_WAITALL(countsend, reqsend, statuses, info)


   ! ---------------------------------------------------------------------
   ! STAGE 6 : Reply with local MG grid status and build reception comms
   ! ---------------------------------------------------------------------
   ! Receive MG mappings from other CPUs back into comm_send2
   countrecv=0; reqrecv=0
   do icpu=1,ncpu
      if(icpu==myid) cycle
      ngrids=nreq2(icpu)
      if(ngrids>0) then
         countrecv=countrecv+1
         call MPI_IRECV(comm_send2(icpu)%igrid,ngrids,MPI_INTEGER,icpu-1, &
            & tag, MPI_COMM_WORLD, reqrecv(countrecv), info)
      end if
   end do

   ! Send local MG mappings to other CPUs from comm_receive
   countsend=0; reqsend=0
   do icpu=1,ncpu
      if(icpu==myid) cycle
      ngrids=recvbuf2(icpu)
      if(ngrids>0) then
         countsend=countsend+1
         call MPI_ISEND(comm_receive2(icpu)%igrid,ngrids,MPI_INTEGER, &
            & icpu-1, tag, MPI_COMM_WORLD, reqsend(countsend), info)
      end if
   end do

   ! Wait for full completion of receives
   call MPI_WAITALL(countrecv, reqrecv, statuses, info)

   ! Count remotely active MG grids, and allocate and fill reception comms
   do icpu=1,ncpu
      if(icpu==myid) cycle
      ! Count requested grids which are MG-active remotely
      newgrids=0
      do i=1,nreq2(icpu)
         if(comm_send2(icpu)%igrid(i)>0) newgrids=newgrids+1
      end do
      ! Allocate and fill reception communicators on the fly
      ngrids=nreq(icpu)+newgrids
      active_mg(icpu,icoarselevel)%ngrid=ngrids
      if(ngrids>0) then
         allocate(active_mg(icpu,icoarselevel)%igrid(1:ngrids))
         allocate(active_mg(icpu,icoarselevel)%u(1:ngrids*twotondim,1:4))
         allocate(active_mg(icpu,icoarselevel)%f(1:ngrids*twotondim,1:1))
         active_mg(icpu,icoarselevel)%igrid=0
         active_mg(icpu,icoarselevel)%u=0.0d0
         active_mg(icpu,icoarselevel)%f=0
      end if
   end do
   
   nreq=0
   do i=1,nreq_tot
      cur_grid=flag2(ngridmax+i)
      cur_cpu =cpu_map(father(cur_grid))
      nreq(cur_cpu)=nreq(cur_cpu)+1
      ! Add to reception comm
      active_mg(cur_cpu,icoarselevel)%igrid(nreq(cur_cpu))=cur_grid
      ! Backup lookup_mg into flag2
      flag2(cur_grid)=lookup_mg(cur_grid)
      ! Update lookup_mg
      lookup_mg(cur_grid)=nreq(cur_cpu)
   end do

   nreq2=0; indx=nreq
   do i=1,nreq_tot2
      cur_grid=flag2(ngridmax+nreq_tot+i)
      cur_cpu =cpu_map(father(cur_grid))
      nreq2(cur_cpu)=nreq2(cur_cpu)+1
      if(comm_send2(cur_cpu)%igrid(nreq2(cur_cpu))>0) then
         indx(cur_cpu)=indx(cur_cpu)+1
         ! Add to reception comm
         active_mg(cur_cpu,icoarselevel)%igrid(indx(cur_cpu))=cur_grid
         ! Backup lookup_mg
         flag2(cur_grid)=-lookup_mg(cur_grid)
         ! Update lookup_mg
         lookup_mg(cur_grid)=indx(cur_cpu)
      end if
   end do

   ! Wait for full completion of sends
   call MPI_WAITALL(countsend, reqsend, statuses, info)


   ! Cleanup
   do icpu=1,ncpu
      if(comm_send (icpu)%ngrid>0) deallocate(comm_send (icpu)%igrid)
      if(comm_send2(icpu)%ngrid>0) deallocate(comm_send2(icpu)%igrid)
      if(comm_receive (icpu)%ngrid>0) deallocate(comm_receive (icpu)%igrid)
      if(comm_receive2(icpu)%ngrid>0) deallocate(comm_receive2(icpu)%igrid)
   end do
#endif

end subroutine build_parent_comms_mg


! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Multigrid level cleanup
! ------------------------------------------------------------------------
subroutine cleanup_mg_level(ilevel)
   use amr_commons
   use pm_commons
   use poisson_commons
   implicit none

   integer, intent(in) :: ilevel

   integer :: igrid, icpu, cur_grid, cur_cpu

   ! ---------------------------------------------------------------------
   ! Cleanup lookup table
   ! ---------------------------------------------------------------------
   do icpu=1,ncpu
      do igrid=1,active_mg(icpu,ilevel)%ngrid
         cur_grid=active_mg(icpu,ilevel)%igrid(igrid)
         cur_cpu=cpu_map(father(cur_grid))
         if(cur_cpu==myid) then
            lookup_mg(cur_grid)=0
         else
            lookup_mg(cur_grid)=-mod(flag2(cur_grid),ngridmax)
         end if
      end do
   end do

   ! ---------------------------------------------------------------------
   ! Deallocate communicators
   ! ---------------------------------------------------------------------
   do icpu=1,ncpu
      if(active_mg(icpu,ilevel)%ngrid>0)then
         deallocate(active_mg(icpu,ilevel)%igrid)
         deallocate(active_mg(icpu,ilevel)%u)
         deallocate(active_mg(icpu,ilevel)%f)
      endif
      active_mg(icpu,ilevel)%ngrid=0
      if(emission_mg(icpu,ilevel)%ngrid>0)then
         deallocate(emission_mg(icpu,ilevel)%igrid)
         deallocate(emission_mg(icpu,ilevel)%u)
         deallocate(emission_mg(icpu,ilevel)%f)
      endif
      emission_mg(icpu,ilevel)%ngrid=0
   end do

end subroutine cleanup_mg_level

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Initialize mask at fine level into f(:,3)
! ------------------------------------------------------------------------
subroutine make_fine_mask(ilevel)
   use ind_accessing
   use amr_commons
   use pm_commons
   use poisson_commons
   implicit none
   integer, intent(in) :: ilevel

   integer  :: ngrid
   integer  :: ind, igrid_mg, icpu, ibound
   integer  :: igrid_amr, icell_amr, iskip_amr
#ifdef MPROFILING
  include 'mpif.h'
  integer::info
  real(kind=8)::tt1,tt2
  tt1=MPI_WTIME(info)
#endif

   ngrid=active(ilevel)%ngrid

   !$acc data present(f,active,reception,boundary) 

   !$acc parallel loop collapse(2) gang vector independent
   do ind=1,twotondim
      do igrid_mg=1,ngrid
         iskip_amr = ncoarse+(ind-1)*ngridmax
         igrid_amr = active(ilevel)%igrid(igrid_mg)
         icell_amr = iskip_amr + igrid_amr
         ! Init mask to 1.0 on active cells :
         f(icell_amr,3) = 1.0d0
      end do
   end do

   !$acc parallel loop collapse(2) gang independent
   do icpu=1,ncpu
   do ind=1,twotondim
      ngrid=reception(icpu,ilevel)%ngrid
      !$acc loop vector
      do igrid_mg=1,ngrid
         iskip_amr = ncoarse+(ind-1)*ngridmax
         igrid_amr = reception(icpu,ilevel)%igrid(igrid_mg)
         icell_amr = iskip_amr + igrid_amr
         ! Init mask to 1.0 on virtual cells :
         f(icell_amr,3) = 1.0d0
      end do
   end do
   end do

   !$acc parallel loop collapse(2) gang independent
   do ibound=1,nboundary
   do ind=1,twotondim
      ngrid=boundary(ibound,ilevel)%ngrid
      !$acc loop vector
      do igrid_mg=1,ngrid
         iskip_amr=ncoarse+(ind-1)*ngridmax
         igrid_amr = boundary(ibound,ilevel)%igrid(igrid_mg)
         icell_amr = iskip_amr + igrid_amr
         ! Init mask to -1.0 on boundary cells :
         f(icell_amr,3) = -1.0d0
      end do
   end do
   end do

   !$acc end data
   
#ifdef MPROFILING
  tt2=MPI_WTIME(info)
  acc_t_make_fine_mask = acc_t_make_fine_mask + (tt2-tt1)
#endif 

end subroutine make_fine_mask

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! Preprocess the fine (AMR) level RHS to account for boundary conditions
!
!  _____#_____
! |     #     |      Cell I is INSIDE active domain (mask > 0)
! |  I  #  O  |      Cell O is OUTSIDE (mask <= 0 or nonexistent cell)
! |_____#_____|      # is the boundary
!       #
!
! phi(I) and phi(O) must BOTH be set at call time, if applicable
! phi(#) is computed from phi(I), phi(O) and the mask values
! If AMR cell O does not exist, phi(O) is computed by interpolation
!
! Sets BC-modified RHS    into f(:,2)
!
! ------------------------------------------------------------------------
subroutine make_fine_bc_rhs(ilevel,icount)
   use amr_commons
   use pm_commons
   use poisson_commons
   use mfbrv
   use mfbrv_subroutines
   use ind_accessing
   implicit none
   integer, intent(in) :: ilevel,icount
   integer, dimension(1:3,1:2,1:8) :: iii, jjj
   real(dp) :: dx, oneoverdx2, phi_b, nb_mask, nb_phi, w
   integer  :: ngrid
   integer  :: ind, igrid_mg, idim, inbor
   integer  :: igrid_amr, icell_amr, iskip_amr
   integer  :: igshift, igrid_nbor_amr, icell_nbor_amr
   integer  :: ifathercell_nbor_amr
   integer  :: nx_loc
   real(dp) :: scale, fourpi

   integer  :: num_blocks
   integer  :: block, igrid_mg_m, ngrid_initial

   ! These 6 variables cannot be used from the module mfbrv. The compiler cannot 
   ! accept it
   integer,dimension(:,:,:),allocatable::ix,iy,iz,iix,iiy,iiz 
#ifdef MPROFILING
  include 'mpif.h'
  integer::info
  real(kind=8)::tt1,tt2
  tt1=MPI_WTIME(info)
#endif  

   ! Set constants
   nx_loc = icoarse_max-icoarse_min+1
   scale  = boxlen/dble(nx_loc)
   fourpi = 4.D0*ACOS(-1.0D0)*scale
   if(cosmo) fourpi = 1.5D0*omega_m*aexp*scale
   dx  = 0.5d0**ilevel
   oneoverdx2 = 1.0d0/(dx*dx)

   iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
   iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

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

   ! For 1 GPU this part crashes because of memory bottleneck. For this reason
   ! we divide in blocks the ngrid and solve every block after the other (seq!)
   num_blocks = 8
   if(active(ilevel)%ngrid.LT.100) num_blocks = 1

   ngrid = ceiling(real(active(ilevel)%ngrid)/real(num_blocks))
   ngrid_initial = active(ilevel)%ngrid

   allocate(phi_int (1:twotondim,1:ngrid,1:nvector,1:twotondim))
   allocate(ind_cell(1:twotondim,1:ngrid,1:nvector))
   allocate(nbors_father_grids(1:twotondim,1:ngrid,1:nvector,1:twotondim))
   allocate(nbors_father_cells(1:twotondim,1:ngrid,1:nvector,1:threetondim))
   allocate(     pos(1:twotondim,1:ngrid,1:nvector), & 
   & ind_grid_father(1:twotondim,1:ngrid,1:nvector), &
   &     ind_grid_ok(1:twotondim,1:ngrid,1:nvector))
   allocate(nbors_father_ok(1:twotondim,1:ngrid,1:nvector,1:threetondim))
   allocate(nbors_grids_ok(1:twotondim,1:ngrid,1:nvector,1:twotondim))
   allocate(ind_grid1(1:twotondim,1:ngrid,1:nvector), &
   &        ind_grid2(1:twotondim,1:ngrid,1:nvector), &
   &        ind_grid3(1:twotondim,1:ngrid,1:nvector))
   allocate(nbors_grids(1:twotondim,1:ngrid,1:nvector,1:twotondim))
   allocate(ix(1:twotondim,1:ngrid,1:nvector), &
   & iy(1:twotondim,1:ngrid,1:nvector), &
   & iz(1:twotondim,1:ngrid,1:nvector), &
   & iix(1:twotondim,1:ngrid,1:nvector), & 
   & iiy(1:twotondim,1:ngrid,1:nvector), & 
   & iiz(1:twotondim,1:ngrid,1:nvector)) 

   !$acc data present(phi,phi_old,f,rho,father,son,nbor,active) &
   !$acc copyin(iii,jjj,bbbb,ccc,lll,mmm,iii_aux,jjj_aux,kkk_aux) &
   !$acc create(phi_int,ind_cell,nbors_father_grids,nbors_father_cells) &
   !$acc create(pos,ind_grid_father,ind_grid_ok,nbors_father_ok,nbors_grids_ok) &
   !$acc create(ind_grid1,ind_grid2,ind_grid3,nbors_grids) &
   !$acc create(ix,iy,iz,iix,iiy,iiz) &
   !$acc copyin(dtold,dtnew)
#if defined(_OPENACC)
   call update_globalvar_dp_to_device(rho,ilevel)
#endif

   do block=1,num_blocks
   ! Loop over cells
   !$acc parallel loop collapse(2) gang vector vector_length(512) independent
   do ind=1,twotondim
      ! Loop over active grids
      do igrid_mg=(block-1)*ngrid+1,merge(ngrid_initial,block*ngrid,block.eq.num_blocks)
         igrid_mg_m = igrid_mg - (block-1)*ngrid
         iskip_amr = ncoarse+(ind-1)*ngridmax
         igrid_amr = active(ilevel)%igrid(igrid_mg)
         icell_amr = iskip_amr + igrid_amr

         ! Init BC-modified RHS to rho - rho_tot :
         f(icell_amr,2) = fourpi*(rho(icell_amr) - rho_tot)
         if(f(icell_amr,3)<=0.0) cycle ! Do not process masked cells

         ! Separate directions 
         do idim=1,ndim
            ! Loop over the 2 neighbors
            do inbor=1,2
               ! Get neighbor grid shift
               igshift = iii(idim,inbor,ind)

               ! Get neighbor grid and its parent cell
               if(igshift==0) then
                  igrid_nbor_amr = igrid_amr
                  ifathercell_nbor_amr = father(igrid_nbor_amr)
               else
                  igrid_nbor_amr = son(nbor(igrid_amr,igshift))
                  ifathercell_nbor_amr = nbor(igrid_amr,igshift)
               end if

               if(igrid_nbor_amr==0) then
                  ! No neighbor: set mask to -1 and interp. phi
                  nb_mask = -1.0d0

                  ! Interpolate from upper level
                  ind_cell(ind,igrid_mg_m,1)=ifathercell_nbor_amr

                  call interpol_phi_mfbrv(1,ilevel,icount,ind,igrid_mg_m, &
                                          ix,iy,iz,iix,iiy,iiz,ngrid)

                  nb_phi = phi_int(ind,igrid_mg_m,1,jjj(idim,inbor,ind))
               else
                  ! Fetch neighbor cell id
                  icell_nbor_amr = igrid_nbor_amr + (ncoarse + (jjj(idim,inbor,ind)-1)*ngridmax)
                  ! Check neighbor cell mask
                  nb_mask = f(icell_nbor_amr,3)
                  if(nb_mask>0) cycle ! Neighbor cell is active too: cycle
                  nb_phi  = phi(icell_nbor_amr)
               end if
               ! phi(#) interpolated with mask:
               w = nb_mask/(nb_mask-f(icell_amr,3)) ! Linear parameter
               phi_b = ((1.0d0-w)*nb_phi + w*phi(icell_amr))

               ! Increment correction for current cell
               f(icell_amr,2) = f(icell_amr,2) - 2.0d0*oneoverdx2*phi_b
            end do
         end do
      end do
   end do
   end do

   !$acc end data

   ! clean up
   deallocate(phi_int)
   deallocate(ind_cell)
   deallocate(nbors_father_grids)
   deallocate(nbors_father_cells)
   deallocate(pos,ind_grid_father,ind_grid_ok)
   deallocate(nbors_father_ok)
   deallocate(nbors_grids_ok)
   deallocate(ind_grid1,ind_grid2,ind_grid3)
   deallocate(nbors_grids)
   deallocate(ix,iy,iz,iix,iiy,iiz)

#ifdef MPROFILING
  tt2=MPI_WTIME(info)
  acc_t_make_fine_bc_rhs = acc_t_make_fine_bc_rhs + (tt2-tt1)
#endif    
   
end subroutine make_fine_bc_rhs

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

! ------------------------------------------------------------------------
! MPI routines for MG communication for CPU boundaries,
! Those are the MG versions of the make_virtual_* AMR routines
! ------------------------------------------------------------------------

subroutine make_virtual_mg_dp(ivar,ilevel)
  use amr_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer,dimension(MPI_STATUS_SIZE,ncpu)::statuses
#endif
  integer::ilevel,ivar,icell
  integer::icpu,i,j,ncache,iskip,step
  integer::countsend,countrecv
  integer::info,tag=101
  integer,dimension(ncpu)::reqsend,reqrecv

#ifndef WITHOUTMPI
  ! Receive all messages
  countrecv=0
  do icpu=1,ncpu
     if(icpu==myid)cycle
     ncache=active_mg(icpu,ilevel)%ngrid
     if(ncache>0) then
       countrecv=countrecv+1
       call MPI_IRECV(active_mg(icpu,ilevel)%u(1,ivar),ncache*twotondim, &
            & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
     end if
  end do

  ! Gather emission array
  do icpu=1,ncpu
     if (emission_mg(icpu,ilevel)%ngrid>0) then
        do j=1,twotondim
           step=(j-1)*emission_mg(icpu,ilevel)%ngrid
           iskip=(j-1)*active_mg(myid,ilevel)%ngrid
           do i=1,emission_mg(icpu,ilevel)%ngrid
              icell=emission_mg(icpu,ilevel)%igrid(i)+iskip
              emission_mg(icpu,ilevel)%u(i+step,1)=active_mg(myid,ilevel)%u(icell,ivar)
           end do
        end do
     end if
  end do

  ! Send all messages
  countsend=0
  do icpu=1,ncpu
     ncache=emission_mg(icpu,ilevel)%ngrid
     if(ncache>0) then
       countsend=countsend+1
       call MPI_ISEND(emission_mg(icpu,ilevel)%u,ncache*twotondim, &
            & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
     end if
  end do

  ! Wait for full completion of receives
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)

#endif

111 format('   Entering make_virtual_mg for level ',I2)

end subroutine make_virtual_mg_dp

! ########################################################################
! ########################################################################

subroutine make_virtual_mg_dp_acc(ivar,ilevel)
  use amr_commons
  use poisson_commons
  use poisson_commons_acc
  use ind_accessing
  implicit none
#ifndef WITHOUTMPI
#ifndef INCLUDEOK
#define INCLUDEOK
  include 'mpif.h'
#endif
  integer,dimension(MPI_STATUS_SIZE,ncpu)::statuses
#endif
  integer::ilevel,ivar,icell
  integer::icpu,i,j,ncache,iskip,step
  integer::countsend,countrecv
  integer::info,tag=101
  integer,dimension(ncpu)::reqsend,reqrecv
  integer  :: ngrid_maxval
#ifdef MPROFILING
#ifndef INCLUDEOK
#define INCLUDEOK
  include 'mpif.h'
#endif
  real(kind=8)::tt1,tt2
  tt1=MPI_WTIME(info)
#endif
#undef INCLUDEOK

#ifndef WITHOUTMPI

  !$acc data present(active_mg_u,emission_mg_u,ngrid_a_mg,ngrid_e_mg,emission_mg_igrid)

  ! Receive all messages
  countrecv=0
  do icpu=1,ncpu
     if(icpu==myid)cycle
     ncache=ngrid_a_mg(icpu,ilevel)
     if(ncache>0) then
       countrecv=countrecv+1
       active_ind=ind_acc_a_mg_u(ilevel,icpu,ivar)
       !$acc host_data use_device(active_mg_u)
       call MPI_IRECV(active_mg_u(active_ind+1),ncache*twotondim, &
            & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
       !$acc end host_data
     end if
  end do

  ngrid_maxval = maxval(ngrid_e_mg(:,ilevel))

  ! Gather emission array
  !$acc parallel loop collapse(3) gang worker vector independent &
  !$acc vector_length(32) num_workers(16)
  do icpu=1,ncpu
  do j=1,twotondim
  do i=1,ngrid_maxval
     if (i.GT.ngrid_e_mg(icpu,ilevel)) cycle
     step=(j-1)*ngrid_e_mg(icpu,ilevel)
     iskip=(j-1)*ngrid_a_mg(myid,ilevel)
     tmp_ind = ind_acc_e_mg_igrid(ilevel,icpu)+i
     icell=emission_mg_igrid(tmp_ind)+iskip
     active_ind  =ind_acc_a_mg_u(ilevel,myid,ivar)+icell
     emission_ind=ind_acc_e_mg_u(ilevel,icpu,1)   +i+step
     emission_mg_u(emission_ind)=active_mg_u(active_ind)
  end do
  end do
  end do

  ! Send all messages
  countsend=0
  do icpu=1,ncpu
     ncache=ngrid_e_mg(icpu,ilevel)
     if(ncache>0) then
       countsend=countsend+1
       emission_ind=ind_acc_e_mg_u(ilevel,icpu,1)
       !$acc host_data use_device(emission_mg_u)
       call MPI_ISEND(emission_mg_u(emission_ind+1),ncache*twotondim, &
            & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
       !$acc end host_data
     end if
  end do

  ! Wait for full completion of receives
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)

  !$acc end data

#endif

#ifdef MPROFILING
  tt2=MPI_WTIME(info)
  acc_t_make_virtual_mg_dp = acc_t_make_virtual_mg_dp + (tt2-tt1)
#endif  

111 format('   Entering make_virtual_mg for level ',I2)

end subroutine make_virtual_mg_dp_acc

! ########################################################################
! ########################################################################

subroutine make_virtual_mg_int(ilevel)
  use amr_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer,dimension(MPI_STATUS_SIZE,ncpu)::statuses
#endif
  integer::ilevel
  integer::icpu,i,j,ncache,iskip,step,icell
  integer::countsend,countrecv
  integer::info,tag=101
  integer,dimension(ncpu)::reqsend,reqrecv

#ifndef WITHOUTMPI
  ! Receive all messages
  countrecv=0
  do icpu=1,ncpu
     if(icpu==myid)cycle
     ncache=active_mg(icpu,ilevel)%ngrid
     if(ncache>0) then
       countrecv=countrecv+1
       call MPI_IRECV(active_mg(icpu,ilevel)%f(1,1),ncache*twotondim, &
            & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
     end if
  end do

  ! Gather emission array
  do icpu=1,ncpu
     if (emission_mg(icpu,ilevel)%ngrid>0) then
        do j=1,twotondim
           step=(j-1)*emission_mg(icpu,ilevel)%ngrid
           iskip=(j-1)*active_mg(myid,ilevel)%ngrid
           do i=1,emission_mg(icpu,ilevel)%ngrid
              icell=emission_mg(icpu,ilevel)%igrid(i)+iskip
              emission_mg(icpu,ilevel)%f(i+step,1)=active_mg(myid,ilevel)%f(icell,1)
           end do
        end do
     end if
  end do

  ! Send all messages
  countsend=0
  do icpu=1,ncpu
     ncache=emission_mg(icpu,ilevel)%ngrid
     if(ncache>0) then
       countsend=countsend+1
       call MPI_ISEND(emission_mg(icpu,ilevel)%f,ncache*twotondim, &
            & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
     end if
  end do

  ! Wait for full completion of receives
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)

#endif

111 format('   Entering make_virtual_mg for level ',I2)

end subroutine make_virtual_mg_int

! ########################################################################
! ########################################################################

subroutine make_reverse_mg_dp(ivar,ilevel)
  use amr_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer,dimension(MPI_STATUS_SIZE,ncpu)::statuses
#endif
  integer::ilevel,ivar,icell
  integer::icpu,i,j,ncache,iskip,step
  integer::countsend,countrecv
  integer::info,tag=101
  integer,dimension(ncpu)::reqsend,reqrecv

#ifndef WITHOUTMPI
  ! Receive all messages
  countrecv=0
  do icpu=1,ncpu
     ncache=emission_mg(icpu,ilevel)%ngrid
     if(ncache>0) then
       countrecv=countrecv+1
       call MPI_IRECV(emission_mg(icpu,ilevel)%u,ncache*twotondim, &
            & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
     end if
  end do

  ! Send all messages
  countsend=0
  do icpu=1,ncpu
     if(icpu==myid)cycle
     ncache=active_mg(icpu,ilevel)%ngrid
     if(ncache>0) then
       countsend=countsend+1
       call MPI_ISEND(active_mg(icpu,ilevel)%u(1,ivar),ncache*twotondim, &
            & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
     end if
  end do

  ! Wait for full completion of receives
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)

  ! Gather emission array
  do icpu=1,ncpu
     if (emission_mg(icpu,ilevel)%ngrid>0) then
        do j=1,twotondim
           step=(j-1)*emission_mg(icpu,ilevel)%ngrid
           iskip=(j-1)*active_mg(myid,ilevel)%ngrid
           do i=1,emission_mg(icpu,ilevel)%ngrid
              icell=emission_mg(icpu,ilevel)%igrid(i)+iskip
              active_mg(myid,ilevel)%u(icell,ivar)=active_mg(myid,ilevel)%u(icell,ivar)+ &
                   & emission_mg(icpu,ilevel)%u(i+step,1)
           end do
        end do
     end if
  end do

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)

#endif

111 format('   Entering make_reverse_mg for level ',I2)

end subroutine make_reverse_mg_dp

! ########################################################################
! ########################################################################

subroutine make_reverse_mg_dp_acc(ivar,ilevel)
  use amr_commons
  use poisson_commons
  use poisson_commons_acc
  use ind_accessing
  implicit none
#ifndef WITHOUTMPI
#ifndef INCLUDEOK
#define INCLUDEOK
  include 'mpif.h'
#endif
  integer,dimension(MPI_STATUS_SIZE,ncpu)::statuses
#endif
  integer::ilevel,ivar,icell
  integer::icpu,i,j,ncache,iskip,step
  integer::countsend,countrecv
  integer::info,tag=101
  integer,dimension(ncpu)::reqsend,reqrecv
#ifdef MPROFILING
#ifndef INCLUDEOK
#define INCLUDEOK
  include 'mpif.h'
#endif
  real(kind=8)::tt1,tt2
  tt1=MPI_WTIME(info)
#endif
#undef INCLUDEOK

#ifndef WITHOUTMPI

  !$acc data present(active_mg_u,emission_mg_u,ngrid_a_mg,ngrid_e_mg,emission_mg_igrid)
  
  ! Receive all messages
  countrecv=0
  do icpu=1,ncpu
     ncache=ngrid_e_mg(icpu,ilevel)
     if(ncache>0) then
       countrecv=countrecv+1
       emission_ind = ind_acc_e_mg_u(ilevel,icpu,1)
       !$acc host_data use_device(emission_mg_u)
       call MPI_IRECV(emission_mg_u(emission_ind+1),ncache*twotondim, &
            & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
       !$acc end host_data
     end if
  end do

  ! Send all messages
  countsend=0
  do icpu=1,ncpu
     if(icpu==myid)cycle
     ncache=ngrid_a_mg(icpu,ilevel)
     if(ncache>0) then
       countsend=countsend+1
       active_ind = ind_acc_a_mg_u(ilevel,icpu,ivar)
       !$acc host_data use_device(active_mg_u)
       call MPI_ISEND(active_mg_u(active_ind+1),ncache*twotondim, &
            & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
       !$acc end host_data
     end if
  end do

  ! Wait for full completion of receives
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)

  ! Scatter reception array
  !$acc parallel
  do icpu=1,ncpu
    if (ngrid_e_mg(icpu,ilevel)>0) then
     !$acc loop collapse(2) worker vector
     do j=1,twotondim
     do i=1,ngrid_e_mg(icpu,ilevel)
        step=(j-1)*ngrid_e_mg(icpu,ilevel)
        iskip=(j-1)*ngrid_a_mg(myid,ilevel)
        emission_ind = ind_acc_e_mg_igrid(ilevel,icpu)+i
        icell=emission_mg_igrid(emission_ind)+iskip
        active_ind  =ind_acc_a_mg_u(ilevel,myid,ivar)+icell
        emission_ind=ind_acc_e_mg_u(ilevel,icpu,1)   +i+step
        active_mg_u(active_ind)=active_mg_u(active_ind)+emission_mg_u(emission_ind)
     end do
     end do
    end if
  end do
  !$acc end parallel

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)

  !$acc end data

#endif

#ifdef MPROFILING
  tt2=MPI_WTIME(info)
  acc_t_make_reverse_mg_dp = acc_t_make_reverse_mg_dp + (tt2-tt1)
#endif

111 format('   Entering make_reverse_mg for level ',I2)

end subroutine make_reverse_mg_dp_acc

! ########################################################################
! ########################################################################

subroutine make_reverse_mg_int(ilevel)
  use amr_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer,dimension(MPI_STATUS_SIZE,ncpu)::statuses
#endif
  integer::ilevel,icell
  integer::icpu,i,j,ncache,iskip,step
  integer::countsend,countrecv
  integer::info,tag=101
  integer,dimension(ncpu)::reqsend,reqrecv

#ifndef WITHOUTMPI
  ! Receive all messages
  countrecv=0
  do icpu=1,ncpu
     ncache=emission_mg(icpu,ilevel)%ngrid
     if(ncache>0) then
       countrecv=countrecv+1
       call MPI_IRECV(emission_mg(icpu,ilevel)%f,ncache*twotondim, &
            & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
     end if
  end do

  ! Send all messages
  countsend=0
  do icpu=1,ncpu
     if(icpu==myid)cycle
     ncache=active_mg(icpu,ilevel)%ngrid
     if(ncache>0) then
       countsend=countsend+1
       call MPI_ISEND(active_mg(icpu,ilevel)%f,ncache*twotondim, &
            & MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
     end if
  end do

  ! Wait for full completion of receives
  call MPI_WAITALL(countrecv,reqrecv,statuses,info)

  ! Gather emission array
  do icpu=1,ncpu
     if (emission_mg(icpu,ilevel)%ngrid>0) then
        do j=1,twotondim
           step=(j-1)*emission_mg(icpu,ilevel)%ngrid
           iskip=(j-1)*active_mg(myid,ilevel)%ngrid
           do i=1,emission_mg(icpu,ilevel)%ngrid
              icell=emission_mg(icpu,ilevel)%igrid(i)+iskip
              active_mg(myid,ilevel)%f(icell,1)=active_mg(myid,ilevel)%f(icell,1)+&
                 & emission_mg(icpu,ilevel)%f(i+step,1)
           end do
        end do
     end if
  end do

  ! Wait for full completion of sends
  call MPI_WAITALL(countsend,reqsend,statuses,info)

#endif

111 format('   Entering make_reverse_mg for level ',I2)

end subroutine make_reverse_mg_int

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

subroutine dump_mg_levels(ilevel,idout)
   use amr_commons
   use poisson_commons
   implicit none
#ifndef WITHOUTMPI
   include 'mpif.h'
#endif
   integer, intent(in) :: idout, ilevel

   character(len=24)  :: cfile
   character(len=5)   :: ccpu='00000'
   character(len=5)   :: cout='00000'

   integer :: i, ngrids, igrid, icpu, idim

   integer,parameter::tag=1119
   integer::dummy_io,info2

   write(ccpu,'(I5.5)') myid
   write(cout,'(I5.5)') idout
   cfile='multigrid_'//cout//'.out'//ccpu

   ! Wait for the token
#ifndef WITHOUTMPI
   if(IOGROUPSIZE>0) then
      if (mod(myid-1,IOGROUPSIZE)/=0) then
         call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
              & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
      end if
   endif
#endif

   open(unit=10,file=cfile,status='unknown',form='formatted')

   write(10,'(I1)') ndim
   write(10,'(I1)') myid
   write(10,'(I1)') ncpu
   write(10,'(I2)') ilevel

   ! Loop over levels
   do i=1,ilevel-1
      ! Active grids
      ngrids=active_mg(myid,i)%ngrid
      write(10,*) ngrids
      do igrid=1,ngrids
         do idim=1,ndim
            write(10,*) xg(active_mg(myid,i)%igrid(igrid),idim)
         end do
      end do

      ! Reception grids
      do icpu=1,ncpu
         if(icpu==myid)cycle
         ngrids=active_mg(icpu,i)%ngrid
         write(10,*) ngrids
         do igrid=1,ngrids
            do idim=1,ndim
               write(10,*) xg(active_mg(icpu,i)%igrid(igrid),idim)
            end do
         end do
      end do

   end do

   close(10)

        ! Send the token
#ifndef WITHOUTMPI
   if(IOGROUPSIZE>0) then
      if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
         dummy_io=1
         call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
              & MPI_COMM_WORLD,info2)
      end if
   endif
#endif

end subroutine dump_mg_levels

! ########################################################################
! ########################################################################
! ########################################################################
! ########################################################################

