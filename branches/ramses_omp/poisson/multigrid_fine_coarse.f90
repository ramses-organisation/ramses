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
!!!$OMP PARALLEL DEFAULT(NONE) REDUCTION(.and.:allmasked) SHARED(active_mg,son,f,cpu_map,lookup_mg) PRIVATE(ind_c_cell,iskip_c_amr,iskip_c_mg,igrid_c_mg,igrid_c_amr,
!!!icell_c_amr,icell_c_mg,ngpmask,igrid_f_amr,ind_f_cell,iskip_f_mg,icell_f_mg,cpu_amr,igrid_f_mg) FIRSTPRIVATE(ncoarse,ngridmax,myid,ifinelevel,icoarselevel,dtwotondim)
!$OMP PARALLEL DEFAULT(private) REDUCTION(.and.:allmasked) SHARED(active_mg,son,f,cpu_map,lookup_mg) FIRSTPRIVATE(ncoarse,ngridmax,myid,ifinelevel,icoarselevel,dtwotondim)
   do ind_c_cell=1,twotondim
      iskip_c_amr=ncoarse+(ind_c_cell-1)*ngridmax
      iskip_c_mg =(ind_c_cell-1)*active_mg(myid,icoarselevel)%ngrid

      ! Loop over coarse grids of myid
!$OMP DO
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
!$OMP END DO NOWAIT
   end do
!$OMP END PARALLEL

end subroutine restrict_mask_coarse

! ------------------------------------------------------------------------
! Mask restriction (bottom-up)
! ------------------------------------------------------------------------

subroutine restrict_mask_coarse_reverse(ifinelevel)
   use amr_commons
   use poisson_commons
   implicit none
   integer, intent(in) :: ifinelevel

   integer :: ind_c_cell, ind_f_cell, cpu_amr
   
   integer :: iskip_c_mg
   integer :: igrid_c_amr, igrid_c_mg
   integer :: icell_c_amr, icell_c_mg

   integer :: iskip_f_mg
   integer :: igrid_f_amr, igrid_f_mg
   integer :: icell_f_mg

   real(dp) :: ngpmask
   real(dp) :: dtwotondim = (twotondim)

   integer :: icoarselevel

   icoarselevel=ifinelevel-1

   ! Loop over fine cells of the myid active comm
   do ind_f_cell=1,twotondim
      iskip_f_mg =(ind_f_cell-1)*active_mg(myid,ifinelevel)%ngrid

      ! Loop over fine grids of myid
      do igrid_f_mg=1,active_mg(myid,ifinelevel)%ngrid
         icell_f_mg=iskip_f_mg+igrid_f_mg
         igrid_f_amr=active_mg(myid,ifinelevel)%igrid(igrid_f_mg)
         ! Get coarse grid AMR index and CPU id
         icell_c_amr=father(igrid_f_amr)
         ind_c_cell=(icell_c_amr-ncoarse-1)/ngridmax+1
         igrid_c_amr=icell_c_amr-ncoarse-(ind_c_cell-1)*ngridmax
         cpu_amr=cpu_map(father(igrid_c_amr))
         ! Convert to MG index, get MG coarse cell id
         igrid_c_mg=lookup_mg(igrid_c_amr)
         iskip_c_mg=(ind_c_cell-1)*active_mg(cpu_amr,icoarselevel)%ngrid
         icell_c_mg=iskip_c_mg+igrid_c_mg
         ! Stack cell volume fraction in coarse cell
         ngpmask=(1d0+active_mg(myid,ifinelevel)%u(icell_f_mg,4))/2d0/dtwotondim
         active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,4)=&
            active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,4)+ngpmask
      end do
   end do

end subroutine restrict_mask_coarse_reverse

! ------------------------------------------------------------------------
! Residual computation
! ------------------------------------------------------------------------

subroutine cmp_residual_mg_coarse(ilevel)
   ! Computes the residual for pure MG levels, and stores it into active_mg(myid,ilevel)%u(:,3)
   use amr_commons
   use poisson_commons
   implicit none
   integer, intent(in) :: ilevel

   integer, dimension(1:3,1:2,1:8) :: iii, jjj

  integer ,dimension(1:nvector)::igrid_amr,igrid_mg,icell_mg,icell_mg_ok
  integer ,dimension(1:nvector,0:twondim)::igrid_nbor_amr,igrid_nbor_mg,cpu_nbor_amr
  integer ,dimension(1:nvector,0:twondim)::igrid_nbor_amr_ok,igrid_nbor_mg_ok,cpu_nbor_amr_ok
  integer ,dimension(1:nvector,1:ndim)::ind_left,ind_right
  real(dp),dimension(1:nvector)::nb_sum,weight

   real(dp) :: dx, oneoverdx2
   integer  :: ngrid, ncache
   integer  :: i, ind, idim, inbor, iskip, igrid, ncell_ok
   integer  :: igshift, icell_nbor_mg, icpu_nbor_amr
   real(dp) :: dtwondim = (twondim)

   ! Set constants
   dx  = 0.5d0**ilevel
   oneoverdx2 = 1.0d0/(dx*dx)

   iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
   iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

  ! Loop over active grids
  ncache=active_mg(myid,ilevel)%ngrid
!!!$OMP PARALLEL DEFAULT(NONE) SHARED(active_mg,son,cpu_map,nbor,lookup_mg) PRIVATE(ind,iskip,igrid,icell_mg,icell_mg_ok,nb_sum,inbor,idim,igshift,igrid_nbor_amr,cpu_nbor_amr,
!!igrid_amr,igrid_mg,ind_left,ind_right,ncell_ok,igrid_nbor_mg,igrid_nbor_amr_ok,cpu_nbor_amr_ok,igrid_nbor_mg_ok,icell_nbor_mg,icpu_nbor_amr) FIRSTPRIVATE(ngrid,ilevel,myid,iii,jjj,oneoverdx2,dtwondim,ncache)
!$OMP PARALLEL DEFAULT(private) SHARED(active_mg,son,cpu_map,nbor,lookup_mg) FIRSTPRIVATE(ngrid,ilevel,myid,iii,jjj,oneoverdx2,dtwondim,ncache)
!$OMP DO SCHEDULE(DYNAMIC) 
  do igrid=1,ncache,nvector
     
     ! Gather nvector grids
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        igrid_amr(i)=active_mg(myid,ilevel)%igrid(igrid+i-1)
        igrid_mg(i)=igrid+i-1
     end do

     ! Gather neighboring grids
     do i=1,ngrid
        igrid_nbor_amr(i,0)=igrid_amr(i)
        cpu_nbor_amr(i,0)=myid
        igrid_nbor_mg(i,0)=lookup_mg(igrid_amr(i))
     end do
     do idim=1,ndim
        do i=1,ngrid
           ind_left (i,idim)=nbor(igrid_amr(i),2*idim-1)
           ind_right(i,idim)=nbor(igrid_amr(i),2*idim  )
           igrid_nbor_amr(i,2*idim-1)=son(ind_left (i,idim))
           igrid_nbor_amr(i,2*idim  )=son(ind_right(i,idim))
           cpu_nbor_amr(i,2*idim-1)=cpu_map(ind_left (i,idim))
           cpu_nbor_amr(i,2*idim  )=cpu_map(ind_right(i,idim))
           igrid_nbor_mg(i,2*idim-1)=lookup_mg(igrid_nbor_amr(i,2*idim-1))
           igrid_nbor_mg(i,2*idim  )=lookup_mg(igrid_nbor_amr(i,2*idim  ))
        end do
     end do

     ! Loop over cells
     do ind=1,twotondim
     
        ! Compute central cell index
        iskip=(ind-1)*ncache
        do i=1,ngrid
           icell_mg(i)=iskip+igrid_mg(i)
        end do

        ! Gather inner cells
        ncell_ok=0
        do i=1,ngrid
           if( .not. btest(active_mg(myid,ilevel)%f(icell_mg(i),1),0))then
              ncell_ok=ncell_ok+1
              icell_mg_ok(ncell_ok)=icell_mg(i)
              igrid_nbor_mg_ok(ncell_ok,0:twondim)=igrid_nbor_mg(i,0:twondim)
              cpu_nbor_amr_ok(ncell_ok,0:twondim)=cpu_nbor_amr(i,0:twondim)
           endif
        end do

        ! Use max-speed "dumb" Gauss-Seidel for "inner" cells
        ! Those cells are active, have all their neighbors active
        ! and all neighbors are in the AMR+MG trees
        nb_sum(1:ncell_ok)=0.
        do inbor=1,2
           do idim=1,ndim
              ! Get neighbor grid shift
              igshift = iii(idim,inbor,ind)
              do i=1,ncell_ok
                 icpu_nbor_amr=cpu_nbor_amr_ok(i,igshift)
                 icell_nbor_mg=igrid_nbor_mg_ok(i,igshift)+((jjj(idim,inbor,ind)-1)*active_mg(icpu_nbor_amr,ilevel)%ngrid)
                 nb_sum(i)=nb_sum(i)+active_mg(icpu_nbor_amr,ilevel)%u(icell_nbor_mg,1)
              end do
           end do
        end do
        do i=1,ncell_ok
           active_mg(myid,ilevel)%u(icell_mg_ok(i),3)= &
                & -oneoverdx2*(nb_sum(i)-dtwondim*active_mg(myid,ilevel)%u(icell_mg_ok(i),1) ) &
                & +active_mg(myid,ilevel)%u(icell_mg_ok(i),2)
        end do

        ! Gather boundary cells
        ncell_ok=0
        do i=1,ngrid
           if(btest(active_mg(myid,ilevel)%f(icell_mg(i),1),0)) then
              if (active_mg(myid,ilevel)%u(icell_mg(i),4)>0.0) then
                 ncell_ok=ncell_ok+1
                 icell_mg_ok(ncell_ok)=icell_mg(i)
                 igrid_nbor_amr_ok(ncell_ok,0:twondim)=igrid_nbor_amr(i,0:twondim)
                 igrid_nbor_mg_ok(ncell_ok,0:twondim)=igrid_nbor_mg(i,0:twondim)
                 cpu_nbor_amr_ok(ncell_ok,0:twondim)=cpu_nbor_amr(i,0:twondim)
              else
                 active_mg(myid,ilevel)%u(icell_mg(i),3)=0.0
              endif
           endif
        end do

        ! Use the finer "solve" Gauss-Seidel near boundaries,
        ! with all necessary checks
        nb_sum(1:ncell_ok)=0.
        do inbor=1,2
           do idim=1,ndim
              ! Get neighbor grid shift
              igshift = iii(idim,inbor,ind)
              do i=1,ncell_ok
                 ! In case neighbor grid does not exist
                 if(igrid_nbor_amr_ok(i,igshift)==0) then
                    ! No neighbor cell, set mask=-1 on nonexistent neighbor cell
                    nb_sum(i)=nb_sum(i)-active_mg(myid,ilevel)%u(icell_mg_ok(i),1)/active_mg(myid,ilevel)%u(icell_mg_ok(i),4)
                 else
                    if(igrid_nbor_mg_ok(i,igshift)<=0)then
                       ! No MG neighbor
                       nb_sum(i)=nb_sum(i)-active_mg(myid,ilevel)%u(icell_mg_ok(i),1)/active_mg(myid,ilevel)%u(icell_mg_ok(i),4)
                    else
                       ! Fetch neighbor cell
                       icpu_nbor_amr=cpu_nbor_amr_ok(i,igshift)
                       icell_nbor_mg=igrid_nbor_mg_ok(i,igshift)+((jjj(idim,inbor,ind)-1)*active_mg(icpu_nbor_amr,ilevel)%ngrid)
                       if(active_mg(icpu_nbor_amr,ilevel)%u(icell_nbor_mg,4)<=0.0) then
                          ! Neighbor cell is masked
                          nb_sum(i)=nb_sum(i)+active_mg(myid,ilevel)%u(icell_mg_ok(i),1)* &
                               & active_mg(icpu_nbor_amr,ilevel)%u(icell_nbor_mg,4) &
                               & /active_mg(myid,ilevel)%u(icell_mg_ok(i),4)
                       else
                          ! Neighbor cell is active, increment neighbor sum
                          nb_sum(i)=nb_sum(i)+active_mg(icpu_nbor_amr,ilevel)%u(icell_nbor_mg,1)
                       end if
                    endif
                 end if
              end do
           end do
        end do
        do i=1,ncell_ok
           active_mg(myid,ilevel)%u(icell_mg_ok(i),3)= &
                & -oneoverdx2*(nb_sum(i)-dtwondim*active_mg(myid,ilevel)%u(icell_mg_ok(i),1) ) &
                & +active_mg(myid,ilevel)%u(icell_mg_ok(i),2)
        end do

     end do
     ! End loop over cells
  end do
!$OMP END PARALLEL

end subroutine cmp_residual_mg_coarse

! ##################################################################
! ##################################################################

subroutine cmp_uvar_norm2_coarse(ivar, ilevel, norm2)
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
  implicit none
  integer, intent(in) :: ilevel
  logical, intent(in) :: safe
  logical, intent(in) :: redstep
  
  integer, dimension(1:3,1:2,1:8) :: iii, jjj
  integer, dimension(1:3,1:4)     :: ired, iblack
  
  integer ,dimension(1:nvector)::igrid_amr,igrid_mg,icell_mg,icell_mg_ok
  integer ,dimension(1:nvector,0:twondim)::igrid_nbor_amr,igrid_nbor_mg,cpu_nbor_amr
  integer ,dimension(1:nvector,0:twondim)::igrid_nbor_amr_ok,igrid_nbor_mg_ok,cpu_nbor_amr_ok
  integer ,dimension(1:nvector,1:ndim)::ind_left,ind_right
  real(dp),dimension(1:nvector)::nb_sum,weight

  real(dp) :: dx2
  integer  :: ngrid, ncache
  integer  :: i, ind, ind0, idim, inbor
  integer  :: iskip, igrid, ncell_ok
  integer  :: igshift, icell_nbor_mg, icpu_nbor_amr
  real(dp) :: dtwondim = (twondim)
  
  ! Set constants
  dx2  = (0.5d0**ilevel)**2
  
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
  
  ! Loop over active grids
  ncache=active_mg(myid,ilevel)%ngrid
!!!$OMP PARALLEL DEFAULT(NONE) SHARED(active_mg,son,cpu_map,nbor,lookup_mg) PRIVATE(ind0,ind,iskip,igrid,icell_mg,icell_mg_ok,nb_sum,inbor,idim,igshift,igrid_nbor_amr,cpu_nbor_amr,
!!igrid_amr,igrid_mg,ind_left,ind_right,ncell_ok,igrid_nbor_mg,igrid_nbor_amr_ok,cpu_nbor_amr_ok,igrid_nbor_mg_ok,icell_nbor_mg,icpu_nbor_amr,weight) FIRSTPRIVATE(ired,iblack,ngrid,ilevel,myid,iii,jjj,dx2,dtwondim,safe,redstep,ncache)
!$OMP PARALLEL DEFAULT(private) SHARED(active_mg,son,cpu_map,nbor,lookup_mg) FIRSTPRIVATE(ired,iblack,ngrid,ilevel,myid,iii,jjj,dx2,dtwondim,safe,redstep,ncache)
!$OMP DO SCHEDULE(DYNAMIC) 
  do igrid=1,ncache,nvector
     
     ! Gather nvector grids
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        igrid_amr(i)=active_mg(myid,ilevel)%igrid(igrid+i-1)
        igrid_mg(i)=igrid+i-1
     end do

     ! Gather neighboring grids
     do i=1,ngrid
        igrid_nbor_amr(i,0)=igrid_amr(i)
        cpu_nbor_amr(i,0)=myid
        igrid_nbor_mg(i,0)=lookup_mg(igrid_amr(i))
     end do
     do idim=1,ndim
        do i=1,ngrid
           ind_left (i,idim)=nbor(igrid_amr(i),2*idim-1)
           ind_right(i,idim)=nbor(igrid_amr(i),2*idim  )
           igrid_nbor_amr(i,2*idim-1)=son(ind_left (i,idim))
           igrid_nbor_amr(i,2*idim  )=son(ind_right(i,idim))
           cpu_nbor_amr(i,2*idim-1)=cpu_map(ind_left (i,idim))
           cpu_nbor_amr(i,2*idim  )=cpu_map(ind_right(i,idim))
           igrid_nbor_mg(i,2*idim-1)=lookup_mg(igrid_nbor_amr(i,2*idim-1))
           igrid_nbor_mg(i,2*idim  )=lookup_mg(igrid_nbor_amr(i,2*idim  ))
        end do
     end do

     ! Only half of the cells for a red or black sweep
     do ind0=1,twotondim/2    
        if(redstep) then
           ind=ired  (ndim,ind0)
        else
           ind=iblack(ndim,ind0)
        end if
     
        ! Compute central cell index
        iskip=(ind-1)*ncache
        do i=1,ngrid
           icell_mg(i)=iskip+igrid_mg(i)
        end do

        ! Gather inner cells
        ncell_ok=0
        do i=1,ngrid
           if( .not. btest(active_mg(myid,ilevel)%f(icell_mg(i),1),0))then
              ncell_ok=ncell_ok+1
              icell_mg_ok(ncell_ok)=icell_mg(i)
              igrid_nbor_mg_ok(ncell_ok,0:twondim)=igrid_nbor_mg(i,0:twondim)
              cpu_nbor_amr_ok(ncell_ok,0:twondim)=cpu_nbor_amr(i,0:twondim)
           endif
        end do

        ! Use max-speed "dumb" Gauss-Seidel for "inner" cells
        ! Those cells are active, have all their neighbors active
        ! and all neighbors are in the AMR+MG trees
        nb_sum(1:ncell_ok)=0.
        do inbor=1,2
           do idim=1,ndim
              ! Get neighbor grid shift
              igshift = iii(idim,inbor,ind)
              do i=1,ncell_ok
                 icpu_nbor_amr=cpu_nbor_amr_ok(i,igshift)
                 icell_nbor_mg=igrid_nbor_mg_ok(i,igshift)+((jjj(idim,inbor,ind)-1)*active_mg(icpu_nbor_amr,ilevel)%ngrid)
                 nb_sum(i)=nb_sum(i)+active_mg(icpu_nbor_amr,ilevel)%u(icell_nbor_mg,1)
              end do
           end do
        end do
        do i=1,ncell_ok
           active_mg(myid,ilevel)%u(icell_mg_ok(i),1)=(nb_sum(i)-dx2*active_mg(myid,ilevel)%u(icell_mg_ok(i),2))/dtwondim
        end do

        ! Gather boundary cells
        ncell_ok=0
        do i=1,ngrid
           if(btest(active_mg(myid,ilevel)%f(icell_mg(i),1),0)) then
              if (active_mg(myid,ilevel)%u(icell_mg(i),4)>0.0) then
                 if (.not. safe .or.  active_mg(myid,ilevel)%u(icell_mg(i),4)>=1.0) then
                    ncell_ok=ncell_ok+1
                    icell_mg_ok(ncell_ok)=icell_mg(i)
                    igrid_nbor_amr_ok(ncell_ok,0:twondim)=igrid_nbor_amr(i,0:twondim)
                    igrid_nbor_mg_ok(ncell_ok,0:twondim)=igrid_nbor_mg(i,0:twondim)
                    cpu_nbor_amr_ok(ncell_ok,0:twondim)=cpu_nbor_amr(i,0:twondim)
                 endif
              endif
           endif
        end do

        ! Use the finer "solve" Gauss-Seidel near boundaries,
        ! with all necessary checks
        nb_sum(1:ncell_ok)=0.
        weight(1:ncell_ok)=0.
        do inbor=1,2
           do idim=1,ndim
              ! Get neighbor grid shift
              igshift = iii(idim,inbor,ind)
              do i=1,ncell_ok
                 ! In case neighbor grid does not exist
                 if(igrid_nbor_amr_ok(i,igshift)==0) then
                    ! No neighbor cell, set mask=-1 on nonexistent neighbor cell
                    weight(i)=weight(i)-1.0d0/active_mg(myid,ilevel)%u(icell_mg_ok(i),4)
                 else
                    if(igrid_nbor_mg_ok(i,igshift)<=0)then
                       ! No MG neighbor
                       weight(i)=weight(i)-1.0d0/active_mg(myid,ilevel)%u(icell_mg_ok(i),4)
                    else
                       ! Fetch neighbor cell
                       icpu_nbor_amr=cpu_nbor_amr_ok(i,igshift)
                       icell_nbor_mg=igrid_nbor_mg_ok(i,igshift)+((jjj(idim,inbor,ind)-1)*active_mg(icpu_nbor_amr,ilevel)%ngrid)
                       if(active_mg(icpu_nbor_amr,ilevel)%u(icell_nbor_mg,4)<=0.0) then
                          ! Neighbor cell is masked
                          weight(i)=weight(i)+ &
                               & active_mg(icpu_nbor_amr,ilevel)%u(icell_nbor_mg,4) &
                               & /active_mg(myid,ilevel)%u(icell_mg_ok(i),4)
                       else
                          ! Neighbor cell is active, increment neighbor sum
                          nb_sum(i)=nb_sum(i)+active_mg(icpu_nbor_amr,ilevel)%u(icell_nbor_mg,1)
                       end if
                    endif
                 end if
              end do
           end do
        end do
        do i=1,ncell_ok
           active_mg(myid,ilevel)%u(icell_mg_ok(i),1)=(nb_sum(i)-dx2*active_mg(myid,ilevel)%u(icell_mg_ok(i),2))/(dtwondim-weight(i))
        end do

     end do
     ! End loop over cells
  end do
!$OMP END PARALLEL
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
   icoarselevel=ifinelevel-1

   ! Loop over fine cells of the myid active comm
!!!$OMP PARALLEL DEFAULT(NONE) SHARED(active_mg,father,cpu_map,lookup_mg) PRIVATE(ind_f_cell,iskip_f_mg,igrid_f_mg,icell_f_mg,igrid_f_amr,icell_c_amr,ind_c_cell,igrid_c_amr,cpu_amr,
!!!igrid_c_mg,iskip_c_mg,icell_c_mg,res) FIRSTPRIVATE(myid,ifinelevel,icoarselevel,ngridmax,ncoarse,dtwotondim)
!$OMP PARALLEL DEFAULT(private) SHARED(active_mg,father,cpu_map,lookup_mg) FIRSTPRIVATE(myid,ifinelevel,icoarselevel,ngridmax,ncoarse,dtwotondim)
   do ind_f_cell=1,twotondim
      iskip_f_mg =(ind_f_cell-1)*active_mg(myid,ifinelevel)%ngrid

      ! Loop over fine grids of myid
!$OMP DO
      do igrid_f_mg=1,active_mg(myid,ifinelevel)%ngrid
         icell_f_mg=iskip_f_mg+igrid_f_mg
         ! Is fine cell masked?
         if(active_mg(myid,ifinelevel)%u(icell_f_mg,4)<=0d0) cycle

         ! Get coarse grid AMR index and CPU id
         igrid_f_amr=active_mg(myid,ifinelevel)%igrid(igrid_f_mg)
         icell_c_amr=father(igrid_f_amr)
         ind_c_cell=(icell_c_amr-ncoarse-1)/ngridmax+1
         igrid_c_amr=icell_c_amr-ncoarse-(ind_c_cell-1)*ngridmax
         cpu_amr=cpu_map(father(igrid_c_amr))

         ! Convert to MG index, get MG coarse cell id
         igrid_c_mg=lookup_mg(igrid_c_amr)
         iskip_c_mg=(ind_c_cell-1)*active_mg(cpu_amr,icoarselevel)%ngrid
         icell_c_mg=iskip_c_mg+igrid_c_mg

         ! Is coarse cell masked?
         if(active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,4)<=0d0) cycle

         ! Stack fine cell residual in coarse cell rhs
         res=active_mg(myid,ifinelevel)%u(icell_f_mg,3)/dtwotondim
         active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,2)=&
            active_mg(cpu_amr,icoarselevel)%u(icell_c_mg,2)+res
      end do
!$OMP END DO
! WILL THIS WORK? IT'S A REDUCTION SUM ON ACTIVE_MG !!! MAYBE USE SMALL BUFFER ARRAY ?
   end do
!$OMP END PARALLEL
end subroutine restrict_residual_coarse_reverse

! ------------------------------------------------------------------------
! Interpolation and correction
! ------------------------------------------------------------------------

subroutine interpolate_and_correct_coarse(ifinelevel)
   use amr_commons
   use poisson_commons
   implicit none
   integer, intent(in) :: ifinelevel

   integer  :: i, ind_father, ind_average, ind_f, iskip_f_amr, ngrid_f, istart, nbatch
   integer  :: icell_c_amr, igrid_c_amr, igrid_c_mg, icell_c_mg, iskip_f_mg, icell_f_mg
   integer  :: icoarselevel, ind_c, cpu_c_amr

   real(dp) :: a, b, c, d, coeff
   real(dp), dimension(1:8)     :: bbb
   integer,  dimension(1:8,1:8) :: ccc

   integer,  dimension(1:nvector)                :: igrid_f_amr, icell_amr, cpu_amr
   integer,  dimension(1:nvector,1:threetondim)  :: nbors_father_cells
   integer,  dimension(1:nvector,1:twotondim)    :: nbors_father_grids
   real(dp), dimension(1:nvector)                :: corr

   ! Local constants
   a = 1.0D0/4.0D0**ndim
   b = 3.0D0*a
   c = 9.0D0*a
   d = 27.D0*a
   icoarselevel=ifinelevel-1

   bbb(:)  =(/a ,b ,b ,c ,b ,c ,c ,d/)

   ccc(:,1)=(/1 ,2 ,4 ,5 ,10,11,13,14/)
   ccc(:,2)=(/3 ,2 ,6 ,5 ,12,11,15,14/)
   ccc(:,3)=(/7 ,8 ,4 ,5 ,16,17,13,14/)
   ccc(:,4)=(/9 ,8 ,6 ,5 ,18,17,15,14/)
   ccc(:,5)=(/19,20,22,23,10,11,13,14/)
   ccc(:,6)=(/21,20,24,23,12,11,15,14/)
   ccc(:,7)=(/25,26,22,23,16,17,13,14/)
   ccc(:,8)=(/27,26,24,23,18,17,15,14/)

   ! Loop over fine grids by vector sweeps
   ngrid_f=active_mg(myid,ifinelevel)%ngrid
!!!!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(DYNAMIC) SHARED(active_mg,father,cpu_map,lookup_mg) PRIVATE(istart,nbatch,i,igrid_f_amr,icell_amr,cpu_amr,nbors_father_cells,nbors_father_grids,ind_f,
!!!iskip_f_amr,iskip_f_mg,corr,ind_average,ind_father,coeff,icell_f_mg,icell_c_amr,ind_c,igrid_c_amr,igrid_c_mg,cpu_c_amr,icell_c_mg) FIRSTPRIVATE(ngrid_f,ifinelevel,bbb,ccc,icoarselevel,ncoarse,ngridmax,myid)
!$OMP PARALLEL DO DEFAULT(private) SCHEDULE(DYNAMIC) SHARED(active_mg,father,cpu_map,lookup_mg) FIRSTPRIVATE(ngrid_f,ifinelevel,bbb,ccc,icoarselevel,ncoarse,ngridmax,myid)
   do istart=1,ngrid_f,nvector

      ! Gather nvector grids
      nbatch=MIN(nvector,ngrid_f-istart+1)
      do i=1,nbatch
         igrid_f_amr(i)=active_mg(myid,ifinelevel)%igrid(istart+i-1)
      end do

      ! Compute father (coarse) cell index
      do i=1,nbatch
         icell_amr(i)=father(igrid_f_amr(i))
         cpu_amr(i)  =cpu_map(icell_amr(i))
      end do

      ! Gather 3x3x3 neighboring parent cells
      call get3cubefather(icell_amr,nbors_father_cells,nbors_father_grids,nbatch,ifinelevel)

      ! Update solution for fine grid cells
      do ind_f=1,twotondim
         iskip_f_amr = ncoarse+(ind_f-1)*ngridmax
         iskip_f_mg  = (ind_f-1)*ngrid_f

         do i=1,nbatch
            ! Compute fine cell indices
            icell_amr(i) = iskip_f_amr + igrid_f_amr(i)
         end do
         corr=0.0d0

         ! Loop over relevant parent cells
         do ind_average=1,twotondim
            ind_father = ccc(ind_average,ind_f)
            coeff      = bbb(ind_average)
            do i=1,nbatch
               icell_f_mg  = iskip_f_mg + istart+i-1
               if(active_mg(cpu_amr(i),ifinelevel)%u(icell_f_mg,4)<=0.0) then
                  corr(i)=0.0d0        ! Fine cell is masked : no correction
                  cycle
               end if
               icell_c_amr = nbors_father_cells(i,ind_father)
               ind_c       = (icell_c_amr-ncoarse)/ngridmax + 1
               igrid_c_amr = icell_c_amr - ncoarse - (ind_c-1)*ngridmax
               igrid_c_mg  = lookup_mg(igrid_c_amr)
               cpu_c_amr   = cpu_map(father(igrid_c_amr))
               if(igrid_c_mg<=0) cycle

               icell_c_mg  = (ind_c-1)*active_mg(cpu_c_amr,icoarselevel)%ngrid + igrid_c_mg
               corr(i)=corr(i)+coeff*active_mg(cpu_c_amr,icoarselevel)%u(icell_c_mg,1)
            end do
         end do

         ! Correct potential
         do i=1,nbatch
            icell_f_mg  = iskip_f_mg + istart+i-1
            active_mg(cpu_amr(i),ifinelevel)%u(icell_f_mg,1) = active_mg(cpu_amr(i),ifinelevel)%u(icell_f_mg,1) + corr(i)
         end do

      end do
      ! End loop over cells

   end do
!$OMP END PARALLEL DO
   ! End loop over grids
end subroutine interpolate_and_correct_coarse


! ------------------------------------------------------------------------
! Flag setting
! ------------------------------------------------------------------------

subroutine set_scan_flag_coarse(ilevel)
  use amr_commons
  use poisson_commons
  implicit none
  
  integer, intent(in) :: ilevel
  
  integer :: i, ind, ngrid, ncache
  integer :: inbor, idim, igshift
  integer :: icpu_nbor_amr, icell_nbor_mg
  
  integer :: igrid, iskip, ncell_ok
  logical ,dimension(1:nvector)::nbor_son_exist
  integer ,dimension(1:nvector)::igrid_amr,igrid_mg,icell_mg,icell_mg_ok
  integer ,dimension(1:nvector,0:twondim)::igrid_nbor_amr,igrid_nbor_mg,cpu_nbor_amr
  integer ,dimension(1:nvector,0:twondim)::igrid_nbor_amr_ok,igrid_nbor_mg_ok,cpu_nbor_amr_ok
  integer ,dimension(1:nvector,1:ndim)::ind_left,ind_right
  integer ,dimension(1:nvector)::scan_flag
  
  integer, dimension(1:3,1:2,1:8) :: iii, jjj
  
  iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)
  
  ncache=active_mg(myid,ilevel)%ngrid
  if(ncache==0) return
  
  ! Loop over cells and set coarse SCAN flag
!!!!$OMP PARALLEL DEFAULT(none) SHARED(active_mg,son,nbor,cpu_map,lookup_mg) PRIVATE(i,ind,ind_left,ind_right,ncell_ok,iskip,igrid,igrid_mg,igrid_amr,icell_mg,icell_mg_ok,scan_flag,idim,
!!!igshift,inbor,igrid_nbor_amr,igrid_nbor_amr_ok,igrid_nbor_mg,igrid_nbor_mg_ok,cpu_nbor_amr,cpu_nbor_amr_ok,icell_nbor_mg,icpu_nbor_amr) FIRSTPRIVATE(myid,ilevel,ngrid,iii,jjj,ncache)
!$OMP PARALLEL DEFAULT(private) SHARED(active_mg,son,nbor,cpu_map,lookup_mg) FIRSTPRIVATE(myid,ilevel,ngrid,iii,jjj,ncache)
!$OMP DO SCHEDULE(DYNAMIC) 
  do igrid=1,ncache,nvector
     
     ! Gather nvector grids
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        igrid_amr(i)=active_mg(myid,ilevel)%igrid(igrid+i-1)
        igrid_mg(i)=igrid+i-1
     end do
     
     ! Gather neighboring grids
     do i=1,ngrid
        igrid_nbor_amr(i,0)=igrid_amr(i)
        cpu_nbor_amr(i,0)=myid
        igrid_nbor_mg(i,0)=lookup_mg(igrid_amr(i))
     end do
     do idim=1,ndim
        do i=1,ngrid
           ind_left (i,idim)=nbor(igrid_amr(i),2*idim-1)
           ind_right(i,idim)=nbor(igrid_amr(i),2*idim  )
           if(son(ind_left (i,idim))/=0.and.son(ind_right(i,idim))/=0)then
              nbor_son_exist(i)=.true.
              igrid_nbor_amr(i,2*idim-1)=son(ind_left (i,idim))
              igrid_nbor_amr(i,2*idim  )=son(ind_right(i,idim))
              cpu_nbor_amr(i,2*idim-1)=cpu_map(ind_left (i,idim))
              cpu_nbor_amr(i,2*idim  )=cpu_map(ind_right(i,idim))
              igrid_nbor_mg(i,2*idim-1)=lookup_mg(igrid_nbor_amr(i,2*idim-1))
              igrid_nbor_mg(i,2*idim  )=lookup_mg(igrid_nbor_amr(i,2*idim  ))
           else
              nbor_son_exist(i)=.false.
           endif
        end do
     end do
     
     ! Loop over cells
     do ind=1,twotondim
        
        ! Compute central cell index
        iskip=(ind-1)*ncache
        do i=1,ngrid
           icell_mg(i)=iskip+igrid_mg(i)
        end do
        
        ! Select flagged cells
        ncell_ok=0
        do i=1,ngrid
           if(active_mg(myid,ilevel)%u(icell_mg(i),4)==1d0.and.nbor_son_exist(i))then
              ncell_ok=ncell_ok+1
              icell_mg_ok(ncell_ok)=icell_mg(i)
              igrid_nbor_amr_ok(ncell_ok,0:twondim)=igrid_nbor_amr(i,0:twondim)
              igrid_nbor_mg_ok(ncell_ok,0:twondim)=igrid_nbor_mg(i,0:twondim)
              cpu_nbor_amr_ok(ncell_ok,0:twondim)=cpu_nbor_amr(i,0:twondim)
           else
              active_mg(myid,ilevel)%f(icell_mg(i),1)=1
           endif
        end do
        
        ! Init flag to 'no scan needed'
        scan_flag(1:ncell_ok)=0
        do inbor=1,2
           do idim=1,ndim
              ! Get neighbor grid shift
              igshift=iii(idim,inbor,ind)
              do i=1,ncell_ok
                 if(igrid_nbor_amr_ok(i,igshift)==0) then
                    scan_flag(i)=1
                 else
                    if(igrid_nbor_mg_ok(i,igshift)<=0)then
                       scan_flag(i)=1
                    else
                       icpu_nbor_amr=cpu_nbor_amr_ok(i,igshift)
                       icell_nbor_mg=igrid_nbor_mg_ok(i,igshift)+((jjj(idim,inbor,ind)-1)*active_mg(icpu_nbor_amr,ilevel)%ngrid)
                       if(active_mg(icpu_nbor_amr,ilevel)%u(icell_nbor_mg,4)<=0.0) then
                          scan_flag(i)=1
                       endif
                    endif
                 endif
              end do
           end do
        end do
        do i=1,ncell_ok
           active_mg(myid,ilevel)%f(icell_mg_ok(i),1)=scan_flag(i)
        end do
        
     end do
     ! End loop over cells
  end do
!$OMP END PARALLEL 
  ! End loop over grids

end subroutine set_scan_flag_coarse
