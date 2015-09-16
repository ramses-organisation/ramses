!************************************************************************!
! This file contains essential modules for the GPU-ported Gravity kernel !
!************************************************************************!

! ########################################################################
! ########################################################################

module poisson_commons_acc
  use amr_parameters
  implicit none
  ! ---------------------------------------------------------------------
  ! Replacement of derived data types 
  ! ---------------------------------------------------------------------
  integer                                 :: lvl_ind, component, cpu_ind
  integer                                 :: active_ind, emission_ind, tmp_ind
  integer, dimension(:,:),allocatable     :: ngrid_a_mg, ngrid_e_mg

  real(kind=8), allocatable, dimension(:) :: active_mg_u
  real(kind=8), allocatable, dimension(:) :: emission_mg_u
  integer,      allocatable, dimension(:) :: active_mg_f

  integer,      allocatable, dimension(:) :: active_mg_igrid
  integer,      allocatable, dimension(:) :: emission_mg_igrid

  ! Constants (avoid copies from CPU to GPU)
  integer, dimension(1:3,1:2,1:8)         :: iii, jjj
  integer, dimension(1:3,1:4)             :: ired, iblack
  real(dp)                                :: a, b, c, d
  real(dp), dimension(1:8)                :: bbb
  integer,  dimension(1:8,1:8)            :: ccc
  integer,dimension(1:8)                  ::iii_aux,jjj_aux,kkk_aux
  integer,dimension(1:27,1:8,1:3)         ::lll,mmm
end module poisson_commons_acc

! ########################################################################
! ########################################################################

module ffv ! f(f)orce_f(f)ine_v(v)ersion
  use amr_parameters
  implicit none
  integer,dimension(:),allocatable       ::ngrid_aux
  integer,dimension(:,:),allocatable     ::ind_grid_aux
  integer ,dimension(:,:),allocatable    ::ind_cell_aux

  integer ,dimension(:,:,:),allocatable  ::ind_left,ind_right
  integer ,dimension(:,:,:),allocatable  ::igridn
  real(dp),dimension(:,:),allocatable    ::phi1,phi2,phi3,phi4
  real(dp),dimension(:,:,:,:),allocatable::phi_left,phi_right

  integer ,dimension(:,:,:,:),allocatable::nbors_father_grids
  integer ,dimension(:,:,:,:),allocatable::nbors_father_cells

  integer,dimension(:,:,:),allocatable   ::pos,ind_grid_father,ind_grid_ok
  integer,dimension(:,:,:,:),allocatable ::nbors_father_ok
  integer,dimension(:,:,:,:),allocatable ::nbors_grids_ok

  integer,dimension(:,:,:),allocatable   ::ind_grid1,ind_grid2,ind_grid3
  integer,dimension(:,:,:,:),allocatable ::nbors_grids

  integer,dimension(1:27,1:8,1:3)        ::lll,mmm
  integer,dimension(1:8)                 ::iii_aux,jjj_aux,kkk_aux
  integer,dimension(1:3,1:4,1:8)         ::ggg,hhh
  real(dp)                               ::aa,bb,cc,dd
  integer,dimension(1:8,1:8)             ::ccc
  real(dp),dimension(1:8)                ::bbbb

end module ffv

! ########################################################################
! ########################################################################

module iacv !i(i)nterpolate_a(a)nd_c(c)orrect_v(v)ersion
  use amr_parameters
  implicit none
  real(dp), dimension(:,:),allocatable ::corr

  integer,dimension(:),allocatable     ::nbatch
  integer ,dimension(:,:),allocatable  ::igrid_f_amr,icell_amr,cpu_amr_aux
  integer ,dimension(:,:,:),allocatable::nbors_father_cells,nbors_father_grids

  integer ,dimension(:,:),allocatable  ::pos,ind_grid_father,ind_grid_ok
  integer ,dimension(:,:,:),allocatable::nbors_father_ok,nbors_grids_ok

  integer ,dimension(:,:),allocatable  ::ind_grid1,ind_grid2,ind_grid3
  integer ,dimension(:,:,:),allocatable::nbors_grids
end module iacv

! ########################################################################
! ########################################################################

module mfbrv ! m(m)ake_f(f)ine_b(b)c_r(r)hs_v(v)ersion
  use amr_parameters
  implicit none
  ! Arrays for vectorized interpol_phi
  real(dp), dimension(:,:,:,:),allocatable :: phi_int
  integer,  dimension(:,:,:)  ,allocatable :: ind_cell

  integer ,dimension(:,:,:,:),allocatable  ::nbors_father_grids
  integer ,dimension(:,:,:,:),allocatable  ::nbors_father_cells

  integer,dimension(:,:,:)  ,allocatable   ::pos,ind_grid_father,ind_grid_ok
  integer,dimension(:,:,:,:),allocatable   ::nbors_father_ok
  integer,dimension(:,:,:,:),allocatable   ::nbors_grids_ok

  integer,dimension(:,:,:)  ,allocatable   ::ind_grid1,ind_grid2,ind_grid3
  integer,dimension(:,:,:,:),allocatable   ::nbors_grids

  integer,dimension(1:8)                   ::iii_aux=(/1,2,1,2,1,2,1,2/)
  integer,dimension(1:8)                   ::jjj_aux=(/3,3,4,4,3,3,4,4/)
  integer,dimension(1:8)                   ::kkk_aux=(/5,5,5,5,6,6,6,6/)
  integer,dimension(1:27,1:8,1:3)          ::lll,mmm
  real(dp)                                 ::aa,bb,cc,dd
  integer,dimension(1:8,1:8)               ::ccc
  real(dp),dimension(1:8)                  ::bbbb
end module mfbrv

! ########################################################################
! ########################################################################

! Index accessing
module ind_accessing
implicit none
contains

! for active_mg_u
function ind_acc_a_mg_u(lvl,cpu,comp)
   use poisson_commons_acc
   implicit none

   integer:: ind_acc_a_mg_u, lvl, cpu, comp

   ind_acc_a_mg_u = (sum(ngrid_a_mg(:,1:lvl))-sum(ngrid_a_mg(:,lvl)))*twotondim*4+&
                       sum(ngrid_a_mg(:,lvl))*twotondim*(comp-1)+&
                       (sum(ngrid_a_mg(1:cpu,lvl))-ngrid_a_mg(cpu,lvl))*twotondim

end function ind_acc_a_mg_u

! for active_mg_f
function ind_acc_a_mg_f(lvl,cpu)
   use poisson_commons_acc
   implicit none

   integer:: ind_acc_a_mg_f, lvl, cpu

   ind_acc_a_mg_f = (sum(ngrid_a_mg(:,1:lvl))-sum(ngrid_a_mg(:,lvl)))*twotondim+&
                       (sum(ngrid_a_mg(1:cpu,lvl))-ngrid_a_mg(cpu,lvl))*twotondim

end function ind_acc_a_mg_f

! for emission_mg
function ind_acc_e_mg_u(lvl,cpu,comp)
   use poisson_commons_acc
   implicit none

   integer:: ind_acc_e_mg_u, lvl, cpu, comp

   ind_acc_e_mg_u = (sum(ngrid_e_mg(:,1:lvl))-sum(ngrid_e_mg(:,lvl)))*twotondim*4+&
                       sum(ngrid_e_mg(:,lvl))*twotondim*(comp-1)+&
                       (sum(ngrid_e_mg(1:cpu,lvl))-ngrid_e_mg(cpu,lvl))*twotondim

end function ind_acc_e_mg_u

! for active_mg_igrid
function ind_acc_a_mg_igrid(lvl,cpu)
   use poisson_commons_acc
   implicit none

   integer:: ind_acc_a_mg_igrid, lvl, cpu

   ind_acc_a_mg_igrid = (sum(ngrid_a_mg(:,1:lvl))-sum(ngrid_a_mg(:,lvl)))+&
                           (sum(ngrid_a_mg(1:cpu,lvl))-ngrid_a_mg(cpu,lvl))

end function ind_acc_a_mg_igrid

! for emission_mg_igrid
function ind_acc_e_mg_igrid(lvl,cpu)
   use poisson_commons_acc
   implicit none

   integer:: ind_acc_e_mg_igrid, lvl, cpu

   ind_acc_e_mg_igrid = (sum(ngrid_e_mg(:,1:lvl))-sum(ngrid_e_mg(:,lvl)))+&
                           (sum(ngrid_e_mg(1:cpu,lvl))-ngrid_e_mg(cpu,lvl))

end function ind_acc_e_mg_igrid
end module ind_accessing

! ########################################################################
! ########################################################################

module iacv_subroutines
contains
subroutine get3cubefather_iacv(ilevel,nbin,ibin,ix,iy,iz,iix,iiy,iiz)
  use amr_commons 
  use iacv
  use poisson_commons_acc
  implicit none
  integer::ilevel,nbin,ibin
  integer ,dimension(1:nbin,1:nvector)::ix,iy,iz,iix,iiy,iiz
  !------------------------------------------------------------------
  ! This subroutine determines the 3^ndim neighboring father cells 
  ! of the input father cell. According to the refinement rule, 
  ! they should be present anytime.
  !------------------------------------------------------------------
  integer::i,j,nxny,i1,j1,k1,ind,iok
  integer::i1min,i1max,j1min,j1max,k1min,k1max,ind_father
  logical::oups

  nxny=nx*ny

  if(ilevel==1)then  ! Easy...

     oups=.false.
     do i=1,nbatch(ibin)
        if(icell_amr(ibin,i)>ncoarse)oups=.true.
     end do
!     if(oups)then
!        write(*,*)'get3cubefather'
!        write(*,*)'oupsssss !'
!        call clean_stop
!     endif

     do i=1,nbatch(ibin)
        iz(ibin,i)=(icell_amr(ibin,i)-1)/nxny
     end do
     do i=1,nbatch(ibin)
        iy(ibin,i)=(icell_amr(ibin,i)-1-iz(ibin,i)*nxny)/nx
     end do
     do i=1,nbatch(ibin)
        ix(ibin,i)=(icell_amr(ibin,i)-1-iy(ibin,i)*nx-iz(ibin,i)*nxny)
     end do

     i1min=0; i1max=0
     if(ndim > 0)i1max=2
     j1min=0; j1max=0
     if(ndim > 1)j1max=2
     k1min=0; k1max=0
     if(ndim > 2)k1max=2

     ! Loop over 3^ndim neighboring father cells
     do k1=k1min,k1max
        iiz=iz
        if(ndim > 2)then
           do i=1,nbatch(ibin)
              iiz(ibin,i)=iz(ibin,i)+k1-1
              if(iiz(ibin,i) < 0   )iiz(ibin,i)=nz-1
              if(iiz(ibin,i) > nz-1)iiz(ibin,i)=0
           end do
        end if
        do j1=j1min,j1max
           iiy=iy
           if(ndim > 1)then
              do i=1,nbatch(ibin)
                 iiy(ibin,i)=iy(ibin,i)+j1-1
                 if(iiy(ibin,i) < 0   )iiy(ibin,i)=ny-1
                 if(iiy(ibin,i) > ny-1)iiy(ibin,i)=0
              end do
           end if
           do i1=i1min,i1max
              iix=ix
              if(ndim > 0)then
                 do i=1,nbatch(ibin)
                    iix(ibin,i)=ix(ibin,i)+i1-1
                    if(iix(ibin,i) < 0   )iix(ibin,i)=nx-1
                    if(iix(ibin,i) > nx-1)iix(ibin,i)=0
                 end do
              end if
              ind_father=1+i1+3*j1+9*k1
              do i=1,nbatch(ibin)
                 nbors_father_cells(ibin,i,ind_father)=1 &
                      & +iix(ibin,i) &
                      & +iiy(ibin,i)*nx &
                      & +iiz(ibin,i)*nxny
              end do
           end do
        end do
     end do

     i1min=0; i1max=0
     if(ndim > 0)i1max=1
     j1min=0; j1max=0
     if(ndim > 1)j1max=1
     k1min=0; k1max=0
     if(ndim > 2)k1max=1

     ! Loop over 2^ndim neighboring father grids
     do k1=k1min,k1max
        iiz=iz
        if(ndim > 2)then
           do i=1,nbatch(ibin)
              iiz(ibin,i)=iz(ibin,i)+2*k1-1
              if(iiz(ibin,i) < 0   )iiz(ibin,i)=nz-1
              if(iiz(ibin,i) > nz-1)iiz(ibin,i)=0
           end do
        end if
        do j1=j1min,j1max
           iiy=iy
           if(ndim > 1)then
              do i=1,nbatch(ibin)
                 iiy(ibin,i)=iy(ibin,i)+2*j1-1
                 if(iiy(ibin,i) < 0   )iiy(ibin,i)=ny-1
                 if(iiy(ibin,i) > ny-1)iiy(ibin,i)=0
              end do
           end if
           do i1=i1min,i1max
              iix=ix
              if(ndim > 0)then
                 do i=1,nbatch(ibin)
                    iix(ibin,i)=ix(ibin,i)+2*i1-1
                    if(iix(ibin,i) < 0   )iix(ibin,i)=nx-1
                    if(iix(ibin,i) > nx-1)iix(ibin,i)=0
                 end do
              end if
              ind_father=1+i1+2*j1+4*k1
              do i=1,nbatch(ibin)
                 nbors_father_grids(ibin,i,ind_father)=1 &
                      & +(iix(ibin,i)/2) &
                      & +(iiy(ibin,i)/2)*(nx/2) &
                      & +(iiz(ibin,i)/2)*(nxny/4)
              end do
           end do
        end do
     end do

  else    ! else, more complicated...
     
     ! Get father cell position in the grid
     do i=1,nbatch(ibin)
        pos(ibin,i)=(icell_amr(ibin,i)-ncoarse-1)/ngridmax+1
     end do
     ! Get father grid
     do i=1,nbatch(ibin)
        ind_grid_father(ibin,i)=icell_amr(ibin,i)-ncoarse-(pos(ibin,i)-1)*ngridmax
     end do

     ! Loop over position
     do ind=1,twotondim

        ! Select father cells that sit at position ind
        iok=0
        do i=1,nbatch(ibin)
           if(pos(ibin,i)==ind)then
              iok=iok+1
              ind_grid_ok(ibin,iok)=ind_grid_father(ibin,i)
           end if
        end do

        if(iok>0)&
        & call get3cubepos_iacv(ind,iok,nbin,ibin)

        ! Store neighboring father cells for selected cells
        do j=1,threetondim
           iok=0
           do i=1,nbatch(ibin)
              if(pos(ibin,i)==ind)then
                 iok=iok+1
                 nbors_father_cells(ibin,i,j)=nbors_father_ok(ibin,iok,j)
              end if
           end do
        end do

        ! Store neighboring father grids for selected cells
        do j=1,twotondim
           iok=0
           do i=1,nbatch(ibin)
              if(pos(ibin,i)==ind)then
                 iok=iok+1
                 nbors_father_grids(ibin,i,j)=nbors_grids_ok(ibin,iok,j)
              end if
           end do
        end do

     end do

  end if

end subroutine get3cubefather_iacv
!##############################################################
!##############################################################
subroutine get3cubepos_iacv(ind,ng,nbin,ibin)
  use amr_commons
  use iacv
  use poisson_commons_acc
  implicit none
  integer::ng,ind,nbin,ibin
  !--------------------------------------------------------------------
  ! This subroutine determines the 3^ndim neighboring father cells 
  ! of the input cell at position ind in grid ind_grid. According to 
  ! the refinements rules and since the input cell is refined, 
  ! they should be present anytime.
  !--------------------------------------------------------------------
  integer::i,j,iskip
  integer::ii,iimin,iimax
  integer::jj,jjmin,jjmax
  integer::kk,kkmin,kkmax
  integer::icell,igrid,inbor

  iimin=0; iimax=0
  if(ndim>0)iimax=1
  jjmin=0; jjmax=0
  if(ndim>1)jjmax=1
  kkmin=0; kkmax=0
  if(ndim>2)kkmax=1

  do kk=kkmin,kkmax
     do i=1,ng
        ind_grid1(ibin,i)=ind_grid_ok(ibin,i)
     end do
     if(kk>0)then
        inbor=kkk_aux(ind)
        do i=1,ng
           ind_grid1(ibin,i)=son(nbor(ind_grid_ok(ibin,i),inbor))
        end do
     end if

     do jj=jjmin,jjmax
        do i=1,ng
           ind_grid2(ibin,i)=ind_grid1(ibin,i)
        end do
        if(jj>0)then
           inbor=jjj_aux(ind)
           do i=1,ng
              ind_grid2(ibin,i)=son(nbor(ind_grid1(ibin,i),inbor))
           end do
        end if
 
        do ii=iimin,iimax
           do i=1,ng
              ind_grid3(ibin,i)=ind_grid2(ibin,i)
           end do
           if(ii>0)then
              inbor=iii_aux(ind)
              do i=1,ng
                 ind_grid3(ibin,i)=son(nbor(ind_grid2(ibin,i),inbor))
              end do
           end if

           inbor=1+ii+2*jj+4*kk
           do i=1,ng
              nbors_grids(ibin,i,inbor)=ind_grid3(ibin,i)
           end do     

        end do
     end do
  end do     

  do j=1,twotondim
     do i=1,ng
        nbors_grids_ok(ibin,i,j)=nbors_grids(ibin,i,j)
     end do
  end do

  do j=1,threetondim
     igrid=lll(j,ind,ndim)
     icell=mmm(j,ind,ndim)
     iskip=ncoarse+(icell-1)*ngridmax
     do i=1,ng
        nbors_father_ok(ibin,i,j)=iskip+nbors_grids(ibin,i,igrid)
     end do
  end do

end subroutine get3cubepos_iacv
end module iacv_subroutines
!###########################################################
!###########################################################

! f(f)orce_f(f)ine_v(v)ersrion
module ffv_subroutines
contains
subroutine interpol_phi_ffv(ind_cell,phi_int,ilevel,icount,ibin,nbin,idim,ix_aux,iy_aux,iz_aux,iix_aux,iiy_aux,iiz_aux)
  use amr_commons
  use poisson_commons, only:phi,phi_old
  use ffv
  implicit none
  integer::ilevel,icount,nbin,idim,ibin
  integer ,dimension(1:nbin,1:nvector,1:ndim)::ind_cell
  real(dp),dimension(1:nbin,1:nvector,1:twotondim,1:ndim)::phi_int
  integer,dimension(1:nbin,1:ndim,1:nvector)::ix_aux,iy_aux,iz_aux,iix_aux,iiy_aux,iiz_aux
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Routine for interpolation at level-boundaries. Interpolation is used for
  ! - boundary conditions for solving poisson equation at fine level
  ! - computing force (gradient_phi) at fine level for cells close to boundary
  ! Interpolation is performed in space (CIC) and - if adaptive timestepping is on -
  ! time (linear extrapolation of the change in phi during the last coarse step 
  ! onto the first fine step)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer::i,ind,indice,ind_average,ind_father
  real(dp)::dx,tfrac
  real(dp)::coeff,add

!  if (icount .ne. 1 .and. icount .ne. 2)then
!     write(*,*)'icount has bad value'
!     call clean_stop
!  endif

  !compute fraction of timesteps for interpolation
  if (dtold(ilevel-1)> 0)then
     !tfrac=0.
     tfrac=1.0*dtnew(ilevel)/dtold(ilevel-1)*(icount-1)
  else
     tfrac=0.
  end if

  ! Mesh siz_auxe at level ilevel
  dx=0.5D0**ilevel

  ! The subroutine and its auxiliaries are in this file (look at the bottom)
  call get3cubefather_ffv(ind_cell,ilevel,ibin,nbin,idim,ix_aux,iy_aux,iz_aux,iix_aux,iiy_aux,iiz_aux)

  ! Third order phi interpolation
  do ind=1,twotondim
     do i=1,ngrid_aux(ibin)
        phi_int(ibin,i,ind,idim)=0d0
     end do
     do ind_average=1,twotondim
        ind_father=ccc(ind_average,ind)
        coeff=bbbb(ind_average)
        do i=1,ngrid_aux(ibin)
           indice=nbors_father_cells(ibin,idim,i,ind_father)
           if (indice==0) then 
!              write(*,*)'no all neighbors present in interpol_phi...'
              add=coeff*(phi(ind_cell(ibin,i,idim))+(phi(ind_cell(ibin,i,idim))-phi_old(ind_cell(ibin,i,idim)))*tfrac)
           else
              add=coeff*(phi(indice)+(phi(indice)-phi_old(indice))*tfrac)
           endif
           phi_int(ibin,i,ind,idim)=phi_int(ibin,i,ind,idim)+add
        end do
     end do
  end do

 end subroutine interpol_phi_ffv
!##############################################################
!##############################################################
subroutine get3cubefather_ffv(ind_cell_father,ilevel,ibin,nbin,idim,ix_aux,iy_aux,iz_aux,iix_aux,iiy_aux,iiz_aux)
  use amr_commons
  use ffv
  implicit none
  integer::ilevel,nbin,idim,ibin
  integer,dimension(1:nbin,1:nvector,1:ndim)::ind_cell_father
  integer,dimension(1:nbin,1:ndim,1:nvector)::ix_aux,iy_aux,iz_aux,iix_aux,iiy_aux,iiz_aux
  !------------------------------------------------------------------
  ! This subroutine determines the 3^ndim neighboring father cells 
  ! of the input father cell. According to the refinement rule, 
  ! they should be present anytime.
  !------------------------------------------------------------------
  integer::i,j,nxny,i1,j1,k1,ind,iok
  integer::i1min,i1max,j1min,j1max,k1min,k1max,ind_father
  logical::oups

  nxny=nx*ny

  if(ilevel==1)then  ! Easy...

     oups=.false.
     do i=1,ngrid_aux(ibin)
        if(ind_cell_father(ibin,i,idim)>ncoarse) oups=.true.
     end do
!     if(oups)then
!        write(*,*)'get3cubefather'
!        write(*,*)'oupsssss !'
!        call clean_stop
!     endif

     do i=1,ngrid_aux(ibin)
        iz_aux(ibin,idim,i)=(ind_cell_father(ibin,i,idim)-1)/nxny
     end do
     do i=1,ngrid_aux(ibin)
        iy_aux(ibin,idim,i)=(ind_cell_father(ibin,i,idim)-1-iz_aux(ibin,idim,i)*nxny)/nx
     end do
     do i=1,ngrid_aux(ibin)
        ix_aux(ibin,idim,i)=(ind_cell_father(ibin,i,idim)-1-iy_aux(ibin,idim,i)*nx-iz_aux(ibin,idim,i)*nxny)
     end do

     i1min=0; i1max=0
     if(ndim > 0)i1max=2
     j1min=0; j1max=0
     if(ndim > 1)j1max=2
     k1min=0; k1max=0
     if(ndim > 2)k1max=2

     ! Loop over 3^ndim neighboring father cells
     do k1=k1min,k1max
        iiz_aux=iz_aux
        if(ndim > 2)then
           do i=1,ngrid_aux(ibin)
              iiz_aux(ibin,idim,i)=iz_aux(ibin,idim,i)+k1-1
              if(iiz_aux(ibin,idim,i) < 0   )iiz_aux(ibin,idim,i)=nz-1
              if(iiz_aux(ibin,idim,i) > nz-1)iiz_aux(ibin,idim,i)=0
           end do
        end if
        do j1=j1min,j1max
           iiy_aux=iy_aux
           if(ndim > 1)then
              do i=1,ngrid_aux(ibin)
                 iiy_aux(ibin,idim,i)=iy_aux(ibin,idim,i)+j1-1
                 if(iiy_aux(ibin,idim,i) < 0   )iiy_aux(ibin,idim,i)=ny-1
                 if(iiy_aux(ibin,idim,i) > ny-1)iiy_aux(ibin,idim,i)=0
              end do
           end if
           do i1=i1min,i1max
              iix_aux=ix_aux
              if(ndim > 0)then
                 do i=1,ngrid_aux(ibin)
                    iix_aux(ibin,idim,i)=ix_aux(ibin,idim,i)+i1-1
                    if(iix_aux(ibin,idim,i) < 0   )iix_aux(ibin,idim,i)=nx-1
                    if(iix_aux(ibin,idim,i) > nx-1)iix_aux(ibin,idim,i)=0
                 end do
              end if
              ind_father=1+i1+3*j1+9*k1
              do i=1,ngrid_aux(ibin)
                 nbors_father_cells(ibin,idim,i,ind_father)=1 &
                      & +iix_aux(ibin,idim,i) &
                      & +iiy_aux(ibin,idim,i)*nx &
                      & +iiz_aux(ibin,idim,i)*nxny
              end do
           end do
        end do
     end do

     i1min=0; i1max=0
     if(ndim > 0)i1max=1
     j1min=0; j1max=0
     if(ndim > 1)j1max=1
     k1min=0; k1max=0
     if(ndim > 2)k1max=1

     ! Loop over 2^ndim neighboring father grids
     do k1=k1min,k1max
        iiz_aux=iz_aux
        if(ndim > 2)then
           do i=1,ngrid_aux(ibin)
              iiz_aux(ibin,idim,i)=iz_aux(ibin,idim,i)+2*k1-1
              if(iiz_aux(ibin,idim,i) < 0   )iiz_aux(ibin,idim,i)=nz-1
              if(iiz_aux(ibin,idim,i) > nz-1)iiz_aux(ibin,idim,i)=0
           end do
        end if
        do j1=j1min,j1max
           iiy_aux=iy_aux
           if(ndim > 1)then
              do i=1,ngrid_aux(ibin)
                 iiy_aux(ibin,idim,i)=iy_aux(ibin,idim,i)+2*j1-1
                 if(iiy_aux(ibin,idim,i) < 0   )iiy_aux(ibin,idim,i)=ny-1
                 if(iiy_aux(ibin,idim,i) > ny-1)iiy_aux(ibin,idim,i)=0
              end do
           end if
           do i1=i1min,i1max
              iix_aux=ix_aux
              if(ndim > 0)then
                 do i=1,ngrid_aux(ibin)
                    iix_aux(ibin,idim,i)=ix_aux(ibin,idim,i)+2*i1-1
                    if(iix_aux(ibin,idim,i) < 0   )iix_aux(ibin,idim,i)=nx-1
                    if(iix_aux(ibin,idim,i) > nx-1)iix_aux(ibin,idim,i)=0
                 end do
              end if
              ind_father=1+i1+2*j1+4*k1
              do i=1,ngrid_aux(ibin)
                 nbors_father_grids(ibin,idim,i,ind_father)=1 &
                      & +(iix_aux(ibin,idim,i)/2) &
                      & +(iiy_aux(ibin,idim,i)/2)*(nx/2) &
                      & +(iiz_aux(ibin,idim,i)/2)*(nxny/4)
              end do
           end do
        end do
     end do

  else    ! else, more complicated...
     
     ! Get father cell position in the grid
     do i=1,ngrid_aux(ibin)
        pos(ibin,idim,i)=(ind_cell_father(ibin,i,idim)-ncoarse-1)/ngridmax+1
     end do
     ! Get father grid
     do i=1,ngrid_aux(ibin)
        ind_grid_father(ibin,idim,i)=ind_cell_father(ibin,i,idim)-ncoarse-(pos(ibin,idim,i)-1)*ngridmax
     end do

     ! Loop over position
     do ind=1,twotondim

        ! Select father cells that sit at position ind
        iok=0
        do i=1,ngrid_aux(ibin)
           if(pos(ibin,idim,i)==ind)then
              iok=iok+1
              ind_grid_ok(ibin,idim,iok)=ind_grid_father(ibin,idim,i)
           end if
        end do

        if(iok>0)&
        & call get3cubepos_ffv(ind,iok,ibin,nbin,idim)

        ! Store neighboring father cells for selected cells
        do j=1,threetondim
           iok=0
           do i=1,ngrid_aux(ibin)
              if(pos(ibin,idim,i)==ind)then
                 iok=iok+1
                 nbors_father_cells(ibin,idim,i,j)=nbors_father_ok(ibin,idim,iok,j)
              end if
           end do
        end do

        ! Store neighboring father grids for selected cells
        do j=1,twotondim
           iok=0
           do i=1,ngrid_aux(ibin)
              if(pos(ibin,idim,i)==ind)then
                 iok=iok+1
                 nbors_father_grids(ibin,idim,i,j)=nbors_grids_ok(ibin,idim,iok,j)
              end if
           end do
        end do

     end do

  end if

end subroutine get3cubefather_ffv
!##############################################################
!##############################################################
subroutine get3cubepos_ffv(ind,iok,ibin,nbin,idim)
  use amr_commons
  use ffv
  implicit none
  integer::iok,ind,nbin,idim,ibin
  !--------------------------------------------------------------------
  ! This subroutine determines the 3^ndim neighboring father cells 
  ! of the input cell at position ind in grid ind_grid. According to 
  ! the refinements rules and since the input cell is refined, 
  ! they should be present anytime.
  !--------------------------------------------------------------------
  integer::i,j,iskip
  integer::ii,iimin,iimax
  integer::jj,jjmin,jjmax
  integer::kk,kkmin,kkmax
  integer::icell,igrid_get,inbor

  iimin=0; iimax=0
  if(ndim>0)iimax=1
  jjmin=0; jjmax=0
  if(ndim>1)jjmax=1
  kkmin=0; kkmax=0
  if(ndim>2)kkmax=1

  do kk=kkmin,kkmax
     do i=1,iok
        ind_grid1(ibin,idim,i)=ind_grid_ok(ibin,idim,i)
     end do
     if(kk>0)then
        inbor=kkk_aux(ind)
        do i=1,iok
           ind_grid1(ibin,idim,i)=son(nbor(ind_grid_ok(ibin,idim,i),inbor))
        end do
     end if

     do jj=jjmin,jjmax
        do i=1,iok
           ind_grid2(ibin,idim,i)=ind_grid1(ibin,idim,i)
        end do
        if(jj>0)then
           inbor=jjj_aux(ind)
           do i=1,iok
              ind_grid2(ibin,idim,i)=son(nbor(ind_grid1(ibin,idim,i),inbor))
           end do
        end if
 
        do ii=iimin,iimax
           do i=1,iok
              ind_grid3(ibin,idim,i)=ind_grid2(ibin,idim,i)
           end do
           if(ii>0)then
              inbor=iii_aux(ind)
              do i=1,iok
                 ind_grid3(ibin,idim,i)=son(nbor(ind_grid2(ibin,idim,i),inbor))
              end do
           end if

           inbor=1+ii+2*jj+4*kk
           do i=1,iok
              nbors_grids(ibin,idim,i,inbor)=ind_grid3(ibin,idim,i)
           end do     

        end do
     end do
  end do     

  do j=1,twotondim
     do i=1,iok
        nbors_grids_ok(ibin,idim,i,j)=nbors_grids(ibin,idim,i,j)
     end do
  end do

  do j=1,threetondim
     igrid_get=lll(j,ind,ndim)
     icell=mmm(j,ind,ndim)
     iskip=ncoarse+(icell-1)*ngridmax
     do i=1,iok
        nbors_father_ok(ibin,idim,i,j)=iskip+nbors_grids(ibin,idim,i,igrid_get)
     end do
  end do

end subroutine get3cubepos_ffv
end module ffv_subroutines

!###########################################################
!###########################################################

! m(m)ake_f(f)ine_b(b)c_r(r)hs_v(v)ersion
module mfbrv_subroutines
contains
subroutine interpol_phi_mfbrv(ncell,ilevel,icount,ind_m,igrid_mg_m, &
                              ix,iy,iz,iix,iiy,iiz,ngrid)
  use amr_commons
  use poisson_commons, only:phi,phi_old
  use mfbrv
  implicit none
  integer::ncell,ilevel,icount,ind_m,igrid_mg_m,idim_m,ngrid
  integer,dimension(1:twotondim,1:ngrid,1:nvector)::ix,iy,iz,iix,iiy,iiz 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Routine for interpolation at level-boundaries. Interpolation is used for
  ! - boundary conditions for solving poisson equation at fine level
  ! - computing force (gradient_phi) at fine level for cells close to boundary
  ! Interpolation is performed in space (CIC) and - if adaptive timestepping is on -
  ! time (linear extrapolation of the change in phi during the last coarse step 
  ! onto the first fine step)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer::i,ind,indice,ind_average,ind_father
  real(dp)::dx,tfrac
  real(dp)::coeff,add

!  if (icount .ne. 1 .and. icount .ne. 2)then
!     write(*,*)'icount has bad value'
!     call clean_stop
!  endif

  !compute fraction of timesteps for interpolation
  if (dtold(ilevel-1)> 0)then
     !tfrac=0.
     tfrac=1.0*dtnew(ilevel)/dtold(ilevel-1)*(icount-1)
  else
     tfrac=0.
  end if

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel
  call get3cubefather_mfbrv(ncell,ilevel,ind_m,igrid_mg_m, &
                            ix,iy,iz,iix,iiy,iiz,ngrid)

  ! Third order phi interpolation
  do ind=1,twotondim
     do i=1,ncell
        phi_int(ind_m,igrid_mg_m,i,ind)=0d0
     end do
     do ind_average=1,twotondim
        ind_father=ccc(ind_average,ind)
        coeff=bbbb(ind_average)
        do i=1,ncell
           indice=nbors_father_cells(ind_m,igrid_mg_m,i,ind_father)
           if (indice==0) then 
!              write(*,*)'no all neighbors present in interpol_phi...'
              add=coeff*(phi(ind_cell(ind_m,igrid_mg_m,i))+(phi(ind_cell(ind_m,igrid_mg_m,i))-phi_old(ind_cell(ind_m,igrid_mg_m,i)))*tfrac)
              !               add=coeff*(-3d0/8d0*dx**2*boxlen*rho(ind_cell(i))+phi(ind_cell(i)))
           else
              add=coeff*(phi(indice)+(phi(indice)-phi_old(indice))*tfrac)
              !               add=coeff*(-3d0/8d0*dx**2*boxlen*rho(indice)+phi(indice)) 
           endif
           phi_int(ind_m,igrid_mg_m,i,ind)=phi_int(ind_m,igrid_mg_m,i,ind)+add
        end do
     end do
  end do

end subroutine interpol_phi_mfbrv
!###########################################################
!###########################################################
subroutine get3cubefather_mfbrv(ncell,ilevel,ind_m,igrid_mg_m, &
                            ix,iy,iz,iix,iiy,iiz,ngrid)
  use amr_commons 
  use mfbrv
  implicit none
  integer::ncell,ilevel,ind_m,igrid_mg_m,idim_m,ngrid
  integer,dimension(1:twotondim,1:ngrid,1:nvector)::ix,iy,iz,iix,iiy,iiz 
  !------------------------------------------------------------------
  ! This subroutine determines the 3^ndim neighboring father cells 
  ! of the input father cell. According to the refinement rule, 
  ! they should be present anytime.
  !------------------------------------------------------------------
  integer::i,j,nxny,i1,j1,k1,ind,iok
  integer::i1min,i1max,j1min,j1max,k1min,k1max,ind_father
  logical::oups

  nxny=nx*ny

  if(ilevel==1)then  ! Easy...

     oups=.false.
     do i=1,ncell
        if(ind_cell(ind_m,igrid_mg_m,i)>ncoarse)oups=.true.
     end do
!     if(oups)then
!        write(*,*)'get3cubefather'
!        write(*,*)'oupsssss !'
!        call clean_stop
!     endif

     do i=1,ncell
        iz(ind_m,igrid_mg_m,i)=(ind_cell(ind_m,igrid_mg_m,i)-1)/nxny
     end do
     do i=1,ncell
        iy(ind_m,igrid_mg_m,i)=(ind_cell(ind_m,igrid_mg_m,i)-1-iz(ind_m,igrid_mg_m,i)*nxny)/nx
     end do
     do i=1,ncell
        ix(ind_m,igrid_mg_m,i)=(ind_cell(ind_m,igrid_mg_m,i)-1-iy(ind_m,igrid_mg_m,i)*nx-iz(ind_m,igrid_mg_m,i)*nxny)
     end do

     i1min=0; i1max=0
     if(ndim > 0)i1max=2
     j1min=0; j1max=0
     if(ndim > 1)j1max=2
     k1min=0; k1max=0
     if(ndim > 2)k1max=2

     ! Loop over 3^ndim neighboring father cells
     do k1=k1min,k1max
        iiz=iz
        if(ndim > 2)then
           do i=1,ncell
              iiz(ind_m,igrid_mg_m,i)=iz(ind_m,igrid_mg_m,i)+k1-1
              if(iiz(ind_m,igrid_mg_m,i) < 0   )iiz(ind_m,igrid_mg_m,i)=nz-1
              if(iiz(ind_m,igrid_mg_m,i) > nz-1)iiz(ind_m,igrid_mg_m,i)=0
           end do
        end if
        do j1=j1min,j1max
           iiy=iy
           if(ndim > 1)then
              do i=1,ncell
                 iiy(ind_m,igrid_mg_m,i)=iy(ind_m,igrid_mg_m,i)+j1-1
                 if(iiy(ind_m,igrid_mg_m,i) < 0   )iiy(ind_m,igrid_mg_m,i)=ny-1
                 if(iiy(ind_m,igrid_mg_m,i) > ny-1)iiy(ind_m,igrid_mg_m,i)=0
              end do
           end if
           do i1=i1min,i1max
              iix=ix
              if(ndim > 0)then
                 do i=1,ncell
                    iix(ind_m,igrid_mg_m,i)=ix(ind_m,igrid_mg_m,i)+i1-1
                    if(iix(ind_m,igrid_mg_m,i) < 0   )iix(ind_m,igrid_mg_m,i)=nx-1
                    if(iix(ind_m,igrid_mg_m,i) > nx-1)iix(ind_m,igrid_mg_m,i)=0
                 end do
              end if
              ind_father=1+i1+3*j1+9*k1
              do i=1,ncell
                 nbors_father_cells(ind_m,igrid_mg_m,i,ind_father)=1 &
                      & +iix(ind_m,igrid_mg_m,i) &
                      & +iiy(ind_m,igrid_mg_m,i)*nx &
                      & +iiz(ind_m,igrid_mg_m,i)*nxny
              end do
           end do
        end do
     end do

     i1min=0; i1max=0
     if(ndim > 0)i1max=1
     j1min=0; j1max=0
     if(ndim > 1)j1max=1
     k1min=0; k1max=0
     if(ndim > 2)k1max=1

     ! Loop over 2^ndim neighboring father grids
     do k1=k1min,k1max
        iiz=iz
        if(ndim > 2)then
           do i=1,ncell
              iiz(ind_m,igrid_mg_m,i)=iz(ind_m,igrid_mg_m,i)+2*k1-1
              if(iiz(ind_m,igrid_mg_m,i) < 0   )iiz(ind_m,igrid_mg_m,i)=nz-1
              if(iiz(ind_m,igrid_mg_m,i) > nz-1)iiz(ind_m,igrid_mg_m,i)=0
           end do
        end if
        do j1=j1min,j1max
           iiy=iy
           if(ndim > 1)then
              do i=1,ncell
                 iiy(ind_m,igrid_mg_m,i)=iy(ind_m,igrid_mg_m,i)+2*j1-1
                 if(iiy(ind_m,igrid_mg_m,i) < 0   )iiy(ind_m,igrid_mg_m,i)=ny-1
                 if(iiy(ind_m,igrid_mg_m,i) > ny-1)iiy(ind_m,igrid_mg_m,i)=0
              end do
           end if
           do i1=i1min,i1max
              iix=ix
              if(ndim > 0)then
                 do i=1,ncell
                    iix(ind_m,igrid_mg_m,i)=ix(ind_m,igrid_mg_m,i)+2*i1-1
                    if(iix(ind_m,igrid_mg_m,i) < 0   )iix(ind_m,igrid_mg_m,i)=nx-1
                    if(iix(ind_m,igrid_mg_m,i) > nx-1)iix(ind_m,igrid_mg_m,i)=0
                 end do
              end if
              ind_father=1+i1+2*j1+4*k1
              do i=1,ncell
                 nbors_father_grids(ind_m,igrid_mg_m,i,ind_father)=1 &
                      & +(iix(ind_m,igrid_mg_m,i)/2) &
                      & +(iiy(ind_m,igrid_mg_m,i)/2)*(nx/2) &
                      & +(iiz(ind_m,igrid_mg_m,i)/2)*(nxny/4)
              end do
           end do
        end do
     end do

  else    ! else, more complicated...
     
     ! Get father cell position in the grid
     do i=1,ncell
        pos(ind_m,igrid_mg_m,i)=(ind_cell(ind_m,igrid_mg_m,i)-ncoarse-1)/ngridmax+1
     end do
     ! Get father grid
     do i=1,ncell
        ind_grid_father(ind_m,igrid_mg_m,i)=ind_cell(ind_m,igrid_mg_m,i)-ncoarse-(pos(ind_m,igrid_mg_m,i)-1)*ngridmax
     end do

     ! Loop over position
     do ind=1,twotondim

        ! Select father cells that sit at position ind
        iok=0
        do i=1,ncell
           if(pos(ind_m,igrid_mg_m,i)==ind)then
              iok=iok+1
              ind_grid_ok(ind_m,igrid_mg_m,iok)=ind_grid_father(ind_m,igrid_mg_m,i)
           end if
        end do

        if(iok>0)&
        & call get3cubepos_mfbrv(ind,iok,ind_m,igrid_mg_m)

        ! Store neighboring father cells for selected cells
        do j=1,threetondim
           iok=0
           do i=1,ncell
              if(pos(ind_m,igrid_mg_m,i)==ind)then
                 iok=iok+1
                 nbors_father_cells(ind_m,igrid_mg_m,i,j)=nbors_father_ok(ind_m,igrid_mg_m,iok,j)
              end if
           end do
        end do

        ! Store neighboring father grids for selected cells
        do j=1,twotondim
           iok=0
           do i=1,ncell
              if(pos(ind_m,igrid_mg_m,i)==ind)then
                 iok=iok+1
                 nbors_father_grids(ind_m,igrid_mg_m,i,j)=nbors_grids_ok(ind_m,igrid_mg_m,iok,j)
              end if
           end do
        end do

     end do

  end if

end subroutine get3cubefather_mfbrv
!##############################################################
!##############################################################
subroutine get3cubepos_mfbrv(ind,ng,ind_m,igrid_mg_m)
  use amr_commons
  use mfbrv
  implicit none
  integer::ng,ind,ind_m,igrid_mg_m
  !--------------------------------------------------------------------
  ! This subroutine determines the 3^ndim neighboring father cells 
  ! of the input cell at position ind in grid ind_grid. According to 
  ! the refinements rules and since the input cell is refined, 
  ! they should be present anytime.
  !--------------------------------------------------------------------
  integer::i,j,iskip
  integer::ii,iimin,iimax
  integer::jj,jjmin,jjmax
  integer::kk,kkmin,kkmax
  integer::icell,igrid,inbor

  iimin=0; iimax=0
  if(ndim>0)iimax=1
  jjmin=0; jjmax=0
  if(ndim>1)jjmax=1
  kkmin=0; kkmax=0
  if(ndim>2)kkmax=1

  do kk=kkmin,kkmax
     do i=1,ng
        ind_grid1(ind_m,igrid_mg_m,i)=ind_grid_ok(ind_m,igrid_mg_m,i)
     end do
     if(kk>0)then
        inbor=kkk_aux(ind)
        do i=1,ng
           ind_grid1(ind_m,igrid_mg_m,i)=son(nbor(ind_grid_ok(ind_m,igrid_mg_m,i),inbor))
        end do
     end if

     do jj=jjmin,jjmax
        do i=1,ng
           ind_grid2(ind_m,igrid_mg_m,i)=ind_grid1(ind_m,igrid_mg_m,i)
        end do
        if(jj>0)then
           inbor=jjj_aux(ind)
           do i=1,ng
              ind_grid2(ind_m,igrid_mg_m,i)=son(nbor(ind_grid1(ind_m,igrid_mg_m,i),inbor))
           end do
        end if
 
        do ii=iimin,iimax
           do i=1,ng
              ind_grid3(ind_m,igrid_mg_m,i)=ind_grid2(ind_m,igrid_mg_m,i)
           end do
           if(ii>0)then
              inbor=iii_aux(ind)
              do i=1,ng
                 ind_grid3(ind_m,igrid_mg_m,i)=son(nbor(ind_grid2(ind_m,igrid_mg_m,i),inbor))
              end do
           end if

           inbor=1+ii+2*jj+4*kk
           do i=1,ng
              nbors_grids(ind_m,igrid_mg_m,i,inbor)=ind_grid3(ind_m,igrid_mg_m,i)
           end do     

        end do
     end do
  end do     

  do j=1,twotondim
     do i=1,ng
        nbors_grids_ok(ind_m,igrid_mg_m,i,j)=nbors_grids(ind_m,igrid_mg_m,i,j)
     end do
  end do

  do j=1,threetondim
     igrid=lll(j,ind,ndim)
     icell=mmm(j,ind,ndim)
     iskip=ncoarse+(icell-1)*ngridmax
     do i=1,ng
        nbors_father_ok(ind_m,igrid_mg_m,i,j)=iskip+nbors_grids(ind_m,igrid_mg_m,i,igrid)
     end do
  end do

end subroutine get3cubepos_mfbrv
end module mfbrv_subroutines
