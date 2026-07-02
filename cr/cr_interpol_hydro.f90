! CR AMR averaging for cruold(:,1:ncrvar): restriction (cr_upload_fine/cr_upl)
! and prolongation (cr_interpol_hydro).
SUBROUTINE cr_upload_fine(ilevel)
! Restriction (averaging down) of the CR variables cruold(:,1:ncrvar).
!------------------------------------------------------------------------
  use amr_commons
  use cr_parameters
  use cr_hydro_commons
  implicit none
  integer::ilevel
  integer::i,ncache,igrid,ngrid,ind,iskip,nsplit,icell
  integer,dimension(1:nvector),save::ind_grid,ind_cell,ind_split
  logical,dimension(1:nvector),save::ok
!------------------------------------------------------------------------
  if(.not.cr_advect)return
  if(ilevel==nlevelmax)return
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

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
           ind_cell(i)=iskip+ind_grid(i)
        end do

        ! Gather split cells
        do i=1,ngrid
           ok(i)=son(ind_cell(i))>0
        end do

        ! Count split cells
        nsplit=0
        do i=1,ngrid
           if(ok(i))nsplit=nsplit+1
        end do

        ! Upload for selected cells
        if(nsplit>0)then
           icell=0
           do i=1,ngrid
              if(ok(i))then
                 icell=icell+1
                 ind_split(icell)=ind_cell(i)
              end if
           end do
           call cr_upl(ind_split,nsplit)
        end if

     end do
     ! End loop over cells

  end do
  ! End loop over grids

111 format('   Entering cr_upload_fine for level',i2)

END SUBROUTINE cr_upload_fine
!################################################################
!################################################################
!################################################################
!################################################################
SUBROUTINE cr_upl(ind_cell,ncell)
! Average each coarse split cell's 2^ndim children into the parent,
! over the CR variables 1:ncrvar.
!------------------------------------------------------------------------
  use amr_commons
  use cr_parameters
  use cr_hydro_commons
  implicit none
  integer::ncell
  integer,dimension(1:nvector)::ind_cell
  integer ::ivar,i,ind_son,iskip_son
  integer ,dimension(1:nvector),save::igrid_son,ind_cell_son
  real(dp),dimension(1:nvector),save::getx
!------------------------------------------------------------------------
  ! Get child oct index
  do i=1,ncell
     igrid_son(i)=son(ind_cell(i))
  end do

  ! Loop over CR variables
  do ivar=1,ncrvar
     getx(1:ncell)=0.0d0
     do ind_son=1,twotondim
        iskip_son=ncoarse+(ind_son-1)*ngridmax
        do i=1,ncell
           ind_cell_son(i)=iskip_son+igrid_son(i)
        end do
        do i=1,ncell
           getx(i)=getx(i)+cruold(ind_cell_son(i),ivar)
        end do
     end do

     ! Scatter result to cells
     do i=1,ncell
        cruold(ind_cell(i),ivar)=getx(i)/dble(twotondim)
     end do

  end do
  ! End loop over variables

END SUBROUTINE cr_upl
!################################################################
!################################################################
!################################################################
!################################################################
SUBROUTINE cr_interpol_hydro(u1,u2,nn)
! Interpolate the CR buffer (E,F over 1:ncrvar) from a coarse father cell to
! its 2^ndim children as cell-centered scalars with the chosen slope limiter.
!------------------------------------------------------------------------
  use amr_commons
  use cr_parameters
  use cr_hydro_commons
  use hydro_parameters, only: interpol_type
  implicit none
  integer::nn
  real(dp),dimension(1:nvector,0:twondim  ,1:ncrvar)::u1
  real(dp),dimension(1:nvector,1:twotondim,1:ncrvar)::u2
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
  do ivar=1,ncrvar

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
!################################################################
!################################################################
!################################################################
!################################################################
