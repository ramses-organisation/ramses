!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_flow
  use amr_commons
  use hydro_commons, ONLY: nvar, uold
  implicit none

  integer::ilevel,ivar

  if(verbose)write(*,*)'Entering init_flow'
  do ilevel=nlevelmax,1,-1
     if(ilevel>=levelmin)call init_flow_fine(ilevel)
     call upload_fine(ilevel)
     do ivar=1,nvar
        call make_virtual_fine_dp(uold(1,ivar),ilevel)
     end do
     if(simple_boundary)call make_boundary_hydro(ilevel)
  end do
  if(verbose)write(*,*)'Complete init_flow'

end subroutine init_flow
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_flow_fine(ilevel)
  use amr_commons
  use hydro_commons
  use cooling_module
  implicit none
  integer::ilevel

  integer::i,icell,igrid,ncache,iskip,ngrid
  integer::ind,idim,ivar,imat,ix,iy,iz,nx_loc
  integer::i1,i2,i3,i1_min,i1_max,i2_min,i2_max,i3_min,i3_max
  integer ,dimension(1:nvector),save::ind_grid,ind_cell

  real(dp)::dx,rr,vx,vy,vz,ek,ei,scale_T2,xx1,xx2,xx3,dx_loc,scale
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector)       ,save::vv
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector,1:nvar),save::uu

  real(dp),allocatable,dimension(:,:,:)::init_array
  real(sp),allocatable,dimension(:,:)  ::init_plane

  logical::error
  character(LEN=80)::filename

  ncache=active(ilevel)%ngrid
  if(ncache==0)return
  if(verbose)write(*,111)ilevel

  ! Mesh size at level ilevel in coarse cell units
  dx=0.5D0**ilevel

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  !-------------------------------------------------------
  ! Compute initial conditions from subroutine condinit
  !-------------------------------------------------------
  ! Local constants
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  dx_loc=dx*scale

  ! Loop over grids by vector sweeps
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
        ! Call initial condition routine
        call condinit(xx,uu,dx_loc,ngrid)
        ! Scatter variables
        do ivar=1,nvar
           do i=1,ngrid
              uold(ind_cell(i),ivar)=uu(i,ivar)
           end do
        end do
     end do
     ! End loop over cells
  end do
  ! End loop over grids

111 format('   Entering init_flow_fine for level ',I2)

end subroutine init_flow_fine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine region_condinit(x,q,f,g,dx,nn)
  use amr_parameters
  use hydro_parameters
  implicit none
  integer ::nn
  real(dp)::dx
  real(dp),dimension(1:nvector,1:npri)::q
  real(dp),dimension(1:nvector,1:nmat)::f,g
  real(dp),dimension(1:nvector,1:ndim)::x

  integer::i,ivar,imat,k
  real(dp)::vol,r,xn,yn,zn,en,ftot,radius,twopi,fourpi

  fourpi=4.0d0*ACOS(-1.0d0)
  twopi=2.0d0*ACOS(-1.0d0)

  ! Set some (tiny) default values in case n_region=0
  q(1:nn,1)=0.0d0
#if NDIM>1
  q(1:nn,2)=0.0d0
#endif
#if NDIM>2
  q(1:nn,3)=0.0d0
#endif
  q(1:nn,ndim+1:ndim+2*nmat)=smallr*smallc**2
  do imat=1,nmat
     f(1:nn,imat)=smallf
     g(1:nn,imat)=smallf
  end do

  ! Loop over initial conditions regions
  do k=1,nregion

     ! For "square" regions only:
     if(region_type(k) .eq. 'square')then

        ! Exponent of choosen norm
        en=exp_region(k)

        do i=1,nn

           ! Compute position in normalized coordinates
           xn=0.0d0; yn=0.0d0; zn=0.0d0
           xn=2.0d0*abs(x(i,1)-x_center(k))/length_x(k)
#if NDIM>1
           yn=2.0d0*abs(x(i,2)-y_center(k))/length_y(k)
#endif
#if NDIM>2
           zn=2.0d0*abs(x(i,3)-z_center(k))/length_z(k)
#endif
           ! Compute cell "radius" relative to region center
           if(exp_region(k)<10)then
              r=(xn**en+yn**en+zn**en)**(1.0/en)
           else
              r=max(xn,yn,zn)
           end if

           ! If cell lies within region,
           ! REPLACE primitive variables by region values
           if(r<1.0)then
              imat=1
              f(i,imat)         = max(f1_region(k),smallf)
              g(i,imat)         = max(d1_region(k),smallr)
              q(i,ndim+imat)    = max(p1_region(k),smallr*smallc**2)
              if(nmat>1)then
                 imat=imat+1
                 f(i,imat)      = max(f2_region(k),smallf)
                 g(i,imat)      = max(d2_region(k),smallr)
                 q(i,ndim+imat) = max(p2_region(k),smallr*smallc**2)
              end if
              if(nmat>2)then
                 imat=imat+1
                 f(i,imat)      = max(f3_region(k),smallf)
                 g(i,imat)      = max(d3_region(k),smallr)
                 q(i,ndim+imat) = max(p3_region(k),smallr*smallc**2)
              end if
              if(nmat>3)then
                 imat=imat+1
                 f(i,imat)      = max(f4_region(k),smallf)
                 g(i,imat)      = max(d4_region(k),smallr)
                 q(i,ndim+imat) = max(p4_region(k),smallr*smallc**2)
              end if
              ! Normalize volume fractions
              ftot=1d-15
              do imat=1,nmat
                 ftot=ftot+f(i,imat)
              end do
              do imat=1,nmat
                 f(i,imat)=f(i,imat)/ftot
              end do
              ! Compute velocities
              q(i,1)=u_region(k)
#if NDIM>1
              q(i,2)=v_region(k)
#endif
#if NDIM>2
              q(i,3)=w_region(k)
#endif
           end if

        end do

     end if

     ! For "point" regions only:
     if(region_type(k) .eq. 'point')then
        radius=x_center(k)
        if(geom>1)then
           radius=max(radius,dx/2.0)
        endif
        ! Volume elements
        vol=dx**ndim
        ! Compute CIC weights relative to region center
        do i=1,nn
           xn=1.0; yn=1.0; zn=1.0
           xn=max(1.0-abs(x(i,1)-radius     )/dx,0.0_dp)
#if NDIM>1
           yn=max(1.0-abs(x(i,2)-y_center(k))/dx,0.0_dp)
#endif
#if NDIM>2
           zn=max(1.0-abs(x(i,3)-z_center(k))/dx,0.0_dp)
#endif
           r=xn*yn*zn
           if(geom==2)r=r/(twopi*x(i,1))
           if(geom==3)r=r/(fourpi*(x(i,1)**2+dx**2/12.0))

           ! If cell lies within CIC cloud,
           ! ADD to primitive variables the region values
           q(i,1)          = q(i,1) + u_region(k)*r
#if NDIM>1
           q(i,2)          = q(i,2) + v_region(k)*r
#endif
#if NDIM>2
           q(i,3)          = q(i,3) + w_region(k)*r
#endif
           do imat=1,nmat
            q(i,ndim+imat) = q(i,ndim+imat) + p_region(k)*r/vol
           end do
        end do
     end if

  end do

  return
end subroutine region_condinit
