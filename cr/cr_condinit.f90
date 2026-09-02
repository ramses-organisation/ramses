!================================================================
!================================================================
!================================================================
!================================================================
subroutine cr_condinit(x,u,dx,nn,ilevel)
  use amr_parameters
  use cr_parameters
  implicit none
  integer ::nn                              ! Number of cells
  integer:: ilevel                          ! Refinement level
  real(dp)::dx                              ! Cell size
  real(dp),dimension(1:nvector,1:ncrvar)::u ! CR conservative variables
  real(dp),dimension(1:nvector,1:ndim  )::x ! Cell center position.
  !================================================================
  ! CR initial conditions; u holds per-group [energy, ndim fluxes].
  ! cr_condinit_kind selects an analytic profile; by default the CR state
  ! comes from the regions declared in &CR_PARAMS.
  !================================================================
  integer::i
  real(dp),dimension(1:nvector)::tmp
  real(dp)::pi
#if NDIM>1
  real(dp)::xx,yy,rr,theta
#endif

  u(1:nn,1:ncrvar)=0d0

  select case(trim(cr_condinit_kind))

  case('jiang_411','jiang_414')              ! Gaussian CR-energy pulse
     ! exp(-40*(x-boxlen/2)^2); flux=4/3*v_gas*E_cr=0 (static gas).
     tmp(1:nn)=(x(1:nn,1)-boxlen*0.5d0)**2
     do i=1,nn
        u(i,Ecr_idx(1))=exp(-40d0*tmp(i))
     end do

  case('jiang_412')                          ! 2D Gaussian CR-energy pulse
     tmp(1:nn)=(x(1:nn,1)-boxlen*0.5d0)**2
#if NDIM>1
     tmp(1:nn)=tmp(1:nn)+(x(1:nn,2)-boxlen*0.5d0)**2
#endif
     do i=1,nn
        u(i,Ecr_idx(1))=exp(-40d0*tmp(i))
     end do

  case('jiang_411_triangular')               ! Triangular CR-energy profile
     ! E0-slope*|x-boxlen/2| (E0=2,slope=1); the streaming wave is driven by
     ! the time-dependent imposed boundary in cr_boundana.
     tmp(1:nn)=(x(1:nn,1)-boxlen*0.5d0)**2
     do i=1,nn
        u(i,Ecr_idx(1))=2d0-1d0*sqrt(tmp(i))
     end do

  case('jiang_421')                          ! Sinusoidal CR-energy profile
     pi=acos(-1d0)
     do i=1,nn
        u(i,Ecr_idx(1))=20d0+10d0*sin(pi*(x(i,1)-boxlen*0.5d0))
     end do

  case('jiang_413','jiang_424')              ! CR floor; CRs enter at the wall
     do i=1,nn
        u(i,Ecr_idx(1))=1d-6
        u(i,Ecr_idx(1)+1:Ecr_idx(1)+ndim)=0d0
     end do

#if NDIM>1
  case('jiang_415','jiang_415_donut')        ! CR arc on a magnetic loop
     ! E_cr=12 on the arc (0.25 box < r < 0.35 box, |theta|<pi/12, xx>0),
     ! 10 elsewhere. atan2 avoids a divide-by-zero FPE at xx=0.
     pi=acos(-1d0)
     do i=1,nn
        xx=x(i,1)-boxlen*0.5d0
        yy=x(i,2)-boxlen*0.5d0
        rr=sqrt(xx**2+yy**2)
        theta=atan2(yy,xx)
        if(rr>0.25d0*boxlen .and. rr<0.35d0*boxlen .and. theta>-pi/12d0 .and. &
             & theta<pi/12d0 .and. xx>0d0)then
           u(i,Ecr_idx(1))=1.2d1
        else
           u(i,Ecr_idx(1))=1.0d1
        endif
     end do
#endif

  case('jiang_422')                          ! Regions, plus advected CR flux
     call cr_region_condinit(x,u,dx,nn,ilevel)
     call cr_flux_from_region_velocity(x,u,dx,nn)

  case DEFAULT                               ! CR regions from &CR_PARAMS
     call cr_region_condinit(x,u,dx,nn,ilevel)

  end select

end subroutine cr_condinit
!================================================================
!================================================================
!================================================================
!================================================================
subroutine cr_flux_from_region_velocity(x,u,dx,nn)
  ! Set the first CR group's x-flux to F = 4/3 u_region(k) E_cr on the same
  ! square CR regions cr_region_condinit filled.
  use amr_parameters
  use cr_parameters
  use hydro_parameters, only: u_region
  implicit none
  integer ::nn
  real(dp)::dx
  real(dp),dimension(1:nvector,1:ncrvar)::u
  real(dp),dimension(1:nvector,1:ndim  )::x
  integer::i,k
  real(dp)::r,xn,yn,zn,en

  do k=1,cr_nregion
     if(cr_region_type(k) .ne. 'square')cycle
     en=cr_exp_region(k)
     do i=1,nn
        xn=0.0d0; yn=0.0d0; zn=0.0d0
        xn=2.0d0*abs(x(i,1)-cr_reg_x_center(k))/cr_reg_length_x(k)
#if NDIM>1
        yn=2.0d0*abs(x(i,2)-cr_reg_y_center(k))/cr_reg_length_y(k)
#endif
#if NDIM>2
        zn=2.0d0*abs(x(i,3)-cr_reg_z_center(k))/cr_reg_length_z(k)
#endif
        if(cr_exp_region(k)<10)then
           r=(xn**en+yn**en+zn**en)**(1.0/en)
        else
           r=max(xn,yn,zn)
        end if
        if(r<1.0)then
           u(i,Ecr_idx(1)+1)=4d0/3d0*u_region(k)*u(i,Ecr_idx(1))
        end if
     end do
  end do

end subroutine cr_flux_from_region_velocity
