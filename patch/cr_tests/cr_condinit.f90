!================================================================
! Jiang & Oh (2018) CR test suite -- CR initial conditions patch: sets only the
! CR conservative variables u(1:nn,1:ncrvar); test chosen by namelist jiang_test.
!================================================================
subroutine cr_condinit(x,u,dx,nn,ilevel)
  use amr_parameters
  use cr_parameters
  implicit none
  integer ::nn                              ! Number of cells
  integer:: ilevel                          ! Refinement level
  real(dp)::dx                              ! Cell size
  real(dp),dimension(1:nvector,1:ncrvar)::u ! CR conservative variables
  real(dp),dimension(1:nvector,1:ndim  )::x ! Cell center position in [0,boxlen]**ndim
  !----------------------------------------------------------------
  ! u(i,1:ncrvar) is the CR conservative vector for group g:
  !   energy at Ecr_idx(g), the ndim fluxes in the slots after it.
  !----------------------------------------------------------------
  integer::i
  real(dp),dimension(1:nvector)::tmp
  real(dp)::pi
#if NDIM>1
  real(dp)::xx,yy,rr,theta
#endif

  ! Default for every group: zero. Test branches below overwrite the
  ! first group's energy (and, where applicable, flux).
  u(1:nn,1:ncrvar)=0d0

  select case(trim(jiang_test))

  case('')
     ! No override: pure default (cr_efloor floor, zero flux).

  case('411','414')                          ! Gaussian CR-energy pulse
     ! exp(-40*(x-boxlen/2)^2) pulse; flux=4/3*v_gas*E_cr=0 (static gas).
     tmp(1:nn)=(x(1:nn,1)-boxlen*0.5d0)**2
     do i=1,nn
        u(i,Ecr_idx(1))=exp(-40d0*tmp(i))
     end do

  case('412')                                ! 2D Gaussian CR-energy pulse
     ! exp(-40*((x-bc)^2+(y-bc)^2)) pulse; flux=0 (static gas, u_region=0).
     tmp(1:nn)=(x(1:nn,1)-boxlen*0.5d0)**2
#if NDIM>1
     tmp(1:nn)=tmp(1:nn)+(x(1:nn,2)-boxlen*0.5d0)**2
#endif
     do i=1,nn
        u(i,Ecr_idx(1))=exp(-40d0*tmp(i))
     end do

  case('411_triangular')                     ! Triangular CR-energy profile
     ! E0-slope*|x-boxlen/2| (E0=2,slope=1); flux=0 (static gas). The streaming
     ! wave is driven by the time-dependent analytic boundary (cr_boundana, bound_type=3).
     tmp(1:nn)=(x(1:nn,1)-boxlen*0.5d0)**2
     do i=1,nn
        u(i,Ecr_idx(1))=2d0-1d0*sqrt(tmp(i))
     end do

  case('421')                                ! Sinusoidal CR-energy profile
     pi=acos(-1d0)
     do i=1,nn
        u(i,Ecr_idx(1))=20d0+10d0*sin(pi*(x(i,1)-boxlen*0.5d0))
     end do

  case('413','424')                          ! CR floor (CR injected at reflexive BC)
     ! Density step is part of the GAS state (region/condinit); here CR only.
     do i=1,nn
        u(i,Ecr_idx(1))=1d-6
        u(i,Ecr_idx(1)+1:Ecr_idx(1)+ndim)=0d0
     end do

#if NDIM>1
  case('415', '415_donut')                                ! 2D magnetic loop: CR-energy arc enhancement
     ! CR energy =12 on the loop arc (0.25 box < r < 0.35 box, theta in
     ! [-pi/12,pi/12], xx>0), =10 elsewhere. atan2 (not atan(yy/xx)) avoids
     ! a divide-by-zero FPE at xx=0.
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

  case('422')                                ! CR energy from cr_region_u; flux from gas velocity
     ! E_cr from cr_region_u (square regions); x-flux = 4/3 u_region(k) E_cr.
     call cr_region_condinit(x,u,dx,nn,ilevel)
     call cr_flux_from_region_velocity(x,u,dx,nn)

  case('tp-nostream','tp-stream-va075', 'tp-stream-va15', '423')                  ! Region-based: CR energy+flux from cr_region_u
     ! Entire CR state comes from cr_region_u (energy at Ecr_idx(1), fluxes after).
     ! tp-*: E_cr contact discontinuity (3|1) across x=5, zero flux.
     ! 423 (2D CR blast/Sedov): uniform E_cr=1d-10 floor, gas blast sweeps it.
     call cr_region_condinit(x,u,dx,nn,ilevel)

  case default
     write(*,*)'cr_condinit: unknown jiang_test = "'//trim(jiang_test)//'"'
     call clean_stop
  end select

end subroutine cr_condinit
!================================================================
!================================================================
!================================================================
!================================================================
! cr_region_condinit lives in core cr/cr_init_flow_fine.f90; defining it here too
! would be a duplicate symbol at link (both are in AMRLIB), so it is omitted here.
!================================================================
!================================================================
!================================================================
!================================================================
subroutine cr_flux_from_region_velocity(x,u,dx,nn)
  ! jiang_test='422': set the first CR group's x-flux F = 4/3 u_region(k) E_cr
  ! (gas region velocity) on the same square regions cr_region_condinit filled.
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
