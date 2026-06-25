!================================================================
!================================================================
!================================================================
!================================================================
subroutine  condinit(x,u,dx,nn)
  use amr_commons
  use hydro_parameters
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar+3)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  logical,save:: first_call = .true.       ! True if this is the first call to condinit

  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:4): d.u,d.v,d.w, U(i,5): E, U(i,6:8): Bleft,
  ! U(i,nvar+1:nvar+3): Bright
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:4):u,v,w, Q(i,5): P, Q(i,6:8): Bleft,
  ! Q(i,nvar+1:nvar+3): Bright, Q(i,9:8+nener): Er (if FLD)
  ! If nvar > 8+nener, remaining variables (9+nener:nvar) are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
#if NENER>0
  integer::irad
#endif
#if NVAR>NHYDRO+NENER
  integer::ivar
#endif
  real(dp),dimension(1:nvector,1:nvar+3),save::q   ! Primitive variables

  select case (condinit_kind)

  case('region')
    ! Call built-in initial condition generator
     call region_condinit(x, q, dx, nn)

  case('orzag_tang')
     call orzag_tang_condinit(x, q, dx, nn)

  case('coeur')
     call coeur_condinit(x, q, dx, nn)

  case('ponomarenko')
      call ponomarenko_condinit(x, q, dx, nn)

  case('collapse')
     call collapse_condinit(x, q, dx, nn)

  case('nimhd_diffusion')
     call nimhd_diffusion_condinit(x, q, dx, nn)

  ! Add here, if you wish, some user-defined initial conditions
  ! ........

  case DEFAULT
     if (myid == 1.and. first_call)  write(*,*) "[condinit] Void or invalid condinit_kind, using default IC"
     call region_condinit(x, q, dx, nn)

  end select

  first_call = .false.

  ! Convert primitive to conservative variables
  ! density -> density
  u(1:nn,1)=q(1:nn,1)
  ! velocity -> momentum
  u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
  u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
  u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
  ! kinetic energy
  u(1:nn,neul)=0.0d0
  u(1:nn,neul)=u(1:nn,neul)+0.5*q(1:nn,1)*q(1:nn,2)**2
  u(1:nn,neul)=u(1:nn,neul)+0.5*q(1:nn,1)*q(1:nn,3)**2
  u(1:nn,neul)=u(1:nn,neul)+0.5*q(1:nn,1)*q(1:nn,4)**2
  ! pressure -> total fluid energy
  select case(condinit_kind)
    case('collapse')
      u(1:nn,neul)=u(1:nn,neul)+q(1:nn,neul)
    case DEFAULT
      u(1:nn,neul)=u(1:nn,neul)+q(1:nn,neul)/(gamma-1.0d0)
  end select
  ! magnetic energy -> total fluid energy
  u(1:nn,neul)=u(1:nn,neul)+0.125d0*(q(1:nn,6)+q(1:nn,nvar+1))**2
  u(1:nn,neul)=u(1:nn,neul)+0.125d0*(q(1:nn,7)+q(1:nn,nvar+2))**2
  u(1:nn,neul)=u(1:nn,neul)+0.125d0*(q(1:nn,8)+q(1:nn,nvar+3))**2
  u(1:nn,6:8)=q(1:nn,6:8)
  u(1:nn,nvar+1:nvar+3)=q(1:nn,nvar+1:nvar+3)
#if USE_FLD==0
#if NENER>0
  ! radiative pressure -> radiative energy
  ! radiative energy -> total fluid energy
  do irad=1,nener
     u(1:nn,nhydro+irad)=q(1:nn,nhydro+irad)/(gamma_rad(irad)-1.0d0)
     u(1:nn,neul)=u(1:nn,neul)+u(1:nn,nhydro+irad)
  enddo
#endif
#if NVAR>NHYDRO+NENER
  ! passive scalars
  do ivar=nhydro+1+nener,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do
#endif
#else
#if NENER>0
  ! non-thermal pressure -> non-thermal energy
  ! non-thermal energy   -> total fluid energy
  do irad=1,nener-ngrp
     u(1:nn,8+irad)=q(1:nn,8+irad)/(gamma_rad(irad)-1.0d0)
     u(1:nn,5)=u(1:nn,5)+u(1:nn,8+irad)
  enddo
 ! Radiative transfer
#if NGRP>0
  ! radiative energy   -> total fluid energy
  do ivar=1,ngrp
     u(1:nn,firstindex_er+ivar)= q(1:nn,firstindex_er+ivar)
     u(1:nn,5)=u(1:nn,5)+ u(1:nn,firstindex_er+ivar)
  enddo
#if USE_M_1==1
  ! radiative flux
  do ivar=1,ndim*ngrp
     do i=1,ncache
        u(1:nn,fisrtindex_fr+ivar)=q(1:nn,firstindex+ivar)
     end do
!      write(ilun)xdp
  end do
#endif
#endif
#endif
#if NPSCAL>0
  ! passive scalars
  do ivar=1,npscal
     u(1:nn,firstindex_pscal+ivar)=q(1:nn,1)*q(1:nn,firstindex_pscal+ivar)
  end do
  ! Internal energy
  select case(condinit_kind)
    case('collapse')
       u(1:nn,nvar)=q(1:nn,5)
    case DEFAULT
      u(1:nn,nvar)=q(1:nn,5)/(gamma-1.0d0)
  end select
#endif
#endif

end subroutine condinit
!================================================================
!================================================================
!================================================================
!================================================================
subroutine orzag_tang_condinit(x,q,dx,nn)
  use amr_parameters
  use hydro_parameters
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar+3)::q ! Primitive variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates Orzag-Tang initial conditions for RAMSES.
  !================================================================
  integer::i
  real(dp)::xc,xr,xl,yl,yr,yc,al,ar,B0,pi

  pi=ACOS(-1.0_dp)
  B0=1.0_dp/sqrt(4.0_dp*pi)
  do i=1,nn

     xl=x(i,1)-0.5_dp*dx
     xr=x(i,1)+0.5_dp*dx
     xc=x(i,1)
     yl=x(i,2)-0.5_dp*dx
     yr=x(i,2)+0.5_dp*dx
     yc=x(i,2)

     q(i,1)=25.0_dp/(36.0_dp*pi)
     q(i,2)=-sin(2.0_dp*pi*yc)
     q(i,3)=+sin(2.0_dp*pi*xc)
     q(i,4)=0.0_dp
     q(i,5)=5.0_dp/(12.0_dp*pi)

     Ar = B0*(cos(4.0_dp*pi*xl)/(4.0_dp*pi)+cos(2.0_dp*pi*yr)/(2.0_dp*pi))
     Al = B0*(cos(4.0_dp*pi*xl)/(4.0_dp*pi)+cos(2.0_dp*pi*yl)/(2.0_dp*pi))
     q(i,6)=(Ar-Al)/dx
     Ar = B0*(cos(4.0_dp*pi*xr)/(4.0_dp*pi)+cos(2.0_dp*pi*yr)/(2.0_dp*pi))
     Al = B0*(cos(4.0_dp*pi*xr)/(4.0_dp*pi)+cos(2.0_dp*pi*yl)/(2.0_dp*pi))
     q(i,nvar+1)=(Ar-Al)/dx
     Ar = B0*(cos(4.0_dp*pi*xr)/(4.0_dp*pi)+cos(2.0_dp*pi*yl)/(2.0_dp*pi))
     Al = B0*(cos(4.0_dp*pi*xl)/(4.0_dp*pi)+cos(2.0_dp*pi*yl)/(2.0_dp*pi))
     q(i,7)=(Al-Ar)/dx
     Ar = B0*(cos(4.0_dp*pi*xr)/(4.0_dp*pi)+cos(2.0_dp*pi*yr)/(2.0_dp*pi))
     Al = B0*(cos(4.0_dp*pi*xl)/(4.0_dp*pi)+cos(2.0_dp*pi*yr)/(2.0_dp*pi))
     q(i,nvar+2)=(Al-Ar)/dx

     q(i,8)=0.0_dp
     q(i,nvar+3)=0.0_dp
  end do

end subroutine orzag_tang_condinit
!================================================================
!================================================================
!================================================================
!================================================================
subroutine coeur_condinit(x,q,dx,nn)
  use amr_parameters
  use hydro_parameters
  use collapse_parameters
  use constants,ONLY:kB,mH,pi,M_sun
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar+3)::q ! Primitive variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates collapsing core initial conditions for RAMSES.
  !================================================================

  integer :: i,id,iu,iv,iw,ip
  real(dp):: x0,y0,z0,rc,rs,xx,yy,zz,r0,d0,B0,p0,omega0,C_s,Temp
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v

 ! Cloud parameters
  mass_c = 1.0d0 ! in solar masses
  delta_rho=0.1
  alpha_dense_core=0.1
  beta_dense_core=0.01
  crit_dense_core=0.1

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  id=1; iu=2; iv=3; iw=4; ip=5
  x0=0.5d0*boxlen
  y0=0.5d0*boxlen
  z0=0.5d0*boxlen
  Temp = 10.0d0
  r0=(alpha_dense_core*2.*6.67d-8*(mass_c*M_sun)*mu_gas*mH/(5.*kB*Temp))/scale_l
  ! cloud density equal to unity
  d0 = 3.0d0*(mass_c*M_sun/scale_l**3/scale_d)/(4.0d0*pi*r0**3.)
  ! cloud rotation
  omega0 = sqrt(beta_dense_core*4.*pi*d0)

  ! sound speed
  C_s = sqrt( Temp / scale_T2 )
  ! cloud pressure
  p0 = C_s**2*d0
  ! vertical magnetic field

  B0 = sqrt(4.*pi/5.)/0.53*(crit_dense_core*d0*r0)  / 2.26d0
  DO i=1,nn
     xx=x(i,1)-x0
     yy=x(i,2)-y0
     zz=x(i,3)-z0
     rc=sqrt(xx**2+yy**2)
     rs=sqrt(xx**2+yy**2+zz**2)

     !Bx component
     q(i,6     ) = 0.
     q(i,nvar+1) = 0.

     !By component
     q(i,7     ) = 0.
     q(i,nvar+2) = 0.

     !Bz component
     q(i,8     ) = B0
     q(i,nvar+3) = B0

     IF(rs .le. r0) THEN
       q(i,id) = d0*(1.0+delta_rho*cos(2.*atan(yy/xx)))!(2.0*(xx/rc)**2-1.0))
       q(i,iu) = omega0 * yy
       q(i,iv) = -omega0 * xx
       q(i,iw) = 0.0
       q(i,ip) = p0
     ELSE
       q(i,id) = d0/100.
       xx = r0 * xx / rc
       yy = r0 * yy / rc
       q(i,iu) = 0.0! omega0 * yy
       q(i,iv) = 0.0!-omega0 * xx
       q(i,iw) = 0.0
       q(i,ip) = p0/100.
     ENDIF
  ENDDO

end subroutine coeur_condinit
!================================================================
!================================================================
!================================================================
!================================================================
! subroutine collapse_condinit(x,q,dx,nn)
!   use amr_commons, only:myid
!   use amr_parameters
!   use hydro_commons
!   use poisson_parameters
!   use constants, only:mH,kB,M_sun,pc2cm
!   implicit none
!   integer ::nn                              ! Number of cells
!   real(dp)::dx                              ! Cell size
!   real(dp),dimension(1:nvector,1:nvar+3)::q ! Primitive variables
!   real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
!   !================================================================
!   ! This routine generates initial conditions of a collapsing core
!   !================================================================
!   integer :: i,j,k,id,iu,iv,iw,ip
!   real(dp):: x0,y0,z0,rc,rs,xx,yy,zz,pi,r0,d0,B0,p0,omega0,mass_c_cu,scale_m
!   integer :: ivar, np
!   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
!   real(dp),dimension(1:3,1:3):: rot_M,rot_invM,rot_tilde
!   real(dp):: theta_mag_radians

!   logical,save:: first=.true.
!   real(dp),dimension(1:3,1:100,1:100,1:100),save::q_idl
!   real(dp),save::vx_tot,vy_tot,vz_tot,vx2_tot,vy2_tot,vz2_tot
!   integer,save:: n_size
!   integer:: ind_i, ind_j, ind_k
!   real(dp),save:: ind,seed1,seed2,seed3,xi,yi,zi,vx,vy,vz
!   real(dp),save:: C_s,v_rms
!   integer, save :: count_vrms
! #if USE_FLD==1
!   real(dp)::radiation_source
! #endif

!   id=1; iu=2; iv=3; iw=4; ip=5
!   x0=0.5*boxlen
!   y0=0.5*boxlen
!   z0=0.5*boxlen
!   pi=acos(-1.0d0)

!   ! Conversion factor from user units to cgs units
!   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
!   scale_m=scale_d*scale_l**3

!   ! cloud mass (warning mass_c should not be changed, because condinit called not just once)
!   ! mass_c    is in solar mass
!   ! mass_c_cu is in code units
!   mass_c_cu = mass_c * (M_sun / scale_m )

!   ! cloud radius
!   r0=(alpha_dense_core*2.*6.67d-8*mass_c_cu*scale_m*mu_gas*mH/(5.*kB*T_eos))/scale_l

!   ! cloud density
!   d0 = 3.0d0*mass_c_cu/(4.0d0*pi*r0**3.)

!   ! cloud rotation ! remember that G=1 in code units
!   omega0 = sqrt(beta_dense_core*4.*pi*d0)

!   ! cloud pressure ! remember that G=1 in code units
!   C_s = sqrt(kB*T_eos/(mu_gas*mH))/scale_v

!   p0 = d0*c_s**2
!   ! vertical magnetic field ! remember that G=1 in code units (and that B as a factor 1/sqrt(4pi) between SI and Gaussian units)
!   B0 = sqrt(4.*pi/5.)/0.53*crit_dense_core*d0*r0

!   ! angle between the rotation axis and the magnetic field
!   theta_mag_radians= theta_mag/180.0d0*pi

  
!   rot_M(1,1:3) = (/cos(theta_mag_radians),0.0d0,-sin(theta_mag_radians)/)
!   rot_M(2,1:3) = (/0.0d0,1.0d0,0.0d0/)
!   rot_M(3,1:3) = (/sin(theta_mag_radians),0.0d0,cos(theta_mag_radians)/)

!   rot_invM(1,1:3) = (/cos(theta_mag_radians),0.0d0,sin(theta_mag_radians)/)
!   rot_invM(2,1:3) = (/0.0d0,1.0d0,0.0d0/)
!   rot_invM(3,1:3) = (/-sin(theta_mag_radians),0.0d0,cos(theta_mag_radians)/)

!   rot_tilde(1,1:3) = (/0.0d0,1.0d0,0.0d0/)
!   rot_tilde(2,1:3) = (/-1.0d0,0.0d0,0.0d0/)
!   rot_tilde(3,1:3) = (/0.0d0,0.0d0,0.0d0/)

!   if(first) then
!     ! sound speed
!     !C_s could be defined equivalently as sqrt( T_eos / (mu_gas*scale_T2) )
    
!     vx_tot=0.d0
!     vy_tot=0.d0
!     vz_tot=0.d0
!     vx2_tot=0.d0
!     vy2_tot=0.d0
!     vz2_tot=0.d0
!     v_rms=0.d0
!     count_vrms=0
!     if(Mach .ne. 0)then
!       if (myid==1) write(*,*) 'Read the file which contains the initial turbulent velocity field'
!       open(20,file='init_turb.data',form='formatted')
!       read(20,*) n_size, ind, seed1,seed2,seed3
!       if(n_size .ne. 100) then
!           write(*,*) 'Unexpected field size'
!           stop
!       endif
!       do k=1,n_size
!         do j=1,n_size
!           do i=1,n_size
!             read(20,*)xi,yi,zi,vx,vy,vz
!             q_idl(1,i,j,k) = vx
!             q_idl(2,i,j,k) = vy
!             q_idl(3,i,j,k) = vz
!             xi = boxlen*((i-0.5)/n_size)-x0
!             yi = boxlen*((j-0.5)/n_size)-y0
!             zi = boxlen*((k-0.5)/n_size)-z0
!             rs=sqrt(xi**2+yi**2+zi**2)

!             IF(rs .le. r0) THEN
!               !print*, vx_tot,vy_tot,vz_tot,vx2_tot,vy2_tot,vz2_tot
!               vx_tot = vx_tot + vx
!               vy_tot = vy_tot + vy
!               vz_tot = vz_tot + vz

!               vx2_tot = vx2_tot + vx**2
!               vy2_tot = vy2_tot + vy**2
!               vz2_tot = vz2_tot + vz**2

!               count_vrms=count_vrms+1
!             end if
!           end do
!         end do
!       end do
!       close(20)
!       v_rms=sqrt((vx2_tot+vy2_tot+vz2_tot)/dble(count_vrms)-((vx_tot+vy_tot+vz_tot)/dble(count_vrms))**2)
!       if (myid == 1) print *, 'v_rms for given seed =',v_rms
!       ! correction factor to have the expected Mach number stored in v_rms
!       v_rms = Mach*C_s/v_rms
!       if (myid == 1) print *, 'correction factor for turbulent field =',v_rms
!    end if

!    if(myid==1)then
!       print*,'alpha_dense_core=',alpha_dense_core
!       print*,'beta_dense_core=',beta_dense_core
!       print*,'Mass=',mass_c,' Msun'
!       print*, 'Mass_cu =', mass_c_cu
!       print*,'d0 (in g/cc)=',d0*scale_d
!       print*,'Turbulent Mach,cs (km/s)=',Mach,C_s*scale_v/1e5
!       print*,'r0,boxlen (in code units)=',r0,boxlen
!       print*,'r0,boxlen (in pc)=',r0*scale_l/pc2cm,boxlen*scale_l/pc2cm
!     endif
!     first = .false.
!   end if



!   DO i=1,nn
!      xx=x(i,1)-x0
!      yy=x(i,2)-y0
!      zz=x(i,3)-z0
!      rc=sqrt(xx**2+yy**2)
!      rs=sqrt(xx**2+yy**2+zz**2)

!      q(i,iu) = 0.0d0
!      q(i,iv) = 0.0d0
!      q(i,iw) = 0.0d0
!      if(Mach .ne. 0)then
!       !initialise the turbulent velocity field
!       !make a zero order interpolation (should be improved)
!       ind_i = int((x(i,1)/boxlen)*n_size)+1
!       ind_j = int((x(i,2)/boxlen)*n_size)+1
!       ind_k = int((x(i,3)/boxlen)*n_size)+1
!       ! safe check
!       if( ind_i .lt. 1 .or. ind_i .gt. n_size) write(*,*) 'ind_i ',ind_i,(x(i,1)/boxlen)*n_size+1,n_size
!       if( ind_j .lt. 1 .or. ind_j .gt. n_size) write(*,*) 'ind_j ',ind_j
!       if( ind_k .lt. 1 .or. ind_k .gt. n_size) write(*,*) 'ind_k ',ind_k
!     end if

!      IF(rs .le. r0) THEN
!         ! if(theta_mag.eq.0.0d0) then
!         !   q(i,id) = d0*(1.0+delta_rho*cos(2.*atan(yy/xx)))!(2.0*(xx/rc)**2-1.0))
!         !   q(i,iu) = omega0 * yy
!         !   q(i,iv) = -omega0 * xx
!         !   q(i,iw) = 0.0
!         ! else
!           q(i,id) = d0*(1.0+delta_rho*cos(2.*atan(yy/(cos(theta_mag_radians)*xx-sin(theta_mag_radians)*zz))))

!           if(Mach .ne. 0)then
!             !print*,omega0*yy,omega0*xx,v_rms,v_rms*(q_idl(1,ind_i,ind_j,ind_k)-vx_tot/dble(count_vrms))
!             q(i,iu) =  v_rms*(q_idl(1,ind_i,ind_j,ind_k)-vx_tot/dble(count_vrms))
!             q(i,iv) =  v_rms*(q_idl(2,ind_i,ind_j,ind_k)-vy_tot/dble(count_vrms))
!             q(i,iw) =  v_rms*(q_idl(3,ind_i,ind_j,ind_k)-vz_tot/dble(count_vrms))
!           end if
!           q(i,iu:iw) = q(i,iu:iw) + matmul(rot_invM,omega0*matmul(rot_tilde,matmul(rot_M,(/xx,yy,zz/))))
!         ! endif
!         q(i,ip) = p0
!         !Bx component
!         q(i,6     ) = 0.0d0
!         q(i,nvar+1) = 0.0d0

!         !By component
!         q(i,7     ) = 0.0d0
!         q(i,nvar+2) = 0.0d0

!         !Bz component
!         q(i,8     ) = B0
!         q(i,nvar+3) = q(i,8     )
!      ELSE
!        q(i,id) = d0/100.
!        xx = r0 * xx / rc
!        yy = r0 * yy / rc

!        if(Mach .ne. 0)then
!         q(i,iu) = v_rms*(q_idl(1,ind_i,ind_j,ind_k)-vx_tot/dble(count_vrms))! omega0 * yy
!         q(i,iv) = v_rms*(q_idl(2,ind_i,ind_j,ind_k)-vy_tot/dble(count_vrms))!-omega0 * xx
!         q(i,iw) = v_rms*(q_idl(3,ind_i,ind_j,ind_k)-vz_tot/dble(count_vrms))
!        end if

!       !  q(i,iu) = 0.0! omega0 * yy
!       !  q(i,iv) = 0.0!-omega0 * xx
!       !  q(i,iw) = 0.0

!         q(i,ip) = p0/100.d0
!         !Bx component
!         q(i,6     ) = 0.0d0
!         q(i,nvar+1) = 0.0d0

!         !By component
!         q(i,7     ) = 0.0d0
!         q(i,nvar+2) = 0.0d0

!         !Bz component
!         q(i,8     ) = B0/100**(2./3.)
!         q(i,nvar+3) = q(i,8     )
!      ENDIF

! #if USE_FLD==1
! #if NGRP>0
!     do ivar=1,ngrp
!         q(i,firstindex_er+ivar) = radiation_source(T_eos,ivar)/(scale_d*scale_v**2)
!     enddo
! #endif     
! #endif

!   ENDDO 

! end subroutine collapse_condinit


!================================================================
!================================================================
!================================================================
!================================================================
subroutine collapse_condinit(x,q,dx,nn)
  use amr_commons, only:myid
  use amr_parameters
  use hydro_commons
  use poisson_parameters
  use constants, only:mH,kB,M_sun,pc2cm
  use cloud_module
  implicit none
  integer ::nn                              ! Number of cells
  real(dp)::dx                              ! Cell size
  real(dp),dimension(1:nvector,1:nvar+3)::q ! Primitive variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates initial conditions of a collapsing core
  !================================================================
  integer :: i,j,k,id,iu,iv,iw,ip
  real(dp):: x0,y0,z0,rc,rs,xx,yy,zz,pi,r0,d0,B0,p0,omega0,mass_c_cu,scale_m
  integer :: ivar, np
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:3,1:3):: rot_M,rot_invM,rot_tilde
  real(dp):: ee,theta_mag_radians

  logical,save:: first=.true.
  real(dp),dimension(1:3,1:100,1:100,1:100),save::q_idl
  real(dp),save::vx_tot,vy_tot,vz_tot,vx2_tot,vy2_tot,vz2_tot
  integer,save:: n_size
  integer:: ind_i, ind_j, ind_k
  real(dp),save:: ind,seed1,seed2,seed3,xi,yi,zi,vx,vy,vz
  real(dp),save:: C_s,v_rms
  integer, save :: count_vrms
#if USE_FLD==1
  real(dp)::radiation_source
#endif

  integer :: igroup
  real(dp):: mass
  real(dp)::Temp
  real(dp),save:: d_c,B_c,zeta
  real(dp),save:: res_int,r_0,omega,cont_ic,mass_total,mass_tot2,min_col_d,max_col_d
  real(dp):: col_d,eli,sph
  real(dp)::xl,yl,zl,bx,by,bz,dxmin
  integer:: ii,jj,kk,nticks
  real(dp)::ener_rot,ener_grav,ener_therm,ener_grav2,ener_turb,dd
  real(dp),dimension(1000):: mass_rad    


  id=1; iu=2; iv=3; iw=4; ip=5
  x0=0.5*boxlen
  y0=0.5*boxlen
  z0=0.5*boxlen
  pi=acos(-1.0d0)

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=scale_d*scale_l**ndim

  ! cloud mass (warning mass_c should not be changed, because condinit called not just once)
  ! mass_c    is in solar mass
  ! mass_c_cu is in code units
  mass_c_cu = mass_c * (M_sun / scale_m )


  if(bb_test) then

  ! cloud radius
  r0=(alpha_dense_core*2.*6.67d-8*mass_c_cu*scale_m*mu_gas*mH/(5.*kB*T_eos))/scale_l
  ! cloud density
  d0 = 3.0d0*mass_c_cu/(4.0d0*pi*r0**3.)

  ! cloud rotation ! remember that G=1 in code units
  omega0 = sqrt(beta_dense_core*4.*pi*d0)

  ! cloud pressure ! remember that G=1 in code units
  p0 = alpha_dense_core*d0*d0*r0*r0*8.*pi/15.

  ! vertical magnetic field ! remember that G=1 in code units (and that B as a factor 1/sqrt(4pi) between SI and Gaussian units)
  B0 = sqrt(4.*pi/5.)/0.53*crit_dense_core*d0*r0
  ! B0 could be defined equivalently as
  !B0 = mass_c_cu*3./sqrt(5.)/0.53*crit_dense_core/r0**2/sqrt(4.*pi)

  ! angle between the rotation axis and the magnetic field
  theta_mag_radians= theta_mag/180.0d0*pi

  rot_M(1,1:3) = (/cos(theta_mag_radians),0.0d0,-sin(theta_mag_radians)/)
  rot_M(2,1:3) = (/0.0d0,1.0d0,0.0d0/)
  rot_M(3,1:3) = (/sin(theta_mag_radians),0.0d0,cos(theta_mag_radians)/)

  rot_invM(1,1:3) = (/cos(theta_mag_radians),0.0d0,sin(theta_mag_radians)/)
  rot_invM(2,1:3) = (/0.0d0,1.0d0,0.0d0/)
  rot_invM(3,1:3) = (/-sin(theta_mag_radians),0.0d0,cos(theta_mag_radians)/)

  rot_tilde(1,1:3) = (/0.0d0,1.0d0,0.0d0/)
  rot_tilde(2,1:3) = (/-1.0d0,0.0d0,0.0d0/)
  rot_tilde(3,1:3) = (/0.0d0,0.0d0,0.0d0/)

  if(first) then
    ! sound speed
    C_s = sqrt(kB*T_eos/(mu_gas*mH))/scale_v
    !C_s could be defined equivalently as sqrt( T_eos / (mu_gas*scale_T2) )

    vx_tot=0.d0
    vy_tot=0.d0
    vz_tot=0.d0
    vx2_tot=0.d0
    vy2_tot=0.d0
    vz2_tot=0.d0
    v_rms=0.d0
    count_vrms=0
    if(Mach .ne. 0)then
      if (myid==1) write(*,*) 'Read the file which contains the initial turbulent velocity field'
      open(20,file='init_turb.data',form='formatted')
      read(20,*) n_size, ind, seed1,seed2,seed3
      if(n_size .ne. 100) then
          write(*,*) 'Unexpected field size'
          stop
      endif
      do k=1,n_size
        do j=1,n_size
          do i=1,n_size
            read(20,*)xi,yi,zi,vx,vy,vz
            q_idl(1,i,j,k) = vx
            q_idl(2,i,j,k) = vy
            q_idl(3,i,j,k) = vz
            xi = boxlen*((i-0.5)/n_size)-x0
            yi = boxlen*((j-0.5)/n_size)-y0
            zi = boxlen*((k-0.5)/n_size)-z0
            rs=sqrt(xi**2+yi**2+zi**2)

            IF(rs .le. r0) THEN
              !print*, vx_tot,vy_tot,vz_tot,vx2_tot,vy2_tot,vz2_tot
              vx_tot = vx_tot + vx
              vy_tot = vy_tot + vy
              vz_tot = vz_tot + vz

              vx2_tot = vx2_tot + vx**2
              vy2_tot = vy2_tot + vy**2
              vz2_tot = vz2_tot + vz**2

              count_vrms=count_vrms+1
            end if
          end do
        end do
      end do
      close(20)
      v_rms=sqrt((vx2_tot+vy2_tot+vz2_tot)/dble(count_vrms)-((vx_tot+vy_tot+vz_tot)/dble(count_vrms))**2)
      if (myid == 1) print *, 'v_rms for given seed =',v_rms
      ! correction factor to have the expected Mach number stored in v_rms
      v_rms = Mach*C_s/v_rms
      if (myid == 1) print *, 'correction factor for turbulent field =',v_rms
   end if

   if(myid==1)then
      print*,'alpha_dense_core=',alpha_dense_core
      print*,'beta_dense_core=',beta_dense_core
      print*,'Mass=',mass_c,' Msun'
      print*,'d0 (in g/cc)=',d0*scale_d
      print*,'Turbulent Mach,cs (km/s)=',Mach,C_s*scale_v/1e5
      print*,'r0,boxlen (in code units)=',r0,boxlen
      print*,'r0,boxlen (in pc)=',r0*scale_l/pc2cm,boxlen*scale_l/pc2cm
    endif
    first = .false.
  end if



  DO i=1,nn
     xx=x(i,1)-x0
     yy=x(i,2)-y0
     zz=x(i,3)-z0
     rc=sqrt(xx**2+yy**2)
     rs=sqrt(xx**2+yy**2+zz**2)

     !Bx component
     q(i,6     ) = 0.
     q(i,nvar+1) = 0.

     !By component
     q(i,7     ) = 0.
     q(i,nvar+2) = 0.

     !Bz component
     if(rc .le. r0)then
      q(i,8     ) = B0
      q(i,nvar+3) = B0
     else
      q(i,8     ) = B0/100.d0**(2./3.)
      q(i,nvar+3) = B0/100.d0**(2./3.)
     end if

     q(i,iu) = 0.
     q(i,iv) = 0.
     q(i,iw) = 0.
     if(Mach .ne. 0)then
      !initialise the turbulent velocity field
      !make a zero order interpolation (should be improved)
      ind_i = int((x(i,1)/boxlen)*n_size)+1
      ind_j = int((x(i,2)/boxlen)*n_size)+1
      ind_k = int((x(i,3)/boxlen)*n_size)+1
      ! safe check
      if( ind_i .lt. 1 .or. ind_i .gt. n_size) write(*,*) 'ind_i ',ind_i,(x(i,1)/boxlen)*n_size+1,n_size
      if( ind_j .lt. 1 .or. ind_j .gt. n_size) write(*,*) 'ind_j ',ind_j
      if( ind_k .lt. 1 .or. ind_k .gt. n_size) write(*,*) 'ind_k ',ind_k
    end if

     IF(rs .le. r0) THEN
        ! if(theta_mag.eq.0.0d0) then
        !   q(i,id) = d0*(1.0+delta_rho*cos(2.*atan(yy/xx)))!(2.0*(xx/rc)**2-1.0))
        !   q(i,iu) = omega0 * yy
        !   q(i,iv) = -omega0 * xx
        !   q(i,iw) = 0.0
        ! else
          q(i,id) = d0*(1.0+delta_rho*cos(2.*atan(yy/(cos(theta_mag_radians)*xx-sin(theta_mag_radians)*zz))))

          if(Mach .ne. 0)then
            !print*,omega0*yy,omega0*xx,v_rms,v_rms*(q_idl(1,ind_i,ind_j,ind_k)-vx_tot/dble(count_vrms))
            q(i,iu) =  v_rms*(q_idl(1,ind_i,ind_j,ind_k)-vx_tot/dble(count_vrms))
            q(i,iv) =  v_rms*(q_idl(2,ind_i,ind_j,ind_k)-vy_tot/dble(count_vrms))
            q(i,iw) =  v_rms*(q_idl(3,ind_i,ind_j,ind_k)-vz_tot/dble(count_vrms))
          end if
          q(i,iu:iw) = q(i,iu:iw) + matmul(rot_invM,omega0*matmul(rot_tilde,matmul(rot_M,(/xx,yy,zz/))))
        ! endif
       q(i,ip) = p0
     ELSE
       q(i,id) = d0/100.
       xx = r0 * xx / rc
       yy = r0 * yy / rc

       if(Mach .ne. 0)then
        q(i,iu) = v_rms*(q_idl(1,ind_i,ind_j,ind_k)-vx_tot/dble(count_vrms))! omega0 * yy
        q(i,iv) = v_rms*(q_idl(2,ind_i,ind_j,ind_k)-vy_tot/dble(count_vrms))!-omega0 * xx
        q(i,iw) = v_rms*(q_idl(3,ind_i,ind_j,ind_k)-vz_tot/dble(count_vrms))
       end if

      !  q(i,iu) = 0.0! omega0 * yy
      !  q(i,iv) = 0.0!-omega0 * xx
      !  q(i,iw) = 0.0

       q(i,ip) = p0/100.
     ENDIF

    if(eos) then
      !call enerint_eos(q(i,1),T_eos,ee)
      !q(i,  ip) = ee
      !q(i,nvar) = ee
    else
      q(i,  ip) = q(i,1)* C_s**2!/(gamma-1.0d0)
      q(i,nvar) = q(i,ip)
    endif
 
#if USE_FLD==1
#if NGRP>0
    do ivar=1,ngrp
        q(i,firstindex_er+ivar) = radiation_source(T_eos,ivar)/(scale_d*scale_v**2)
    enddo
#endif     
#endif

  ENDDO 

else ! bb_test = .false.

    !do various things which needs to be done only one time
     if( first ) then 
        id=1; iu=2; iv=3; iw=4; ip=5
        pi=acos(-1.0d0)

        if(myid==1) write(*,*) '** ENTER  in condinit **'



        !calculate the mass in code units (Msolar / Mparticle / pc^3
        !    mass_c = mass_c * (Msun / (scale_d * scale_l**3) )
        !    done in calc_boxlen


        if(myid ==1) write(*,*) 'cloud mass (code units) ',mass_c_cu

        !calculate the sound speed
        C_s = sqrt( T_eos / scale_T2 )


        if(myid == 1)  write(*,*) 'T_0 (K) ', T_eos
        if(myid == 1)  write(*,*)  'C_s (code units) ', C_s

        !cont_ic is the density contrast between the edge of the cloud and the intercloud medium
        cont_ic = 10.

        !calculate  zeta=r_ext/r_0
        zeta = sqrt(cont - 1.)


        !calculate an integral used to compute the cloud radius 
        res_int=0.
        do i=1,1000
           res_int = res_int + log(1.+(zeta/1000.*i)**2) * zeta/1000.
           mass_rad(i) = i*zeta/1000. * log(1+(zeta/1000.*i)**2) - res_int
        enddo
        res_int = zeta*log(1.+zeta**2) - res_int


        !now we determine the central density and the external cloud radius
        !we have mass = 2 pi rho_c r_0^2 z_0 * res_int
        !which results from the integration of rho = dc/(1.+(x^2+y^2)/r_O^2+z^2/z_0^2)
        !for (x^2+y^2)/r_O^2+z^2/z_0^2 < zeta
        !we also have ff_sct = sqrt(3. pi / 32 / G / d_c) C_s / (r_0 ) 
        !which just state the ratio of freefall time over sound crossing time 
        !from these 2 formula, rho_c and r_0 are found to be:

        !ph 01/09 new definition entails r_0 instead of r_0 * zeta, the external radius
        r_0 = mass_c_cu / (2.*pi*rap*res_int) * (ff_sct)**2 / (3.*pi/32.) / C_s**2

        if (myid ==1) write(*,*) 'inner radius (pc) ',r_0

        d_c = mass_c_cu / (2.*pi*rap*res_int) / r_0**3

        if(myid ==1) write(*,*) 'central density ',d_c

        ener_therm = 3./2.*mass_c_cu*C_s**2
        ener_grav  = 3./5.*(mass_c_cu**2)/(r_0*zeta)
        ener_grav2=0.
        do i=1,1000
           ener_grav2 = ener_grav2 + (i*zeta/1000.) / (1.+(zeta/1000.*i)**2) * zeta/1000. * mass_rad(i)
        enddo
        ener_grav2 = ener_grav2 * 8.*(pi**2)*(d_c**2)*(r_0**5)



        !angular velocity
        omega = ff_rt * 2.*pi * sqrt( 32.*d_c/3./pi)    
        if (myid==1) write(*,*)'Angular velocity Omega = ',omega/(2.0*pi)/scale_t, 'Hz'
        !central value of magnetic field
        !remember magnetic variable is B/sqrt(4pi)
        if(.not.uniform_bmag)then
           !ph 01/09 new definition entails r_0 instead of r_0 * zeta, the external radius
           B_c = ff_act * sqrt( 32./3./pi) * d_c * r_0 

           mass_sph = d_c / cont * (boxlen*(0.5**levelmin))**3

           !the smallest initial column density
           min_col_d = boxlen * d_c / cont / cont_ic 

           !the largest initial column density
           !obtained by integrating the density distribution through the box
           max_col_d = r_0*d_c*atan(zeta) + (boxlen -2.*r_0*zeta) * d_c / cont / cont_ic 

           if (myid==1) write(*,*) 'valeur du champ magnetique central non normalise B_c', B_c
           if (myid==1) write(*,*) 'valeur du champ magnetique a l exterieur ', B_c*min_col_d/max_col_d

           !calculate the value of mu the mass to flux over critical mass to flux ratio
           !from Mouschovias & Spitzer 1979 M/phi)_crit = 1/(3pi) * sqrt(5/G) * 0.53
           !since B(r)=B_c * sig(r)/sig(0), phi = B_c * mass_c / sig(0)
           !thus mass_c / phi = sig(0) / B_c 
           !taking into account the fact that B_c = champ mag / sqrt(4 pi)
           ! we have in code units mu = sig(0) / (B_c*sqrt(4 pi)) / (sqrt(5)/(3 pi) * 0.53)
           if (myid ==1) write(*,*) 'the mass to flux over critical mass to flux ratio in the case of a spheroidal cloud (not correct if rap ne 1)'
           if (myid ==1) then
              if (B_c.lt.1.e-10) then
                 write(*,*) 'mu= +Infinity'
              else
                 write(*,*) 'mu= ',max_col_d / (B_c*sqrt(4.*pi)) / (sqrt(5.)/(3.*pi) * 0.53)
              endif
           end if
           !note here we make the approximation that max_col_d is equal to the column density through the cloud which is note exactly 
           !the case since the column density of the external medium is also taken into account
        else
           !bc switch to uniform magnetic field if keyword uniform_bmag 
           !calculate the value of mu the mass to flux over critical mass to flux ratio at the core boundary
           !from Mouschovias & Spitzer 1979 M/phi)_crit = 1/(3pi) * sqrt(5/G) * 0.53

           B_c =mass_c_cu*crit_dense_core/(pi*(r_0*zeta)**2)/ (sqrt(5.)/(3.*pi) * 0.53)/sqrt(4.*pi)

           mass_sph = d_c / cont * (boxlen*(0.5**levelmin))**3

           !the smallest initial column density
           min_col_d = boxlen * d_c / cont / cont_ic 

           !the largest initial column density
           !obtained by integrating the density distribution through the box
           max_col_d = r_0*d_c*atan(zeta) + (boxlen -2.*r_0*zeta) * d_c / cont / cont_ic 

           if (myid ==1) write(*,*) 'the mass to flux over critical mass to flux ratio in the case of a spheroidal cloud (not correct if rap ne 1)'
           if (myid ==1) then
              if (B_c.lt.1.e-10) then
                 write(*,*) 'mu= +Infinity'
              else
                 write(*,*) 'mu= ',1.0/crit_dense_core
              endif
           end if

        end if
        !now read the turbulent velocity field used as initial condition
        v_rms=0.
        mass_total=0.
        mass_tot2 =0.
        ener_turb=0.
        ener_rot=0.
        n_size=0
        vx_tot=0.d0
        vy_tot=0.d0
        vz_tot=0.d0
        vx2_tot=0.d0
        vy2_tot=0.d0
        vz2_tot=0.d0

        if(ff_vct .ne. 0)then
           if( myid ==1) write(*,*) 'Read the file which contains the initial turbulent velocity field'
           open(20,file='init_turb.data',form='formatted')
           read(20,*) n_size, ind, seed1,seed2,seed3

           if(n_size .ne. 100) then 
              write(*,*) 'Unexpected field size'
              stop
           endif

           do k=1,n_size
              do j=1,n_size
                 do i=1,n_size
                    read(20,*)xi,yi,zi,vx,vy,vz
                    q_idl(1,i,j,k) = vx
                    q_idl(2,i,j,k) = vy
                    q_idl(3,i,j,k) = vz

                    xi = boxlen*((i-0.5)/n_size-0.5)
                    yi = boxlen*((j-0.5)/n_size-0.5)
                    zi = boxlen*((k-0.5)/n_size-0.5)
                    eli =  (xi/r_0)**2+(yi/r_0)**2+(zi/(r_0*rap))**2

                    if( eli .lt. zeta**2) then

                       vx_tot = vx_tot + d_c/(1.+eli)*vx
                       vy_tot = vy_tot + d_c/(1.+eli)*vy
                       vz_tot = vz_tot + d_c/(1.+eli)*vz

                       vx2_tot = vx2_tot + d_c/(1.+eli)*vx**2
                       vy2_tot = vy2_tot + d_c/(1.+eli)*vy**2
                       vz2_tot = vz2_tot + d_c/(1.+eli)*vz**2

                       ener_turb = ener_turb + d_c/(1.+eli)*(vx**2+vy**2+vz**2)
                       mass_total = mass_total +  d_c / (1.+eli)
                    endif
                 enddo
                 !       eli = (yi/r_0)**2 + (zi/r_0/rap)**2
                 !        if( eli .lt. zeta**2) then
                 !          col_d = r_0*d_c/sqrt(1.+eli)*atan( sqrt( (zeta**2-eli)/(1.+eli) ) )
                 !          mass_tot2 = mass_tot2 + col_d
                 !        endif 
              enddo
           enddo
           close(20)

           vx_tot = vx_tot / mass_total
           vy_tot = vy_tot / mass_total
           vz_tot = vz_tot / mass_total

           vx2_tot = vx2_tot / mass_total
           vy2_tot = vy2_tot / mass_total
           vz2_tot = vz2_tot / mass_total

           v_rms = sqrt( vx2_tot-vx_tot**2 + vy2_tot-vy_tot**2 + vz2_tot-vz_tot**2 ) 

           mass_total = mass_total*(boxlen/n_size)**3
           if (myid ==1) write(*,*) 'We verify the calculation for the mass. The 2 following values must be very close:'
           if (myid ==1) write(*,*) 'mass_total, mass_c_cu ',mass_total, mass_c_cu !,mass_tot2

        end if

        n_size=100
        do k=1,n_size
           do j=1,n_size
              do i=1,n_size

                 xi = boxlen*((i-0.5)/n_size-0.5)
                 yi = boxlen*((j-0.5)/n_size-0.5)
                 zi = boxlen*((k-0.5)/n_size-0.5)
                 eli =  (xi/r_0)**2+(yi/r_0)**2+(zi/(r_0*rap))**2

                 if( eli .lt. zeta**2) then
                    ener_rot = ener_rot + d_c/(1.+eli) * omega**2 * (yi**2 + zi**2)
                 endif
              enddo
           enddo
        enddo


        ener_rot  = 0.5 * ener_rot*(boxlen/n_size)**3
        ener_turb = 0.5 * ener_turb*(boxlen/n_size)**3

        !estimate of the thermal over gravitational energy 
        if (myid == 1) write(*,*) 'estimate (uniform density is assumed) of the ratio of thermal over gravitational energy'
        if (myid == 1) write(*,*)  ener_therm / ener_grav

        if (myid == 1) write(*,*) 'good estimate of the ratio of thermal over gravitational energy'
        if (myid == 1) write(*,*)  ener_therm / ener_grav2

        if (myid .eq. 1) write(*,*) 'estimate of the rotational over gravitational energy ratio'
        if (myid .eq. 1) write(*,*) 'ener_rot/ener_grav2 ', ener_rot / ener_grav2, ener_rot , ener_grav2


        !calculate now the coefficient by which the turbulence velocity needs
        !to be multiplied 
        if(v_rms .ne.0)then
           if (myid .eq. 1) write(*,*) 'vrms non norm ',v_rms

           !ph 01/09 new definition entails r_0 instead of r_0 * zeta, the external radius
           v_rms = ff_vct * sqrt(32.*d_c/3./pi)*r_0 / v_rms

           if (myid .eq. 1) write(*,*) 'vrms mult ',v_rms


           !estimate of the turbulent over gravitational energy ratio
           if (myid .eq. 1) write(*,*) 'estimate of the turbulent over gravitational energy ratio'
           if (myid .eq. 1) write(*,*) 'ener_turb/ener_grav2 ', ener_turb*(v_rms**2) / ener_grav2

        end if

100     format(i5,4e12.5)
101     format(6e12.5)
102     format(i5)

        if (myid ==1)  write(*,*) 'Reading achieved'
        first = .false.
     endif

     Temp = T_eos

     DO i=1,nn

        x(i,1) = x(i,1) - 0.5*boxlen
        x(i,2) = x(i,2) - 0.5*boxlen
        x(i,3) = x(i,3) - 0.5*boxlen

        q(i,2) = 0.0d0
        q(i,3) = 0.0d0
        q(i,4) = 0.0d0

        !initialise the density field
        eli =  (x(i,1)/r_0)**2+(x(i,2)/r_0)**2+(x(i,3)/(r_0*rap))**2

        if( eli .gt. zeta**2) then 
           q(i,1) = d_c / cont / cont_ic
           if(eos)then
              call enerint_eos(q(i,1),Temp,ee)
              q(i,5) = ee
              !if the cloud is in pressure equilibrium with the surrounding medium
              !remove this line if the IC gas is isothermal as well 
              !        q(i,5) = q(i,5) * cont_ic 
              q(i,nvar) = ee
           else
              q(i,5) = q(i,1) * C_s**2/(gamma-1.0d0)
              q(i,nvar) = q(i,5)
           endif
        else
           q(i,1) = d_c / (1.+eli)
           if(eos)then
              call enerint_eos(q(i,1),Temp,ee)
              q(i,5   ) = ee
              q(i,nvar) = ee
           else
              q(i,5) = q(i,1) * C_s**2/(gamma-1.0d0)
              q(i,nvar) = q(i,5)
           endif
        endif

#if NGRP>0
        do igroup=1,ngrp
           q(i,firstindex_er+igroup) = radiation_source(Temp,igroup)/(scale_d*scale_v**2)
        enddo
#endif

        if(v_rms .ne. 0)then
           !initialise the turbulent velocity field
           !make a zero order interpolation (should be improved)
           ind_i = int((x(i,1)/boxlen+0.5)*n_size)+1
           ind_j = int((x(i,2)/boxlen+0.5)*n_size)+1
           ind_k = int((x(i,3)/boxlen+0.5)*n_size)+1


           if( ind_i .lt. 1 .or. ind_i .gt. n_size) write(*,*) 'ind_i ',ind_i,boxlen,x(i,1),n_size
           if( ind_j .lt. 1 .or. ind_j .gt. n_size) write(*,*) 'ind_j ',ind_j
           if( ind_k .lt. 1 .or. ind_k .gt. n_size) write(*,*) 'ind_k ',ind_k

           q(i,2) = v_rms*(q_idl(1,ind_i,ind_j,ind_k)-vx_tot)
           q(i,3) = v_rms*(q_idl(2,ind_i,ind_j,ind_k)-vy_tot)
           q(i,4) = v_rms*(q_idl(3,ind_i,ind_j,ind_k)-vz_tot)

        end if

        !add  rotation. x cos(theta_mag_radians) + y sin(theta_mag_radians) is the rotation axis
        sph =  (x(i,1)/r_0)**2+(x(i,2)/r_0)**2+(x(i,3)/(r_0))**2
        if( sph .lt. (zeta*rap)**2 ) then 

           !to check these formulae one can verify that those arrays are perpendicular
           !with (cos(theta_mag_radians),sin(theta_mag_radians),0) and that the norm of the vectorial product of
           ! (cos(theta_mag_radians),sin(theta_mag_radians),0) by the above arrays is equal to the distance
           !  (x sin(thet)-y cos(thet))^2 + z^2 
           q(i,2) = q(i,2) - (omega*x(i,3)*sin(theta_mag_radians))
           q(i,3) = q(i,3) + (omega*x(i,3)*cos(theta_mag_radians))
           q(i,4) = q(i,4) + (omega*(x(i,1)*sin(theta_mag_radians)-x(i,2)*cos(theta_mag_radians)))

        endif


     ENDDO

     if(uniform_bmag) then

        DO i=1,nn
           q(i,6     ) = B_c
           q(i,nvar+1) = q(i,6)
           
           !By component
           q(i,7     ) = 0.
           q(i,nvar+2) = 0.
           
           !Bz component
           q(i,8     ) = 0.
           q(i,nvar+3) = 0.

        end DO
           
     else

        dxmin=boxlen*0.5d0**(levelmin+3)
        
        if( dx .lt. dxmin) then 
           write(*,*) 'dxmin too large'
           write(*,*) 'dx ',dx/boxlen
           write(*,*) 'dxmin ',dxmin/boxlen
           stop
        endif

        nticks=dx/dxmin



        DO i=1,nn
           q(i,6)=0.
           
           xl=x(i,1)-0.5*dx
           yl=x(i,2)-0.5*dx
           zl=x(i,3)-0.5*dx
           
           
           !the magnetic field in cells must be subdivided in order to insure that the magnetic
           !flux is the same in coarse and refined grids
           DO jj=1,nticks
              DO kk=1,nticks
                 
                 yy=yl+(dble(jj)-0.5d0)*dxmin
                 zz=zl+(dble(kk)-0.5d0)*dxmin

                 !this formula comes from the integration of the density distribution along x
                 eli = (yy/r_0)**2 + (zz/r_0/rap)**2
                 if( eli .lt. zeta**2) then
                    col_d = r_0*d_c/sqrt(1.+eli)*atan( sqrt( (zeta**2-eli)/(1.+eli) ) )
                    col_d = max(col_d,min_col_d)
                 else 
                    col_d = min_col_d
                 endif
                 
                 !Bx component
                 if(uniform_bmag)then
                    q(i,6     ) = q(i,6) + B_c
                 else
                    q(i,6     ) = q(i,6) + B_c * col_d / max_col_d
                 end if
                 q(i,nvar+1) = q(i,6)

                 !By component
                 q(i,7     ) = 0.
                 q(i,nvar+2) = 0.
                 
                 !Bz component
                 q(i,8     ) = 0.
                 q(i,nvar+3) = 0.

              ENDDO
           ENDDO

           q(i,6:8)           = q(i,6:8)           / dble(nticks)**2

           !new version rotates the rotation velocity 
           !rotates the magnetic field of an angle theta
           !       bx=q(i,6)
           !       by=q(i,7)
           !       q(i,6) =  bx*cos(theta_mag_radians) + by*sin(theta_mag_radians)
           !       q(i,7) =  bx*sin(theta_mag_radians) - by*cos(theta_mag_radians)
           
           q(i,nvar+1:nvar+3) = q(i,6:8)
        ENDDO
     end if

endif

end subroutine collapse_condinit

!================================================================
!================================================================
!================================================================
!================================================================
subroutine nimhd_diffusion_condinit(x,q,dx,nn)
  use amr_parameters
  use hydro_parameters
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar_all)::q ! Primitive variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for the NIMHD diffusion test
  !================================================================
  integer::i
  real(dp) :: xx,yy,zz

  ! Call built-in initial condition generator
  call region_condinit(x,q,dx,nn)
  
  do i = 1,nn
     xx = x(i,1) - 0.5_dp*(1.0_dp + dx)
     yy = x(i,2) - 0.5_dp*(1.0_dp + dx)
     zz = x(i,3) - 0.5_dp*(1.0_dp + dx)
     if((abs(xx) < 0.9_dp*dx) .and. (abs(zz) < 0.9_dp*dx))then
       q(i,     7) = 1.0_dp
       q(i,nvar+2) = 1.0_dp
     endif
  enddo

end subroutine nimhd_diffusion_condinit
!================================================================
!================================================================
!================================================================
!================================================================
subroutine velana(x,v,dx,t,ncell)
  use amr_parameters
  use hydro_parameters
  implicit none
  integer ::ncell                         ! Size of input arrays
  real(dp)::dx                            ! Cell size
  real(dp)::t                             ! Current time
  real(dp),dimension(1:nvector,1:3)::v    ! Velocity field
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine computes the user defined velocity fields.
  ! x(i,1:ndim) are cell center position in [0,boxlen] (user units).
  ! v(i,1:3) is the imposed 3-velocity in user units.
  !================================================================
  integer::i
  real(dp)::xx,yy=0.,zz=0.,vx,vy,vz,aa,twopi
!!$  real(dp)::rr,tt,omega

  ! Add here, if you wish, some user-defined initial conditions
  aa=1.0+0.*t
  twopi=2d0*ACOS(-1d0)
  do i=1,ncell

     xx=x(i,1)+0.*dx
#if NDIM > 1
     yy=x(i,2)
#endif
#if NDIM > 2
     zz=x(i,3)
#endif
     ! ABC
     vx=aa*(cos(twopi*yy)+sin(twopi*zz))
     vy=aa*(sin(twopi*xx)+cos(twopi*zz))
     vz=aa*(cos(twopi*xx)+sin(twopi*yy))

!!$     ! 1D advection test
!!$     vx=1.0_dp
!!$     vy=0.0_dp
!!$     vz=0.0_dp

!!$     ! Ponomarenko
!!$     xx=xx-boxlen/2.0
!!$     yy=yy-boxlen/2.0
!!$     rr=sqrt(xx**2+yy**2)
!!$     if(yy>0)then
!!$        tt=acos(xx/rr)
!!$     else
!!$        tt=-acos(xx/rr)+twopi
!!$     endif
!!$     if(rr<1.0)then
!!$        omega=0.609711
!!$        vz=0.792624
!!$     else
!!$        omega=0.0
!!$        vz=0.0
!!$     endif
!!$     vx=-sin(tt)*rr*omega
!!$     vy=+cos(tt)*rr*omega

     v(i,1)=vx
#if NDIM > 1
     v(i,2)=vy
#endif
#if NDIM > 2
     v(i,3)=vz
#endif
  end do


end subroutine velana
!================================================================
!================================================================
!================================================================
!================================================================
subroutine ponomarenko_condinit(x,u,dx,nn)
  use amr_parameters
  use hydro_parameters
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar+3)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:4): d.u,d.v,d.w, U(i,5): E, U(i,6:8): Bleft,
  ! U(i,nvar+1:nvar+3): Bright
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:4):u,v,w, Q(i,5): P, Q(i,6:8): Bleft,
  ! Q(i,nvar+1:nvar+3): Bright
  ! If nvar > 8, remaining variables (9:nvar) are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
  integer::ivar,i
  real(dp),dimension(1:nvector,1:nvar+3),save::q   ! Primitive variables
  real(dp)::xc,xr,xl,yl,yr,yc,al,ar,r0,a0,twopi,rr,ss,tt

  ! density
  q(1:nn,1)=1.0

  ! velocity
  q(1:nn,2:4)=0.0
  R0=1.0
  do i=1,nn
     xc=x(i,1)-boxlen/2.0
     yc=x(i,2)-boxlen/2.0
     rr = sqrt(xc**2+yc**2)
     if(rr < R0)q(i,4)=1.0
  end do

  ! pressure
  q(1:nn,5)=1.0

  ! magnetic field
  call mag_screw(x,q,dx,nn)

  ! Convert primitive to conservative variables
  ! density -> density
  u(1:nn,1)=q(1:nn,1)
  ! velocity -> momentum
  u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
  u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
  u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
  ! kinetic energy
  u(1:nn,5)=0.0d0
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,2)**2
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,3)**2
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,4)**2
  ! pressure -> total fluid energy
  u(1:nn,5)=u(1:nn,5)+q(1:nn,5)/(gamma-1.0d0)
  ! magnetic energy -> total fluid energy
  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,6)+q(1:nn,nvar+1))**2
  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,7)+q(1:nn,nvar+2))**2
  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,8)+q(1:nn,nvar+3))**2
  u(1:nn,6:8)=q(1:nn,6:8)
  u(1:nn,nvar+1:nvar+3)=q(1:nn,nvar+1:nvar+3)
  ! passive scalars
#if NDIM > 8
  do ivar=9,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do
#endif

end subroutine ponomarenko_condinit
!================================================================
!================================================================
!================================================================
!================================================================
subroutine mag_screw(x,q,dx,nn)
  use amr_parameters
  use hydro_parameters, only: nvar
  use const
  implicit none

  real(dp),dimension(1:nvector,1:ndim)::x    ! Cell center position
  real(dp),dimension(1:nvector,1:nvar+3)::q  ! Primitive variables
  integer::nn                                ! Number of cells
  real(dp)::dx                               ! Cell size
  real(dp)::xx,yy,zz,A_screw
  integer::i,it,nticks
  real(dp)::dxmin,zfl,zceil,rtick,x_edge,y_edge
  real(dp),dimension(1:2,1:2,1:2)::A_mag     ! x/y,left/right,up/down
  real(dp)::dmax,dnfw

  ! Ponomarenko
  dxmin=boxlen*0.5d0**nlevelmax
  nticks=nint(dx/dxmin)

  do i=1,nn
    ! box-centered coordinates
    xx=x(i,1) - boxlen * 0.5d0
    yy=x(i,2) - boxlen * 0.5d0
    zz=x(i,3) - boxlen * 0.5d0

    zfl   = zz - 0.5d0*dx
    zceil = zz + 0.5d0*dx

    A_mag = 0.0

    ! A_(x,l)
    x_edge = xx + 0.5*(dxmin-dx)
    y_edge = yy - 0.5*dx
    do it=1,nticks
      ! floor
      A_mag(1,1,1) = A_mag(1,1,1) + A_screw(x_edge,y_edge,zfl,1)
      ! ceiling
      A_mag(1,1,2) = A_mag(1,1,2) + A_screw(x_edge,y_edge,zceil,1)
      ! advance along x
      x_edge=x_edge + dxmin
    end do

    ! A_(x,r)
    x_edge = xx + 0.5*(dxmin-dx)
    y_edge = yy + 0.5*dx
    do it=1,nticks
      ! floor
      A_mag(1,2,1) = A_mag(1,2,1) + A_screw(x_edge,y_edge,zfl,1)
      ! ceiling
      A_mag(1,2,2) = A_mag(1,2,2) + A_screw(x_edge,y_edge,zceil,1)
      ! advance along x
      x_edge=x_edge + dxmin
    end do

    ! A_(y,l)
    x_edge = xx - 0.5*dx
    y_edge = yy + 0.5*(dxmin-dx)
    do it=1,nticks
      ! floor
      A_mag(2,1,1) = A_mag(2,1,1) + A_screw(x_edge,y_edge,zfl,2)
      ! ceiling
      A_mag(2,1,2) = A_mag(2,1,2) + A_screw(x_edge,y_edge,zceil,2)
      ! advance along y
      y_edge=y_edge + dxmin
    end do

    ! A_(y,r)
    x_edge = xx + 0.5*dx
    y_edge = yy + 0.5*(dxmin-dx)
    do it=1,nticks
      ! floor
      A_mag(2,2,1) = A_mag(2,2,1) + A_screw(x_edge,y_edge,zfl,2)
      ! ceiling
      A_mag(2,2,2) = A_mag(2,2,2) + A_screw(x_edge,y_edge,zceil,2)
      ! advance along y
      y_edge=y_edge + dxmin
    end do

    ! average value
    A_mag = A_mag / DBLE(nticks)

    ! B left
    ! X direction
    q(i,6)= - (A_mag(2,1,2)-A_mag(2,1,1)) / dx ! [-d/dz A_y]
    ! Y direction
    q(i,7)=   (A_mag(1,1,2)-A_mag(1,1,1)) / dx ! [ d/dz A_x]
    ! Z direction  [d/dx A_y - d/dy A_x]
    q(i,8)=   ( A_mag(2,2,1)-A_mag(2,1,1) - A_mag(1,2,1)+A_mag(1,1,1) ) / dx

    ! B right
    ! X direction
    q(i,nvar+1)= - (A_mag(2,2,2)-A_mag(2,2,1)) / dx ! [-d/dz A_y]
    ! Y direction
    q(i,nvar+2)=   (A_mag(1,2,2)-A_mag(1,2,1)) / dx ! [ d/dz A_x]
    ! Z direction  [d/dx A_y - d/dy A_x]
    q(i,nvar+3)=   ( A_mag(2,2,2)-A_mag(2,1,2) - A_mag(1,2,2)+A_mag(1,1,2) ) / dx
  end do

end subroutine mag_screw
!================================================================
!================================================================
!================================================================
!================================================================
function A_screw(xx,yy,zz,dir)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none
  real(kind=8)::xx,yy,zz,A_screw
  real(kind=8)::A,rr,tt,twopi
  integer::dir

  twopi=2d0*acos(-1d0)

  rr=sqrt(xx**2+yy**2)
  if(yy>0)then
     tt=acos(xx/rr)
  else
     tt=-acos(xx/rr)+twopi
  endif
  A_screw=0d0
  if(rr<2.0)then
     A_screw=1d-1*(2d0-rr)*rr*cos(tt)*sin(twopi*(zz+boxlen/2)/boxlen)
     if(dir==1)A_screw=A_screw*cos(tt)
     if(dir==2)A_screw=A_screw*sin(tt)
  endif

  return
end function A_screw
!================================================================
!================================================================
!================================================================
!================================================================
subroutine velana_ponomarenko(x,v,dx,t,ncell)
  use amr_parameters
  use hydro_parameters
  implicit none
  integer ::ncell                         ! Size of input arrays
  real(dp)::dx                            ! Cell size
  real(dp)::t                             ! Current time
  real(dp),dimension(1:nvector,1:3)::v    ! Velocity field
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine computes the user defined velocity fields.
  ! x(i,1:ndim) are cell center position in [0,boxlen] (user units).
  ! v(i,1:3) is the imposed 3-velocity in user units.
  !================================================================
  integer::i
  real(dp)::xx,yy,zz,vx,vy,vz,rr,tt,omega,aa,twopi

  ! Add here, if you wish, some user-defined initial conditions
  aa=1.0
  twopi=2d0*ACOS(-1d0)
  do i=1,ncell

     xx=x(i,1)
#if NDIM > 1
     yy=x(i,2)
#endif
#if NDIM > 2
     zz=x(i,3)
#endif
     ! Ponomarenko
     xx=xx-boxlen/2.0
     yy=yy-boxlen/2.0
     rr=sqrt(xx**2+yy**2)
     if(yy>0)then
        tt=acos(xx/rr)
     else
        tt=-acos(xx/rr)+twopi
     endif
     if(rr<1.0)then
        omega=0.609711
        vz=0.792624
     else
        omega=0.0
        vz=0.0
     endif
     vx=-sin(tt)*rr*omega
     vy=+cos(tt)*rr*omega

     v(i,1)=vx
#if NDIM > 1
     v(i,2)=vy
#endif
#if NDIM > 2
     v(i,3)=vz
#endif
  end do


end subroutine velana_ponomarenko
!================================================================
!================================================================
!================================================================
!================================================================
subroutine velana_nimhd_diffusion(x,v,dx,t,ncell)
  use amr_parameters
  use hydro_parameters
  implicit none
  integer ::ncell                         ! Size of input arrays
  real(dp)::dx                            ! Cell size
  real(dp)::t                             ! Current time
  real(dp),dimension(1:nvector,1:3)::v    ! Velocity field
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine computes the user defined velocity fields.
  ! x(i,1:ndim) are cell center position in [0,boxlen] (user units).
  ! v(i,1:3) is the imposed 3-velocity in user units.
  !================================================================
  integer::i
  real(dp)::aa,twopi

  ! Add here, if you wish, some user-defined initial conditions
  aa=1.0
  twopi=2d0*ACOS(-1d0)
  v=0.0d0

end subroutine velana_nimhd_diffusion