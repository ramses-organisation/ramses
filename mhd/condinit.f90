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
  ! Q(i,nvar+1:nvar+3): Bright
  ! If nvar > 8, remaining variables (9:nvar) are treated as passive
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

  case('ponomarenko')
      call ponomarenko_condinit(x, q, dx, nn)

  case('collapse')
     call collapse_condinit(x, q, dx, nn)

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
  u(1:nn,neul)=u(1:nn,neul)+q(1:nn,neul)/(gamma-1.0d0)
  ! magnetic energy -> total fluid energy
  u(1:nn,neul)=u(1:nn,neul)+0.125d0*(q(1:nn,6)+q(1:nn,nvar+1))**2
  u(1:nn,neul)=u(1:nn,neul)+0.125d0*(q(1:nn,7)+q(1:nn,nvar+2))**2
  u(1:nn,neul)=u(1:nn,neul)+0.125d0*(q(1:nn,8)+q(1:nn,nvar+3))**2
  u(1:nn,6:8)=q(1:nn,6:8)
  u(1:nn,nvar+1:nvar+3)=q(1:nn,nvar+1:nvar+3)
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
subroutine collapse_condinit(x,q,dx,nn)
  use amr_commons, only:myid
  use amr_parameters
  use hydro_commons
  use poisson_parameters
  use constants, only:mH,kB,M_sun,pc2cm
  implicit none
  integer ::nn                              ! Number of cells
  real(dp)::dx                              ! Cell size
  real(dp),dimension(1:nvector,1:nvar+3)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:ndim+1): d.u,d.v,d.w and U(i,ndim+2): E.
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:ndim+1):u,v,w and Q(i,ndim+2): P.
  ! If nvar >= ndim+3, remaining variables are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
  integer :: i,j,k,id,iu,iv,iw,ip
  real(dp):: x0,y0,z0,rc,rs,xx,yy,zz,pi,r0,d0,B0,p0,omega0,mass_c_cu,scale_m
  integer :: ivar, np
  real(dp),dimension(1:nvector,1:nvar+3)::q   ! Primitive variables
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:3,1:3):: rot_M,rot_invM,rot_tilde
  real(dp):: theta_mag_radians


  logical,save:: first=.true.
  real(dp),dimension(1:3,1:100,1:100,1:100),save::q_idl
  real(dp),save::vx_tot,vy_tot,vz_tot,vx2_tot,vy2_tot,vz2_tot
  integer,save:: n_size
  integer:: ind_i, ind_j, ind_k
  real(dp),save:: ind,seed1,seed2,seed3,xi,yi,zi,vx,vy,vz
  real(dp),save:: C_s,v_rms
  integer, save :: count_vrms

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
     q(i,8     ) = B0
     q(i,nvar+3) = B0

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
  ENDDO

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
  !kinetic + magnetic energy
  u(1:nn,5)=u(1:nn,5)+0.125*(q(1:nn,6)+q(1:nn,nvar+1))**2
  u(1:nn,5)=u(1:nn,5)+0.125*(q(1:nn,7)+q(1:nn,nvar+2))**2
  u(1:nn,5)=u(1:nn,5)+0.125*(q(1:nn,8)+q(1:nn,nvar+3))**2
  ! pressure -> total fluid energy
  u(1:nn,5)=u(1:nn,5)+q(1:nn,5)/(gamma-1.0d0)
  ! magnetic field
  u(1:nn,6:8)=q(1:nn,6:8)
  u(1:nn,nvar+1:nvar+3)=q(1:nn,nvar+1:nvar+3)
  ! passive scalars
  do ivar=9,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do

end subroutine collapse_condinit
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
