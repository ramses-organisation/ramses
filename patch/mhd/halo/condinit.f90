!================================================================
!================================================================
!================================================================
!================================================================
module halo_parameters!{{{
  use amr_commons

  character(len=16)::mag_topology='constant' ! magnetic topology:
                                             ! 'constant' (along x-axis)
                                             ! 'toroidal'
                                             ! 'dipole'
                                             ! 'quadrupole'

end module halo_parameters!}}}

module halo_commons

  use halo_parameters

end module halo_commons

subroutine read_halo_params! {{{
  use halo_parameters
  implicit none

#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  logical::nml_ok=.true.
  character(LEN=80)::infile

  !--------------------------------------------------
  ! Local variables
  !--------------------------------------------------
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2

  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------
  namelist/halo_params/ mag_topology


  CALL getarg(1,infile)
  open(1,file=infile)
  rewind(1)
  read(1,NML=halo_params,END=106)
  goto 107
  106 write(*,*)' You need to set up namelist &HALO_PARAMS in parameter file'
  call clean_stop
  107 continue
  close(1)

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  !--------------------------------------------------------------
  ! This section deals with the halo initial conditions
  !--------------------------------------------------------------

  ! Magnetic Topology
  select case (mag_topology)
    case ('constant')
      if(myid==1) write(*,*) "Chosen magnetic topology: 'constant'"
    case ('toroidal')
      if(myid==1) write(*,*) "Chosen magnetic topology: 'toroidal'"
    case ('dipole')
      if(myid==1) write(*,*) "Chosen magnetic topology: 'dipole'"
    case ('quadrupole')
      if(myid==1) write(*,*) "Chosen magnetic topology: 'quadrupole'"
    case default
      if(myid==1) write(*,*) "Magnetic field topology not recognised."
      nml_ok=.false.
  end select

  if(.not. nml_ok)then
    if(myid==1)write(*,*)'Too many errors in the namelist'
    if(myid==1)write(*,*)'Aborting...'
    call clean_stop
  end if

end subroutine read_halo_params
! }}}


subroutine condinit(x,u,dx,nn)
  use amr_parameters
  use hydro_parameters
  use poisson_parameters, ONLY: gravity_params
  use halo_commons
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
  integer::ivar,i,it,nticks
  real(dp),dimension(1:nvector,1:nvar+3),save::q   ! Primitive variables
  real(dp)::rx,ry,rz,rr,x0,y0,z0,eps,rcut,lr,lrmax,lrcut
  real(dp)::pi,v200,spin,c,fgas,dmax,mmax,Bnfw,Bz,dxmin
  real(dp)::fcut,vcirc,vrot,pnfw,dnfw,mnfw,fc,cosphi
  real(dp)::Axll,Axlr,Axrl,Axrr
  real(dp)::Ayll,Aylr,Ayrl,Ayrr
  real(dp)::xl,xr,yl,yr,zl,zr,xc,yc,zc
  real(dp)::xxl,xxr,yyl,yyr,zzl,zzr,xx,yy
  logical, save:: init_nml=.false.

  ! Add here, if you wish, some user-defined initial conditions

  ! Read user-defined parameters from the namelist
  if (.not. init_nml) then
    call read_halo_params
    init_nml = .true.
  end if

  eps=2.*boxlen*0.5**nlevelmax
  x0=boxlen/2.0
  y0=boxlen/2.0
  z0=boxlen/2.0
  rcut=boxlen/2.0*0.75

  v200=gravity_params(1)
  c   =gravity_params(2)
  spin=gravity_params(3)
  fgas=gravity_params(4)
  dmax=1d0/eps/(1d0+eps)**2
  pi  =acos(-1d0)

  mmax=4d0*pi*(log(1d0+rcut)-rcut/(1d0+rcut))
  vcirc=sqrt(4.*pi*(log(1d0+c)-c/(1d0+c))/c)
  fc=c*(0.5-0.5/(1d0+c)**2-log(1d0+c)/(1d0+c))/(log(1d0+c)-c/(1d0+c))**2
  Bnfw=B_ave*vcirc ! Scaled magnetic field

  dxmin=boxlen*0.5d0**nlevelmax
  nticks=max(ceiling(dx/dxmin),1)

  lrmax=log(10*boxlen)
  lrcut=log(rcut)
  do i=1,nn

     xl=x(i,1)-0.5*dx-x0
     xr=x(i,1)+0.5*dx-x0
     xc=x(i,1)-x0
     yl=x(i,2)-0.5*dx-y0
     yr=x(i,2)+0.5*dx-y0
     yc=x(i,2)-y0
     zl=x(i,3)-0.5*dx-z0
     zr=x(i,3)+0.5*dx-z0
     zc=x(i,3)-z0

     rr=sqrt(xc**2+yc**2+zc**2)
     fcut=exp(-(rr/rcut)**50)
     dnfw=1d0/rr/(1d0+rr)**2
     dnfw=MIN(dmax,dnfw)
     dnfw=MAX(dnfw*fcut,1d-6)
     mnfw=(log(1d0+rr)-rr/(1d0+rr))/(log(1d0+c)-c/(1d0+c))
     vrot=2d0*sqrt(2d0/fc)*spin*vcirc*c*mnfw/rr*fcut
     cosphi=sqrt(xc**2+yc**2)/rr
     lr=MIN(log(rr),lrcut)
     pnfw=romberg(lr,lrcut,1d-4)+mmax/MAX(rr,rcut)*1d-6

     ! Euler variables
     q(i,1)=fgas*dnfw
     q(i,2)=-vrot/rr*yc
     q(i,3)=+vrot/rr*xc
     q(i,4)=0.0
     q(i,5)=fgas*pnfw

     ! Passive scalar
     if(metal)q(i,9)=0.

  enddo

  ! Magnetic Topology
  select case (mag_topology)
    case ('constant')
      call mag_constant(q,nn,Bnfw)
    case ('toroidal')
      call mag_toroidal(x,q,dx,nn,Bnfw)
    case ('dipole')
      call mag_dipole(x,q,dx,nn,Bnfw)
    case ('quadrupole')
      call mag_quadrupole(x,q,dx,nn,Bnfw)
  end select

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
  ! Magnetic field components
  u(1:nn,6:8)=q(1:nn,6:8)
  u(1:nn,nvar+1:nvar+3)=q(1:nn,nvar+1:nvar+3)
  ! passive scalars
  do ivar=9,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do

contains
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function fffy(rint)
    implicit none
    !      Computes the integrand
    real(dp)::fffy
    real(dp)::rint,rrr,M,rho,vphi
    rrr=exp(rint)
    rho=1d0/rrr/(1d0+rrr)**2
    rho=MIN(rho,dmax)
    rho=MAX(rho*exp(-(rrr/rcut)**50),1d-6)
    M=4d0*pi*(log(1d0+rrr)-rrr/(1d0+rrr))

    vphi=2d0*sqrt(2d0/fc)*spin*vcirc*c/rrr &
         & *exp(-(rrr/rcut)**50) &
         & *(log(1d0+rrr)-rrr/(1d0+rrr))/(log(1d0+c)-c/(1d0+c))
    vphi=vphi*cosphi

    fffy=M/rrr*rho-vphi**2*rho
    return
  end function fffy
  !cccccccccccccccccccccccccccccccccccccccccccccccccccc
  function romberg(a,b,tol)
    implicit none
    real(dp)::romberg
    !
    !     Romberg returns the integral from a to b of f(x)dx using Romberg
    !     integration. The method converges provided that f(x) is continuous
    !     in (a,b). The function f must be double precision and must be
    !     declared external in the calling routine.
    !     tol indicates the desired relative accuracy in the integral.
    !
    integer::maxiter=16,maxj=5
    real(dp),dimension(100):: g
    real(dp)::a,b,tol,fourj
    real(dp)::h,error,gmax,g0,g1
    integer::nint,i,j,k,jmax

    h=0.5d0*(b-a)
    gmax=h*(fffy(a)+fffy(b))
    g(1)=gmax
    nint=1
    error=1.0d20
    i=0
10  i=i+1
    if(.not.  (i>maxiter.or.(i>5.and.abs(error)<tol)))then
       !  Calculate next trapezoidal rule approximation to integral.

       g0=0.0d0
       do k=1,nint
          g0=g0+fffy(a+(k+k-1)*h)
       end do
       g0=0.5d0*g(1)+h*g0
       h=0.5d0*h
       nint=nint+nint
       jmax=min(i,maxj)
       fourj=1.0d0

       do j=1,jmax
          ! Use Richardson extrapolation.
          fourj=4.0d0*fourj
          g1=g0+(g0-g(j))/(fourj-1.0d0)
          g(j)=g0
          g0=g1
       enddo
       if (abs(g0).gt.tol) then
          error=1.0d0-gmax/g0
       else
          error=gmax
       end if
       gmax=g0
       g(jmax+1)=g0
       go to 10
    end if
    romberg=g0

    if (i>maxiter.and.abs(error)>tol) &
         &    write(*,*) 'Romberg failed to converge; integral, error=', &
         &    romberg,error



    return
  end function romberg

end subroutine condinit

subroutine mag_constant(q,nn,B_0)
  use amr_parameters
  use hydro_parameters, only: nvar
  implicit none

  real(dp),dimension(1:nvector,1:nvar+3)::q  ! Primitive variables
  integer::nn                                ! Number of cells
  real(dp)::B_0                              ! Default B-strength

  ! B-left
  q(1:nn,6) = B_0
  q(1:nn,7) = 0.0
  q(1:nn,8) = 0.0

  ! B-right
  q(1:nn,nvar+1) = B_0
  q(1:nn,nvar+2) = 0.0
  q(1:nn,nvar+3) = 0.0

end subroutine mag_constant

subroutine mag_toroidal(x,q,dx,nn,B_0)
  use amr_parameters
  use hydro_parameters, only: nvar
  use const
  implicit none

  real(dp),dimension(1:nvector,1:ndim)::x    ! Cell center position
  real(dp),dimension(1:nvector,1:nvar+3)::q  ! Primitive variables
  integer::nn                                ! Number of cells
  real(dp)::dx                               ! Cell size
  real(dp)::B_0                              ! Default B-strength
  real(dp)::xx,yy,zz
  integer::i,it,nticks
  real(dp)::dxmin,Al,Ar,ztick,rc,rs,Bscale
  real(dp)::dmax,dnfw

  ! Toroidal field, prop. to rho**(2/3)

  dxmin=boxlen*0.5d0**nlevelmax
  nticks=nint(dx/dxmin)

  dmax=1d0/(2.0d0*dxmin)/(1d0+(2.0d0*dxmin))**2

  do i=1,nn
    ! box-centered coordinates
    xx=x(i,1) - boxlen * 0.5d0
    yy=x(i,2) - boxlen * 0.5d0
    zz=x(i,3) - boxlen * 0.5d0

    Bscale = B_0/DBLE(nticks)/dx

    ! B left
    ! X direction
    Al = 0.0
    Ar = 0.0
    ztick = zz + 0.5d0*(dxmin-dx)
    do it=1,nticks
      ! left
      rs = SQRT( (xx - 0.5d0*dx)**2 + (yy - 0.5d0*dx)**2 + ztick**2 )
      dnfw=1d0/rs/(1d0+rs)**2
      dnfw=MIN(dmax,dnfw)
      rc = SQRT( (xx - 0.5d0*dx)**2 + (yy - 0.5d0*dx)**2 )
      Al= Al + rc * dnfw**two3rd

      ! right
      rs = SQRT( (xx - 0.5d0*dx)**2 + (yy + 0.5d0*dx)**2 + ztick**2 )
      dnfw=1d0/rs/(1d0+rs)**2
      dnfw=MIN(dmax,dnfw)
      rc = SQRT( (xx - 0.5d0*dx)**2 + (yy + 0.5d0*dx)**2 )
      Ar= Ar + rc * dnfw**two3rd

      ! advance in z direction
      ztick=ztick + dxmin
    end do

    q(i,6)=Bscale * (Ar-Al)

    ! B left
    ! Y direction
    Al = 0.0
    Ar = 0.0
    ztick = zz + 0.5d0*(dxmin-dx)
    do it=1,nticks
      ! left
      rs = SQRT( (xx - 0.5d0*dx)**2 + (yy - 0.5d0*dx)**2 + ztick**2 )
      dnfw=1d0/rs/(1d0+rs)**2
      dnfw=MIN(dmax,dnfw)
      rc = SQRT( (xx - 0.5d0*dx)**2 + (yy - 0.5d0*dx)**2 )
      Al= Al + rc * dnfw**two3rd

      ! right
      rs = SQRT( (xx + 0.5d0*dx)**2 + (yy - 0.5d0*dx)**2 + ztick**2 )
      dnfw=1d0/rs/(1d0+rs)**2
      dnfw=MIN(dmax,dnfw)
      rc = SQRT( (xx + 0.5d0*dx)**2 + (yy - 0.5d0*dx)**2 )
      Ar= Ar + rc * dnfw**two3rd

      ! advance in z direction
      ztick=ztick + dxmin
    end do

    q(i,7)=Bscale * (Al-Ar)

    ! B left
    ! Z direction
    q(i,8)=0.0

    ! B right
    ! X direction
    Al = 0.0
    Ar = 0.0
    ztick = zz + 0.5d0*(dxmin-dx)
    do it=1,nticks
      ! left
      rs = SQRT( (xx + 0.5d0*dx)**2 + (yy - 0.5d0*dx)**2 + ztick**2 )
      dnfw=1d0/rs/(1d0+rs)**2
      dnfw=MIN(dmax,dnfw)
      rc = SQRT( (xx + 0.5d0*dx)**2 + (yy - 0.5d0*dx)**2 )
      Al= Al + rc * dnfw**two3rd

      ! right
      rs = SQRT( (xx + 0.5d0*dx)**2 + (yy + 0.5d0*dx)**2 + ztick**2 )
      dnfw=1d0/rs/(1d0+rs)**2
      dnfw=MIN(dmax,dnfw)
      rc = SQRT( (xx + 0.5d0*dx)**2 + (yy + 0.5d0*dx)**2 )
      Ar= Ar + rc * dnfw**two3rd

      ! advance in z direction
      ztick=ztick + dxmin
    end do

    q(i,nvar+1)=Bscale * (Ar-Al)

    ! B right
    ! Y direction
    Al = 0.0
    Ar = 0.0
    ztick = zz + 0.5d0*(dxmin-dx)
    do it=1,nticks
      ! left
      rs = SQRT( (xx - 0.5d0*dx)**2 + (yy + 0.5d0*dx)**2 + ztick**2 )
      dnfw=1d0/rs/(1d0+rs)**2
      dnfw=MIN(dmax,dnfw)
      rc = SQRT( (xx - 0.5d0*dx)**2 + (yy + 0.5d0*dx)**2 )
      Al= Al + rc * dnfw**two3rd

      ! right
      rs = SQRT( (xx + 0.5d0*dx)**2 + (yy + 0.5d0*dx)**2 + ztick**2 )
      dnfw=1d0/rs/(1d0+rs)**2
      dnfw=MIN(dmax,dnfw)
      rc = SQRT( (xx + 0.5d0*dx)**2 + (yy + 0.5d0*dx)**2 )
      Ar= Ar + rc * dnfw**two3rd

      ! advance in z direction
      ztick=ztick + dxmin
    end do

    q(i,nvar+2)=Bscale * (Al-Ar)

    ! B right
    ! Z direction
    q(i,nvar+3)=0.0
  end do

end subroutine mag_toroidal

subroutine mag_dipole(x,q,dx,nn,B_0)
  use amr_parameters
  use hydro_parameters, only: nvar
  use const
  implicit none

  real(dp),dimension(1:nvector,1:ndim)::x    ! Cell center position
  real(dp),dimension(1:nvector,1:nvar+3)::q  ! Primitive variables
  integer::nn                                ! Number of cells
  real(dp)::dx                               ! Cell size
  real(dp)::B_0                              ! Default B-strength
  real(dp)::xx,yy,zz
  integer::i,it,nticks
  real(dp)::dxmin,zfl,zceil,rtick,x_edge,y_edge
  real(dp),dimension(1:2,1:2,1:2)::A_mag     ! x/y,left/right,up/down
  real(dp)::dmax,dnfw

  ! Dipole field, prop. to rho**(2/3)

  dxmin=boxlen*0.5d0**nlevelmax
  nticks=nint(dx/dxmin)

  dmax=1d0/(2.0d0*dxmin)/(1d0+(2.0d0*dxmin))**2

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
      rtick = SQRT( x_edge**2 + y_edge**2 + zfl**2 )
      dnfw=1d0/rtick/(1d0+rtick)**2
      dnfw=MIN(dmax,dnfw)
      A_mag(1,1,1) = A_mag(1,1,1) - dnfw**two3rd * y_edge

      ! ceiling
      rtick = SQRT( x_edge**2 + y_edge**2 + zceil**2 )
      dnfw=1d0/rtick/(1d0+rtick)**2
      dnfw=MIN(dmax,dnfw)
      A_mag(1,1,2) = A_mag(1,1,2) - dnfw**two3rd * y_edge

      ! advance along x
      x_edge=x_edge + dxmin
    end do

    ! A_(x,r)
    x_edge = xx + 0.5*(dxmin-dx)
    y_edge = yy + 0.5*dx
    do it=1,nticks
      ! floor
      rtick = SQRT( x_edge**2 + y_edge**2 + zfl**2 )
      dnfw=1d0/rtick/(1d0+rtick)**2
      dnfw=MIN(dmax,dnfw)
      A_mag(1,2,1) = A_mag(1,2,1) - dnfw**two3rd * y_edge

      ! ceiling
      rtick = SQRT( x_edge**2 + y_edge**2 + zceil**2 )
      dnfw=1d0/rtick/(1d0+rtick)**2
      dnfw=MIN(dmax,dnfw)
      A_mag(1,2,2) = A_mag(1,2,2) - dnfw**two3rd * y_edge

      ! advance along x
      x_edge=x_edge + dxmin
    end do

    ! A_(y,l)
    x_edge = xx - 0.5*dx
    y_edge = yy + 0.5*(dxmin-dx)
    do it=1,nticks
      ! floor
      rtick = SQRT( x_edge**2 + y_edge**2 + zfl**2 )
      dnfw=1d0/rtick/(1d0+rtick)**2
      dnfw=MIN(dmax,dnfw)
      A_mag(2,1,1) = A_mag(2,1,1) + dnfw**two3rd * x_edge

      ! ceiling
      rtick = SQRT( x_edge**2 + y_edge**2 + zceil**2 )
      dnfw=1d0/rtick/(1d0+rtick)**2
      dnfw=MIN(dmax,dnfw)
      A_mag(2,1,2) = A_mag(2,1,2) + dnfw**two3rd * x_edge

      ! advance along y
      y_edge=y_edge + dxmin
    end do

    ! A_(y,r)
    x_edge = xx + 0.5*dx
    y_edge = yy + 0.5*(dxmin-dx)
    do it=1,nticks
      ! floor
      rtick = SQRT( x_edge**2 + y_edge**2 + zfl**2 )
      dnfw=1d0/rtick/(1d0+rtick)**2
      dnfw=MIN(dmax,dnfw)
      A_mag(2,2,1) = A_mag(2,2,1) + dnfw**two3rd * x_edge

      ! ceiling
      rtick = SQRT( x_edge**2 + y_edge**2 + zceil**2 )
      dnfw=1d0/rtick/(1d0+rtick)**2
      dnfw=MIN(dmax,dnfw)
      A_mag(2,2,2) = A_mag(2,2,2) + dnfw**two3rd * x_edge

      ! advance along y
      y_edge=y_edge + dxmin
    end do

    ! average value
    A_mag = A_mag / DBLE(nticks)

    ! B left
    ! X direction
    q(i,6)= - B_0 * (A_mag(2,1,2)-A_mag(2,1,1)) / dx ! [-d/dz A_y]
    ! Y direction
    q(i,7)=   B_0 * (A_mag(1,1,2)-A_mag(1,1,1)) / dx ! [ d/dz A_x]
    ! Z direction  [d/dx A_y - d/dy A_x]
    q(i,8)=   B_0 * ( A_mag(2,2,1)-A_mag(2,1,1) - A_mag(1,2,1)+A_mag(1,1,1) ) / dx

    ! B right
    ! X direction
    q(i,nvar+1)= - B_0 * (A_mag(2,2,2)-A_mag(2,2,1)) / dx ! [-d/dz A_y]
    ! Y direction
    q(i,nvar+2)=   B_0 * (A_mag(1,2,2)-A_mag(1,2,1)) / dx ! [ d/dz A_x]
    ! Z direction  [d/dx A_y - d/dy A_x]
    q(i,nvar+3)=   B_0 * ( A_mag(2,2,2)-A_mag(2,1,2) - A_mag(1,2,2)+A_mag(1,1,2) ) / dx
  end do

end subroutine mag_dipole

subroutine mag_quadrupole(x,q,dx,nn,B_0)
  use amr_parameters
  use hydro_parameters, only: nvar
  use const
  implicit none
  real(dp),dimension(1:nvector,1:ndim)::x    ! Cell center position
  real(dp),dimension(1:nvector,1:nvar+3)::q  ! Primitive variables
  integer::nn                                ! Number of cells
  real(dp)::dx                               ! Cell size
  real(dp)::B_0                              ! Default B-strength
  real(dp)::xx,yy,zz
  integer::i,it,nticks
  real(dp)::dxmin,zfl,zceil,rc,rs,x_edge,y_edge
  real(dp),dimension(1:2,1:2,1:2)::A_mag     ! x/y,left/right,up/down
  real(dp)::dmax,dnfw

  ! Quadrupole field, prop. to rho**(2/3)

  dxmin=boxlen*0.5d0**nlevelmax
  nticks=nint(dx/dxmin)

  dmax=1d0/(2.0d0*dxmin)/(1d0+(2.0d0*dxmin))**2

  do i=1,nn
    ! box-centered coordinates
    xx=x(i,1) - boxlen * 0.5d0
    yy=x(i,2) - boxlen * 0.5d0
    zz=x(i,3) - boxlen * 0.5d0

    zfl   = zz - 0.5*dx
    zceil = zz + 0.5*dx

    A_mag = 0.0

    ! A_(x,l)
    x_edge = xx + 0.5*(dxmin-dx)
    y_edge = yy - 0.5*dx
    do it=1,nticks
      rc = SQRT( x_edge**2 + y_edge**2 )

      ! floor
      rs = SQRT( x_edge**2 + y_edge**2 + zfl**2 )
      dnfw=1d0/rs/(1d0+rs)**2
      dnfw=MIN(dmax,dnfw)

      A_mag(1,1,1) = A_mag(1,1,1) - zfl * dnfw**two3rd * y_edge / rc

      ! ceiling
      rs = SQRT( x_edge**2 + y_edge**2 + zceil**2 )
      dnfw=1d0/rs/(1d0+rs)**2
      dnfw=MIN(dmax,dnfw)
      A_mag(1,1,2) = A_mag(1,1,2) - zceil * dnfw**two3rd * y_edge / rc

      ! advance along x
      x_edge=x_edge + dxmin
    end do

    ! A_(x,r)
    x_edge = xx + 0.5*(dxmin-dx)
    y_edge = yy + 0.5*dx
    do it=1,nticks
      rc = SQRT( x_edge**2 + y_edge**2 )

      ! floor
      rs = SQRT( x_edge**2 + y_edge**2 + zfl**2 )
      dnfw=1d0/rs/(1d0+rs)**2
      dnfw=MIN(dmax,dnfw)
      A_mag(1,2,1) = A_mag(1,2,1) - zfl * dnfw**two3rd * y_edge / rc

      ! ceiling
      rs = SQRT( x_edge**2 + y_edge**2 + zceil**2 )
      dnfw=1d0/rs/(1d0+rs)**2
      dnfw=MIN(dmax,dnfw)
      A_mag(1,2,2) = A_mag(1,2,2) - zceil * dnfw**two3rd * y_edge / rc

      ! advance along x
      x_edge=x_edge + dxmin
    end do

    ! A_(y,l)
    x_edge = xx - 0.5*dx
    y_edge = yy + 0.5*(dxmin-dx)
    do it=1,nticks
      rc = SQRT( x_edge**2 + y_edge**2 )

      ! floor
      rs = SQRT( x_edge**2 + y_edge**2 + zfl**2 )
      dnfw=1d0/rs/(1d0+rs)**2
      dnfw=MIN(dmax,dnfw)
      A_mag(2,1,1) = A_mag(2,1,1) + zfl * dnfw**two3rd * x_edge / rc

      ! ceiling
      rs = SQRT( x_edge**2 + y_edge**2 + zceil**2 )
      dnfw=1d0/rs/(1d0+rs)**2
      dnfw=MIN(dmax,dnfw)
      A_mag(2,1,2) = A_mag(2,1,2) + zceil * dnfw**two3rd * x_edge / rc

      ! advance along y
      y_edge=y_edge + dxmin
    end do

    ! A_(y,r)
    x_edge = xx + 0.5*dx
    y_edge = yy + 0.5*(dxmin-dx)
    do it=1,nticks
      rc = SQRT( x_edge**2 + y_edge**2 )

      ! floor
      rs = SQRT( x_edge**2 + y_edge**2 + zfl**2 )
      dnfw=1d0/rs/(1d0+rs)**2
      dnfw=MIN(dmax,dnfw)
      A_mag(2,2,1) = A_mag(2,2,1) + zfl * dnfw**two3rd * x_edge / rc

      ! ceiling
      rs = SQRT( x_edge**2 + y_edge**2 + zceil**2 )
      dnfw=1d0/rs/(1d0+rs)**2
      dnfw=MIN(dmax,dnfw)
      A_mag(2,2,2) = A_mag(2,2,2) + zceil * dnfw**two3rd * x_edge / rc

      ! advance along y
      y_edge=y_edge + dxmin
    end do

    ! average value
    A_mag = A_mag / DBLE(nticks)

    ! B left
    ! X direction
    q(i,6)= - B_0 * (A_mag(2,1,2)-A_mag(2,1,1)) / dx ! [-d/dz A_y]
    ! Y direction
    q(i,7)=   B_0 * (A_mag(1,1,2)-A_mag(1,1,1)) / dx ! [ d/dz A_x]
    ! Z direction  [d/dx A_y - d/dy A_x]
    q(i,8)=   B_0 * ( A_mag(2,2,1)-A_mag(2,1,1) - A_mag(1,2,1)+A_mag(1,1,1) ) / dx

    ! B right
    ! X direction
    q(i,nvar+1)= - B_0 * (A_mag(2,2,2)-A_mag(2,2,1)) / dx ! [-d/dz A_y]
    ! Y direction
    q(i,nvar+2)=   B_0 * (A_mag(1,2,2)-A_mag(1,2,1)) / dx ! [ d/dz A_x]
    ! Z direction  [d/dx A_y - d/dy A_x]
    q(i,nvar+3)=   B_0 * ( A_mag(2,2,2)-A_mag(2,1,2) - A_mag(1,2,2)+A_mag(1,1,2) ) / dx
  end do

end subroutine mag_quadrupole

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
  real(dp)::xx,yy,zz,vx,vy,vz,rr,tt,omega,aa,twopi

  ! Add here, if you wish, some user-defined initial conditions
  aa=1.0
  twopi=2d0*ACOS(-1d0)
  do i=1,ncell

     xx=x(i,1)
     yy=x(i,2)
     zz=x(i,3)

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
     v(i,2)=vy
     v(i,3)=vz

  end do


end subroutine velana
