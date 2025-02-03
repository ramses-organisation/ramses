!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine rho_ana(x,d,dx,ncell)
  use amr_parameters
  use hydro_parameters
  use poisson_parameters
  use constants
  implicit none
  integer ::ncell                         ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector)       ::d ! Density
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates analytical Poisson source term.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! d(i) is the density field in user units.
  !================================================================
  integer::i
  real(dp)::gmass,dmass,emass,xmass,ymass,zmass,rr,rx,ry,rz,dd
  real(dp):: a1,a2,z0,a1_rho,a2_rho
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2

  select case(gravity_type)
  ! Remark: no obvious equivalent for gravana's gravity_type==1
  case(-2)
    ! point mass
    emass=dx
    emass=gravity_params(1) ! Softening length
    xmass=gravity_params(2) ! Point mass coordinates
    ymass=gravity_params(3)
    zmass=gravity_params(4)
    dmass=1.0d0/(emass*(1.0d0+emass)**2)

    do i=1,ncell
       rx=0.0d0; ry=0.0d0; rz=0.0d0
       rx=x(i,1)-xmass
#if NDIM>1
       ry=x(i,2)-ymass
#endif
#if NDIM>2
       rz=x(i,3)-zmass
#endif
       rr=sqrt(rx**2+ry**2+rz**2)
       dd=1.0d0/(rr*(1.0d0+rr)**2)
       d(i)=MIN(dd,dmass)
  end do

  case(-3)
    ! vertical galactic gravitational field
    ! Kuijken & Gilmore 1989 taken from Joung & MacLow (2006)
    a1 = gravity_params(1) ! Star potential coefficient in kpc Myr-2
    a2 = gravity_params(2) ! DM potential coefficient in Myr-2
    z0 = gravity_params(3) ! Scale height in pc
    ! standard values are: a1 = 1.42d-3, a2 = 5.49d-4, z0 = 0.18d3 pc

    ! The gravitational field is given by
    ! g = -a1 z / sqrt(z^2+z0^2) - a2 z
    ! rho = [(a1 / z0) ( (z/z0)^2 + 1)^(-3/2) + a2] / (4piG)

    ! convert to code units
    call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
    a1_rho = (a1*1d3 / z0) / (Myr2sec**2) / (4d0*pi*factG_in_cgs) / scale_d
    a2_rho = a2/(Myr2sec)**2 / (4d0*pi*factG_in_cgs) / scale_d
    z0=z0*pc2cm / scale_l

    do i=1,ncell
       ! the last dimension is vertical (1D -> x, 2D -> y, 3D -> z)
       rz=x(i,ndim)-0.5*boxlen
       d(i)= a1_rho / (1.+(rz/z0)**2)**(1.5) + a2_rho
    end do
  end select

end subroutine rho_ana
