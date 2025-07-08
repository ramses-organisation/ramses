! metal_yields_module.f90
module metal_yields_module
  use amr_parameters, only: dp
  implicit none

  private  ! everything is private by default
  public :: initialize_portinari_yields, get_portinari_ejecta_mass

  integer, parameter::portinati_N_metal = 5
  integer, parameter::portinati_N_mass = 14

  real(dp):: portinati_metal(portinati_N_metal)
  real(dp):: portinati_mass(portinati_N_mass)
  real(dp):: portinati_lifetimes(portinati_N_metal,portinati_N_mass)
  ! Yields are loaded as H, He, C, N, O, Ne, Mg, Si, S, Fe, Ca
  ! The last one is out of order...
  real(dp):: portinati_yields(portinati_N_metal,portinati_N_mass,11)

CONTAINS

SUBROUTINE initialize_portinari_yields()
  use amr_commons, only: myid
  implicit none

  integer:: i, j, k, unit_num, ios

  if(myid.eq.1) write(*,*) 'Initializing Portinati 1998 Yields'

  ! Load yields file
  open(newunit=unit_num, file='./data/yields/Portinari_yields.tsv', status='old', action='read', iostat=ios)
  if (ios /= 0) then
      write(*,*) 'Error: Could not open Portinari yields file'
      return
  end if

  do i=1,portinati_N_metal ! loop over metallicities

     ! Read in the metallicity
     read(unit_num, *, iostat=ios) portinati_metal(i) 

     do j=1,portinati_N_mass ! Loop over stellar masses

        ! Read in mass, stellar lifetime and 
        read(unit_num, *, iostat=ios) portinati_mass(j), portinati_lifetimes(i,j), portinati_yields(i,j,:)

     end do ! end loop over stellar masses

  end do ! end loop over metallicities

END SUBROUTINE initialize_portinari_yields

FUNCTION get_portinari_ejecta_mass(metallicity, mass, species) result(fq)
  implicit none
  real(dp), intent(in) :: metallicity, mass
  integer, intent(in) :: species
  real(dp):: xq, yq
  real(dp) :: fq
  integer :: i, j, element_idx
  real(dp) :: x1, x2, y1, y2
  real(dp) :: f11, f21, f12, f22
  real(dp) :: dx, dy
  real(dp) :: rescale_factor

  rescale_factor = 1.0 ! Rescale factor for metallicity and mass

  ! Enforce bounds
  xq = MIN(MAX(metallicity,portinati_N_metal(1)),portinati_N_metal(portinati_N_metal))
  yq = MIN(MAX(mass,portinati_N_mass(1)),portinati_N_mass(portinati_N_mass))

  ! Downweight the yields if we are below lower limits --> do not extrapolate in the case where we are above!
  if (metallicity.lt.portinati_N_metal(1)) then
     rescale_factor = rescale_factor * (metallicity / portinati_N_metal(1))
  endif

  if (mass.lt.portinati_N_mass(1)) then
     rescale_factor = rescale_factor * (mass / portinati_N_mass(1))
  endif

  ! Find i such that x(i) <= xq <= x(i+1)
  do i = 1, portinati_N_metal-1
    if (xq.ge.portinati_N_metal(i) .and. xq.le.portinati_N_metal(i+1)) exit
  end do

  ! Find j such that y(j) <= yq <= y(j+1)
  do j = 1, portinati_N_mass-1
    if (yq.ge.portinati_N_mass(j) .and. yq.le.portinati_N_mass(j+1)) exit
  end do

  ! Rescale factors from Appendix A3 of https://articles.adsabs.harvard.edu/pdf/2009MNRAS.399..574W
  select case (species)
    case (1) ! Hydrogen
      element_idx = 1
    case (2) ! Helium
      element_idx = 2
    case (6) ! Carbon
      element_idx = 3
      rescale_factor = rescale_factor * 0.5d0
    case (7) ! Nitrogen
      element_idx = 4
    case (8) ! Oxygen
      element_idx = 5
    case (10) ! Neon
      element_idx = 6
    case (12) ! Magnesium
      element_idx = 7
      rescale_factor = rescale_factor * 2.d0
    case (14) ! Silicon
      element_idx = 8
    case (16) ! Sulfur
      element_idx = 9
    case (20) ! Calcium
      element_idx = 11
    case (26) ! Iron
      element_idx = 10
      rescale_factor = rescale_factor * 0.5d0
    case default
      write(*,*) "Element not available in portinari yields", species
      stop
  end select

  ! Extract corner values
  x1 = portinati_N_metal(i);   x2 = portinati_N_metal(i+1)
  y1 = portinati_N_mass(j);   y2 = portinati_N_mass(j+1)

  f11 = portinati_yields(i, j, element_idx)
  f21 = portinati_yields(i+1, j, element_idx)
  f12 = portinati_yields(i, j+1, element_idx)
  f22 = portinati_yields(i+1, j+1, element_idx)

  dx = x2 - x1
  dy = y2 - y1

  ! Bilinear interpolation
  fq = (1.0 / (dx * dy)) * (
        f11 * (x2 - xq) * (y2 - yq) +
        f21 * (xq - x1) * (y2 - yq) +
        f12 * (x2 - xq) * (yq - y1) +
        f22 * (xq - x1) * (yq - y1)
       ) * rescale_factor
END FUNCTION get_portinari_ejecta_mass

end module metal_yields_module