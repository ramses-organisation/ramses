! cosmic_ray_ionization.f90
module cosmic_ray_ionization_module
  use amr_parameters, only: dp
  implicit none

  private  ! everything is private by default
  public :: cosmic_ray_ionization_rates, cosmic_ray_ionization_rates_induced_UV, cosmic_ray_ionization_rates_induced_UV_heat
  public :: initialize_cr_rates, secondary_cr_rates

  real(dp) :: cosmic_ray_ionization_rates(27,27)
  real(dp) :: cosmic_ray_ionization_rates_induced_UV(27)
  real(dp) :: cosmic_ray_ionization_rates_induced_UV_heat(27)

  ! Hydrogen
  real(dp), parameter :: cosmic_ray_ionization_rates_hydrogen(27)  = (/ 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
  
  ! Helium
  real(dp), parameter :: cosmic_ray_ionization_rates_helium(27)    = (/ 1.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
  
  ! Carbon
  real(dp), parameter :: cosmic_ray_ionization_rates_carbon(27)    = (/ 3.83, 1.6638695, 0.8312262, 0.4541473, 0.06937006, 0.027755102, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
  
  ! Nitrogen
  real(dp), parameter :: cosmic_ray_ionization_rates_nitrogen(27)  = (/ 4.52, 1.8237852, 0.8574613, 0.5265334, 0.30504152, 0.04926644, 0.020386748, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
  
  ! Oxygen
  real(dp), parameter :: cosmic_ray_ionization_rates_oxygen(27)    = (/ 5.637, 1.9201249, 0.98636097, 0.5261845, 0.36398333, 0.21907325, 0.03679156, 0.015607069, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
  
  ! Neon
  real(dp), parameter :: cosmic_ray_ionization_rates_neon(27)      = (/ 4.8251953, 2.2843935, 1.2726027, 0.6982404, 0.43035796, 0.2581805, 0.204098, 0.12864962, 0.022742474, 0.009985316, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
  
  ! Magnesium
  real(dp), parameter :: cosmic_ray_ionization_rates_magnesium(27) = (/ 9.016716, 5.037889, 1.345602, 0.86565775, 0.5749281, 0.36411676, 0.24167651, 0.15331018, 0.13045996, 0.084565096, 0.015437003, 0.006928171, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
  
  ! Silicon
  real(dp), parameter :: cosmic_ray_ionization_rates_silicon(27)   = (/ 6.37329, 2.4664342, 2.5934594, 1.8422631, 0.64957255, 0.46263242, 0.3302073, 0.22406998, 0.15485357, 0.10161172, 0.09062556, 0.059769776, 0.011156686, 0.0050879163, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
  
  ! Sulfur
  real(dp), parameter :: cosmic_ray_ionization_rates_sulfur(27)    = (/ 7.5164065, 2.8824131, 1.5510099, 0.8586138, 1.2691398, 0.9713803, 0.38636714, 0.28946567, 0.21488377, 0.15199806, 0.107721165, 0.07223445, 0.066572525, 0.044464245, 0.008436725, 0.003892387, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /)
  
  ! Iron
  real(dp), parameter :: cosmic_ray_ionization_rates_iron(27)      = (/ 3.192021, 5.842174, 4.724112, 2.5858068, 1.7483323, 1.2082134, 0.85902447, 0.62585855, 0.46490052, 0.3626024, 0.28070357, 0.20538796, 0.15056798, 0.10395802, 0.21419027, 0.18064684, 0.086154655, 0.07006637, 0.056015324, 0.042974997, 0.032204125, 0.02267751, 0.022580078, 0.015563778, 0.0030807566, 0.0014658332, 0.0 /)

CONTAINS

SUBROUTINE initialize_cr_rates()
  implicit none
  integer :: i,j

  ! Initialize
  do i = 1, 27
     do j = 1, 27
        cosmic_ray_ionization_rates(i,j) = 0.d0
        cosmic_ray_ionization_rates_induced_UV_heat(i) = 0.d0
     end do
  end do

  ! Now populate
  do i = 1, 27
     cosmic_ray_ionization_rates(1,i)  = cosmic_ray_ionization_rates_hydrogen(i)
     cosmic_ray_ionization_rates(2,i)  = cosmic_ray_ionization_rates_helium(i)
     cosmic_ray_ionization_rates(6,i)  = cosmic_ray_ionization_rates_carbon(i)
     cosmic_ray_ionization_rates(7,i)  = cosmic_ray_ionization_rates_nitrogen(i)
     cosmic_ray_ionization_rates(8,i)  = cosmic_ray_ionization_rates_oxygen(i)
     cosmic_ray_ionization_rates(10,i) = cosmic_ray_ionization_rates_neon(i)
     cosmic_ray_ionization_rates(12,i) = cosmic_ray_ionization_rates_magnesium(i)
     cosmic_ray_ionization_rates(14,i) = cosmic_ray_ionization_rates_silicon(i)
     cosmic_ray_ionization_rates(16,i) = cosmic_ray_ionization_rates_sulfur(i)
     cosmic_ray_ionization_rates(26,i) = cosmic_ray_ionization_rates_iron(i)
  end do

  ! Cosmic ray data for ionization from induced UV emission
  ! Data from: https://home.strw.leidenuniv.nl/~ewine/photo/cosmic_ray_rates.html
  ! Note that this assumes 10-16 s-1 H2-1 for the primary H2 cosmic ray dissociation rate
  ! Initialize to 0
  do i = 1, 27
     cosmic_ray_ionization_rates_induced_UV(i) = 0.d0
  end do

  cosmic_ray_ionization_rates_induced_UV(1)  = 0.d0 !4.08d-16 ! Hydrogen
  cosmic_ray_ionization_rates_induced_UV(2)  = 0.00d+00 ! Helium  --> ground state too high
  cosmic_ray_ionization_rates_induced_UV(6)  = 2.60d-14 ! Carbon
  cosmic_ray_ionization_rates_induced_UV(7)  = 7.34d-17 ! Nitrogen
  cosmic_ray_ionization_rates_induced_UV(8)  = 2.70d-16 ! Oxygen
  cosmic_ray_ionization_rates_induced_UV(10) = 0.00d+00 ! Neon  --> ground state too high
  cosmic_ray_ionization_rates_induced_UV(12) = 1.12d-14 ! Magnesium
  cosmic_ray_ionization_rates_induced_UV(14) = 4.16d-13 ! Silicon
  cosmic_ray_ionization_rates_induced_UV(16) = 7.91d-14 ! Sulfur
  cosmic_ray_ionization_rates_induced_UV(26) = 4.81d-14 ! Iron

  cosmic_ray_ionization_rates_induced_UV_heat(1)  = 0.d0 ! Hydrogen
  cosmic_ray_ionization_rates_induced_UV_heat(2)  = 0.00d0 ! Helium  --> ground state too high
  cosmic_ray_ionization_rates_induced_UV_heat(6)  = 0.96d0 ! Carbon
  cosmic_ray_ionization_rates_induced_UV_heat(7)  = 1.17d0 ! Nitrogen
  cosmic_ray_ionization_rates_induced_UV_heat(8)  = 1.72d0 ! Oxygen
  cosmic_ray_ionization_rates_induced_UV_heat(10) = 0.00d0 ! Neon  --> ground state too high
  cosmic_ray_ionization_rates_induced_UV_heat(12) = 2.55d0 ! Magnesium
  cosmic_ray_ionization_rates_induced_UV_heat(14) = 2.08d0 ! Silicon
  cosmic_ray_ionization_rates_induced_UV_heat(16) = 1.32d0 ! Sulfur
  cosmic_ray_ionization_rates_induced_UV_heat(26) = 2.32d0 ! Iron

END SUBROUTINE initialize_cr_rates

FUNCTION secondary_cr_rates(xe) result(phi_s)
    ! Secondary CR ionization rate
    implicit none
    real(dp), intent(in) :: xe
    real(dp) :: phi_s

    phi_s = (1.d0 - (xe / 1.2d0)) * (0.670d0 / (1.d0 + (xe / 0.05d0)));
END FUNCTION secondary_cr_rates

end module cosmic_ray_ionization_module