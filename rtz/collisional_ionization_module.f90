! collisional_ionization_module.f90
module collisional_ionization_module
  use amr_parameters, only: dp
  implicit none

  private  ! everything is private by default
  public :: collisional_ionization

  ! Carbon
  real(dp), parameter :: dE_carbon(6)  = [11.3d0, 24.4d0, 47.9d0, 64.5d0, 392.1d0, 490.0d0]
  real(dp), parameter :: A_carbon(6)   = [0.685d-7, 0.186d-7, 0.635d-8, 0.150d-8, 0.299d-9, 0.123d-9]
  real(dp), parameter :: X_carbon(6)   = [0.193d0, 0.286d0, 0.427d0, 0.416d0, 0.666d0, 0.620d0]
  real(dp), parameter :: K_carbon(6)   = [0.25d0, 0.24d0, 0.21d0, 0.13d0, 0.02d0, 0.16d0]
  real(dp), parameter :: P_carbon(6)   = [0.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0]

  ! Oxygen
  real(dp), parameter :: dE_oxygen(8)  = [13.6d0, 35.1d0, 54.9d0, 77.4d0, 113.9d0, 138.1d0, 739.3d0, 871.5d0]
  real(dp), parameter :: A_oxygen(8)   = [0.359d-7, 0.139d-7, 0.931d-8, 0.102d-7, 0.219d-8, 0.195d-8, 0.212d-9, 0.521d-10]
  real(dp), parameter :: X_oxygen(8)   = [0.073d0, 0.212d0, 0.270d0, 0.614d0, 0.630d0, 0.360d0, 0.396d0, 0.629d0]
  real(dp), parameter :: K_oxygen(8)   = [0.34d0, 0.22d0, 0.27d0, 0.27d0, 0.17d0, 0.54d0, 0.35d0, 0.16d0]
  real(dp), parameter :: P_oxygen(8)   = [0.d0, 1.d0, 1.d0, 0.d0, 1.d0, 0.d0, 0.d0, 1.d0]

  ! Nitrogen
  real(dp), parameter :: dE_nitrogen(7) = [14.5d0, 29.6d0, 47.5d0, 77.5d0, 97.9d0, 552.1d0, 667.0d0]
  real(dp), parameter :: A_nitrogen(7)  = [0.482d-7, 0.298d-7, 0.810d-8, 0.371d-8, 0.151d-8, 0.371d-9, 0.777d-10]
  real(dp), parameter :: X_nitrogen(7)  = [0.0652d0, 0.310d0, 0.350d0, 0.549d0, 0.0167d0, 0.546d0, 0.624d0]
  real(dp), parameter :: K_nitrogen(7)  = [0.42d0, 0.30d0, 0.24d0, 0.18d0, 0.74d0, 0.29d0, 0.16d0]
  real(dp), parameter :: P_nitrogen(7)  = [0.d0, 0.d0, 1.d0, 1.d0, 0.d0, 0.d0, 1.d0]

  ! Neon
  real(dp), parameter :: dE_neon(10) = [21.6d0, 41.0d0, 63.5d0, 97.1d0, 126.2d0, 157.9d0, 207.3d0, 239.1d0, 1196.0d0, 1360.6d0]
  real(dp), parameter :: A_neon(10)  = [0.150d-7, 0.198d-7, 0.703d-8, 0.424d-8, 0.279d-8, 0.345d-8, 0.956d-9, 0.473d-9, 0.392d-10, 0.277d-10]
  real(dp), parameter :: X_neon(10)  = [0.0329d0, 0.295d0, 0.0677d0, 0.0482d0, 0.305d0, 0.581d0, 0.749d0, 0.992d0, 0.262d0, 0.661d0]
  real(dp), parameter :: K_neon(10)  = [0.43d0, 0.20d0, 0.39d0, 0.58d0, 0.25d0, 0.28d0, 0.14d0, 0.04d0, 0.20d0, 0.13d0]
  real(dp), parameter :: P_neon(10)  = [1.d0, 0.d0, 1.d0, 1.d0, 1.d0, 0.d0, 1.d0, 1.d0, 1.d0, 1.d0]

  ! Magnesium
  real(dp), parameter :: dE_magnesium(12) = [7.6d0, 15.2d0, 80.1d0, 109.3d0, 141.3d0, 186.5d0, 224.9d0, 266.0d0, 328.2d0, 367.5d0, 1761.8d0, 1962.7d0]
  real(dp), parameter :: A_magnesium(12)  = [0.621d-6, 0.192d-7, 0.556d-8, 0.435d-8, 0.710d-8, 0.170d-8, 0.122d-8, 0.220d-8, 0.486d-9, 0.235d-9, 0.206d-10, 0.175d-10]
  real(dp), parameter :: X_magnesium(12)  = [0.592d0, 0.0027d0, 0.107d0, 0.159d0, 0.658d0, 0.242d0, 0.343d0, 0.897d0, 0.751d0, 1.030d0, 0.196d0, 0.835d0]
  real(dp), parameter :: K_magnesium(12)  = [0.39d0, 0.85d0, 0.30d0, 0.31d0, 0.25d0, 0.28d0, 0.23d0, 0.22d0, 0.14d0, 0.10d0, 0.25d0, 0.11d0]
  real(dp), parameter :: P_magnesium(12)  = [0.d0, 0.d0, 1.d0, 1.d0, 0.d0, 1.d0, 1.d0, 0.d0, 1.d0, 1.d0, 1.d0, 1.d0]

  ! Silicon
  real(dp), parameter :: dE_silicon(14) = [8.2d0, 16.4d0, 33.5d0, 54.0d0, 166.8d0, 205.3d0, 246.5d0, 303.5d0, 351.1d0, 401.4d0, 476.4d0, 523.5d0, 2437.7d0, 2673.2d0]
  real(dp), parameter :: A_silicon(14)  = [0.188d-6, 0.643d-7, 0.201d-7, 0.494d-8, 0.176d-8, 0.174d-8, 0.123d-8, 0.827d-9, 0.601d-9, 0.465d-9, 0.263d-9, 0.118d-9, 0.336d-10, 0.119d-10]
  real(dp), parameter :: X_silicon(14)  = [0.376d0, 0.632d0, 0.473d0, 0.172d0, 0.102d0, 0.180d0, 0.518d0, 0.239d0, 0.305d0, 0.666d0, 0.666d0, 0.734d0, 0.336d0, 0.989d0]
  real(dp), parameter :: K_silicon(14)  = [0.25d0, 0.20d0, 0.22d0, 0.23d0, 0.31d0, 0.29d0, 0.07d0, 0.28d0, 0.25d0, 0.04d0, 0.16d0, 0.16d0, 0.37d0, 0.08d0]
  real(dp), parameter :: P_silicon(14)  = [1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 0.d0, 1.d0]

  ! Sulfur data
  real(dp), parameter :: dE_sulfur(16) = [ 10.4d0, 23.3d0, 34.8d0, 47.3d0, 72.6d0, 88.1d0, &
      280.9d0, 328.2d0, 379.1d0, 447.1d0, 504.8d0, 564.7d0, 651.6d0, 707.2d0, 3223.9d0, 3494.2d0 ]
  real(dp), parameter :: A_sulfur(16) = [ 0.549d-7, 0.681d-7, 0.214d-7, 0.166d-7, 0.612d-8, 0.133d-8, &
      0.493d-8, 0.873d-9, 0.135d-8, 0.459d-9, 0.349d-9, 0.523d-9, 0.259d-9, 0.750d-10, 0.267d-10, 0.632d-11 ]
  real(dp), parameter :: X_sulfur(16) = [ 0.100d0, 0.693d0, 0.353d0, 1.030d0, 0.580d0, 0.0688d0, &
      1.130d0, 0.193d0, 0.431d0, 0.242d0, 0.305d0, 0.428d0, 0.854d0, 0.734d0, 0.572d0, 0.585d0 ]
  real(dp), parameter :: K_sulfur(16) = [ 0.25d0, 0.21d0, 0.24d0, 0.14d0, 0.19d0, 0.35d0, &
      0.16d0, 0.28d0, 0.32d0, 0.28d0, 0.25d0, 0.35d0, 0.12d0, 0.16d0, 0.28d0, 0.17d0 ]
  real(dp), parameter :: P_sulfur(16) = [ 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, &
      0.d0, 1.d0, 0.d0, 1.d0, 1.d0, 0.d0, 0.d0, 1.d0, 0.d0, 1.d0 ]

  ! Iron data
  real(dp), parameter :: dE_iron(26) = [ 7.9d0, 16.2d0, 30.6d0, 54.8d0, 75.0d0, 99.0d0, 125.0d0, 151.1d0, &
      233.6d0, 262.1d0, 290.0d0, 331.0d0, 361.0d0, 392.0d0, 457.0d0, 489.3d0, 1262.0d0, 1360.0d0, &
      1470.0d0, 1582.0d0, 1690.0d0, 1800.0d0, 1960.0d0, 2046.0d0, 8828.0d0, 9277.7d0 ]
  real(dp), parameter :: A_iron(26) = [ 0.252d-6, 0.221d-7, 0.410d-7, 0.353d-7, 0.104d-7, 0.123d-7, 0.947d-8, 0.471d-8, &
      0.302d-8, 0.234d-8, 0.176d-8, 0.114d-8, 0.866d-9, 0.661d-9, 0.441d-9, 0.118d-9, 0.361d-9, 0.245d-9, &
      0.187d-9, 0.133d-9, 0.784d-10, 0.890d-10, 0.229d-10, 0.112d-10, 0.246d-11, 0.979d-12 ]
  real(dp), parameter :: X_iron(26) = [ 0.701d0, 0.033d0, 0.366d0, 0.243d0, 0.285d0, 0.411d0, 0.458d0, 0.280d0, &
      0.697d0, 0.764d0, 0.805d0, 0.773d0, 0.805d0, 0.762d0, 0.698d0, 0.211d0, 1.160d0, 0.978d0, 0.988d0, &
      1.030d0, 0.848d0, 1.200d0, 0.936d0, 0.034d0, 1.020d0, 0.664d0 ]
  real(dp), parameter :: K_iron(26) = [ 0.25d0, 0.45d0, 0.17d0, 0.39d0, 0.17d0, 0.21d0, 0.21d0, 0.28d0, &
      0.15d0, 0.14d0, 0.14d0, 0.15d0, 0.14d0, 0.14d0, 0.16d0, 0.15d0, 0.09d0, 0.13d0, 0.14d0, 0.12d0, &
      0.14d0, 0.35d0, 0.12d0, 0.81d0, 0.02d0, 0.14d0 ]
  real(dp), parameter :: P_iron(26) = [ 0.d0, 1.d0, 0.d0, 0.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, &
      1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 1.d0, 0.d0, 1.d0, 0.d0, 1.d0, 1.d0 ]

CONTAINS

FUNCTION coll_ion(T, dE, A, X, K, P) result(rate)
    implicit none
    real(dp), intent(in) :: T, dE, A, X, K, P
    real(dp) :: rate
    real(dp) :: U

    ! Eqn 1 of Voronov 1997
    U = dE / (T * 8.61732814974056D-05)
    rate = A * (1.D0 + P * sqrt(U)) * U**K * exp(-U) / (X + U)
end FUNCTION coll_ion

FUNCTION collisional_ionization(T, ion, element_idx) result(rate)
    implicit none
    real(dp), intent(in) :: T
    integer, intent(in) :: ion, element_idx
    real(dp) :: T5, f
    real(dp) :: rate

    rate = 0.D0

    select case (element_idx)
      case (1) ! Hydrogen
        T5 = T / 1D5
        f = 1.D0 + sqrt(T5)
        rate = 5.85D-11 * (sqrt(T) / f) * exp(-157809.1D0 / T)

      case (2) ! Helium
        T5 = T / 1.D5
        f = 1.D0 + sqrt(T5)
          select case (ion)
            case (1) ! HeI -> HeII
              rate = 2.38D-11 * (sqrt(T) / f) * exp(-285335.4D0 / T)
            case (2) ! HeII -> HeIII
              rate = 5.68D-12 * (sqrt(T) / f) * exp(-631515.0D0 / T)
          end select

      case (6) ! Carbon
        rate = coll_ion(T, dE_carbon(ion), A_carbon(ion), X_carbon(ion), K_carbon(ion), P_carbon(ion))

      case (7) ! Nitrogen
        rate = coll_ion(T, dE_nitrogen(ion), A_nitrogen(ion), X_nitrogen(ion), K_nitrogen(ion), P_nitrogen(ion))

      case (8) ! Oxygen
        rate = coll_ion(T, dE_oxygen(ion), A_oxygen(ion), X_oxygen(ion), K_oxygen(ion), P_oxygen(ion))

      case (10) ! Neon
        rate = coll_ion(T, dE_neon(ion), A_neon(ion), X_neon(ion), K_neon(ion), P_neon(ion))

      case (12) ! Magnesium
        rate = coll_ion(T, dE_magnesium(ion), A_magnesium(ion), X_magnesium(ion), K_magnesium(ion), P_magnesium(ion))

      case (14) ! Silicon
        rate = coll_ion(T, dE_silicon(ion), A_silicon(ion), X_silicon(ion), K_silicon(ion), P_silicon(ion))

      case (16) ! Sulfur
        rate = coll_ion(T, dE_sulfur(ion), A_sulfur(ion), X_sulfur(ion), K_sulfur(ion), P_sulfur(ion))

      case (26) ! Iron
        rate = coll_ion(T, dE_iron(ion), A_iron(ion), X_iron(ion), K_iron(ion), P_iron(ion))

    end select

END FUNCTION collisional_ionization

end module collisional_ionization_module