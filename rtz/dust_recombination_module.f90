! dust_recombination_module.f90
module dust_recombination_module
  use amr_parameters, only: dp
  implicit none

  private  ! everything is private by default
  public :: dust_recombination

  real(dp), dimension(7,26) :: dust_rec_coefs = reshape( &
  [ 12.25d0, 8.074d-6, 1.378d0, 5.087d2, 1.586d-2, 0.4723d0, 1.102d-5, &   ! Hydrogen
    5.572d0, 3.185d-7, 1.512d0, 5.115d3, 3.903d-7, 0.4956d0, 5.494d-7, &   ! Helium
    0.000d0, 0.000d0,  0.000d0, 0.000d0, 0.000d0,  0.0000d0, 0.000d0, &    ! NA
    0.000d0, 0.000d0,  0.000d0, 0.000d0, 0.000d0,  0.0000d0, 0.000d0, &    ! NA
    0.000d0, 0.000d0,  0.000d0, 0.000d0, 0.000d0,  0.0000d0, 0.000d0, &    ! NA
    45.58d0, 6.089d-3, 1.128d0, 4.331d2, 4.845d-2, 0.8120d0, 1.333d-4, &   ! Carbon
    0.000d0, 0.000d0,  0.000d0, 0.000d0, 0.000d0,  0.0000d0, 0.000d0, &    ! NA
    0.000d0, 0.000d0,  0.000d0, 0.000d0, 0.000d0,  0.0000d0, 0.000d0, &    ! NA
    0.000d0, 0.000d0,  0.000d0, 0.000d0, 0.000d0,  0.0000d0, 0.000d0, &    ! NA
    0.000d0, 0.000d0,  0.000d0, 0.000d0, 0.000d0,  0.0000d0, 0.000d0, &    ! NA
    0.000d0, 0.000d0,  0.000d0, 0.000d0, 0.000d0,  0.0000d0, 0.000d0, &    ! NA
    2.510d0, 8.116d-8, 1.864d0, 6.170d4, 2.169d-6, 0.9605d0, 7.232d-5, &   ! Magnesium
    0.000d0, 0.000d0,  0.000d0, 0.000d0, 0.000d0,  0.0000d0, 0.000d0, &    ! NA
    2.166d0, 5.678d-8, 1.874d0, 4.375d4, 1.635d-6, 0.8964d0, 7.538d-5, &   ! Silicon
    0.000d0, 0.000d0,  0.000d0, 0.000d0, 0.000d0,  0.0000d0, 0.000d0, &    ! NA
    3.064d0, 7.769d-5, 1.319d0, 1.087d2, 3.475d-1, 0.4790d0, 4.689d-2, &   ! Sulfur
    0.000d0, 0.000d0,  0.000d0, 0.000d0, 0.000d0,  0.0000d0, 0.000d0, &    ! NA
    0.000d0, 0.000d0,  0.000d0, 0.000d0, 0.000d0,  0.0000d0, 0.000d0, &    ! NA
    0.000d0, 0.000d0,  0.000d0, 0.000d0, 0.000d0,  0.0000d0, 0.000d0, &    ! NA
    0.000d0, 0.000d0,  0.000d0, 0.000d0, 0.000d0,  0.0000d0, 0.000d0, &    ! NA
    0.000d0, 0.000d0,  0.000d0, 0.000d0, 0.000d0,  0.0000d0, 0.000d0, &    ! NA
    0.000d0, 0.000d0,  0.000d0, 0.000d0, 0.000d0,  0.0000d0, 0.000d0, &    ! NA
    0.000d0, 0.000d0,  0.000d0, 0.000d0, 0.000d0,  0.0000d0, 0.000d0, &    ! NA
    0.000d0, 0.000d0,  0.000d0, 0.000d0, 0.000d0,  0.0000d0, 0.000d0, &    ! NA
    0.000d0, 0.000d0,  0.000d0, 0.000d0, 0.000d0,  0.0000d0, 0.000d0, &    ! NA
    1.701d0, 9.554d-8, 1.851d0, 5.763d4, 4.116d-8, 0.9456d0, 2.198d-5 ], & ! Iron
    shape=[7,26])

CONTAINS

FUNCTION dust_recombination(ion, nelem, T, G, ne) result(rate)
  implicit none

  integer, intent(in):: ion, nelem
  real(dp), intent(in):: T, G, ne
  real(dp)::rate
  real(dp)::dr_sf, row_sum, phi
  real(dp)::a1, a2, a3
  integer::i

  dr_sf = 1.d0;

  ! initialize rate to 0.0
  rate = 0.d0

  ! No dust recombination except for the first excited state.
  ! Maybe this will change later...
  if (ion.ne.2) then
     return
  end if

  if (T.gt.1d4) then
     ! No dust recombination at high temperatures
     return
  else if (T.gt.1d3) then
     ! Scale down if greater than 1.d3
     dr_sf = exp(-1.d0 * T / 1.d3) / exp(-1.d0)
  end if

  if (T.lt.10.0) then
     ! No dust recombination at very low temperatures
     return 
  end if

  ! First check to make sure that all elements are not zero
  row_sum = 0.0
  do i = 1, 7
     row_sum = row_sum + abs(dust_rec_coefs(nelem,i))
  end do

  ! In this case there is nothing to compute
  if (row_sum.le.0.0) then
     return
  end if

  ! Extra fac on the denominator to avoid divide by zero
  phi = (G + 1d-8) * sqrt(T) / (ne + 1d-10) ! units K^1/2 cm^3

  a1 = dust_rec_coefs(2,nelem) * (phi**dust_rec_coefs(3,nelem))
  a2 = dust_rec_coefs(4,nelem) * (T**dust_rec_coefs(5,nelem))
  a3 = (-1.d0 * dust_rec_coefs(6,nelem)) - (dust_rec_coefs(7,nelem) * log(T))

  rate = 1.d-14 * dust_rec_coefs(1,nelem)
  rate = rate / (1.d0 + (a1 * (1.d0 + (a2 * (phi**a3)))))
  rate = rate * dr_sf

  ! Rescale the dust recombination rates due to a PAH normalization issue
  ! Zubko assumes 3.3d-5 C in PAH / H atom
  ! Weingartner & Draine assume 6d-5 C in PAH / H atom
  rate = rate * (3.3d-5 / 6.0d-5)

END FUNCTION dust_recombination

end module dust_recombination_module
