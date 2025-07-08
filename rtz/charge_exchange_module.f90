! charge_exchange_module.f90
module charge_exchange_module
  use amr_parameters, only: dp
  implicit none

  private  ! everything is private by default
  public :: load_ct_rates, charge_transfer_recombination, charge_transfer_ionization

  real(dp):: CTRecomb(6,4,30)
  real(dp):: CTIon(7,4,30)

CONTAINS

SUBROUTINE load_ct_rates()
  !Load the charge transfer ionization
  !and recombination rates from file
  use amr_commons, only: myid
  implicit none

  integer:: i, j, k, unit_num, ios

  if(myid.eq.1) write(*,*) 'Initializing charge exchange rates'

  ! Load charge exchange file
  open(newunit=unit_num, file='./data/charge_transfer/ct_ionization.dat', status='old', action='read', iostat=ios)
  if (ios /= 0) then
      write(*,*) 'Error: Could not open CT Ionization file'
      return
  end if

  !Reading data from the file into the 3D array
  do i = 3, 30
     do j = 1, 3
         read(unit_num, *, iostat=ios) CTIon(:,j,i)
        if (ios /= 0) exit
     end do
     if (ios /= 0) exit
  end do
  close(unit_num)

  ! Zero out helium
  do j = 1, 3
     do k = 1, 7
        CTIon(k,j,2) = 0.d0
     end do
  end do

  ! Load Recombination
  open(newunit=unit_num, file='./data/charge_transfer/ct_recombination.dat', status='old', action='read', iostat=ios)
  if (ios /= 0) then
      write(*,*) 'Error: Could not open CT Recombination file'
      return
  end if

  !Reading data from the file into the 3D array
  do i = 2, 30
     do j = 1, 4
         read(unit_num, *, iostat=ios) CTRecomb(:,j,i)
        if (ios /= 0) exit
     end do
     if (ios /= 0) exit
  end do
  close(unit_num)

END SUBROUTINE load_ct_rates

FUNCTION charge_transfer_recombination(ion, nelem, T) result(rate)
  ! ion is stage of ionization, 2 for the ion going to the atom
  ! nelem is atomic number of element, 2 up to 30
  ! Example:  O+ + H => O + H+ is HCTRecom(2,8,1e4)
  ! Note that temperature is in linear scale
  implicit none

  integer, intent(in):: ion, nelem
  real(dp), intent(in):: T
  real(dp):: rate
  real(dp):: a_op, b_op, c_op, d_op, e_op, f_op
  real(dp):: logT, tused
  integer:: ipIon

  rate = 0.0
  ! No recombination on the ground state
  if (ion .eq. 1) then
     return
  end if

  ipIon = ion - 1
  !use statistical charge transfer for ion > 5
  if (ion .gt. 5) then
     rate = 1.92d-9 * real(ipIon, dp)
     return
  end if

  !Make sure T is between temp. boundaries
  tused = max(min(T,CTRecomb(6,ipIon,nelem)),CTRecomb(5,ipIon,nelem))
  tused = tused * 1d-4
  tused = max(tused,1d-10) ! harley added to prevent zero temperature

  ! The interpolation equation
  rate = CTRecomb(1,ipIon,nelem) * 1d-9 * (tused**CTRecomb(2,ipIon,nelem)) * (1.d0 + CTRecomb(3,ipIon,nelem) * exp(CTRecomb(4,ipIon,nelem)*tused) )

END FUNCTION charge_transfer_recombination

FUNCTION charge_transfer_ionization(ion, nelem, T) result(rate)
  ! ion is stage of ionization, 1 for atom
  ! nelem is atomic number of element, 2 up to 30
  ! Example:  O + H+ => O+ + H is HCTIon(1,8,1e4)
  ! Note that temperature is in linear scale
  implicit none
  
  integer, intent(in):: ion, nelem
  real(dp), intent(in):: T
  real(dp):: rate
  real(dp):: a, b, c
  real(dp):: a_o, b_o, c_o, d_o, e_o, f_o, g_o
  real(dp):: logT, tused
  integer:: ipIon

  rate = 0.d0

  ipIon = ion
  if (ipIon.gt.2) then 
     rate = 0.d0
     return
  end if

  ! Make sure T is between temp. boundaries
  tused = max(min(T,CTIon(6,ipIon,nelem)),CTIon(5,ipIon,nelem))
  tused = tused * 1d-4
  tused = max(tused,1d-10) ! harley added to prevent zero temperature

  ! the interpolation equation
  rate = CTIon(1,ipIon,nelem) * 1d-9 * (tused**CTIon(2,ipIon,nelem)) * (1.d0 + CTIon(3,ipIon,nelem) * exp(CTIon(4,ipIon,nelem)*tused) ) * exp(-1.d0 * CTIon(7,ipIon,nelem)/tused)

END FUNCTION charge_transfer_ionization

end module charge_exchange_module