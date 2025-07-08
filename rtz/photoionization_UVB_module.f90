! photoionization_UVB_module.f90
module photoionization_UVB_module
  use amr_parameters, only: dp
  implicit none

  private  ! everything is private by default
  public :: HM12_UVB_z
  public :: load_UVB_data, update_UVB

  integer, parameter :: N_UVB_POINTS = 60

  ! Initialize the UVB arrays
  real(dp):: HM12_UVB_redshifts(N_UVB_POINTS)
  real(dp):: HM12_UVB_hydrogen(N_UVB_POINTS,1,2)
  real(dp):: HM12_UVB_helium(N_UVB_POINTS,2,2)
  real(dp):: HM12_UVB_carbon(N_UVB_POINTS,6,2) 
  real(dp):: HM12_UVB_nitrogen(N_UVB_POINTS,7,2)
  real(dp):: HM12_UVB_oxygen(N_UVB_POINTS,8,2)
  real(dp):: HM12_UVB_neon(N_UVB_POINTS,10,2)
  real(dp):: HM12_UVB_magnesium(N_UVB_POINTS,12,2)
  real(dp):: HM12_UVB_silicon(N_UVB_POINTS,14,2)
  real(dp):: HM12_UVB_sulfur(N_UVB_POINTS,16,2)
  real(dp):: HM12_UVB_iron(N_UVB_POINTS,26,2)

  ! Array that holds the UVB for each ion at a given redshift
  real(dp):: HM12_UVB_z(27,27,2)

CONTAINS

subroutine load_UVB_data()
    use amr_commons, only: myid
    implicit none
    
    integer :: unit_num, ios, i, j
    
    if(myid.eq.1) write(*,*) 'Initializing UV background data'
    
    ! Load redshifts
    open(newunit=unit_num, file='./data/HM12/redshifts.dat', status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(*,*) 'Error: Could not open Redshifts UVB file'
        return
    end if
    
    do i = 1, N_UVB_POINTS
        read(unit_num, *, iostat=ios) HM12_UVB_redshifts(i)
        if (ios /= 0) exit
    end do
    close(unit_num)
    
    ! Hydrogen photoionization
    open(newunit=unit_num, file='./data/HM12/hydrogen_pi.dat', status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(*,*) 'Error: Could not open Hydrogen UVB file'
        return
    end if
    
    do i = 1, N_UVB_POINTS
        read(unit_num, *, iostat=ios) HM12_UVB_hydrogen(i, :, 1)
        if (ios /= 0) exit
    end do
    close(unit_num)
    
    ! Hydrogen photoheating
    open(newunit=unit_num, file='./data/HM12/hydrogen_ph.dat', status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(*,*) 'Error: Could not open Hydrogen UVB heat file'
        return
    end if
    
    do i = 1, N_UVB_POINTS
        read(unit_num, *, iostat=ios) HM12_UVB_hydrogen(i, :, 2)
        if (ios /= 0) exit
    end do
    close(unit_num)
    
    ! Helium photoionization
    open(newunit=unit_num, file='./data/HM12/helium_pi.dat', status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(*,*) 'Error: Could not open Helium UVB file'
        return
    end if
    
    do i = 1, N_UVB_POINTS
        read(unit_num, *, iostat=ios) HM12_UVB_helium(i, :, 1)
        if (ios /= 0) exit
    end do
    close(unit_num)
    
    ! Helium photoheating
    open(newunit=unit_num, file='./data/HM12/helium_ph.dat', status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(*,*) 'Error: Could not open Helium UVB heat file'
        return
    end if
    
    do i = 1, N_UVB_POINTS
        read(unit_num, *, iostat=ios) HM12_UVB_helium(i, :, 2)
        if (ios /= 0) exit
    end do
    close(unit_num)
    
    ! Carbon photoionization
    open(newunit=unit_num, file='./data/HM12/carbon_pi.dat', status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(*,*) 'Error: Could not open Carbon UVB file'
        return
    end if
    
    do i = 1, N_UVB_POINTS
        read(unit_num, *, iostat=ios) HM12_UVB_carbon(i, :, 1)
        if (ios /= 0) exit
    end do
    close(unit_num)
    
    ! Carbon photoheating
    open(newunit=unit_num, file='./data/HM12/carbon_ph.dat', status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(*,*) 'Error: Could not open Carbon UVB heat file'
        return
    end if
    
    do i = 1, N_UVB_POINTS
        read(unit_num, *, iostat=ios) HM12_UVB_carbon(i, :, 2)
        if (ios /= 0) exit
    end do
    close(unit_num)
    
    ! Nitrogen photoionization
    open(newunit=unit_num, file='./data/HM12/nitrogen_pi.dat', status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(*,*) 'Error: Could not open Nitrogen UVB file'
        return
    end if
    
    do i = 1, N_UVB_POINTS
        read(unit_num, *, iostat=ios) HM12_UVB_nitrogen(i, :, 1)
        if (ios /= 0) exit
    end do
    close(unit_num)
    
    ! Nitrogen photoheating
    open(newunit=unit_num, file='./data/HM12/nitrogen_ph.dat', status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(*,*) 'Error: Could not open Nitrogen UVB heat file'
        return
    end if
    
    do i = 1, N_UVB_POINTS
        read(unit_num, *, iostat=ios) HM12_UVB_nitrogen(i, :, 2)
        if (ios /= 0) exit
    end do
    close(unit_num)
    
    ! Oxygen photoionization
    open(newunit=unit_num, file='./data/HM12/oxygen_pi.dat', status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(*,*) 'Error: Could not open Oxygen UVB file'
        return
    end if
    
    do i = 1, N_UVB_POINTS
        read(unit_num, *, iostat=ios) HM12_UVB_oxygen(i, :, 1)
        if (ios /= 0) exit
    end do
    close(unit_num)
    
    ! Oxygen photoheating
    open(newunit=unit_num, file='./data/HM12/oxygen_ph.dat', status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(*,*) 'Error: Could not open Oxygen UVB heat file'
        return
    end if
    
    do i = 1, N_UVB_POINTS
        read(unit_num, *, iostat=ios) HM12_UVB_oxygen(i, :, 2)
        if (ios /= 0) exit
    end do
    close(unit_num)
    
    ! Neon photoionization
    open(newunit=unit_num, file='./data/HM12/neon_pi.dat', status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(*,*) 'Error: Could not open Neon UVB file'
        return
    end if
    
    do i = 1, N_UVB_POINTS
        read(unit_num, *, iostat=ios) HM12_UVB_neon(i, :, 1)
        if (ios /= 0) exit
    end do
    close(unit_num)
    
    ! Neon photoheating
    open(newunit=unit_num, file='./data/HM12/neon_ph.dat', status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(*,*) 'Error: Could not open Neon UVB heat file'
        return
    end if
    
    do i = 1, N_UVB_POINTS
        read(unit_num, *, iostat=ios) HM12_UVB_neon(i, :, 2)
        if (ios /= 0) exit
    end do
    close(unit_num)
    
    ! Magnesium photoionization
    open(newunit=unit_num, file='./data/HM12/magnesium_pi.dat', status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(*,*) 'Error: Could not open Magnesium UVB file'
        return
    end if
    
    do i = 1, N_UVB_POINTS
        read(unit_num, *, iostat=ios) HM12_UVB_magnesium(i, :, 1)
        if (ios /= 0) exit
    end do
    close(unit_num)
    
    ! Magnesium photoheating
    open(newunit=unit_num, file='./data/HM12/magnesium_ph.dat', status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(*,*) 'Error: Could not open Magnesium UVB heat file'
        return
    end if
    
    do i = 1, N_UVB_POINTS
        read(unit_num, *, iostat=ios) HM12_UVB_magnesium(i, :, 2)
        if (ios /= 0) exit
    end do
    close(unit_num)
    
    ! Silicon photoionization
    open(newunit=unit_num, file='./data/HM12/silicon_pi.dat', status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(*,*) 'Error: Could not open Silicon UVB file'
        return
    end if
    
    do i = 1, N_UVB_POINTS
        read(unit_num, *, iostat=ios) HM12_UVB_silicon(i, :, 1)
        if (ios /= 0) exit
    end do
    close(unit_num)
    
    ! Silicon photoheating
    open(newunit=unit_num, file='./data/HM12/silicon_ph.dat', status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(*,*) 'Error: Could not open Silicon UVB heat file'
        return
    end if
    
    do i = 1, N_UVB_POINTS
        read(unit_num, *, iostat=ios) HM12_UVB_silicon(i, :, 2)
        if (ios /= 0) exit
    end do
    close(unit_num)
    
    ! Sulfur photoionization
    open(newunit=unit_num, file='./data/HM12/sulfur_pi.dat', status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(*,*) 'Error: Could not open Sulfur UVB file'
        return
    end if
    
    do i = 1, N_UVB_POINTS
        read(unit_num, *, iostat=ios) HM12_UVB_sulfur(i, :, 1)
        if (ios /= 0) exit
    end do
    close(unit_num)
    
    ! Sulfur photoheating
    open(newunit=unit_num, file='./data/HM12/sulfur_ph.dat', status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(*,*) 'Error: Could not open Sulfur UVB heat file'
        return
    end if
    
    do i = 1, N_UVB_POINTS
        read(unit_num, *, iostat=ios) HM12_UVB_sulfur(i, :, 2)
        if (ios /= 0) exit
    end do
    close(unit_num)
    
    ! Iron photoionization
    open(newunit=unit_num, file='./data/HM12/iron_pi.dat', status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(*,*) 'Error: Could not open Iron UVB file'
        return
    end if
    
    do i = 1, N_UVB_POINTS
        read(unit_num, *, iostat=ios) HM12_UVB_iron(i, :, 1)
        if (ios /= 0) exit
    end do
    close(unit_num)
    
    ! Iron photoheating
    open(newunit=unit_num, file='./data/HM12/iron_ph.dat', status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(*,*) 'Error: Could not open Iron UVB heat file'
        return
    end if
    
    do i = 1, N_UVB_POINTS
        read(unit_num, *, iostat=ios) HM12_UVB_iron(i, :, 2)
        if (ios /= 0) exit
    end do
    close(unit_num)
    
end subroutine load_UVB_data

subroutine update_UVB(redshift)
    !interpolates the UVB for all ions
    !to a particular redshift
    implicit none

    integer:: i, j, k
    real(dp),intent(in):: redshift
    real(dp):: scale_low, scale_high

    ! No interpolation beyond redshift bounds --> set the UVB to zero
    if (redshift .lt. HM12_UVB_redshifts(1) .or. redshift .ge. HM12_UVB_redshifts(N_UVB_POINTS)) then
       do i = 1,27
          do j = 1,27
             do k = 1,2
                HM12_UVB_z(i,j,k) = 0.d0
             end do
          end do
       end do
    else
       ! Loop through the array to find the interval in which x_interp lies
       do i = 1,N_UVB_POINTS-1
          if (redshift.ge.HM12_UVB_redshifts(i) .and. redshift.lt.HM12_UVB_redshifts(i+1)) then
             ! Perform linear interpolation
             scale_high = (redshift - HM12_UVB_redshifts(i)) / (HM12_UVB_redshifts(i+1) - HM12_UVB_redshifts(i))
             scale_low  = 1.d0 - scale_high

             do j = 1,27
               if (j.lt.2) then
                  HM12_UVB_z(1,j,1)  = (scale_low * HM12_UVB_hydrogen(i,j,1))  + (scale_high * HM12_UVB_hydrogen(i+1,j,1))
                  HM12_UVB_z(1,j,2)  = (scale_low * HM12_UVB_hydrogen(i,j,2))  + (scale_high * HM12_UVB_hydrogen(i+1,j,2))
               end if
               if (j.lt.3) then
                  HM12_UVB_z(2,j,1)  = (scale_low * HM12_UVB_helium(i,j,1))  + (scale_high * HM12_UVB_helium(i+1,j,1))
                  HM12_UVB_z(2,j,2)  = (scale_low * HM12_UVB_helium(i,j,2))  + (scale_high * HM12_UVB_helium(i+1,j,2))               
               end if
               if (j.lt.7) then
                  HM12_UVB_z(6,j,1)  = (scale_low * HM12_UVB_carbon(i,j,1))  + (scale_high * HM12_UVB_carbon(i+1,j,1))
                  HM12_UVB_z(6,j,2)  = (scale_low * HM12_UVB_carbon(i,j,2))  + (scale_high * HM12_UVB_carbon(i+1,j,2)) 
               end if
               if (j.lt.8) then
                  HM12_UVB_z(7,j,1)  = (scale_low * HM12_UVB_nitrogen(i,j,1))  + (scale_high * HM12_UVB_nitrogen(i+1,j,1))
                  HM12_UVB_z(7,j,2)  = (scale_low * HM12_UVB_nitrogen(i,j,2))  + (scale_high * HM12_UVB_nitrogen(i+1,j,2))                
               end if
               if (j.lt.9) then
                  HM12_UVB_z(8,j,1)  = (scale_low * HM12_UVB_oxygen(i,j,1))  + (scale_high * HM12_UVB_oxygen(i+1,j,1))
                  HM12_UVB_z(8,j,2)  = (scale_low * HM12_UVB_oxygen(i,j,2))  + (scale_high * HM12_UVB_oxygen(i+1,j,2))                
               end if
               if (j.lt.11) then
                  HM12_UVB_z(10,j,1)  = (scale_low * HM12_UVB_neon(i,j,1))  + (scale_high * HM12_UVB_neon(i+1,j,1))
                  HM12_UVB_z(10,j,2)  = (scale_low * HM12_UVB_neon(i,j,2))  + (scale_high * HM12_UVB_neon(i+1,j,2))                
               end if
               if (j.lt.13) then
                  HM12_UVB_z(12,j,1)  = (scale_low * HM12_UVB_magnesium(i,j,1))  + (scale_high * HM12_UVB_magnesium(i+1,j,1))
                  HM12_UVB_z(12,j,2)  = (scale_low * HM12_UVB_magnesium(i,j,2))  + (scale_high * HM12_UVB_magnesium(i+1,j,2))
               end if
               if (j.lt.15) then
                  HM12_UVB_z(14,j,1)  = (scale_low * HM12_UVB_silicon(i,j,1))  + (scale_high * HM12_UVB_silicon(i+1,j,1))
                  HM12_UVB_z(14,j,2)  = (scale_low * HM12_UVB_silicon(i,j,2))  + (scale_high * HM12_UVB_silicon(i+1,j,2))
               end if
               if (j.lt.17) then
                  HM12_UVB_z(16,j,1)  = (scale_low * HM12_UVB_sulfur(i,j,1))  + (scale_high * HM12_UVB_sulfur(i+1,j,1))
                  HM12_UVB_z(16,j,2)  = (scale_low * HM12_UVB_sulfur(i,j,2))  + (scale_high * HM12_UVB_sulfur(i+1,j,2))
               end if
               if (j.lt.7) then
                  HM12_UVB_z(26,j,1)  = (scale_low * HM12_UVB_iron(i,j,1))  + (scale_high * HM12_UVB_iron(i+1,j,1))
                  HM12_UVB_z(26,j,2)  = (scale_low * HM12_UVB_iron(i,j,2))  + (scale_high * HM12_UVB_iron(i+1,j,2))
               end if
             end do               
          end if
       end do
    end if
end subroutine update_UVB

end module photoionization_UVB_module