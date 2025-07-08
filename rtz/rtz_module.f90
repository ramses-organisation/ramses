module rtz_module
  use amr_parameters, only: dp
  implicit none

  private
  public:: elements, n_elements, initialize_elements

  type Element
      integer(KIND=4) :: atomic_number
      integer(KIND=4) :: n_ions
      integer(KIND=4) :: n_mol
      real(dp)   :: atomic_mass
      real(dp)   :: z_solar
      real(dp)   :: G0_photo_rate
      real(dp)   :: depletion
      character(LEN=20):: element_name
      character(LEN=2):: symbol
  end type Element

  ! RTZ STUFF
  integer, parameter :: n_elements = 27
  type(Element) :: elements(n_elements)

CONTAINS

SUBROUTINE initialize_elements()
   ! Initializes the atomic data we need for the RTZ module
   implicit none
   integer::i
   
   ! Initialize everything to zero
   do i=1,n_elements
      elements(i)%atomic_number = -1
      elements(i)%n_ions = -1
      elements(i)%atomic_mass = -1.0
      elements(i)%z_solar = -1.0
      elements(i)%G0_photo_rate = 0.0
      elements(i)%n_mol = 0
      elements(i)%depletion = 1.0
      elements(i)%element_name = "NANANANANANANANANANA"
      elements(i)%symbol = "XX"
   enddo

   ! Element 1: Hydrogen
   elements(1)%atomic_number = 1
   elements(1)%atomic_mass = 1.008
   elements(1)%z_solar = 1.0
   elements(1)%G0_photo_rate = 0.0 ! No subionizing PI
   elements(1)%n_ions = 2
#if N_H2 > 0
      elements(1)%n_mol = N_H2
#else
      elements(1)%n_mol = 0
#endif
   elements(1)%depletion = 1.0
   elements(1)%element_name = "HYDROGEN"
   elements(1)%symbol = "H"

#if N_HELIUM_IONS > 0
   ! Element 2: Helium
   elements(2)%atomic_number = 2
   elements(2)%atomic_mass = 4.0026
   elements(2)%z_solar = 8.51E-02
   elements(2)%G0_photo_rate = 0. ! No subionizing PI
   elements(2)%n_ions = N_HELIUM_IONS
   elements(2)%depletion = 1.0
   elements(2)%element_name = "HELIUM"
   elements(2)%symbol = "He"
#endif

#if N_CARBON_IONS > 0
   ! Element 6: Carbon
   elements(6)%atomic_number = 6
   elements(6)%atomic_mass = 12.0107
   elements(6)%z_solar = 2.69E-04
   elements(6)%G0_photo_rate = 3.39E-10
   elements(6)%n_ions = N_CARBON_IONS
   elements(6)%depletion = 0.5
   elements(6)%element_name = "CARBON"
   elements(6)%symbol = "C"
#endif
    
#if N_NITROGEN_IONS > 0
   ! Element 7: Nitrogen
   elements(7)%atomic_number = 7
   elements(7)%atomic_mass = 14.0067
   elements(7)%z_solar = 6.76E-05
   elements(7)%G0_photo_rate = 0.0 ! No subionizing PI
   elements(7)%n_ions = N_NITROGEN_IONS
   elements(7)%depletion = 0.6
   elements(7)%element_name = "NITROGEN"
   elements(7)%symbol = "N"
#endif 

#if N_OXYGEN_IONS > 0
   ! Element 8: Oxygen
   elements(8)%atomic_number = 8
   elements(8)%atomic_mass = 15.9994
   elements(8)%z_solar = 4.90E-04
   elements(8)%G0_photo_rate = 0.0 ! No subionizing PI
   elements(8)%n_ions = N_OXYGEN_IONS
   elements(8)%depletion = 0.73
   elements(8)%element_name = "OXYGEN"
   elements(8)%symbol = "O"
#endif

#if N_NEON_IONS > 0
   ! Element 10: Neon
   elements(10)%atomic_number = 10
   elements(10)%atomic_mass = 20.1797
   elements(10)%z_solar = 8.51E-05
   elements(10)%G0_photo_rate = 0.0 ! No subionizing PI
   elements(10)%n_ions = N_NEON_IONS
   elements(10)%depletion = 1.0
   elements(10)%element_name = "NEON"
   elements(10)%symbol = "Ne"
#endif

#if N_MAGNESIUM_IONS > 0
   ! Element 12: Magnesium
   elements(12)%atomic_number = 12
   elements(12)%atomic_mass = 24.305
   elements(12)%z_solar = 3.98E-05
   elements(12)%G0_photo_rate = 6.59E-11 
   elements(12)%n_ions = N_MAGNESIUM_IONS
   elements(12)%depletion = 0.16
   elements(12)%element_name = "MAGNESIUM"
   elements(12)%symbol = "Mg"
#endif

#if N_SILICON_IONS
   ! Element 14: Silicon
   elements(14)%atomic_number = 14
   elements(14)%atomic_mass = 28.0855
   elements(14)%z_solar = 3.24E-05
   elements(14)%G0_photo_rate = 4.47E-09
   elements(14)%n_ions = N_SILICON_IONS
   elements(14)%depletion = 0.1
   elements(14)%element_name = "SILICON"
   elements(14)%symbol = "Si"
#endif

#if N_SULFUR_IONS > 0
   ! Element 16: Sulfur
   elements(16)%atomic_number = 16
   elements(16)%atomic_mass = 32.065
   elements(16)%z_solar = 1.32E-05
   elements(16)%G0_photo_rate = 1.13E-09
   elements(16)%n_ions = N_SULFUR_IONS
   elements(16)%depletion = 1.0
   elements(16)%element_name = "SULFUR"
   elements(16)%symbol = "S"
#endif

#if N_IRON_IONS > 0
   ! Element 26: Iron
   elements(26)%atomic_number = 26
   elements(26)%atomic_mass = 55.854
   elements(26)%z_solar = 3.16E-05
   elements(26)%G0_photo_rate = 4.71E-10
   elements(26)%n_ions = N_IRON_IONS
   elements(26)%depletion = 0.01
   elements(26)%element_name = "IRON"
   elements(26)%symbol = "Fe"
#endif

END SUBROUTINE initialize_elements

end module rtz_module