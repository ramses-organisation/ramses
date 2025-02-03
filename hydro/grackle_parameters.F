      module grackle_parameters

      USE ISO_C_BINDING
      implicit None
#include "grackle.def"

      TYPE (grackle_units) :: my_grackle_units
      TYPE (grackle_field_data) :: my_grackle_fields
      TYPE (grackle_chemistry_data) :: my_grackle_data

      R_PREC, TARGET :: gr_density(NVECTOR), gr_energy(NVECTOR),
     &     gr_x_velocity(NVECTOR), gr_y_velocity(NVECTOR),
     &     gr_z_velocity(NVECTOR),
     &     gr_HI_density(NVECTOR), gr_HII_density(NVECTOR),
     &     gr_HM_density(NVECTOR),
     &     gr_HeI_density(NVECTOR), gr_HeII_density(NVECTOR),
     &     gr_HeIII_density(NVECTOR),
     &     gr_H2I_density(NVECTOR), gr_H2II_density(NVECTOR),
     &     gr_DI_density(NVECTOR), gr_DII_density(NVECTOR),
     &     gr_HDI_density(NVECTOR),
     &     gr_e_density(NVECTOR), gr_metal_density(NVECTOR),
     &     gr_volumetric_heating_rate(NVECTOR),
     &     gr_specific_heating_rate(NVECTOR),
     &     gr_RT_HI_ionization_rate(NVECTOR),
     &     gr_RT_HeI_ionization_rate(NVECTOR),
     &     gr_RT_HeII_ionization_rate(NVECTOR),
     &     gr_RT_H2_dissociation_rate(NVECTOR),
     &     gr_RT_heating_rate(NVECTOR)

      INTEGER, TARGET :: gr_rank,gr_dimension(3),gr_start(3),gr_end(3)
      INTEGER :: iresult

      INTEGER :: use_grackle = 1
      INTEGER :: grackle_with_radiative_cooling = 1
      INTEGER :: grackle_primordial_chemistry = 0
      INTEGER :: grackle_metal_cooling = 0
      INTEGER :: grackle_UVbackground = 0
      INTEGER :: grackle_cmb_temperature_floor = 1
      INTEGER :: grackle_h2_on_dust = 0
      INTEGER :: grackle_photoelectric_heating = 0
      INTEGER :: grackle_use_volumetric_heating_rate = 0
      INTEGER :: grackle_use_specific_heating_rate = 0
      INTEGER :: grackle_three_body_rate = 0
      INTEGER :: grackle_cie_cooling = 0
      INTEGER :: grackle_h2_optical_depth_approximation = 0
      INTEGER :: grackle_ih2co = 1
      INTEGER :: grackle_ipiht = 1
      INTEGER :: grackle_NumberOfTemperatureBins = 600
      INTEGER :: grackle_CaseBRecombination = 0
      INTEGER :: grackle_Compton_xray_heating = 0
      INTEGER :: grackle_LWbackground_sawtooth_suppression = 0
      INTEGER :: grackle_NumberOfDustTemperatureBins = 250
      INTEGER :: grackle_use_radiative_transfer = 0
      INTEGER :: grackle_radiative_transfer_coupled_rate_solver = 0
      INTEGER :: grackle_radiative_transfer_intermediate_step = 0
      INTEGER :: grackle_radiative_transfer_hydrogen_only = 0
      INTEGER :: grackle_self_shielding_method = 0
      REAL    :: grackle_Gamma = 5.d0/3.d0
      REAL    :: grackle_photoelectric_heating_rate = 8.5D-26
      REAL    :: grackle_HydrogenFractionByMass = 0.76d0
      REAL    :: grackle_DeuteriumToHydrogenRatio = 2.0d0*3.4d-5
      REAL    :: grackle_SolarMetalFractionByMass = 0.01295d0
      REAL    :: grackle_TemperatureStart = 1.0d0
      REAL    :: grackle_TemperatureEnd = 1.0D9
      REAL    :: grackle_DustTemperatureStart = 1.0d0
      REAL    :: grackle_DustTemperatureEnd = 1500.0d0
      REAL    :: grackle_LWbackground_intensity = 0.0d0
      REAL    :: grackle_UVbackground_redshift_on = 7.0d0
      REAL    :: grackle_UVbackground_redshift_off = 0.0d0
      REAL    :: grackle_UVbackground_redshift_fullon = 6.0d0
      REAL    :: grackle_UVbackground_redshift_drop = 0.0d0
      REAL    :: grackle_cloudy_electron_fraction_factor = 9.153959D-3
      CHARACTER(LEN=256), TARGET :: grackle_data_file

      end module grackle_parameters
