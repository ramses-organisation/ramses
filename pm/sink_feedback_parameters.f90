module sink_feedback_parameters
  use amr_parameters,only:dp

  ! STELLAR_PARAMS namelist: parameters for spawning stellar particles

  integer:: nstellarmax                   ! maximum number of stellar objects
  real(dp):: imf_index, imf_low, imf_high ! power-law IMF model: PDF index (dN/dM), lower and higher mass bounds (Msun)
  ! Stellar lifetime model: t(M) = lt_t0 * exp(lt_a * (log(lt_m0 / M))**lt_b)
  ! default: Woosley et al 2002
  real(dp):: lt_t0=3.26515722d0   !Myr
  real(dp):: lt_m0=148.158972d0   !Msun
  real(dp):: lt_a=0.23840797d0
  real(dp):: lt_b=2.20522946d0

  character(LEN=100)::stellar_strategy='local' ! local: create stellar particles from each sink
                                               ! global: create when the total mass in sinks exceeds stellar_msink_th
  real(dp):: stellar_msink_th                  ! sink mass threshold for stellar object creation (Msun)

  ! Allow users to pre-set stellar mass selection for physics comparison runs, etc
  ! Every time mstellar is added to, instead of a random value, use mstellarini
  integer,parameter::nstellarini=100
  real(dp),dimension(nstellarini)::mstellarini ! List of stellar masses to use

  ! STELLAR_PARAMS namelist: SN feedback parameters

  logical::sn_feedback_sink = .false. !SN feedback emanates from the sink
  logical::sn_direct = .false.        ! explode immediately instead of after lifetime

  real(dp):: sn_e_ref=1.d51      ! SN energy for forcing by sinks [erg]
  real(dp):: sn_p_ref=4.d43      ! SN momentum [g cm/s] for 10 H/cc (Iffrig and Hennebelle 2015)

  real(dp):: Tsat=1d99    ! maximum temperature in SN remnants
  real(dp):: Vsat=1d99    ! maximum velocity in SN remnants
  real(dp):: sn_r_sat=0d0 ! minimum radius for SN remnant

  real(dp):: Vdisp=1d0    ! dispersion velocity of the stellar objects [km/s]
                          ! determines how far SN can explode from the sink

  logical::stellar_info=.true.  ! write stellar particles to log file

  ! STELLAR_PARAMS namelist: HII feedback parameters

  ! Stellar ionizing flux model: S(M) = stf_K * (M / stf_m0)**stf_a / (1 + (M / stf_m0)**stf_b)**stf_c
  ! This is a fit from Vacca et al. 1996
  ! Corresponding routine : vaccafit
  real(dp)::stf_K=9.634642584812752d48 ! s**(-1) then normalised in code units in read_stellar
  real(dp)::stf_m0=2.728098824280431d1 ! Msun then normalised in code units in read_stellar
  real(dp)::stf_a=6.840015602892084d0
  real(dp)::stf_b=4.353614230584390d0
  real(dp)::stf_c=1.142166657042991d0

  real(dp):: hii_t=0 !fiducial HII region lifetime [yr?], it is normalised in code units in read_stellar
  integer:: feedback_photon_group=-1 ! index of the photon group where to put the radiation

  ! commons

  ! stellar object arrays
  integer:: nstellar = 0 ! current number of stellar objects
  real(dp), allocatable, dimension(:):: mstellar, tstellar, ltstellar ! mass, birth time, life time
  integer, allocatable, dimension(:):: id_stellar                     !the id  of the sink to which it belongs

end module sink_feedback_parameters
