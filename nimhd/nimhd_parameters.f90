module nimhd_parameters
  use amr_parameters

  logical :: nambipolar = .false. ! flag to activate ambipolar diffusion
  logical :: nmagdiffu  = .false. ! flag to activate magnetic  diffusion
  logical :: nimhdheating_in_flux = .true. ! flag to activate the non-ideal mhd energy fluxes
  logical :: nimhdheating_source_term = .false. ! flag to activate the non-ideal mhd heating as a source term

  logical :: use_nonideal_mhd = .false.     ! true if any of the non-ideal MHD effects is used

  ! Fixed resistivity coefficients (in code units)
  real(dp):: gammaAD=1.0d0            ! ambipolar diffusion coefficient (beta = 1/(gammaAD*rho))
  real(dp):: etaMD=1d0                ! Ohmic (magnetic) diffusion coefficient

  ! timestep regulation
  real(dp):: coefad = 0.1d0    ! CFL condition for ambipolar diffusion
  real(dp):: coefohm = 0.05d0  ! CFL condition for ohmic dissipation
  real(dp):: coefalfven = 1d-10  ! threshold where you cap the timestep. Maximal ratio between AD timestep and ideal MHD timestep. Change name!
  ! Meme si ca n'a rien a voir avec alfven : c'est le coefficient de seuil. Par defaut, on ne seuille pas.
  real(dp):: coefdtohm = 1d-10  ! same as coefalfven but for Ohmic diff.
  logical :: nminitimestep = .false. ! flag to activate timestep reduction hack TODO change name


end module nimhd_parameters
