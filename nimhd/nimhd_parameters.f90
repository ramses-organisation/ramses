module nimhd_parameters
  use amr_parameters
  
  ! TC: Why not just use the value directly?
  integer:: nxx=1
  integer:: nyy=2
  integer:: nzz=3

  logical :: nambipolar = .false. ! flag to activate ambipolar diffusion
  logical :: nmagdiffu  = .false. ! flag to activate magnetic  diffusion
  logical :: use_nonideal_mhd     ! true if any of the non-ideal MHD effects is used

  ! Resistivity parameters
  integer :: resistivity_method=0     ! How to determine the resistivity
                                      ! resistivity_method=0  -> use fixed value
                                      ! resistivity_method=1  -> use analytical model
                                      ! resistivity_method=2  -> use tabulated values from resnh.dat
  ! Mellon & Li 2009 (?) or Hennebelle & Teyssier 2007
  ! WARNING this value is in CGS. The connection with user units
  ! is made in function gammaadbis
  real(dp):: gammaAD=1.0d0
  real(dp):: rho_threshold=1d-10     ! safeguard for the ambipolar flux in high density contrast cases (in code units)
  ! useful to restrict ambipolar diff in low rho regions
  ! rename rho_min_AD or something
  real(dp):: etaMD=1d0                ! fixed magnetic diffusion coefficient

  real(dp), parameter:: H2_fraction = 0.844d0 !remove ! H2 fraction in number of particules (equals 0.73 in mass)
! WARNING !! Think to change xmolaire if proportion are changed
  
  ! timestep regulation
  real(dp):: coefad = 0.1d0    ! CFL condition for ambipolar diffusion
  real(dp):: coefohm = 0.05d0  ! CFL condition for ohmic dissipation
  real(dp):: coefalfven = 1d-10  ! threshold where you cap the timestep. Maximal ratio between AD timestep and ideal MHD timestep. Change name!
  ! Meme si ca n'a rien a voir avec alfven : c'est le coefficient de seuil. Par defaut, on ne seuille pas.
  real(dp):: coefdtohm = 1d-10  ! same as coefalfven but for Ohmic diff.
  logical :: nminitimestep = .false. ! flag to activate timestep reduction hack TODO change name


end module nimhd_parameters
