module nimhd_parameters
  use amr_parameters
  
  ! what are these indices useful for? Can't we just use these directy?
  integer:: nxx=1
  integer:: nyy=2
  integer:: nzz=3

  logical :: nambipolar = .false. ! Ambipolar diffusion ?
  logical :: nmagdiffu  = .false. ! Magnetic  diffusion ?
  logical :: use_nonideal_mhd     ! true if any of the effects is used

! Mellon & Li 2009 (?) or Hennebelle & Teyssier 2007
! WARNING this value is in CGS. The connection with user units
! is made in function gammaadbis in umsucl
! rework to choose from fixed value, model or table resistivity (see read params)
  real(dp):: gammaAD=1.0d0
  real(dp):: rho_threshold=1d-10     ! safeguard for the ambipolar flux in high density contrast cases (in code units)
  ! useful to restrict ambipolar diff in low rho regions
  ! rename rho_min_AD or something

  ! Resistivity parameters
  ! resistivity_method=0  -> use fixed value
  ! resistivity_method=1  -> use analytical model
  ! resistivity_method=2  -> use tabulated values

  integer :: resistivity_method=0    ! How to determine the resistivity
  integer :: resistivity_table_ndim=0       ! number of variables to extrapolate in table: rho, T, Xi
  character(len=80)::res_table_name='marchand2016_table.dat' !filename of the table

  integer :: use_x2d=0               ! use abundances
  integer :: use_x3d=0               ! use abundance table with rho, T, Xi
  integer :: use_res=0               ! use resistivities

! magnetic diffusion coefficient see function etaohmdiss in umsucl
  real(dp):: etaMD=1d0

  real(dp), parameter:: H2_fraction = 0.844d0 !remove ! H2 fraction in number of particules (equals 0.73 in mass)
! WARNING !! Think to change xmolaire if proportion are changed
  !real(dp):: xmneutre=(H2_fraction*2d0 +(1d0-H2_fraction)*4.)*1.667d-24 !remove

! Mellon & Li 2009 (?) or Hennebelle & Teyssier 2007
! WARNING this value is in CGS. The connection with user units
! is made in function densionbis in umsucl
  real(dp):: coefionis=3d-16 !remove ! coefionis*sqrt(n_H)=n_i , empirical value from Shu book 2, p. 363
  real(dp):: coefad = 0.1d0 !CFL conditions
  real(dp):: coefohm = 0.05d0 !CFL conditions
  integer :: nminitimestep = 0 ! flag to activate timestep reduction hack
  real(dp):: coefalfven = 1d-10  ! threshold where you cap the timestep. Maximal ratio between AD timestep and ideal MHD timestep. Change name!
  ! Meme si ca n'a rien a voir avec alfven : c'est le coefficient de seuil. Par defaut, on ne seuille pas.
  real(dp):: coefdtohm = 1d-10  ! same as coefalfven but for Ohmic diff.
  ! TC:why have different courant conditions for different nimhd effects?
  real(dp):: default_ionisrate=1d-17
  ! choose better names
  
end module nimhd_parameters
