module cr_parameters
  ! Cosmic-ray module compile-time configuration, namelist parameters and
  ! global bookkeeping. Leaf module: may only use amr_parameters.
  ! Parameter names and defaults are kept verbatim from ramses_cral
  ! feat/CR_tests (mhd/hydro_parameters.f90) so cral namelists port directly.
  use amr_parameters

  ! Compile-time guards: this file is only compiled when NCR>0
#ifndef NCR
#error "cr/ sources require -DNCR=<ngroups>: set NCR in bin/Makefile"
#endif
#ifndef CRFLX
#error "the CR module is two-moment only: NCR>0 requires CRFLX=1 in bin/Makefile"
#endif
#ifndef SOLVERmhd
#error "the CR module requires SOLVER=mhd (closure and source terms need B)"
#endif

  integer,parameter::ncr=NCR                ! Number of CR groups
  integer,parameter::ncrvars=ncr*(ndim+1)   ! Per group: E_cr + ndim fluxes
  ! CR base index in cruold/crunew. Group g: energy at iCRu+(ndim+1)*(g-1),
  ! fluxes in the ndim slots after it. Never reference nvar in CR indexing.
  integer,parameter::iCRu=1

  ! --- &cr_params namelist ---------------------------------------------
  ! Master switch & scheme
  logical::cr_advect=.false.                ! Master CR transport switch
  logical::cr_HLLE=.true.                   ! HLLE Riemann solver for CR
  logical::cr_use_minmod=.false.            ! Minmod slope limiter
  logical::isotropic_pressure=.true.        ! .true.=P1 closure, .false.=M1
  logical::reduced_CR_flux_correction=.false. ! Rescale superluminal fluxes
  logical::cr_interpolation=.true.          ! Interpolate CR vars on AMR
  ! Physics
  real(dp),dimension(1:ncr)::gamma_cr=4d0/3d0 ! CR adiabatic index
  logical::gradpcr_mom=.true.               ! gradP_cr form of gas back-reaction
  real(dp)::cr_smallr_decouple=1d4          ! Decouple CRs where rho<smallr*this
  real(dp)::smallcr=1d-30                   ! CR energy floor
  ! Timestep / reduced light speed
  real(dp)::cr_c_fraction=1d0               ! Reduced light-speed fraction
  integer::cr_nsubcycle=1                   ! Max CR subcycles per MHD step
  logical::cr_varvmax=.false.               ! Adapt cr_vmax to gas/Alfven speed
  real(dp)::cr_varvmax_fudge=10d0           ! Fudge factor for adaptive cr_vmax
  logical::cr_varvmax_vdvs=.false.          ! Include diff/stream speeds in vmax
  ! Scattering / streaming source terms (used from phase 2)
  real(dp),dimension(1:ncr)::Dcr=1.0d29     ! Diffusion coefficient [cm^2/s]
  real(dp)::DCRmax=0d0                      ! Max diffusion coefficient [cm^2/s]
  real(dp),dimension(1:ncr)::Dcr_perp_factor=1d-6 ! Perpendicular suppression
  logical::mom_streaming_diffusion=.false.  ! Streaming term in sigma
  logical::mom_streaming_heating=.false.    ! Streaming heating of gas
  real(dp)::v_alfven=0d0                    ! Imposed Alfven speed (tests)
  real(dp)::cr_f_taucell=1d0                ! Cell optical-depth stability factor
  ! Cooling: Fitz Axen et al. 2024, as implemented in cral cr_cooling_fine.f90
  ! (used from phase 4)
  logical::cr_cooling=.false.               ! CR collisional cooling
  real(dp)::zeta_cr=7.51d-16                ! Coulomb loss rate coefficient
  real(dp)::ne=1d-3                         ! Free electrons per H nucleus
  real(dp)::fneut=0.875d0                   ! Neutral gas fraction
  ! Boundaries / tests
  real(dp)::cr_bound_floor=-1d0             ! >=0: E_cr in reflexive boundaries
  character(LEN=32)::jiang_test=''          ! Jiang & Oh test IC/BC dispatch
  ! Refinement
  real(dp),dimension(1:ncrvars)::err_grad_crmom=-1d0 ! CR gradient refinement
  ! Output
  logical::cr_legacy_output=.false.         ! .true.: CR columns in hydro files

  ! --- Derived / bookkeeping (not in the namelist) ----------------------
  real(dp),dimension(1:MAXLEVEL)::cr_vmax=0d0 ! Reduced light speed, code units
  real(dp)::cr_vgas_max=0d0                 ! Max gas speed (adaptive cr_vmax)
  real(dp)::cr_va_max=0d0                   ! Max Alfven speed (adaptive cr_vmax)
  real(dp)::c_cu=0d0                        ! Light speed in code units
  real(dp),dimension(1:ncr)::DCR_code=0d0   ! Dcr in code units
  real(dp)::DCRmax_code=0d0                 ! DCRmax in code units
  real(dp)::ecrs_tot=0d0                    ! Total CR energy (log diagnostic)

end module cr_parameters
