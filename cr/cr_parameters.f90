module cr_parameters
  use amr_parameters

  ! Compile-time guards: this file is only compiled when CRPHYS is enabled
#ifndef CRPHYS
#error "cr/ sources require -DCRPHYS: set CRPHYS=1 in bin/Makefile"
#endif
#ifndef NCR_GROUPS
#error "cr/ sources require -DNCR_GROUPS=<ngroups>: set NCR_GROUPS in bin/Makefile"
#endif
#ifndef SOLVERmhd
#error "the CR module requires SOLVER=mhd (closure and source terms need B)"
#endif

  integer,parameter::ncr_groups=NCR_GROUPS         ! Number of CR groups
  integer,parameter::ncrvar=ncr_groups*(ndim+1)   ! Per group: E_cr + ndim fluxes
  ! Energy index of each CR group in cruold/crunew (mirrors rt iGroups);
  ! group g energy at Ecr_idx(g), fluxes at Ecr_idx(g)+1:Ecr_idx(g)+ndim.
  ! Filled in read_cr_params (cr_init.f90), like rt_init fills iGroups.
  integer,dimension(1:ncr_groups)::Ecr_idx=1

  ! --- &cr_params namelist ---------------------------------------------
  ! Master switch & scheme
  logical::cr_advect=.false.                ! Master CR transport switch
  logical::cr_HLLE=.true.                   ! HLLE Riemann solver for CR
  logical::cr_isotropic_pressure=.true.        ! .true.=P1 closure, .false.=M1
  logical::cr_flux_correction=.false. ! Rescale superluminal fluxes
  ! Physics
  real(dp),dimension(1:ncr_groups)::gamma_cr=4d0/3d0 ! CR adiabatic index
  logical::cr_feedback=.true.               ! master CR->gas back-reaction switch (gradPcr momentum + energy)
  real(dp)::cr_smallr_decouple=1d4          ! Decouple CRs where rho<smallr*this
  real(dp)::cr_efloor=1d-30                   ! CR energy floor
  ! Timestep / reduced light speed
  real(dp)::cr_c_fraction=1d0               ! Reduced light-speed fraction
  integer::cr_nsubcycle=1                   ! Max CR subcycles per MHD step
  logical::cr_varvmax=.false.               ! Adapt cr_vmax to gas/Alfven speed
  real(dp)::cr_varvmax_fudge=10d0           ! Fudge factor for adaptive cr_vmax
  logical::cr_varvmax_vdvs=.false.          ! Include diff/stream speeds in vmax
  ! Scattering / streaming source terms (used from phase 2)
  real(dp),dimension(1:ncr_groups)::Dcr=1.0d29     ! Diffusion coefficient [cm^2/s]
  real(dp)::DCRmax=1d30                     ! Max CR streaming diffusion coeff [cm^2/s] (0 would disable streaming diffusion via 1/DCRmax_code=Inf)
  real(dp),dimension(1:ncr_groups)::Dcr_perp_factor=1d-6 ! Perpendicular suppression
  logical::cr_streaming_diffusion=.false.  ! Streaming term in sigma
  logical::cr_streaming_heating=.false.    ! Streaming heating of gas
  real(dp)::cr_v_alfven=0d0                    ! Imposed Alfven speed (tests)
  real(dp)::cr_f_taucell=1d0                ! Cell optical-depth stability factor
  ! Cooling: Fitz Axen et al. 2024 (used from phase 4)
  logical::cr_cooling=.false.               ! CR collisional cooling
  real(dp)::zeta_cr=7.51d-16                ! Coulomb loss rate coefficient
  real(dp)::cr_ne=1d-3                         ! Free electrons per H nucleus
  real(dp)::cr_fneut=0.875d0                   ! Neutral gas fraction
  ! Boundaries / tests
  character(LEN=32)::jiang_test=''          ! Jiang & Oh test IC/BC dispatch
  ! Region-based CR initial conditions (per &init_params region geometry).
  ! cr_region_u(k,1)=E_cr, cr_region_u(k,2:ncrvar)=CR flux in region k.
  real(dp),dimension(1:MAXREGION,1:ncrvar)::cr_region_u=0d0
  ! Per-boundary CR state for imposed (bound_type=3) boundaries: cr_boundary_u(b,1)=E_cr,
  ! cr_boundary_u(b,2:ncrvar)=CR flux on boundary region b. Applied in cr_boundana
  ! for the tp_* two-pressure shock tests.
  real(dp),dimension(1:MAXBOUND,1:ncrvar)::cr_boundary_u=0d0
  ! CR-owned IC region geometry: lets CR regions differ from gas regions;
  ! cr_reg_group selects the target group (inert when NCR_GROUPS=1).
  integer::cr_nregion=0
  character(LEN=10),dimension(1:MAXREGION)::cr_region_type='square'
  real(dp),dimension(1:MAXREGION)::cr_reg_x_center=0d0
  real(dp),dimension(1:MAXREGION)::cr_reg_y_center=0d0
  real(dp),dimension(1:MAXREGION)::cr_reg_z_center=0d0
  real(dp),dimension(1:MAXREGION)::cr_reg_length_x=1d10
  real(dp),dimension(1:MAXREGION)::cr_reg_length_y=1d10
  real(dp),dimension(1:MAXREGION)::cr_reg_length_z=1d10
  real(dp),dimension(1:MAXREGION)::cr_exp_region=2d0
  integer,dimension(1:MAXREGION)::cr_reg_group=1
  ! Refinement
  real(dp),dimension(1:ncrvar)::err_grad_cr=-1d0 ! CR gradient refinement
  ! Output

  ! --- Derived / bookkeeping (not in the namelist) ----------------------
  real(dp),dimension(1:MAXLEVEL)::cr_vmax=0d0 ! Reduced light speed, code units
  real(dp)::cr_va_max=0d0                   ! Max Alfven speed (adaptive cr_vmax)
  real(dp)::cr_c_code=0d0                        ! Light speed in code units
  real(dp),dimension(1:ncr_groups)::DCR_code=0d0   ! Dcr in code units
  real(dp)::DCRmax_code=0d0                 ! DCRmax in code units
  real(dp)::cr_smalld=1d-25                  ! Floor on DCR_code
  real(dp)::ecrs_tot=0d0                    ! Total CR energy (log diagnostic)
  ! Per-cell CR-energy gather buffer for cmpdt's CR-pressure term: filled by
  ! courant_fine, read by cmpdt. Module-level so cmpdt takes no extra argument.
  real(dp),dimension(1:nvector,1:ncr_groups)::cr_egather=0d0

end module cr_parameters
