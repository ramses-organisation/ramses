!################################################################
!################################################################
subroutine read_cr_params(nml_ok)
  ! Read the &CR_PARAMS namelist block (cosmic-ray runtime parameters)
  ! from the parameter file already open on unit 1, and sanity-check it.
  use amr_commons
  use cr_parameters
  implicit none
  logical::nml_ok

  namelist/cr_params/cr_advect,cr_HLLE,cr_use_minmod,isotropic_pressure &
       & ,reduced_CR_flux_correction,cr_interpolation,gamma_cr,gradpcr_mom &
       & ,cr_smallr_decouple,smallcr,cr_c_fraction,cr_nsubcycle &
       & ,cr_varvmax,cr_varvmax_fudge,cr_varvmax_vdvs &
       & ,Dcr,DCRmax,Dcr_perp_factor,mom_streaming_diffusion &
       & ,mom_streaming_heating,v_alfven,cr_f_taucell &
       & ,cr_cooling,zeta_cr,ne,fneut &
       & ,cr_bound_floor,jiang_test,err_grad_crmom,cr_legacy_output &
       & ,crmom_region

  if(myid==1)write(*,'(A50)')"Reading cr_params namelist ..."
  rewind(1)
  read(1,NML=cr_params,END=101)
101 continue

  ! Sanity checks
  if(cr_nsubcycle<1)then
     if(myid==1)write(*,*)'Error in cr_params: cr_nsubcycle must be >= 1'
     nml_ok=.false.
  endif
  if(cr_c_fraction<=0d0 .or. cr_c_fraction>1d0)then
     if(myid==1)write(*,*)'Error in cr_params: cr_c_fraction must be in (0,1]'
     nml_ok=.false.
  endif

end subroutine read_cr_params
