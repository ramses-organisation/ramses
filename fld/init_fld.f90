subroutine init_fld
   use amr_commons, only:ncoarse
   use amr_parameters, only:twotondim,ngridmax,ndim
   use hydro_parameters, only:ngrp
   use fld_commons
   use radiation_parameters
   implicit none
   !------------------------------
   ! Initialize variables for FLD
   !------------------------------
   integer::ncell

   ! allocate global arrays
   ncell=ncoarse+twotondim*ngridmax

   if(fld)then
      allocate(rad_flux(1:ncell,1:nvar_bicg))
      allocate(frad(1:ncell,1:ndim))
      allocate(in_sink(1:ncell))
      rad_flux=0.0d0; frad=0.0d0; in_sink=.false.
   endif

   ! Variables for BICG scheme
   ! 1 : r
   ! 2 : p
   ! 3 : r*
   ! 4 : M-1
   ! 5 : 
   ! 6 : z and Ap
   ! 7 : p*
   ! 8 : p*A
   ! 9 : z*
   allocate(kappaR_bicg(1:ncell,1:ngrp))
   ! if FLD: matrix of size ngrpxngrp (because matrix only on Eg)
   ! if  M1: matrix of size (1+nrad)x(1+nrad) (on T,Eg,Fg)
   allocate(var_bicg(1:ncell,1:nvar_bicg,1:10+2*ndim))
   allocate(precond_bicg(1:ncell,1:nvar_bicg,1:nvar_bicg))
   if(store_matrix) then
      allocate(mat_residual_glob(1:ncell,1:nvar_bicg,1:nvar_bicg),residual_glob(1:ncell,1:nvar_bicg))
      allocate(coeff_glob_left(1:ncell,1:nvar_bicg,1:nvar_bicg,1:ndim),coeff_glob_right(1:ncell,1:nvar_bicg,1:nvar_bicg,1:ndim))
   else
      allocate(mat_residual_glob(1,1:nvar_bicg,1:nvar_bicg),residual_glob(1,1:nvar_bicg))
      allocate(coeff_glob_left(1,1:nvar_bicg,1:nvar_bicg,1:ndim),coeff_glob_right(1,1:nvar_bicg,1:nvar_bicg,1:ndim))
   endif
   kappar_bicg=0.0d0;var_bicg=0.0d0;precond_bicg=0.0d0
   mat_residual_glob=0.0d0;residual_glob=0.0d0
   coeff_glob_left=0.0d0;coeff_glob_right=0.0d0

end subroutine init_fld