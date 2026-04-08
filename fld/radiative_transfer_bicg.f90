subroutine diffusion_cg (ilevel,Nsub)
  use amr_commons,only:myid,numbtot,active,son,ncpu,reception,dtnew,ncoarse,nstep
  use amr_parameters, only : verbose, ndim
  use cloud_module, only: rt_feedback
  use hydro_commons
  use radiation_parameters
  use const
  use units_commons
#ifdef RT
  use pm_commons !hybrid RT, needs sink quantities                                                  
  use rt_hydro_commons !hybrid RT
  use rt_cooling_module !hybrid RT
#endif
  use mpi_mod
  implicit none
  !=========================================================
  ! Iterative solver with Stabilised Bi-Conjugate Gradient method
  ! to solve A x = b
  !  i   : cell index
  !  irad: radiative variable index (from 1 to ngrp if FLD, from 1 to (1+ndim)*ngrp if M1)
  !
  !   r         : stored in var_bicg(i,irad,1)
  !   p         : stored in var_bicg(i,irad,2)
  !   v         : stored in var_bicg(i,irad,3)
  !   K^{-1}    : stored in var_bicg(i,irad,4)
  !   y         : stored in var_bicg(i,irad,5)
  !   z         : stored in var_bicg(i,irad,6)
  !   s         : stored in var_bicg(i,irad,7)
  !   t         : stored in var_bicg(i,irad,8) 
  !   rbar_0    : stored in var_bicg(i,irad,9)
  !   K^{-1}*t  : stored in var_bicg(i,irad,10)
  !
  !  radflux (cell_left ,idim=1)   : stored in var_bicg(i,irad,11)
  !  radflux (cell_right,idim=1)   : stored in var_bicg(i,irad,12)
  !  radflux (cell_left ,idim=2)   : stored in var_bicg(i,irad,13)
  !  radflux (cell_right,idim=2)   : stored in var_bicg(i,irad,14)
  !  radflux (cell_left ,idim=3)   : stored in var_bicg(i,irad,15)
  !  radflux (cell_right,idim=3)   : stored in var_bicg(i,irad,16)
  !
  !  new radiative energy at time n+1 : stored in unew(i,irad)
  !      radiative energy at time n   : stored in uold(i,irad)
  !
  !Tgas (iteration) : stored in unew(i,nvar)
  !Tgas (old)       : stored in uold(i,nvar)
  !
  !
  !=========================================================
  integer,intent(IN)::ilevel,Nsub
  complex*16 :: final_sum
  real(dp)::error,error_ini,epsilon
  real(dp)::Cv,told,rho,dt_exp,wdtB,wdtE,Tr,Trold,cal_Teg
  real(dp)::r2,rhs_norm1,r3
  real(dp)::temp,density,planck_ana,rosseland_ana
  integer::i,ind,iter,iskip,itermax,icpu,igroup,igrp,irad,jrad,ivar
  integer::this,nleaf_tot
  real(dp)::radiation_source,deriv_radiation_source,rhs,lhs

  real(dp)::rho_bicg_new,rho_bicg_old,alpha_bicg,omega_bicg,beta_bicg

  real(dp)::max_loc
  integer ::im1 !hybrid RT
#ifndef WITHOUTMPI
  integer::info,nleaf_all
  real(dp)::max_loc_all
#endif

  logical::exist_leaf_cell=.true.,debug_energy=.false.

  integer::nx_loc
  real(dp)::scale,dx,dx_loc

    real(dp)::ambi_heating,ohm_heating,nimhd_heating,protostellar_heating
!#if NIMHD==1
!  real(dp)::bcell2,bx,by,bz,jsquare,jx,jy,jz,etaohmdiss,betaad,ionisrate
!#endif

  real(dp)::scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l

#ifdef RT
  real(dp)::scale_Np,scale_Fp !hybrid RT
  call rt_units(scale_Np,scale_Fp)
#endif

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  
  if(myid==1 .and. (mod(nstep,ncontrol)==0)) write(*,*) 'entering radiative transfer for level ',ilevel

  if(bicg_to_cg)then
     block_diagonal_precond_bicg=.false.
     i_rho  = 6
     i_beta = 6
     i_y    = 2
     i_pAp  = 2
     i_s    = 1
  else
     block_diagonal_precond_bicg=.true.
     i_rho  = 9
     i_beta = 1
     i_y    = 5
     i_pAp  = 9
     i_s    = 7
  endif

  if(verbose)write(*,111)
  if(numbtot(1,ilevel)==0)return

  ! Rescaling factors
  ! Mesh size at level ilevel
  dx=half**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  allocate(liste_ind (1:twotondim*active(ilevel)%ngrid))

  nb_ind = 0

  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,active(ilevel)%ngrid
        if(son(active(ilevel)%igrid(i)+iskip) == 0)then
           nb_ind = nb_ind+1 
           liste_ind(nb_ind) = active(ilevel)%igrid(i)+iskip
        end if
     end do
  end do

  do irad=1,nvar_bicg
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do icpu=1,ncpu
           do i=1,reception(icpu,ilevel)%ngrid
              rad_flux(reception(icpu,ilevel)%igrid(i)+iskip,irad)=zero
           end do
           do i=1,reception(icpu,ilevel-1)%ngrid
              rad_flux(reception(icpu,ilevel-1)%igrid(i)+iskip,irad)=zero
           end do
        end do
     end do
  enddo

  if(debug_energy)then
  write(*,*) 'At the beginning of uplmde - uold(5-9)-unew(5-9)'
  do i=1,nb_ind
     this = liste_ind(i)
     write(*,'(12(ES15.6))') uold(this,5),uold(this,nvar),uold(this,9),uold(this,10),unew(this,5),unew(this,nvar),unew(this,9),unew(this,10)
  enddo
  read(*,*)
  endif

  !===================================================================
  ! Begin of subcycles....
  !===================================================================
  dt_exp = dtnew(ilevel)
  dt_imp = dtnew(ilevel)

  if (nb_ind == 0)then
     !print*,'No leaf-cell - myid=',myid
     exist_leaf_cell=.false.
  end if

  nleaf_tot=nb_ind
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(nb_ind,nleaf_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  nleaf_tot=nleaf_all
#endif
     
  if(nleaf_tot .eq. 0)then
     !write(*,*)'No leaf cells at level',ilevel,'. Exiting BiCG'
     deallocate(liste_ind)
     return
  end if

  do i=1,nb_ind
     this = liste_ind(i)

     var_bicg(this,:,:)=zero
     if(block_diagonal_precond_bicg) precond_bicg(this,:,:)=zero

     do irad=1,nvar_bicg
        unew(this,ind_bicg(irad))=zero
     enddo
     unew(this,nvar+1)=zero
#if USE_FLD==1
     kappaR_bicg(this,:)=zero
#endif
  end do

  ! Set constants
  epsilon = epsilon_diff

  !===================================================================
  ! Compute gas temperature stored in uold(i,nvar) and in unew(i,nvar)
  !===================================================================
  call cmp_energy(1)

  if(debug_energy)then
  write(*,*) 'After cmp_energy(1) - uold(5-9)-unew(5-9)'
  do i=1,nb_ind
     this = liste_ind(i)
     write(*,'(12(ES15.6))') uold(this,5),uold(this,nvar),uold(this,9),uold(this,10),unew(this,5),unew(this,nvar),unew(this,9),unew(this,10)
  enddo
  read(*,*)
  endif

#if USE_FLD==1
  do i=1,nb_ind
     this = liste_ind(i)

     density = scale_d * max(uold(this,1),smallr)
     temp = uold(this,nvar)*Tr_floor

     ! Compute Rosseland opacity (Compute kappa*rho)
     do igroup=1,ngrp
        Tr = cal_Teg(uold(this,firstindex_er+igroup)*scale_E0,igroup)
        if(rt_feedback)then
           kappaR_bicg(this,igroup)= rosseland_ana(density,temp,temp,igroup,in_sink(this)) / scale_kappa
        else
           kappaR_bicg(this,igroup)= rosseland_ana(density,temp,tr,igroup,in_sink(this)) / scale_kappa
        end if
        if(kappaR_bicg(this,igroup)*dx_loc .lt. min_optical_depth) kappaR_bicg(this,igroup)=min_optical_depth/dx_loc
     enddo
  end do
#endif

  ! Update boundaries
#if USE_FLD==1
  call make_virtual_fine_dp(uold(1,nvar),ilevel)
  call make_virtual_fine_dp(unew(1,nvar),ilevel)
  do igrp=1,ngrp
     call make_virtual_fine_dp(kappaR_bicg(1,igrp),ilevel)
  enddo
#endif
  call make_virtual_fine_dp(unew(1,nvar+1),ilevel)

  call make_virtual_fine_dp(uold(1,5),ilevel)
  call make_virtual_fine_dp(unew(1,5),ilevel)
  do irad=1,nvar_bicg
     call make_virtual_fine_dp(uold(1,ind_bicg(irad)),ilevel)
     call make_virtual_fine_dp(unew(1,ind_bicg(irad)),ilevel)

     do ivar=1,10+2*ndim
        call make_virtual_fine_dp(var_bicg(:,irad,ivar),ilevel)
     enddo

     if(block_diagonal_precond_bicg) then
        do ivar=1,nvar_bicg
           call make_virtual_fine_dp(precond_bicg(:,irad,ivar),ilevel)
        enddo
     endif
  enddo

  call make_boundary_diffusion_tot(ilevel)

  !===========================================
  ! Compute the matrix and vector coefficients
  !===========================================
#if USE_FLD==1
  call cmp_matrix_and_vector_coeff_fld(ilevel)
#endif
#if USE_M_1==1
  call cmp_matrix_and_vector_coeff_m1(ilevel)
#endif
  !==============================================
  ! Update preconditionner M=1/diag(A) boundaries
  !==============================================
  do irad=1,nvar_bicg
     if(block_diagonal_precond_bicg) then
        do i=1,nvar_bicg
           call make_virtual_fine_dp(precond_bicg(:,irad,i),ilevel)
        enddo
     else
        call make_virtual_fine_dp(var_bicg(:,irad,4),ilevel)
     endif
  enddo

!!$  write(*,*) 'debug matrix - vect'
!!$  do i=1,nb_ind
!!$     this = liste_ind(i)
!!$     do irad=1,nvar_bicg
!!$        write(*,'(3(3(ES11.3),2x),5x,ES11.3)') (coeff_glob_left(this,irad,jrad,1),jrad=1,nvar_bicg),(mat_residual_glob(this,irad,jrad),jrad=1,nvar_bicg),(coeff_glob_right(this,irad,jrad,1),jrad=1,nvar_bicg),residual_glob(this,irad)
!!$     enddo
!!$     write(*,*)
!!$  enddo
!!$  read(*,*)



  !==================================================================
  ! Compute r1 = b1 - A1x1 and store it into var_bicg(1:ncell,irad,1)
  !==================================================================
  call cmp_matrix_vector_product(ilevel,1)
  if(bicg_to_cg)then
     do irad=1,nvar_bicg
        do i=1,nb_ind
           this = liste_ind(i)
           var_bicg(this,irad,2) = var_bicg(this,irad,1)
        enddo
     enddo
  endif
  do irad=1,nvar_bicg
     call make_virtual_fine_dp(var_bicg(:,irad,1),ilevel)
     call make_virtual_fine_dp(var_bicg(:,irad,2),ilevel)
  enddo

  if(.not.bicg_to_cg) then
     !=========================================================================
     ! BiCGSTAB: Compute rbar_0 = r1 and store it into var_bicg(1:ncell,irad,9)
     !=========================================================================
     do irad=1,nvar_bicg
        do i=1,nb_ind
           this = liste_ind(i)
           var_bicg(this,irad,9) = var_bicg(this,irad,1)
        enddo
     enddo
  endif

    rho_bicg_old = one
    rho_bicg_new = one
  alpha_bicg     = one
  omega_bicg     = one
  if(bicg_to_cg) omega_bicg=zero

  !================================================================
  ! All     : Set v0 = 0 and store it into var_bicg(1:ncell,irad,3)
  ! BiCGSTAB: Set p0 = 0 and store it into var_bicg(1:ncell,irad,2)
  !================================================================
  do irad=1,nvar_bicg
     do i=1,nb_ind
        this = liste_ind(i)
        if(.not.bicg_to_cg) var_bicg(this,irad,2) = zero
        var_bicg(this,irad,3) = zero
     enddo
  enddo

  !=============================
  ! Compute right-hand side norm
  !=============================
  call dot_product_tot(var_bicg(:,:,1),var_bicg(:,:,1),rhs_norm1,final_sum)
 
  !============================================================
  ! Compute z_0 = K^{-1} r and store it into var_bicg(i,irad,6)
  !============================================================
  if(bicg_to_cg)then!neilneil
     do irad=1,nvar_bicg
        do i=1,nb_ind
           this = liste_ind(i)
           if(block_diagonal_precond_bicg) then
              var_bicg(this,irad,6)=zero
              do jrad=1,nvar_bicg
                 var_bicg(this,irad,6) = var_bicg(this,irad,6) + precond_bicg(this,irad,jrad) * var_bicg(this,jrad,1)
              enddo
           else
              var_bicg(this,irad,6) = var_bicg(this,irad,4) * var_bicg(this,irad,1)
           endif
        end do
     enddo
  endif !neilneil

  !====================
  ! MAIN ITERATION LOOP
  !====================   

  iter=0; itermax=5000

  error_ini=sqrt(rhs_norm1)
  error=error_ini

  max_loc=2.*epsilon

!  do while(error_ini.ne.zero .and. error/error_ini>epsilon .and.iter<itermax .and. error_ini .gt. 1.0e-12_dp)
!  do while(error/error_ini>epsilon .and.iter<itermax .and. error_cg_loc .gt. epsilon)! .and. error_ini/norm_er .gt. 1.0d-15)
  do while((max_loc>epsilon .or. error/error_ini>epsilon).and.iter<itermax)

     iter=iter+1

     !=========================================
     ! BiCGSTAB: Compute rho_bicg_new = rbar0.r
     ! BiCG2CG : Compute rho_bicg_new = r.z
     !=========================================
     call dot_product_tot(var_bicg(:,:,i_rho),var_bicg(:,:,1),r2,final_sum)
     rho_bicg_new = r2 ! real(final_sum)

     !================================================================================
     ! BiCGSTAB: Compute beta_bicg = rho_bicg_new/rho_bicg_old * alpha_bicg/omega_bicg
     ! BiCG2CG : Compute beta_bicg = rho_bicg_new/rho_bicg_old = (r.z)/(rold.zold)
     !================================================================================
     if(bicg_to_cg) then
        if(iter==1) then
           beta_bicg = zero
        else
           beta_bicg = rho_bicg_new/rho_bicg_old
        endif
     else
        beta_bicg = rho_bicg_new/rho_bicg_old * alpha_bicg/omega_bicg
     endif

     !=====================================================================
     ! BiCGSTAB: Recurrence on p = r + beta_bicg*p - omega_bicg*beta_bicg*v 
     ! BiCG2CG : Recurrence on p = z + beta_bicg*p 
     !=====================================================================
     call cX_plus_Y_to_Z_tot (beta_bicg,var_bicg(:,:,2),var_bicg(:,:,i_beta),var_bicg(:,:,2))
     if(.not.bicg_to_cg)then
        call cX_plus_Y_to_Z_tot (-omega_bicg*beta_bicg,var_bicg(:,:,3),var_bicg(:,:,2),var_bicg(:,:,2))
     endif

     call make_boundary_diffusion_tot(ilevel)
     do irad=1,nvar_bicg
        call make_virtual_fine_dp(var_bicg(:,irad,2),ilevel)
     enddo

     if(.not.bicg_to_cg)then
        !====================================================================
        ! BiCGSTAB: Compute y = K^{-1} p and store it into var_bicg(i,irad,5)
        !====================================================================
        do irad=1,nvar_bicg
           do i=1,nb_ind
              this = liste_ind(i)
              
              if(block_diagonal_precond_bicg) then
                 var_bicg(this,irad,5)=zero
                 do jrad=1,nvar_bicg
                    var_bicg(this,irad,5) = var_bicg(this,irad,5) + precond_bicg(this,irad,jrad) * var_bicg(this,jrad,2)
                 enddo
              else
                 var_bicg(this,irad,5) = var_bicg(this,irad,4) * var_bicg(this,irad,2)
              endif
           end do
        enddo
        ! Update boundaries
        call make_boundary_diffusion_tot(ilevel)
        do irad=1,nvar_bicg
           call make_virtual_fine_dp(var_bicg(:,irad,5),ilevel)
        enddo
     endif

     !===============================================================
     ! BiCGSTAB: Compute v = A y and store it into var_bicg(i,irad,3)
     ! BiCG2CG : Compute v = A p and store it into var_bicg(i,irad,3)
     !===============================================================
     call cmp_matrix_vector_product(ilevel,2)

     do irad=1,nvar_bicg
        call make_virtual_fine_dp(var_bicg(:,irad,3),ilevel)
     enddo

     !==========================
     ! BiCGSTAB: Compute rbar0.v
     ! BiCG2CG : Compute p.Ap
     !==========================
     call dot_product_tot(var_bicg(:,:,i_pAp),var_bicg(:,:,3),r2,final_sum)

     !===========================================================================
     ! BiCGSTAB: Compute scalar alpha_bicg = rho_bicg_new (=rbar0.r) / (rbar_0,v)
     ! BiCG2CG : Compute scalar alpha_bicg = rho_bicg_new (=r.z) / p.Ap
     !===========================================================================
     if(r2.eq.zero) then
        alpha_bicg = zero
     else
        alpha_bicg = rho_bicg_new / r2  ! real(final_sum)
     endif

     !===============================================================================
     ! BiCGSTAB: Recurrence on s = r - alpha_bicg*v   and store it in var_bicg(:,:,7)
     ! BiCG2CG : Recurrence on r = r - alpha_bicg*A.p and store it in var_bicg(:,:,1)
     !===============================================================================
     call cX_plus_Y_to_Z_tot (-alpha_bicg,var_bicg(:,:,3),var_bicg(:,:,1),var_bicg(:,:,i_s))

     !====================================================================
     ! BiCGSTAB: Compute z = K^{-1} s and store it into var_bicg(i,irad,6)
     ! BiCG2CG : Compute z = K^{-1} r and store it into var_bicg(i,irad,6)
     !====================================================================
     do irad=1,nvar_bicg
        do i=1,nb_ind
           this = liste_ind(i)

           if(block_diagonal_precond_bicg) then
              var_bicg(this,irad,6)=zero
              do jrad=1,nvar_bicg
                 var_bicg(this,irad,6) = var_bicg(this,irad,6) + precond_bicg(this,irad,jrad) * var_bicg(this,jrad,i_s)
              enddo
           else
              var_bicg(this,irad,6) = var_bicg(this,irad,4) * var_bicg(this,irad,i_s)
           endif
        end do
     enddo

     if(.not.bicg_to_cg)then

        ! Update boundaries
        call make_boundary_diffusion_tot(ilevel)
        do irad=1,nvar_bicg
           call make_virtual_fine_dp(var_bicg(:,irad,6),ilevel)
        enddo

        !===============================================================
        ! BiCGSTAB: Compute t = A z and store it into var_bicg(i,irad,8)
        !===============================================================
        call cmp_matrix_vector_product(ilevel,6)

        !=====================================================================================
        ! BiCGSTAB: Compute K^{-1} t to compute omega_bicg and store it in var_bicg(i,irad,10)
        !=====================================================================================
        do irad=1,nvar_bicg
           do i=1,nb_ind
              this = liste_ind(i)

              if(block_diagonal_precond_bicg) then
                 var_bicg(this,irad,10)=zero
                 do jrad=1,nvar_bicg
                    var_bicg(this,irad,10) = var_bicg(this,irad,10) + precond_bicg(this,irad,jrad) * var_bicg(this,jrad,8)
                 enddo
              else
                 var_bicg(this,irad,10) = var_bicg(this,irad,4) * var_bicg(this,irad,8)
              endif
           end do
        enddo

        !=============================================================================
        ! BiCGSTAB: Compute omega_bicg = (K^{-1} t , K^{-1} s) / (K^{-1} t , K^{-1} t)
        !=============================================================================
        call dot_product_tot(var_bicg(:,:,10),var_bicg(:,:, 6),r2,final_sum)
        call dot_product_tot(var_bicg(:,:,10),var_bicg(:,:,10),r3,final_sum)

        if(r3.eq.zero) then
           omega_bicg = zero
        else
           omega_bicg = r2 / r3
        endif

     else

        omega_bicg = zero

     endif
     
     !===========================
     ! Compute maximum variations
     !===========================
     max_loc=zero
     do irad=1,nvar_bicg
        do i=1,nb_ind
           this = liste_ind(i)
           if(uold(this,ind_bicg(irad)).ne.zero)then
              max_loc=max(max_loc,abs((alpha_bicg*var_bicg(this,irad,i_y)+&
                      omega_bicg*var_bicg(this,irad,6))/uold(this,ind_bicg(irad))))
           endif
        enddo
     enddo
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(max_loc,max_loc_all,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
     max_loc = max_loc_all
#endif

     !=======================================================
     ! BiCGSTAB: Recurrence on x = x + alpha*y + omega_bicg*z
     ! BiCG2CG : Recurrence on x = x + alpha*p
     !=======================================================
     do irad=1,nvar_bicg
        call cX_plus_Y_to_Z (alpha_bicg,var_bicg(:,irad,i_y),unew(:,ind_bicg(irad)),unew(:,ind_bicg(irad)))
        if(.not.bicg_to_cg) call cX_plus_Y_to_Z (omega_bicg,var_bicg(:,irad,6),unew(:,ind_bicg(irad)),unew(:,ind_bicg(irad)))
     enddo

     !=============================================
     ! BiCGSTAB: Recurrence on r = s - omega_bicg*t
     !=============================================
     if(.not.bicg_to_cg)then
        call cX_plus_Y_to_Z_tot (-omega_bicg,var_bicg(:,:,8),var_bicg(:,:,i_s),var_bicg(:,:,1))
     endif

     !===================================
     ! Update rho_bicg_old = rho_bicg_new
     !===================================
     rho_bicg_old = rho_bicg_new

     !===================================
     ! Compute right-hand side norm
     !===================================
     call dot_product_tot(var_bicg(:,:,1),var_bicg(:,:,1),rhs_norm1,final_sum)

     error=SQRT(rhs_norm1)

     if(verbose) then
        if (error_ini.ne.zero) then
           write(*,112)iter,error,error/error_ini,max_loc,error_ini
        else
           write(*,112)iter,error
        endif
     endif

  end do
  ! End main iteration loop

  if(iter >= itermax)then
     if(myid==1)write(*,*)'Radiative transfer failed to converge after ',iter,' iterations'
     call clean_stop
  end if

  !====================================
  ! Copie des flux
  !====================================
  call cmp_matrix_vector_product(ilevel,4)

  niter=niter+iter

#if USE_FLD==1
  !====================================
  ! Update gas temperature
  !====================================
  do i=1,nb_ind

     rho = uold(liste_ind(i),1)
     Told= uold(liste_ind(i),nvar) * Tr_floor
     Cv = unew(liste_ind(i),nvar+1)

     ambi_heating=zero
     ohm_heating=zero
     nimhd_heating=zero
     
!#if NIMHD==1
!     if(radiative_nimhdheating_in_cg)then
!        bx=0.5d0*(uold(liste_ind(i),6)+uold(liste_ind(i),nvar+1))
!        by=0.5d0*(uold(liste_ind(i),7)+uold(liste_ind(i),nvar+2))
!        bz=0.5d0*(uold(liste_ind(i),8)+uold(liste_ind(i),nvar+3))
!        bcell2=(bx**2+by**2+bz**2)
!        jx=uold(liste_ind(i),nvar-3)
!        jy=uold(liste_ind(i),nvar-2)
!        jz=uold(liste_ind(i),nvar-1)
!        jsquare=(jx**2+jy**2+jz**2)
!        ionisrate=default_ionisrate
!        
!        if(nmagdiffu .eq. 1 .or. nmagdiffu2 .eq. 1 )ohm_heating=jsquare*etaohmdiss(rho,bcell2,Told,ionisrate)*dt_imp
!        
!        if(nambipolar .eq. 1 .or. nambipolar2 .eq.1 )then
!           ambi_heating = (jy*bz-jz*by)**2+(jz*bx-jx*bz)**2+(jx*by-jy*bx)**2
!           ambi_heating = ambi_heating * betaad(rho,bcell2,Told,ionisrate)*dt_imp
!        endif
!        nimhd_heating=ambi_heating+ohm_heating
!     end if
!#endif  
     
     rhs=zero
     lhs=zero
     do igrp=1,ngrp
        Trold = cal_Teg(uold(liste_ind(i),firstindex_er+igrp)*scale_E0,igrp)

        wdtB = C_cal*dt_imp*planck_ana(rho*scale_d,Told,Told ,igrp,in_sink(liste_ind(i)))/scale_kappa
        if(rt_feedback)then
           wdtE = C_cal*dt_imp*planck_ana(rho*scale_d,Told,Told,igrp,in_sink(liste_ind(i)))/scale_kappa
        else
           wdtE = C_cal*dt_imp*planck_ana(rho*scale_d,Told,Trold,igrp,in_sink(liste_ind(i)))/scale_kappa
        end if
!!$        rhs=rhs-P_cal*wdt*(radiation_source(Told,igrp)/scale_E0-Told*deriv_radiation_source(Told,igrp)/scale_E0 &
!!$             & -unew(liste_ind(i),firstindex_er+igrp))
        rhs=rhs-P_cal*wdtB*(radiation_source(Told,igrp)/scale_E0-Told*deriv_radiation_source(Told,igrp)/scale_E0) &
             & + P_cal*wdtE*unew(liste_ind(i),firstindex_er+igrp)

        lhs=lhs+P_cal*wdtB*deriv_radiation_source(Told,igrp)/scale_E0
     enddo

     protostellar_heating = 0.0d0 !hybrid RT
#ifdef RT
     if(rt_protostar_m1 .and. rt_advect .and. nsink >0 .and. Teff_sink(1) > 0.0d0)then !heating term by M1 photons absorbed
        do im1=1,ngroups !!! Loop over M1 photon groups
           protostellar_heating = rtuold(liste_ind(i),iGroups(im1))*scale_Np*group_egy(im1)*ev2erg*planck_ana(rho*scale_d,Told,Teff_sink(1),1,in_sink(liste_ind(i)))*rt_c_cgs(ilevel)*z_ave / (scale_d*scale_v**2)
           end do
     end if
#endif
     protostellar_heating = protostellar_heating * dt_imp * scale_t

     unew(liste_ind(i),nvar) = (cv*Told+rhs+protostellar_heating+nimhd_heating)/(cv+lhs) / Tr_floor

  end do
#endif

  if(debug_energy)then
  write(*,*) 'After iterations - uold(5-9)-unew(5-9)'
  do i=1,nb_ind
     this = liste_ind(i)
     write(*,'(12(ES15.6))') uold(this,5),uold(this,nvar),uold(this,9),uold(this,10),unew(this,5),unew(this,nvar),unew(this,9),unew(this,10)
  enddo
  read(*,*)
  endif

  if(myid==1 .and. (mod(nstep,ncontrol)==0)) then
     if(bicg_to_cg) then 
        if(error_ini.ne.zero) then
           write(*,117)ilevel,iter,error/error_ini,max_loc
        else
           write(*,*)' CG :',iter, 'error_ini=',error_ini
        endif
     else
        if(error_ini.ne.zero) then
           write(*,118)ilevel,iter,error/error_ini,max_loc
        else
           write(*,*)' BiCGSTAB :',iter, 'error_ini=',error_ini
        endif
     endif
     write(*,*)'niter tot=',niter
     if(error_ini.ne.zero) then
        write(*,115)ilevel,iter,error,error/error_ini
     else
        write(*,115)ilevel,iter,error
     endif
  endif

  call make_boundary_diffusion_tot(ilevel)

  !====================
  ! Update energy value
  !====================
  if(static) then
     do i=1,nb_ind
        do irad = 1,nvar_bicg
           uold(liste_ind(i),ind_bicg(irad)) = unew(liste_ind(i),ind_bicg(irad))*norm_bicg(irad)
        enddo
     enddo
  else
     call cmp_energy(2)
  end if

  if(debug_energy)then
  write(*,*) 'after cmp_energy 2'
  do i=1,nb_ind
     this = liste_ind(i)
     write(*,'(12(ES15.6))') uold(this,5),uold(this,nvar),uold(this,9),uold(this,10),unew(this,5),unew(this,nvar),unew(this,9),unew(this,10)
  enddo
  read(*,*)
  endif

  ! Update boundaries
  do irad=1,nvar_trad
     call make_virtual_fine_dp(uold(1,ind_trad(irad)),ilevel)
  enddo
  do irad=1,nvar_bicg
     if(ilevel .gt. levelmin)call make_virtual_reverse_dp(rad_flux(1,irad),ilevel-1)
     call make_virtual_reverse_dp(rad_flux(1,irad),ilevel)
  enddo
  call make_virtual_fine_dp(uold(1,5),ilevel)

111 format('   Entering diffusion_cg')
112 format('   ==> Step=',i5,' Error=',2(1pe10.3,1x),e23.15,es18.5)
115 format('   ==> Level=',i5,' Step=',i5,' Error=',2(1pe10.3,1x))
117 format('   ==> Level=',i5,' Iteration CG=',i5,' Error L2=',(1pe10.3,1x),' Error Linf=',(1pe10.3,1x))
118 format('   ==> Level=',i5,' Iteration BiCGSTAB=',i5,' Error L2=',(1pe10.3,1x),' Error Linf=',(1pe10.3,1x))

  deallocate(liste_ind)

contains

  !###########################################################
  !###########################################################

  subroutine cX_plus_Y_to_Z (cste,vectX,vectY,vectZ) ! vectZ = cste*vectX+vectY
    implicit none
    real(dp),dimension(1:ncoarse+twotondim*ngridmax),intent(IN)::vectX,vectY
    real(dp),intent(IN)::cste
    real(dp),dimension(1:ncoarse+twotondim*ngridmax),intent(OUT)::vectZ


    do i=1,nb_ind
       vectZ(liste_ind(i)) = vectY(liste_ind(i)) + cste*vectX(liste_ind(i)) 
    end do

  end subroutine cX_plus_Y_to_Z

  !###########################################################
  !###########################################################

  subroutine dot_product_tot(fact1,fact2,dot_pdt,local_sum) ! dot_pdt = sum(fact1*fact2)
    implicit none
    real(dp),dimension(1:ncoarse+twotondim*ngridmax,1:nvar_bicg),intent(IN)::fact1,fact2
    real(dp),intent(OUT)::dot_pdt
    complex*16,intent(OUT)::local_sum

#ifndef WITHOUTMPI
    real(dp)::dot_pdt_all
#endif
    complex*16 ::global_sum
    integer::this

    dot_pdt=zero
    local_sum = cmplx(zero,zero)
    global_sum = cmplx(zero,zero)

    do irad=1,nvar_bicg
       do i=1,nb_ind
          this = liste_ind(i)
          !call DDPDD (cmplx(fact1(this,irad)*fact2(this,irad), zero,dp), local_sum, 1, itype)
          dot_pdt = dot_pdt + fact1(this,irad)*fact2(this,irad)
       end do
    enddo

    ! Compute global norms
#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(dot_pdt,dot_pdt_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    dot_pdt   = dot_pdt_all
!!$  	call MPI_ALLREDUCE(local_sum,global_sum,1,MPI_COMPLEX,MPI_SUMDD,MPI_COMM_WORLD,info)
!!$	local_sum = global_sum
#endif

  end subroutine dot_product_tot

  !###########################################################
  !###########################################################

  subroutine cX_plus_Y_to_Z_tot (cste,vectX,vectY,vectZ) ! vectZ = cste*vectX+vectY
    implicit none
    real(dp),dimension(1:ncoarse+twotondim*ngridmax,1:nvar_bicg),intent(IN)::vectX,vectY
    real(dp),intent(IN)::cste
    real(dp),dimension(1:ncoarse+twotondim*ngridmax,1:nvar_bicg),intent(OUT)::vectZ

    do irad=1,nvar_bicg
       do i=1,nb_ind
          vectZ(liste_ind(i),irad) = vectY(liste_ind(i),irad) + cste*vectX(liste_ind(i),irad) 
       end do
    enddo

  end subroutine cX_plus_Y_to_Z_tot


end subroutine diffusion_cg

!###########################################################
!###########################################################
!###########################################################
!###########################################################

subroutine cmp_matrix_and_vector_coeff_fld(ilevel)
  !------------------------------------------------
  ! This routine computes the matrix A and vector b
  !------------------------------------------------

  use amr_commons,only:active,ncoarse,nbor,son,myid
  use amr_parameters, only : ndim
  use hydro_commons
  use radiation_parameters
  use const
  use units_commons
  implicit none

  integer,intent(IN)::ilevel
  integer :: igroup,irad

  integer , dimension(1:nvector,1:2*ndim),save:: nbor_ilevel
  integer , dimension(1:nvector,1:ndim),save::   cell_left , cell_right , big_left, big_right
  integer ,dimension(1:nvector,0:2*ndim),save::  igridn
  integer ,dimension(1:nvector),save ::          ind_cell , ind_grid

!!$  real(dp),dimension(1:nvector  ,1:  nvar_bicg),save:: C_g,C_d

  integer :: i,idim,ind,igrid,ngrid,ncache,iskip,igrp,nx_loc
  integer :: supG,sub,supD

  real(dp)::dx,dx_loc,surf_loc,vol_loc,scale

#if NGRP>1
  ! variables used by the LAPACK inversion routines
  integer, parameter                        :: nwork = 256
  integer                                   :: info2
  integer                                   :: lda,lwork
  integer, dimension(      nvar_bicg)       :: ipiv
  integer, dimension(nwork*nvar_bicg)       :: work
  real(dp),dimension(1:nvar_bicg,1:nvar_bicg) ::inv
#endif

  real(dp),dimension(nvar_bicg,nvar_bicg)::coeff_left,coeff_right,mat_residual
  real(dp),dimension(nvar_bicg          )::residual
  
  ! Mesh size at level ilevel
  dx=half**ilevel

  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  surf_loc = dx_loc**(ndim-1)
  vol_loc  = dx_loc**ndim

  ! **************************** LOOP OVER CELLS ********************************** !

  ! Loop over myid grids by vector sweeps
  ncache = active(ilevel)%ngrid
  do igrid=1,ncache,nvector

     ! Gather nvector grids
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i) = active(ilevel)%igrid(igrid+i-1)
        igridn(i,0) = ind_grid(i)
     end do

     do idim=1,ndim
        do i=1,ngrid
           big_left (i,idim)  = nbor(ind_grid(i),2*idim-1)
           big_right(i,idim)  = nbor(ind_grid(i),2*idim  )
           igridn(i,2*idim-1) = son(big_left (i,idim))
           igridn(i,2*idim  ) = son(big_right(i,idim))
        end do
     end do

     ! Loop over cells
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do

        ! Determine the two2ndim and the direction of the grid of neighboors (-1,0,1)
        do idim = 1,ndim

           if (modulo((ind-1)/2**(idim-1),2)==0)then
              supG = (idim-1)*2+1               !direction of left nbor grid
              supD = 0                          !direction of right nbor grid
              sub = ind + 2**(idim-1)           ! position of nbor in its own grid
           else
              supG = 0                          !direction of left nbor grid
              supD = (idim-1)*2+2               !direction of right nbor grid
              sub = ind - 2**(idim-1)           !position of nbor in its own grid
           end if

           sub = ncoarse + (sub-1)*ngridmax     !nbor index offset from its own grid

           do i=1,ngrid

              ! Getting neighboors relative level (-1,0,1)

              if(son(ind_cell(i)) == 0 )then

                 if(igridn(i,supG)>0)then

                    cell_left(i,idim) = igridn(i,supG)+ sub
                    if(son(cell_left(i,idim))>0)then ! Left nbor more refined than me
                       nbor_ilevel(i,2*idim-1) = 1
                    else                             ! Left nbor as refined as me
                       nbor_ilevel(i,2*idim-1) = 0
                    end if

                 else                                ! Left nbor less refined than me

                    nbor_ilevel(i,2*idim-1) = -1
                    cell_left(i,idim)    = big_left(i,idim)
                 end if

                 if(igridn(i,supD)>0)then

                    cell_right(i,idim) = igridn(i,supD)+ sub
                    if(son(cell_right(i,idim))>0)then ! Right nbor more refined than me
                       nbor_ilevel(i,2*idim) = 1
                    else                              ! Right nbor as refined as me
                       nbor_ilevel(i,2*idim) = 0
                    end if

                 else                                 ! Right nbor less refined than me

                    nbor_ilevel(i,2*idim) = -1
                    cell_right(i,idim) = big_right(i,idim)
                 end if

              end if
           end do

        end do !ndim
        
        do i=1,ngrid
           if(son(ind_cell(i)) == 0 )then

              call compute_residual_in_cell(ind_cell(i),ilevel,vol_loc,residual,mat_residual)

              do igroup=1,ngrp
                 do igrp=1,ngrp
                    if(store_matrix) mat_residual_glob(ind_cell(i),igroup,igrp) = mat_residual(igroup,igrp)
                    if(block_diagonal_precond_bicg.or.igroup==igrp) then
                       precond_bicg(ind_cell(i),igroup,igrp) = mat_residual(igroup,igrp)
                    endif
                 enddo
                 if(store_matrix) residual_glob(ind_cell(i),igroup) = residual(igroup)
              enddo

           endif
        enddo

        ! Compute off-diagonal terms
        do idim = 1,ndim

           do i=1,ngrid
              if(son(ind_cell(i)) == 0 )then

                 call compute_coeff_left_right_in_cell(ind_cell(i),idim,cell_left(i,idim),cell_right(i,idim),nbor_ilevel(i,1:2*ndim),dx_loc,coeff_left,coeff_right)

                 do igroup=1,ngrp

                    if(store_matrix)then
                       coeff_glob_left (ind_cell(i),igroup,igroup,idim)=coeff_left(igroup,igroup)
                       coeff_glob_right(ind_cell(i),igroup,igroup,idim)=coeff_right(igroup,igroup)
                    endif

                    precond_bicg(ind_cell(i),igroup,igroup) = precond_bicg(ind_cell(i),igroup,igroup) + (coeff_left(igroup,igroup) + coeff_right(igroup,igroup))*alpha_imp

                 enddo

              end if
           end do

        enddo !ndim

        ! Compute preconditionning matrix                                                               
        do i=1,ngrid
           if(son(ind_cell(i)) == 0 )then
#if NGRP>1
              if(block_diagonal_precond_bicg) then
                 inv = precond_bicg(ind_cell(i),1:nvar_bicg,1:nvar_bicg)
                 lda = nvar_bicg ; lwork = nwork*nvar_bicg
                 
                 ! Invert the (nvar_bicg x nvar_bicg) matrix using LAPACK routines                      
                 !
                 ! DGETRF computes an LU factorization of a general M-by-N matrix A                     
                 ! using partial pivoting with row interchanges                                         
                 call dgetrf(nvar_bicg,nvar_bicg,inv,lda,ipiv,info2)
                 
                 ! DGETRI computes the inverse of a matrix using the LU factorization                   
                 ! computed by DGETRF                                                                   
                 call dgetri(nvar_bicg,inv,lda,ipiv,work,lwork,info2)
                 
                 precond_bicg(ind_cell(i),1:nvar_bicg,1:nvar_bicg)=inv
              else
#endif
                 do irad=1,nvar_bicg
                    var_bicg(ind_cell(i),irad,4) = one/precond_bicg(ind_cell(i),irad,irad)
                 enddo
#if NGRP>1
              endif
#endif
           end if
        end do !ngrid
        
     end do ! twotodim
  end do ! ncache

  return

end subroutine cmp_matrix_and_vector_coeff_fld

!###########################################################
!###########################################################
!###########################################################
!###########################################################

#if USE_M_1==1
subroutine cmp_matrix_and_vector_coeff_m1(ilevel)
  !------------------------------------------------------------------
  ! This routine computes the matrix A to vect_in and create vect_out
  !------------------------------------------------------------------
  use amr_commons,only:active,ncoarse,nbor,son,myid
  use amr_parameters, only : ndim
  use hydro_commons
  use radiation_parameters
  use const
  use units_commons
  use constants, only:pi
  implicit none

  integer,intent(IN)::ilevel
  integer :: igroup,irad,jrad

  integer , dimension(1:nvector,1:2*ndim),save:: nbor_ilevel
  integer , dimension(1:nvector,1:ndim),save::   cell_left , cell_right , big_left, big_right
  integer ,dimension(1:nvector,0:2*ndim),save::  igridn
  integer ,dimension(1:nvector),save ::          ind_cell , ind_grid

  real(dp),dimension(1:nvector  ,1:  nvar_bicg),save:: residu
  real(dp),dimension(1:3,1:3,1:nvar_trad+1+nvar_bicg) :: var_rad_subset
  real(dp),dimension(1:ngrp):: deriv

  ! M1 variables:
  !
  ! var_rad_subset(i,idim,irad)
  !   - dimension 1 : 1 -> 3   : left neighbour, i, right neighbour
  !   - dimension 2 : 1 -> ndim: direction of interest for flux computations
  !   - dimension 3 : 1 -> nvar_trad + 1+nvar_bicg :       1                          is Told 
  !                                                        2 :  ngrp+1                is Er
  !                                                    ngrp+2:2*ngrp+1                is Frx
  !                                                  2*ngrp+2:3*ngrp+1                is Fry
  !            IN 3D                                 3*ngrp+2:4*ngrp+1                is Frz
  !                                               nvar_trad+1                         is rho
  !                                               nvar_trad+2:nvar_trad+1+nvar_bicg   is rad_flux
  real(dp), dimension(       3,  3       ) :: Dedd,Dedd_dE
  real(dp), dimension(       3,  3,3     ) :: Dedd_dF
  real(dp), dimension(     1:3,1:3,2,3   ) :: DeddP,DeddM
  real(dp), dimension(             1:ngrp) :: flux_F_tot
  real(dp), dimension(      1:ndim,1:ngrp) :: flux_P_tot
  real(dp), dimension(           2,1:3   ) :: eps,lm,lp,lmp
  real(dp), dimension(             2,ndim) :: ErayM, ErayP, dEray
  real(dp), dimension(        ndim,2,ndim) :: FrayM, FrayP, dFray
  real(dp), dimension(             2,ndim) :: ffM, ffP
  real(dp), dimension(     1:3,1:3,2,ndim) :: PrayM,PrayP
  real(dp), dimension(             2,ndim) :: ap,am
  real(dp), dimension(         3,3       ) :: Dedd_temp,Dedd_dE_temp
  real(dp), dimension(         3,3,3     ) :: Dedd_dF_temp
  real(dp), dimension(             1:3   ) :: Fr_temp
  real(dp), dimension(3,3,       2,1:ndim) :: Dedd_ndim,Dedd_dE_ndim
  real(dp), dimension(3,3,3,     2,1:ndim) :: Dedd_dF_ndim
  real(dp), dimension(3,3,       2,1:ndim) :: DeddM_dE, DeddP_dE
  real(dp), dimension(3,3,3,     2,1:ndim) :: DeddM_dF, DeddP_dF
  real(dp), dimension(           2,1:ndim) :: flux_F
  real(dp), dimension(    1:ndim,2,1:ndim) :: flux_P
  real(dp)                                 :: conv
  real(dp), dimension(                  2) :: signe
  real(dp)                                 :: xx,v1,v2,v3,v4,r,kappa_rhoM,kappa_rhoP
  real(dp)                                 :: ff1,ff2,thetaM,thetaP,interpol_valp,fM,fP
  real(dp)                                 :: kp,kE,kF,ks,kFM,kFP,ksM,ksP,scattering_ana
  integer                                  :: iface,idim,jdim,kdim,index_e
  logical                                  :: cal_valp,chi_uniforme,derive_chi

  integer :: i,ind,igrid,ngrid,ncache,iskip,igrp,nx_loc
  integer :: supG,sub,supD

  real(dp),dimension(       1:ngrp) :: Erold
  real(dp),dimension(1:ndim,1:ngrp) :: Frold
  real(dp),dimension(1:nvar_bicg  ) :: radflux
  real(dp)::rho,Told,Told_norm,radiation_source,deriv_radiation_source,cv
  real(dp)::dx,dx_loc,surf_loc,vol_loc,scale
  real(dp):: rosseland_ana,planck_ana

  ! variables used by the LAPACK inversion routines
  integer, parameter                        :: nwork = 256
  integer                                   :: info2
  integer                                   :: lda,lwork
  integer, dimension(      nvar_bicg)       :: ipiv
  integer, dimension(nwork*nvar_bicg)       :: work
  real(dp),dimension(1:nvar_bicg,1:nvar_bicg) ::inv

  real(dp)::scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Mesh size at level ilevel
  dx=half**ilevel

  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  surf_loc = dx_loc**(ndim-1)
  vol_loc  = dx_loc**ndim

  ! **************************** LOOP OVER CELLS ********************************** !

  residu = zero

  ! Loop over myid grids by vector sweeps
  ncache = active(ilevel)%ngrid
  do igrid=1,ncache,nvector

     ! Gather nvector grids
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i) = active(ilevel)%igrid(igrid+i-1)
        igridn(i,0) = ind_grid(i)
     end do

     do idim=1,ndim
        do i=1,ngrid
           big_left (i,idim)  = nbor(ind_grid(i),2*idim-1)
           big_right(i,idim)  = nbor(ind_grid(i),2*idim  )
           igridn(i,2*idim-1) = son(big_left (i,idim))
           igridn(i,2*idim  ) = son(big_right(i,idim))
        end do
     end do

     ! Loop over cells
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do

        ! Determine the two2ndim and the direction of the grid of neighboors (-1,0,1)
        do idim = 1,ndim

           if (modulo((ind-1)/2**(idim-1),2)==0)then
              supG = (idim-1)*2+1   		!direction of left nbor grid
              supD = 0              		!direction of right nbor grid
              sub = ind + 2**(idim-1)           ! position of nbor in its own grid
           else
              supG = 0              		!direction of left nbor grid
              supD = (idim-1)*2+2   		!direction of right nbor grid
              sub = ind - 2**(idim-1)           !position of nbor in its own grid
           end if

           sub = ncoarse + (sub-1)*ngridmax     !nbor index offset from its own grid

           do i=1,ngrid

              ! Getting neighboors relative level (-1,0,1)

              if(son(ind_cell(i)) == 0 )then

                 if(igridn(i,supG)>0)then

                    cell_left(i,idim) = igridn(i,supG)+ sub
                    if(son(cell_left(i,idim))>0)then			! Left nbor more refined than me
                       nbor_ilevel(i,2*idim-1) = 1
                    else						! Left nbor as refined as me
                       nbor_ilevel(i,2*idim-1) = 0
                    end if

                 else							! Left nbor less refined than me

                    nbor_ilevel(i,2*idim-1) = -1
                    cell_left(i,idim)    = big_left(i,idim)
                 end if

              end if
           end do

           do i=1,ngrid
              if(son(ind_cell(i)) == 0 )then

                 if(igridn(i,supD)>0)then

                    cell_right(i,idim) = igridn(i,supD)+ sub
                    if(son(cell_right(i,idim))>0)then			! Right nbor more refined than me
                       nbor_ilevel(i,2*idim) = 1
                    else						! Right nbor as refined as me
                       nbor_ilevel(i,2*idim) = 0
                    end if

                 else							! Right nbor less refined than me

                    nbor_ilevel(i,2*idim) = -1
                    cell_right(i,idim) = big_right(i,idim)
                 end if

              end if
           end do

        end do !ndim
        
        do i=1,ngrid
           if(son(ind_cell(i)) == 0 )then

              do idim = 1,ndim

                 do irad=1,nvar_trad
                    var_rad_subset(1,idim,irad) = uold(cell_left (i,idim),ind_trad(irad))
                    var_rad_subset(2,idim,irad) = uold(ind_cell  (i     ),ind_trad(irad))
                    var_rad_subset(3,idim,irad) = uold(cell_right(i,idim),ind_trad(irad))
                 enddo
                 var_rad_subset(1,idim,nvar_trad+1) = uold(cell_left (i,idim),1)
                 var_rad_subset(2,idim,nvar_trad+1) = uold(ind_cell  (i     ),1)
                 var_rad_subset(3,idim,nvar_trad+1) = uold(cell_right(i,idim),1)

                 select case (nbor_ilevel(i,2*idim-1))
                 case (1,-1)
                    do irad=1,nvar_trad
                       var_rad_subset(1,idim,irad) = var_rad_subset(1,idim,irad)/norm_trad(irad)
                    enddo
                 end select

                 select case (nbor_ilevel(i,2*idim))
                 case (1,-1)
                    do irad=1,nvar_trad
                       var_rad_subset(3,idim,irad) = var_rad_subset(3,idim,irad)/norm_trad(irad)
                    enddo
                 end select

                 do irad = 1,nvar_bicg
                    !var_rad_subset(1,idim,nvar_trad+1+irad) = rad_flux(cell_left (i,idim),irad)
                    var_rad_subset(2,idim,nvar_trad+1+irad) = rad_flux(ind_cell  (i     ),irad)
                    !var_rad_subset(3,idim,nvar_trad+1+irad) = rad_flux(cell_right(i,idim),irad)
                 enddo

              enddo


              !==============================================================================

              rho       = var_rad_subset(2,1,nvar_trad+1)
              Told_norm = var_rad_subset(2,1,1)
              Told      = Told_norm * Tr_floor  
              Cv        = unew(ind_cell(i),nvar+1)

              do igrp=1,ngrp
                 Erold(igrp)=var_rad_subset(2,1,1+igrp)
                 ! Store the deriv_radiation_source to save computational time
                 deriv(igrp)=deriv_radiation_source(Told,igrp)
              enddo
              do irad = 1,nvar_bicg
                 radflux(irad)=var_rad_subset(2,1,nvar_trad+1+irad)
              enddo
              
              do idim=1,ndim
                 do igrp=1,ngrp
                    Frold(idim,igrp)=var_rad_subset(2,1,igrp+ngrp*idim+1)
                 enddo
              enddo

              cal_valp=.true.
              chi_uniforme=.false.
              derive_chi=.false.

              coeff_glob_left  (ind_cell(i),:,:,:) = zero
              coeff_glob_right (ind_cell(i),:,:,:) = zero
              residual_glob    (ind_cell(i),:  ) = zero
              mat_residual_glob(ind_cell(i),:,:) = zero
              do irad = 1,nvar_bicg
                 mat_residual_glob(ind_cell(i),irad,irad) = one
              enddo

              conv = one / Cv *P_cal/Tr_floor

              residual_glob(ind_cell(i),1) = Told_norm

              do igroup=1,ngrp

                 ! Fill matrix and vector
                 !******************************************************************
                 index_e = 1 + igroup
                 !index_f=index_e+idim*ngrp !to be done inside a idim loop


                 ! pour flux rentrants ou sortants
                 signe(1) = -one ; signe(2) = one

     
                 PrayM = zero ; PrayP = zero

                 do idim = 1,ndim
                    do iface = 1,2

           
                       ErayM(iface,idim) = var_rad_subset(iface  ,idim,index_e) !+ half*dx_loc*pvar_rad_subset(1,idim,idim,igroup)
                       ErayP(iface,idim) = var_rad_subset(iface+1,idim,index_e) !- half*dx_loc*pvar_rad_subset(1,idim,idim,igroup)
                       if(iface == 1)then
                          dEray(iface,idim) = ErayP(iface,idim) - var_rad_subset(2,idim,index_e)
                       else
                          dEray(iface,idim) = ErayM(iface,idim) - var_rad_subset(2,idim,index_e)
                       endif

                       fM = zero
                       fP = zero
                       do jdim = 1,ndim
                          FrayM(jdim,iface,idim) = var_rad_subset(iface  ,idim,index_e+jdim*ngrp) !+ half*dx_loc*  !!!!pvar_rad_subset(1,idim,idim,igroup)
                          FrayP(jdim,iface,idim) = var_rad_subset(iface+1,idim,index_e+jdim*ngrp) !- half*dx_loc*  !!!!pvar_rad_subset(1,idim,idim,igroup)
                          fM = fM + (FrayM(jdim,iface,idim)/ErayM(iface,idim))**2
                          fP = fP + (FrayP(jdim,iface,idim)/ErayP(iface,idim))**2
                       enddo
                       fM = sqrt(fM)
                       if(fM > one)then
                          r = one / fM
                          do jdim = 1,ndim
                             FrayM(jdim,iface,idim) = FrayM(jdim,iface,idim) * r
                          enddo
                       endif
                       fP = sqrt(fP)
                       if(fP > one)then
                          r = one / fP
                          do jdim = 1,ndim
                             FrayP(jdim,iface,idim) = FrayP(jdim,iface,idim) * r
                          enddo
                       endif

                       ffM(iface,idim) = zero ; ffP(iface,idim) = zero
                       do jdim = 1,ndim
                          if(iface.eq.1) then
                             dFray(jdim,iface,idim) = FrayP(jdim,iface,idim)- var_rad_subset(2,idim,index_e+jdim*ngrp)
                          else
                             dFray(jdim,iface,idim) = FrayM(jdim,iface,idim)- var_rad_subset(2,idim,index_e+jdim*ngrp)
                          endif
                          ffM(iface,idim) = ffM(iface,idim) + FrayM(jdim,iface,idim)**2
                          ffP(iface,idim) = ffP(iface,idim) + FrayP(jdim,iface,idim)**2
                       enddo
                       ffM(iface,idim) = sqrt(ffM(iface,idim))/ErayM(iface,idim)
                       ffP(iface,idim) = sqrt(ffP(iface,idim))/ErayP(iface,idim)

                       ! Calculate theta incidence angle to interpolate eigenvalues
                       ! theta needs to be calculated here before ffM is changed if it is greater than 1
                       if(ffM(iface,idim).gt.1.0e-05_dp) then
                          xx = FrayM(idim,iface,idim)/ErayM(iface,idim)/ffM(iface,idim)
                          if(abs(xx).gt.one) then
                             thetaM = half*pi*(one-sign(one,xx))
                          else
                             thetaM = acos(xx)
                          endif
                       else
                          thetaM = zero
                       endif
                       if(ffP(iface,idim).gt.1.0e-05_dp) then
                          xx = FrayP(idim,iface,idim)/ErayP(iface,idim)/ffP(iface,idim)
                          if(abs(xx).gt.one) then
                             thetaP = half*pi*(one-sign(one,xx))
                          else
                             thetaP = acos(xx)
                          endif
                       else
                          thetaP = zero
                       endif

                       ! If ffM > 1 then set it to 1 - WARNING: NOT VERY SATISFACTORY!!
                       if (ffM(iface,idim).gt. one) ffM(iface,idim)= one
                       if (ffP(iface,idim).gt. one) ffP(iface,idim)= one

                       ! Compute Eddington tensor and its derivatives at the cell centre and interfaces
                       if (Chi_uniforme) then 
                          Fr_temp=zero
                          do jdim = 1,ndim
                             Fr_temp(jdim) = var_rad_subset(iface+1,idim,index_e+jdim*ngrp)
                          enddo
                          call cal_Dedd(var_rad_subset(iface+1,idim,index_e),Fr_temp,Dedd_temp,Dedd_dE_temp,Dedd_dF_temp)
                          DeddP   (:,:  ,iface,idim) = Dedd_temp
                          DeddP_dE(:,:  ,iface,idim) = Dedd_dE_temp
                          DeddP_dF(:,:,:,iface,idim) = Dedd_dF_temp

                          Fr_temp=zero
                          do jdim = 1,ndim
                             Fr_temp(jdim) = var_rad_subset(iface,idim,index_e+jdim*ngrp)
                          enddo
                          call cal_Dedd(var_rad_subset(iface,idim,index_e),Fr_temp,Dedd_temp,Dedd_dE_temp,Dedd_dF_temp)
                          DeddM   (:,:  ,iface,idim) = Dedd_temp
                          DeddM_dE(:,:  ,iface,idim) = Dedd_dE_temp
                          DeddM_dF(:,:,:,iface,idim) = Dedd_dF_temp

                       else

                          Fr_temp=zero
                          Fr_temp(1:ndim) = FrayP(1:ndim,iface,idim)
                          call cal_Dedd(ErayP(iface,idim),Fr_temp,Dedd_temp,Dedd_dE_temp,Dedd_dF_temp)
                          DeddP   (:,:  ,iface,idim) = Dedd_temp
                          DeddP_dE(:,:  ,iface,idim) = Dedd_dE_temp
                          DeddP_dF(:,:,:,iface,idim) = Dedd_dF_temp

                          Fr_temp=zero
                          Fr_temp(1:ndim) = FrayM(1:ndim,iface,idim)
                          call cal_Dedd(ErayM(iface,idim),Fr_temp,Dedd_temp,Dedd_dE_temp,Dedd_dF_temp)
                          DeddM   (:,:  ,iface,idim) = Dedd_temp
                          DeddM_dE(:,:  ,iface,idim) = Dedd_dE_temp
                          DeddM_dF(:,:,:,iface,idim) = Dedd_dF_temp
                          
                       endif

                       do jdim = 1,ndim
                          do kdim = 1,ndim
                             PrayM(jdim,kdim,iface,idim) = DeddM(jdim,kdim,iface,idim)*ErayM(iface,idim)
                             PrayP(jdim,kdim,iface,idim) = DeddP(jdim,kdim,iface,idim)*ErayP(iface,idim)
                          enddo
                       enddo


                       ! Compute epsilon for asymptotic preserving scheme
                       if (verbose)  write(*,*) '     -- compute epsilon --'

                       ff1 = ffM(iface,idim)
                       ff2 = ffP(iface,idim)

                       print*,'WARNING: M1 needs UPDATE in function scattering_ana, planck_ana,rosseland_ana for Tr '
                       stop

                       kFM= rosseland_ana(var_rad_subset(iface  ,idim,nvar_trad+1)*scale_d,var_rad_subset(iface  ,idim,1)*Tr_floor,var_rad_subset(iface  ,idim,1)*Tr_floor,igroup)/scale_kappa
                       kFP= rosseland_ana(var_rad_subset(iface+1,idim,nvar_trad+1)*scale_d,var_rad_subset(iface+1,idim,1)*Tr_floor,var_rad_subset(iface  ,idim,1)*Tr_floor,igroup)/scale_kappa
                       ksM=scattering_ana(var_rad_subset(iface  ,idim,nvar_trad+1)*scale_d,var_rad_subset(iface  ,idim,1)*Tr_floor,var_rad_subset(iface  ,idim,1)*Tr_floor,igroup)/scale_kappa
                       ksP=scattering_ana(var_rad_subset(iface+1,idim,nvar_trad+1)*scale_d,var_rad_subset(iface+1,idim,1)*Tr_floor,var_rad_subset(iface  ,idim,1)*Tr_floor,igroup)/scale_kappa

                       kappa_rhoM = (kFM+ksM)
                       kappa_rhoP = (kFP+ksP)

                       if( (kappa_rhoM .ne.zero) .and. (kappa_rhoM .ne.zero) ) then
                          eps(iface,idim) = max( one/(kappa_rhoM*dx_loc) , one/(kappa_rhoP*dx_loc) ) 
                          eps(iface,idim) = max(eps(iface,idim),two*abs(ff1),two*abs(ff2)) 
                          eps(iface,idim) = min(one,eps(iface,idim))
                       else
                          eps(iface,idim) = one
                       endif

                       ! Compute eigenvalues at the cell interfaces
                       if (verbose)  write(*,*) '     -- compute eigenvalues --'

                       if (cal_valp) then

                          select case(irad_trans_model)

                          case(irad_trans_model_p1) ! 'P1'

                             am(iface,idim)=min(-valp_min,-one/sqrt(three))
                             ap(iface,idim)=max(+valp_min, one/sqrt(three))

                          case(irad_trans_model_m1) ! 'M1'

                             am(iface,idim)=min(-valp_min,interpol_valp(ffM(iface,idim),thetaM,eps(iface,idim),1), &
                                                          interpol_valp(ffP(iface,idim),thetaP,eps(iface,idim),1))
                             ap(iface,idim)=max(+valp_min,interpol_valp(ffM(iface,idim),thetaM,eps(iface,idim),4), &
                                                          interpol_valp(ffP(iface,idim),thetaP,eps(iface,idim),4))

                          end select

                       else
                          am(iface,idim) = -one
                          ap(iface,idim) =  one
                       endif


                       if(ap(iface,idim).ne.am(iface,idim)) then
                          lm (iface,idim) = am(iface,idim)               /(ap(iface,idim)-am(iface,idim))
                          lp (iface,idim) =                ap(iface,idim)/(ap(iface,idim)-am(iface,idim))
                          lmp(iface,idim) = am(iface,idim)*ap(iface,idim)/(ap(iface,idim)-am(iface,idim))
                       else ! then ap=am=0.
                          lm (iface,idim) = -half
                          lp (iface,idim) =  half
                          lmp(iface,idim) =  zero
                       endif
                       
                    enddo
                 enddo

                 if (verbose)  write(*,*) '     -- compute Dedd --'

                 Fr_temp=zero
                 do jdim = 1,ndim
                    Fr_temp(jdim) = var_rad_subset(2,1,index_e+jdim*ngrp)
                 enddo
                 call cal_Dedd(var_rad_subset(2,1,index_e),Fr_temp,Dedd,Dedd_dE,Dedd_dF)

                 flux_F = zero
                 flux_P = zero
                 do idim = 1,ndim

                    v1 = dFray(idim,1,idim)
                    v2 = dEray(     1,idim)
                    v3 = dFray(idim,2,idim)
                    v4 = dEray(     2,idim)

                    flux_F(1,idim) = (+       lp (1,idim)            *FrayM(idim,1,idim) &
                                      -       lm (1,idim)            *v1                 &
                                      -       lmp(1,idim)*eps(1,idim)*ErayM(     1,idim) &
                                      +       lmp(1,idim)*eps(1,idim)*v2 )*surf_loc*C_cal
                    flux_F(2,idim) = (-       lm (2,idim)            *FrayP(idim,2,idim) &
                                      +       lp (2,idim)            *v3                 &
                                      +       lmp(2,idim)*eps(2,idim)*ErayP(     2,idim) &
                                      -       lmp(2,idim)*eps(2,idim)*v4 )*surf_loc*C_cal

                 enddo

                 do jdim = 1,ndim
        
       
                    do idim=1,ndim

                       v1 = dFray(jdim     ,1,idim)
                       v2 = PrayP(jdim,idim,1,idim) - Dedd(jdim,idim)*var_rad_subset(2,idim,index_e)
                       v3 = dFray(jdim     ,2,idim)
                       v4 = PrayM(jdim,idim,2,idim) - Dedd(jdim,idim)*var_rad_subset(2,idim,index_e)

                       flux_P(jdim,1,idim) = ( ( lp (1,idim)*eps(1,idim)**2+(one-eps(1,idim)**2)*half)* PrayM(jdim,idim,1,idim) &
                                              -  lmp(1,idim)*eps(1,idim)                              * FrayM(jdim     ,1,idim) &
                                              +  lmp(1,idim)*eps(1,idim)                              * v1                      &
                                              +(-lm (1,idim)*eps(1,idim)**2+(one-eps(1,idim)**2)*half)* v2)                     &
                                              *surf_loc*C_cal

                       flux_P(jdim,2,idim) = ( (-lm (2,idim)*eps(2,idim)**2+(one-eps(2,idim)**2)*half)* PrayP(jdim,idim,2,idim) &
                                              +  lmp(2,idim)*eps(2,idim)                              * FrayP(jdim     ,2,idim) &
                                              -  lmp(2,idim)*eps(2,idim)                              * v3                      &
                                              +( lp (2,idim)*eps(2,idim)**2+(one-eps(2,idim)**2)*half)* v4)                     &
                                              *surf_loc*C_cal
                    enddo
                 enddo

                 if (verbose)  write(*,*) '     -- compute fluxtot --'

                 flux_F_tot(  igroup) = zero
                 flux_P_tot(:,igroup) = zero
                 do idim=1,ndim
                    do iface=1,2
                       flux_F_tot(igroup) = flux_F_tot(igroup) + signe(iface)*dt_imp*flux_F(iface,idim)
                       do jdim=1,ndim
                          flux_P_tot(jdim,igroup) = flux_P_tot(jdim,igroup) + signe(iface)*dt_imp*flux_P(jdim,iface,idim)
                       enddo
                    enddo
                 enddo

                 flux_F_tot(igroup) = flux_F_tot(igroup)/vol_loc
                 do idim=1,ndim
                    flux_P_tot(idim,igroup) = flux_P_tot(idim,igroup)/vol_loc
                 enddo

                 ! Store flux_F and flux_P in var_rad_subset
                 do idim=1,ndim

                    var_rad_subset(1,idim,nrad+3) = zero
                    var_rad_subset(3,idim,nrad+3) = zero

                    var_rad_subset(1,idim,nrad+3+igroup) = flux_F(1,idim) * dt_imp/vol_loc
                    var_rad_subset(3,idim,nrad+3+igroup) = flux_F(2,idim) * dt_imp/vol_loc

                    do jdim=1,ndim
                       var_rad_subset(1,idim,nrad+3+igroup+jdim*ngrp) = flux_P(jdim,1,idim) * dt_imp/vol_loc
                       var_rad_subset(3,idim,nrad+3+igroup+jdim*ngrp) = flux_P(jdim,2,idim) * dt_imp/vol_loc
                    enddo

                 enddo

                 ! Compute Eddington tensors
                 do idim=1,ndim

                    Fr_temp=zero
                    do jdim = 1,ndim
                       Fr_temp(jdim) = var_rad_subset(1,idim,index_e+jdim*ngrp)
                    enddo
                    call cal_Dedd(var_rad_subset(1,idim,index_e),Fr_temp,Dedd_temp,Dedd_dE_temp,Dedd_dF_temp)
                    Dedd_ndim   (:,:  ,1,idim)=Dedd_temp
                    Dedd_dE_ndim(:,:  ,1,idim)=Dedd_dE_temp
                    Dedd_dF_ndim(:,:,:,1,idim)=Dedd_dF_temp
                    
                    Fr_temp=zero
                    do jdim = 1,ndim
                       Fr_temp(jdim) = var_rad_subset(3,idim,index_e+jdim*ngrp)
                    enddo
                    call cal_Dedd(var_rad_subset(3,idim,index_e),Fr_temp,Dedd_temp,Dedd_dE_temp,Dedd_dF_temp)
                    Dedd_ndim   (:,:  ,2,idim)=Dedd_temp
                    Dedd_dE_ndim(:,:  ,2,idim)=Dedd_dE_temp
                    Dedd_dF_ndim(:,:,:,2,idim)=Dedd_dF_temp

                 enddo

                 ! Compute matrix index for radiative energy
                 ! The matrix is stored in the following way: (for 2 groups)
                 ! [T,E(1),E(2),Fx(1),Fx(2),Fy(1),Fy(2),Fz(1),Fz(2)]

                 if (verbose)  write(*,*) '     -- compute matrix --'

                 ! Fill matrix and vector
                 !******************************************************************

                 ! Er equation without the source term
                 do idim = 1,ndim
                    mat_residual_glob(ind_cell(i),index_e,index_e          ) = mat_residual_glob(ind_cell(i),index_e,index_e) - dt_imp/vol_loc*C_cal* &
                                (lmp(1,idim)*eps(1,idim)*surf_loc+lmp(2,idim)*eps(2,idim)*surf_loc)
                    mat_residual_glob(ind_cell(i),index_e,index_e+idim*ngrp) =                                                 dt_imp/vol_loc*C_cal * &
                                (lm (1,idim)            *surf_loc+lp (2,idim)            *surf_loc)
                 enddo

                 residual_glob(ind_cell(i),index_e) =  Erold(igroup) !- flux_F_tot(igroup)
    
                 ! Tgaz equation (using the coefficients from the Er equation without the source term)
                 mat_residual_glob(ind_cell(i),1,index_e) = mat_residual_glob(ind_cell(i),index_e,index_e) * conv
                 do idim = 1,ndim
                    mat_residual_glob(ind_cell(i),1,index_e+idim*ngrp) = mat_residual_glob(ind_cell(i),index_e,index_e+idim*ngrp)* conv
                 enddo

                 residual_glob(ind_cell(i),1) = residual_glob(ind_cell(i),1) + residual_glob(ind_cell(i),index_e)*conv

                 ! Finally include the energy source term: c (kp aT^4-kE E)

                 ! B(      Tn+1              )=B(Tn)+B'(Tn)*(      Tn+1-      Tn)=      B'(Tn)*Tn+1 + B(Tn)-      B'(Tn)*Tn
                 ! B(alpha*Tn+1+(1-alpha)*Tn))=B(Tn)+B'(Tn)*(alpha*Tn+1-alpha*Tn)=alpha*B'(Tn)*Tn+1 + B(Tn)-alpha*B'(Tn)*Tn

                 kp=    planck_ana(rho*scale_d,Told,Told,igroup)/scale_kappa
                 kE=    planck_ana(rho*scale_d,Told,Told,igroup)/scale_kappa
                 kF= rosseland_ana(rho*scale_d,Told,Told,igroup)/scale_kappa
                 ks=scattering_ana(rho*scale_d,Told,Told,igroup)/scale_kappa

                 mat_residual_glob(ind_cell(i),index_e,1      ) =                                                - dt_imp*kp*C_cal &
                                                                  *deriv(igroup)/scale_E0 *Tr_floor
                 mat_residual_glob(ind_cell(i),index_e,index_e) = mat_residual_glob(ind_cell(i),index_e,index_e) + dt_imp*kE*C_cal

                 residual_glob    (ind_cell(i),index_e        ) = residual_glob    (ind_cell(i),index_e        ) + dt_imp*kp*C_cal &
                                                                  *(radiation_source(Told,igroup)-deriv(igroup)*Told)/scale_E0
                                      !- dt_imp*kE*C_cal*Erold(igroup)
                                     
                 ! Fr equation using only dtFr + c2(Div . Pr) = 0 terms
                 do idim = 1,ndim

                    do jdim = 1,ndim
                       mat_residual_glob(ind_cell(i),index_e+idim*ngrp,index_e          ) = mat_residual_glob(ind_cell(i),index_e+idim*ngrp,index_e) &
                            + C_cal*dt_imp/vol_loc*( ( lp(2,jdim)*eps(2,jdim)**2+(one-eps(2,jdim)**2)*half)*surf_loc &
                                                   - (-lm(1,jdim)*eps(1,jdim)**2+(one-eps(1,jdim)**2)*half)*surf_loc ) * Dedd(idim,jdim)
                       mat_residual_glob(ind_cell(i),index_e+idim*ngrp,index_e+idim*ngrp) = mat_residual_glob(ind_cell(i),index_e+idim*ngrp,index_e+idim*ngrp) &
                            - C_cal*dt_imp/vol_loc*(  lmp(1,jdim)*eps(1,jdim)*surf_loc + lmp(2,jdim)*eps(2,jdim)*surf_loc )
                    enddo

                    mat_residual_glob(ind_cell(i),index_e+idim*ngrp,index_e+idim*ngrp) = mat_residual_glob(ind_cell(i),index_e+idim*ngrp,index_e+idim*ngrp) + C_cal *dt_imp*(kF+ks)

                    residual_glob(ind_cell(i),index_e+idim*ngrp) = Frold(idim,igroup) !- flux_P_tot(idim,igroup)
                 enddo

                 if (verbose)  write(*,*) '     -- compute deriv matrix --'

                 do idim = 1,ndim
                    do jdim = 1,ndim
                       ! The Dedd derivatives
                       mat_residual_glob(ind_cell(i),index_e+idim*ngrp,index_e     ) = mat_residual_glob(ind_cell(i),index_e+idim*ngrp,index_e     ) &
                            + C_cal*dt_imp/vol_loc*( ( lp(2,jdim)*eps(2,jdim)**2+(one-eps(2,jdim)**2)*half)*surf_loc &
                                                   - (-lm(1,jdim)*eps(1,jdim)**2+(one-eps(1,jdim)**2)*half)*surf_loc ) * Dedd_dE(idim,jdim     )  ! d/dEray
                       do kdim=1,ndim
                          mat_residual_glob(ind_cell(i),index_e+idim*ngrp,index_e+kdim*ngrp) = mat_residual_glob(ind_cell(i),index_e+idim*ngrp,index_e+kdim*ngrp) &
                               + C_cal*dt_imp/vol_loc*( ( lp(2,jdim)*eps(2,jdim)**2+(one-eps(2,jdim)**2)*half)*surf_loc &
                                                      - (-lm(1,jdim)*eps(1,jdim)**2+(one-eps(1,jdim)**2)*half)*surf_loc ) * Dedd_dF(idim,jdim,kdim)  ! d/dFray
                       enddo
                    enddo
                 enddo

                 if (verbose)  write(*,*) '     -- compute left and right coefficients --'

                 ! non-diagonal terms

                 do idim = 1,ndim

                    ! Left coefficient ---

                    ! Radiative energy equation
                    coeff_glob_left(ind_cell(i),index_e,index_e          ,idim) =  dt_imp/vol_loc*C_cal*lmp(1,idim)*eps(1,idim)*surf_loc
                    coeff_glob_left(ind_cell(i),index_e,index_e+idim*ngrp,idim) = -dt_imp/vol_loc*C_cal* lp(1,idim)            *surf_loc

                    ! Temperature equation
                    coeff_glob_left(ind_cell(i),1,index_e          ,idim) = coeff_glob_left(ind_cell(i),index_e,index_e          ,idim)*conv
                    coeff_glob_left(ind_cell(i),1,index_e+idim*ngrp,idim) = coeff_glob_left(ind_cell(i),index_e,index_e+idim*ngrp,idim)*conv
                    
                    ! Radiative flux equation
                    do jdim=1,ndim
                       coeff_glob_left(ind_cell(i),index_e+jdim*ngrp,index_e,idim) = -dt_imp/vol_loc*   C_cal*( lp(1,idim)*eps(1,idim)**2+(one-eps(1,idim)**2)*half ) * &
                            (DeddM(jdim,idim,1,idim)+Dedd_dE_ndim(jdim,idim,1,idim))*surf_loc
                       do kdim=1,ndim
                          coeff_glob_left(ind_cell(i),index_e+jdim*ngrp,index_e+kdim*ngrp,idim) = -dt_imp/vol_loc* C_cal*( lp(1,idim)*eps(1,idim)**2+(one-eps(1,idim)**2)*half ) * &
                               Dedd_dF_ndim(jdim,idim,kdim,1,idim)*surf_loc
                       enddo
                       coeff_glob_left(ind_cell(i),index_e+jdim*ngrp,index_e+jdim*ngrp,idim) = coeff_glob_left(ind_cell(i),index_e+jdim*ngrp,index_e+jdim*ngrp,idim) + dt_imp/vol_loc*C_cal*lmp(1,idim)*eps(1,idim)*surf_loc
                    enddo

                    ! Right coefficient ---

                    ! Radiative energy equation
                    coeff_glob_right(ind_cell(i),index_e,index_e          ,idim) =  dt_imp/vol_loc*C_cal*lmp(2,idim)*eps(2,idim)*surf_loc
                    coeff_glob_right(ind_cell(i),index_e,index_e+idim*ngrp,idim) = -dt_imp/vol_loc*C_cal* lm(2,idim)            *surf_loc

                    ! Temperature equation
                    coeff_glob_right(ind_cell(i),1,index_e          ,idim) = coeff_glob_right(ind_cell(i),index_e,index_e          ,idim)*conv
                    coeff_glob_right(ind_cell(i),1,index_e+idim*ngrp,idim) = coeff_glob_right(ind_cell(i),index_e,index_e+idim*ngrp,idim)*conv
                    
                    ! Radiative flux equation
                    do jdim=1,ndim
                       coeff_glob_right(ind_cell(i),index_e+jdim*ngrp,index_e,idim) = dt_imp/vol_loc*   C_cal*(-lm(2,idim)*eps(2,idim)**2+(one-eps(2,idim)**2)*half ) * &
                            (DeddP(jdim,idim,2,idim)+Dedd_dE_ndim(jdim,idim,2,idim))*surf_loc
                       do kdim=1,ndim
                          coeff_glob_right(ind_cell(i),index_e+jdim*ngrp,index_e+kdim*ngrp,idim) = dt_imp/vol_loc* C_cal*(-lm(2,idim)*eps(2,idim)**2+(one-eps(2,idim)**2)*half ) * &
                               Dedd_dF_ndim(jdim,idim,kdim,2,idim)*surf_loc
                       enddo
                       coeff_glob_right(ind_cell(i),index_e+jdim*ngrp,index_e+jdim*ngrp,idim) = coeff_glob_right(ind_cell(i),index_e+jdim*ngrp,index_e+jdim*ngrp,idim) + dt_imp/vol_loc*C_cal*lmp(2,idim)*eps(2,idim)*surf_loc
                    enddo

                 enddo ! end idim loop
                 
              enddo ! end multigroup loop
              
              if(verbose) write(*,*) 'End of multigroup loop'

              !==============================================================================

              do irad=1,nvar_bicg
                 do jrad=1,nvar_bicg
                    if(block_diagonal_precond_bicg.or.irad==jrad) then
                       precond_bicg(ind_cell(i),irad,jrad) = mat_residual_glob(ind_cell(i),irad,jrad)
                    endif
                    coeff_glob_left (ind_cell(i),irad,jrad,:) = -coeff_glob_left (ind_cell(i),irad,jrad,:)
                    coeff_glob_right(ind_cell(i),irad,jrad,:) = -coeff_glob_right(ind_cell(i),irad,jrad,:)
                 enddo
              enddo

           endif
        enddo

        ! Compute preconditionning matrix                                                               
        do i=1,ngrid
           if(son(ind_cell(i)) == 0 )then
!#if NRAD>1 !warning NRAD changed in NGRP in Makefile
              if(block_diagonal_precond_bicg) then
                 inv = precond_bicg(ind_cell(i),1:nvar_bicg,1:nvar_bicg)
                 lda = nvar_bicg ; lwork = nwork*nvar_bicg
                 
                 ! Invert the (nvar_bicg x nvar_bicg) matrix using LAPACK routines                      
                 !                                                                                      
                 ! DGETRF computes an LU factorization of a general M-by-N matrix A                     
                 ! using partial pivoting with row interchanges                                         
                 call dgetrf(nvar_bicg,nvar_bicg,inv,lda,ipiv,info2)
                 
                 ! DGETRI computes the inverse of a matrix using the LU factorization                   
                 ! computed by DGETRF                                                                   
                 call dgetri(nvar_bicg,inv,lda,ipiv,work,lwork,info2)
                 
                 precond_bicg(ind_cell(i),1:nvar_bicg,1:nvar_bicg)=inv
              else
!#endif
                 do irad=1,nvar_bicg
                    var_bicg(ind_cell(i),irad,4) = one/precond_bicg(ind_cell(i),irad,irad)
                 enddo
!#if NRAD>1 !warning NRAD changed in NGRP in Makefile
              endif
!#endif
           end if
        end do !ngrid
        
     end do ! twotodim
  end do ! ncache

  return

end subroutine cmp_matrix_and_vector_coeff_m1
#endif

!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmp_matrix_vector_product(ilevel,compute)
  !------------------------------------------------------------------
  ! This routine computes 
  ! compute = 1 : residual           	  return B - Ax
  ! compute = 2 : Product                 return  A.p
  ! compute = 3 : Preconditionner         return diag(A) or block_diag(A)
  ! compute = 4 : Compute flux in rad_flux
  !
  ! For BICG
  ! compute = 6 : product                 return  A.p
  !------------------------------------------------------------------

  use amr_commons,only:active,ncoarse,nbor,son,myid
  use amr_parameters, only : ndim
  use hydro_commons
  use radiation_parameters
  use hydro_parameters,only:nvar_bicg
  use const
  use units_commons
  implicit none

  integer,intent(IN)::compute,ilevel
  integer :: irad,jrad

  integer , dimension(1:nvector,1:2*ndim),save:: nbor_ilevel
  integer , dimension(1:nvector,1:ndim),save::   cell_left , cell_right , big_left, big_right
  integer ,dimension(1:nvector,0:2*ndim),save::  igridn
  integer ,dimension(1:nvector),save ::          ind_cell , ind_grid

  real(dp),dimension(1:nvector  ,1:  nvar_bicg),save:: residu
  real(dp),dimension(1:nvector  ,1:  nvar_bicg),save:: phi_g,phi_c,phi_d,val_g,val_d

  integer :: i,idim,ind,igrid,ngrid,ncache,iskip,nx_loc
  integer :: supG,sub,supD

  real(dp)::dx,dx_loc,surf_loc,vol_loc,scale

  integer :: ind_res

  ! Mesh size at level ilevel
  dx=half**ilevel

  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  surf_loc = dx_loc**(ndim-1)
  vol_loc  = dx_loc**ndim

  ! **************************** LOOP OVER CELLS ********************************** !

  residu = zero

  ! Loop over myid grids by vector sweeps
  ncache = active(ilevel)%ngrid
  do igrid=1,ncache,nvector
  
     ! Gather nvector grids
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i) = active(ilevel)%igrid(igrid+i-1)
        igridn(i,0) = ind_grid(i)
     end do

     do idim=1,ndim
        do i=1,ngrid
           big_left (i,idim)  = nbor(ind_grid(i),2*idim-1)
           big_right(i,idim)  = nbor(ind_grid(i),2*idim  )
           igridn(i,2*idim-1) = son(big_left (i,idim))
           igridn(i,2*idim  ) = son(big_right(i,idim))
        end do
     end do

     ! Loop over cells
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do

        ! Determine the two2ndim and the direction of the grid of neighboors (-1,0,1)
        do idim = 1,ndim

           if (modulo((ind-1)/2**(idim-1),2)==0)then
              supG = (idim-1)*2+1     !direction of left nbor grid
              supD = 0                !direction of right nbor grid
              sub = ind + 2**(idim-1) !position of nbor in its own grid
           else
              supG = 0                !direction of left nbor grid
              supD = (idim-1)*2+2     !direction of right nbor grid
              sub = ind - 2**(idim-1) !position of nbor in its own grid
           end if

           sub = ncoarse + (sub-1)*ngridmax !nbor index offset from its own grid

           do i=1,ngrid

              ! Getting neighboors relative level (-1,0,1)

              if(son(ind_cell(i)) == 0 )then

                 if(igridn(i,supG)>0)then

                    cell_left(i,idim) = igridn(i,supG)+ sub
                    if(son(cell_left(i,idim))>0)then ! Left nbor more refined than me
                       nbor_ilevel(i,2*idim-1) = 1
                    else                             ! Left nbor as refined as me
                       nbor_ilevel(i,2*idim-1) = 0
                    end if

                 else                                ! Left nbor less refined than me

                    nbor_ilevel(i,2*idim-1) = -1
                    cell_left(i,idim)    = big_left(i,idim)
                 end if

                 if(igridn(i,supD)>0)then

                    cell_right(i,idim) = igridn(i,supD)+ sub
                    if(son(cell_right(i,idim))>0)then ! Right nbor more refined than me
                       nbor_ilevel(i,2*idim) = 1
                    else                              ! Right nbor as refined as me
                       nbor_ilevel(i,2*idim) = 0
                    end if

                 else                                 ! Right nbor less refined than me

                    nbor_ilevel(i,2*idim) = -1
                    cell_right(i,idim) = big_right(i,idim)
                 end if

              end if
           end do

        end do !ndim
        
        do i=1,ngrid
           if(son(ind_cell(i)) == 0 )then

              if(.not.store_matrix)then
                 call compute_residual_in_cell(ind_cell(i),ilevel,vol_loc,residual_glob(1,:),mat_residual_glob(1,:,:))
                 ind_res = 1
              else
                 ind_res = ind_cell(i)
              endif

              select case (compute)

              case (1) ! residu = b - Ix
                 
                 do irad = 1,nvar_bicg
                    residu(i,irad) = residual_glob(ind_res,irad)
                    do jrad = 1,nvar_bicg
                       residu(i,irad) = residu(i,irad) - mat_residual_glob(ind_res,irad,jrad)*uold(ind_cell(i),ind_bicg(jrad))
                    enddo
                 enddo

              case (2) ! residu = Ix

                 residu(i,1:nvar_bicg)=zero
                 do irad=1,nvar_bicg
                    do jrad=1,nvar_bicg
                       residu(i,irad)=residu(i,irad)+mat_residual_glob(ind_res,irad,jrad)*var_bicg(ind_cell(i),jrad,i_y)
                    enddo
                 enddo

              case (4) ! reinitialize rad_flux for this level
                 rad_flux(ind_cell(i),:) = zero

!neil
!!$                 do idim=1,ndim
!!$                    do irad = 1,nvar_bicg
!!$                       var_bicg(ind_cell(i),irad,11+(idim-1)*2) = var_rad_subset(1,idim,nrad+2+irad)
!!$                       var_bicg(ind_cell(i),irad,12+(idim-1)*2) = var_rad_subset(3,idim,nrad+2+irad)
!!$                    enddo
!!$                 enddo
!neil

              case (6) ! residu = Ix

                 residu(i,1:nvar_bicg)=zero
                 do irad=1,nvar_bicg
                    do jrad=1,nvar_bicg
                       residu(i,irad)=residu(i,irad)+mat_residual_glob(ind_res,irad,jrad)*var_bicg(ind_cell(i),jrad,6)
                    enddo
                 enddo

              end select

           endif
        enddo

        do idim = 1,ndim

           select case (compute)! Getting val_g and val_d
           case(1)
              do i=1,ngrid
                 if(son(ind_cell(i)) == 0 )then
                    do irad=1,nvar_bicg
                       val_g(i,irad) = uold(cell_left (i,idim),ind_bicg(irad))
                       val_d(i,irad) = uold(cell_right(i,idim),ind_bicg(irad))
                    enddo
                 end if
              end do

           case(2)
              do i=1,ngrid
                 if(son(ind_cell(i)) == 0 )then
                    do irad=1,nvar_bicg
                       val_g(i,irad) = var_bicg(cell_left (i,idim),irad,i_y)
                       val_d(i,irad) = var_bicg(cell_right(i,idim),irad,i_y)
                    enddo
                 end if
              end do

           case(4)
              do i=1,ngrid
                 if(son(ind_cell(i)) == 0 )then
                    do irad=1,nvar_bicg
                       val_g(i,irad) = unew(cell_left (i,idim),ind_bicg(irad))
                       val_d(i,irad) = unew(cell_right(i,idim),ind_bicg(irad))
                    enddo
                 end if
              end do

           case(6)
              do i=1,ngrid
                 if(son(ind_cell(i)) == 0 )then
                    do irad=1,nvar_bicg
                       val_g(i,irad) = var_bicg(cell_left (i,idim),irad,6)
                       val_d(i,irad) = var_bicg(cell_right(i,idim),irad,6)
                    enddo
                 end if
              end do

           end select

           do i=1,ngrid
              if(son(ind_cell(i)) == 0 )then

                 select case (nbor_ilevel(i,2*idim-1)) ! Gather main characteristics of left neighbour

                 case (1)
                    do irad=1,nvar_bicg
                       if (compute==2  .or. compute==6) then
                          val_g(i,irad) = zero
                       else
                          val_g(i,irad) = uold(cell_left(i,idim),ind_bicg(irad))/norm_bicg(irad)
                       endif
#if USE_FLD==1
                       phi_g (i,irad) = uold(cell_left(i,idim),ind_bicg(irad))/norm_bicg(irad)
#endif
                    enddo
                 case (0)
                    do irad=1,nvar_bicg
                       phi_g (i,irad)       = uold(cell_left(i,idim),ind_bicg(irad))
                    enddo
                 case (-1)
                    do irad=1,nvar_bicg
                       if (compute==2  .or. compute==6) then
                          val_g(i,irad) = zero
                       else
                          val_g(i,irad) = uold(cell_left(i,idim),ind_bicg(irad))/norm_bicg(irad)
                       endif
#if USE_FLD==1
                       phi_g (i,irad) = uold(cell_left(i,idim),ind_bicg(irad))/norm_bicg(irad)
#endif
                    enddo
                 end select

                 select case (nbor_ilevel(i,2*idim)) ! Gather main characteristics of right neighbour
                 case (1)
                    do irad=1,nvar_bicg
                       if (compute==2  .or. compute==6) then
                          val_d(i,irad) = zero
                       else
                          val_d(i,irad) = uold(cell_right(i,idim),ind_bicg(irad))/norm_bicg(irad)
                       endif
#if USE_FLD==1
                       phi_d (i,irad) = uold(cell_right(i,idim),ind_bicg(irad))/norm_bicg(irad)
#endif
                    enddo
                 case (0)
                    do irad=1,nvar_bicg
                       phi_d (i,irad)  = uold(cell_right(i,idim),ind_bicg(irad))
                    enddo
                 case (-1)
                    do irad=1,nvar_bicg
                       if (compute==2  .or. compute==6) then
                          val_d(i,irad) = zero
                       else
                          val_d(i,irad) = uold(cell_right(i,idim),ind_bicg(irad))/norm_bicg(irad)
                       endif
#if USE_FLD==1
                       phi_d (i,irad) = uold(cell_right(i,idim),ind_bicg(irad))/norm_bicg(irad)
#endif
                    enddo
                 end select
              end if
           end do

           if (compute ==4)then ! Computing and saving flux to the coarser ilevel

              do i=1,ngrid
                 if(son(ind_cell(i)) == 0)then

                    if(.not.store_matrix)then
                       call compute_coeff_left_right_in_cell(ind_cell(i),idim,cell_left(i,idim),cell_right(i,idim),nbor_ilevel(i,1:2*ndim),dx_loc,coeff_glob_left(1,:,:,idim),coeff_glob_right(1,:,:,idim))
                       ind_res=1
                    else
                       ind_res=ind_cell(i)
                    endif

                    if( nbor_ilevel(i,2*idim-1) == -1)then
                       do irad=1,nvar_bicg
#if USE_FLD==1
                          phi_c(i,irad) = uold(ind_cell(i),firstindex_er+irad)
                          rad_flux(cell_left(i,idim),irad)  = rad_flux(cell_left(i,idim),irad)  + &
                               & coeff_glob_left(ind_res,irad,irad,idim)*( alpha_imp * (unew(ind_cell(i),ind_bicg(irad)) - val_g(i,irad)) + (one-alpha_imp)*(phi_c(i,irad) - phi_g(i,irad))) 
#endif
#if USE_M_1==1
                          !if M1, C_g(i,irad) is just 0.5 , 1.0 , 1.5
!                          rad_flux(cell_left (i,idim),irad) =  rad_flux(cell_left(i,idim),irad)  + C_g(i,irad)*var_bicg(ind_cell(i),irad,10+(idim-1)*2)
                          rad_flux(cell_left (i,idim),irad) = 0. ! rad_flux(cell_left(i,idim),irad)  + C_g(i,irad)*(- a1d_glob(ind_cell(i),irad,          jrad)*val_g(i,jrad)+precond_bicg(ind_cell(i),irad,jrad)*unew(ind_cell(i),ind_bicg(irad)))
#endif
                       enddo
                    end if

                    if( nbor_ilevel(i,2*idim)   == -1 )then
                       do irad=1,nvar_bicg
#if USE_FLD==1
                          phi_c(i,irad) = uold(ind_cell(i),firstindex_er+irad)
                          rad_flux(cell_right(i,idim),irad) = rad_flux(cell_right(i,idim),irad) + &
                               & coeff_glob_right(ind_res,irad,irad,idim)*( alpha_imp * (unew(ind_cell(i),ind_bicg(irad)) - val_d(i,irad)) + (one-alpha_imp)*(phi_c(i,irad) - phi_d(i,irad)))
#endif
#if USE_M_1==1
                          !if M1, C_d(i,irad) is just 0.5 , 1.0 , 1.5
!                          rad_flux(cell_right (i,idim),irad) =  rad_flux(cell_right(i,idim),irad)  + C_d(i,irad)*var_bicg(ind_cell(i),irad,11+(idim-1)*2)
                          rad_flux(cell_left (i,idim),irad) = 0.! rad_flux(cell_left(i,idim),irad)  + C_d(i,irad)*(- a1d_glob(ind_cell(i),irad,nvar_bicg+jrad)*val_d(i,jrad)+precond_bicg(ind_cell(i),irad,jrad)*unew(ind_cell(i),ind_bicg(irad)))
#endif
                       enddo
                    end if

                 end if
              end do
           end if


           do i=1,ngrid
              if(son(ind_cell(i)) == 0 )then

                 if(.not.store_matrix)then
                    call compute_coeff_left_right_in_cell(ind_cell(i),idim,cell_left(i,idim),cell_right(i,idim),nbor_ilevel(i,1:2*ndim),dx_loc,coeff_glob_left(1,:,:,idim),coeff_glob_right(1,:,:,idim))
                    ind_res=1
                 else
                    ind_res=ind_cell(i)
                 endif
                 
                 select case (compute)
                    
                 case (1) ! compute b-Ax from b-Ix by adding intern flux

                    do irad=1,nvar_bicg
                       do jrad=1,nvar_bicg
                          residu(i,irad) = residu(i,irad) + coeff_glob_left(ind_res,irad,jrad,idim)*val_g(i,jrad) + coeff_glob_right(ind_res,irad,jrad,idim)*val_d(i,jrad)
                       enddo
#if USE_FLD==1
                       residu(i,irad) = residu(i,irad) - (coeff_glob_left(ind_res,irad,irad,idim)+coeff_glob_right(ind_res,irad,irad,idim))* uold(ind_cell(i),ind_bicg(irad))
#endif
                    enddo

                 case (2) ! compute Ap from Ip by adding intern flux

                    do irad=1,nvar_bicg
                       do jrad=1,nvar_bicg
                          residu(i,irad) = residu(i,irad) - (coeff_glob_left(ind_res,irad,jrad,idim)*val_g(i,jrad) + coeff_glob_right(ind_res,irad,jrad,idim)*val_d(i,jrad))*alpha_imp
                       enddo
#if USE_FLD==1
                       residu(i,irad) = residu(i,irad) + (coeff_glob_left(ind_res,irad,irad,idim)+coeff_glob_right(ind_res,irad,irad,idim))* var_bicg(ind_cell(i),irad,i_y)*alpha_imp
#endif
                    enddo

                 case (6) ! compute Ap* from Ip by adding intern flux

                    do irad=1,nvar_bicg
                       do jrad=1,nvar_bicg
                          residu(i,irad) = residu(i,irad) - (coeff_glob_left(ind_res,irad,jrad,idim)*val_g(i,jrad) + coeff_glob_right(ind_res,irad,jrad,idim)*val_d(i,jrad))*alpha_imp
                       enddo
#if USE_FLD==1
                       residu(i,irad) = residu(i,irad) + (coeff_glob_left(ind_res,irad,irad,idim)+coeff_glob_right(ind_res,irad,irad,idim))* var_bicg(ind_cell(i),irad,6)*alpha_imp
#endif
                    enddo

                 end select
              end if
           end do


        end do !ndim


        select case (compute)
           ! get the result out

        case (1)
           do i=1,ngrid
              if(son(ind_cell(i)) == 0 )then
                 do irad=1,nvar_bicg
                    var_bicg(ind_cell(i),irad,1) = residu(i,irad)
                 enddo
              end if
           end do

        case (2)
           do i=1,ngrid
              if(son(ind_cell(i)) == 0 )then
                 do irad=1,nvar_bicg
                    var_bicg(ind_cell(i),irad,3) = residu(i,irad)
                 enddo
              end if
           end do

        case (6)
           do i=1,ngrid
              if(son(ind_cell(i)) == 0 )then
                 do irad=1,nvar_bicg
                    var_bicg(ind_cell(i),irad,8) = residu(i,irad)
                 enddo
              end if
           end do

        end select

     end do ! twotodim
  end do ! ncache

  return

end subroutine cmp_matrix_vector_product

!###########################################################
!###########################################################
!###########################################################
!###########################################################

function lambda_fld(R)
  use radiation_parameters
  use const
  implicit none
  real(dp)::R,lambda_fld

  lambda_fld = one/three
  if(i_fld_limiter==i_fld_limiter_levermore) lambda_fld =(2.0d0+r)/(6.0d0+2.0d0*R+R**2)! (one/tanh(R)-one/R) / R
  if(i_fld_limiter==i_fld_limiter_minerbo) then 
     if(R .le. three/two) then
        lambda_fld = two/(three+sqrt(nine+12.0_dp*R*R))
     else
        lambda_fld = one/(one + R + sqrt(one+two*R))
     end if
  end if
  return 
end function lambda_fld

!################################################################
!################################################################
!################################################################ 
!################################################################

subroutine cmp_energy(Etype)
  use hydro_commons
  use radiation_parameters
  use const
  use units_commons
  implicit none
  integer,intent(in) :: Etype ! Etype=1 : beginning ; Etype=2 : end
  integer ::i,idim,this,ivar,igroup,irad
  real(dp)::usquare,Cv,eps,ekin,emag,rho,erad_loc
  real(dp)::tp_loc,cmp_temp
  
  real(dp)::sum_dust
#if NDUST>0  
  integer::idust
#endif
  do i=1,nb_ind
     this = liste_ind(i)
     rho   = uold(this,1)

     ! Compute total kinetic energy
     usquare=zero
     do idim=1,ndim
        usquare=usquare+(uold(this,idim+1)/uold(this,1))**2
     end do
     ekin  = rho*usquare*half

     ! Compute total magnetic energy
     emag = zero
     do ivar=1,3
        emag = emag + ((uold(this,5+ivar)+uold(this,nvar+ivar))**2)/eight
     end do

     if(Etype==1)then
        ! Compute total non-thermal+radiative energy
        erad_loc = zero
        do igroup=1,nener
           erad_loc = erad_loc + uold(this,8+igroup)
        enddo
        
        eps = uold(this,5)-ekin-emag-erad_loc
        if(energy_fix)eps = uold(this,nvar) ! use energy fix for collapse
        
        Tp_loc = cmp_temp(this)
        Cv = eps/Tp_loc

        unew(this,nvar+1) = Cv
        uold(this,nvar  ) = Tp_loc

        do irad=1,nvar_trad
           uold(this,ind_trad(irad))=uold(this,ind_trad(irad))/norm_trad(irad)
           if(is_radiative_energy(irad)) uold(this,ind_trad(irad)) = max(uold(this,ind_trad(irad)),eray_min/scale_E0)
           unew(this,ind_trad(irad))=uold(this,ind_trad(irad))
        enddo

     elseif(Etype==2)then

        unew(this,nvar)=unew(this,nvar)*unew(this,nvar+1)

        do irad=1,nvar_trad
           if(is_radiative_energy(irad)) unew(this,ind_trad(irad)) = max(unew(this,ind_trad(irad)),eray_min/scale_E0)
           unew(this,ind_trad(irad))=unew(this,ind_trad(irad))*norm_trad(irad)
 !          unew(this,ind_trad(irad))=uold(this,ind_trad(irad))*norm_trad(irad)
           uold(this,ind_trad(irad))=unew(this,ind_trad(irad))
        enddo

        eps = unew(this,nvar)
        uold(this,5) = eps + ekin + emag
        do igroup=1,nener
           uold(this,5) = uold(this,5) + uold(this,8+igroup)
        enddo

     end if
  end do



end subroutine cmp_energy

!################################################################
!################################################################
!################################################################ 
!################################################################
function cmp_temp(this)
  use hydro_commons
!  use constants,only:kB,mH
  use radiation_parameters
  use const
  use units_commons
  implicit none
  integer,intent(in) ::this
  integer ::idim,ivar,igrp,ht
  real(dp)::usquare,eps,ekin,emag,rho,erad_loc
  real(dp)::cmp_temp
  real(dp) :: sum_dust
#if NDUST>0
  integer :: idust
#endif  
  rho   = uold(this,1)
!!$  Cv    = rho*kB/(mu_gas*mH*(gamma-one))/scale_v**2

  ! Compute total kinetic energy
  usquare=zero
  do idim=1,ndim
     usquare=usquare+(uold(this,idim+1)/uold(this,1))**2
  end do
  ekin  = rho*usquare*half

  ! Compute total magnetic energy
  emag = zero
  do ivar=1,3
     emag = emag + ((uold(this,5+ivar)+uold(this,nvar+ivar))**2)/eight
  end do

  ! Compute total non-thermal+radiative energy
  erad_loc  = zero
  do igrp=1,nener
     erad_loc = erad_loc + uold(this,8+igrp) 
  enddo
  eps = uold(this,5)-ekin-emag-erad_loc
  if(energy_fix)eps = uold(this,nvar) ! use energy fix for collapse


  sum_dust =0.0d0
#if NDUST>0
  do idust = 1, ndust
     sum_dust = sum_dust + uold(this,firstindex_ndust+idust)/uold(this,1)
  end do
#endif
  
  call temperature_eos((1.0d0-sum_dust)*rho,eps,cmp_temp,ht)

  return

end function cmp_temp
!################################################################
!################################################################
!################################################################ 
!################################################################
function nu_surf(Er1,Er2,ind1,ind2,dx)
  use hydro_commons
  use const
  implicit none
  integer ::ind1,ind2
  real(dp),INTENT(IN)::Er2,Er1,dx
  real(dp)::nu_surf,nu_harmo,nu_ari

  nu_ari=(Er2+Er1)*half

  nu_harmo=max(Er2*Er1/nu_ari,four/(three*dx))
  nu_surf = nu_ari

  nu_surf=min(nu_harmo,nu_ari)

  return 
end function nu_surf

!###########################################################
!###########################################################
!###########################################################
!###########################################################

subroutine compute_residual_in_cell(i,ilevel,vol_loc,residual,mat_residual)

  use hydro_parameters,only:nvar
  use hydro_commons
  use radiation_parameters
  use cloud_module, only : rt_feedback
  use units_commons
  use const
  use constants, only: eV2erg
#ifdef RT
  use pm_commons !hybrid RT, needs sink quantities
  use rt_hydro_commons !hybrid RT
  use rt_cooling_module !hybrid RT
#endif

  implicit none
  integer,intent(in)::i,ilevel
  real(dp),intent(in)::vol_loc
  real(dp),dimension(nvar_bicg,nvar_bicg),intent(out)::mat_residual
  real(dp),dimension(nvar_bicg          ),intent(out)::residual

  real(dp)::rho,Told_norm,Told,cv,lhs,rhs,planck_ana,radiation_source,deriv_radiation_source,cal_Teg,Trold
  integer::igrp,igroup,im1
  real(dp),dimension(ngrp)::wdtB,wdtE,source,deriv
  real(dp)::ambi_heating,ohm_heating,nimhd_heating,protostellar_heating
!#if NIMHD==1
!  real(dp)::bcell2,bx,by,bz,jsquare,jx,jy,jz,etaohmdiss,betaad,ionisrate
!#endif
#ifdef RT
  real(dp)::scale_Np,scale_Fp !hybrid RT
#endif

  real(dp)::scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

#ifdef RT
  call rt_units(scale_Np,scale_Fp)
#endif

  rho       = uold(i,1          )
  Told_norm = uold(i,ind_trad(1))
  Told      = Told_norm * Tr_floor
  Cv        = unew(i,nvar+1)

  ambi_heating=zero
  ohm_heating=zero
  nimhd_heating=zero

!#if NIMHD==1
!  if(radiative_nimhdheating_in_cg)then
!     bx=0.5d0*(uold(i,6)+uold(i,nvar+1))
!     by=0.5d0*(uold(i,7)+uold(i,nvar+2))
!     bz=0.5d0*(uold(i,8)+uold(i,nvar+3))
!     bcell2=(bx**2+by**2+bz**2)
!     jx=uold(i,nvar-3)
!     jy=uold(i,nvar-2)
!     jz=uold(i,nvar-1)
!     jsquare=(jx**2+jy**2+jz**2)
!     ionisrate=default_ionisrate
!
!     if(nmagdiffu .eq. 1 .or. nmagdiffu2 .eq. 1 )ohm_heating=jsquare*etaohmdiss(rho,bcell2,Told,ionisrate)*dt_imp
!     
!     if(nambipolar .eq. 1 .or. nambipolar2 .eq.1 )then
!        ambi_heating = (jy*bz-jz*by)**2+(jz*bx-jx*bz)**2+(jx*by-jy*bx)**2
!        ambi_heating = ambi_heating * betaad(rho,bcell2,Told,ionisrate)*dt_imp
!     endif
!     nimhd_heating=ambi_heating+ohm_heating
!  end if
!#endif  

  lhs=zero
  rhs=zero
  do igrp=1,ngrp
     ! Store radiation_source, deriv_radiation_source and planck opacity to save cpu time
     Trold=cal_Teg(uold(i,firstindex_er+igrp)*scale_E0,igrp)
     source(igrp)=radiation_source(Told,igrp)
     deriv(igrp)=deriv_radiation_source(Told,igrp)
     wdtB(igrp) = C_cal*dt_imp*planck_ana(rho*scale_d,Told,Told ,igrp,in_sink(i))/scale_kappa
     if(rt_feedback)then
        wdtE(igrp) = C_cal*dt_imp*planck_ana(rho*scale_d,Told,Told,igrp,in_sink(i))/scale_kappa
     else
        wdtE(igrp) = C_cal*dt_imp*planck_ana(rho*scale_d,Told,Trold,igrp,in_sink(i))/scale_kappa
     end if
     lhs=lhs+P_cal*wdtB(igrp)*deriv(igrp)/scale_E0
     rhs=rhs-P_cal*wdtB(igrp)*(source(igrp)/scale_E0-Told*deriv(igrp)/scale_E0)
  enddo

  mat_residual(:,:) = zero
  residual    (:  ) = zero
  do igroup=1,ngrp
     mat_residual(igroup,igroup) =  (one+wdtE(igroup))*vol_loc

     ! Terms of coupling radiative groups
     do igrp=1,ngrp
        mat_residual(igroup,igrp) = mat_residual(igroup,igrp) - wdtB(igroup)*(deriv(igroup)*P_cal*wdtE(igrp)/scale_E0/(cv+lhs))*vol_loc
     enddo
     protostellar_heating = 0.0d0
#ifdef RT
     if(rt_protostar_m1 .and. rt_advect .and. nsink >0 .and. Teff_sink(1) > 0.0d0)then !heating by M1 photons into CvT
       do im1=1,ngroups
           protostellar_heating = rtuold(i,iGroups(im1))*scale_Np*group_egy(im1)*ev2erg*planck_ana(rho*scale_d,Told,Teff_sink(1),1,in_sink(i))*rt_c_cgs(ilevel)*z_ave / (scale_d*scale_v**2)
        end do
     end if
#endif

     residual(igroup) = uold(i,firstindex_er+igroup)*vol_loc  &
          & + vol_loc*wdtB(igroup)*(source(igroup)/scale_E0-Told*deriv(igroup)/scale_E0) &
          & + vol_loc*wdtB(igroup)*deriv(igroup)/scale_E0*protostellar_heating*dt_imp*scale_t/(cv+lhs) &
          & + vol_loc*wdtB(igroup)*deriv(igroup)/scale_E0*(cv*Told+rhs+nimhd_heating)/(cv+lhs)
  enddo

  return

end subroutine compute_residual_in_cell

!###########################################################
!###########################################################
!###########################################################
!###########################################################

subroutine compute_coeff_left_right_in_cell(i,idim,cell_left,cell_right,nbor_ilevel,dx_loc,coeff_left,coeff_right)

  use hydro_parameters,only:nvar
  use hydro_commons
  use radiation_parameters
  use units_commons
  use const
  use cloud_module, only : rt_feedback !hybrid RT
  
  implicit none
  integer,intent(in)::i,idim,cell_left,cell_right
  integer,dimension(2*ndim),intent(in)::nbor_ilevel
  real(dp),intent(in)::dx_loc
  real(dp),dimension(nvar_bicg,nvar_bicg),intent(out)::coeff_left,coeff_right

  real(dp)::rho,Told,Trold,cal_Teg,cmp_temp,rosseland_ana,lambda,lambda_fld,R,nu_surf,surf_loc
  integer::igroup,irad
  real(dp),dimension(nvar_bicg)::C_g,C_d,phi_g,phi_c,phi_d,nu_g,nu_c,nu_d

  real(dp)::scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  surf_loc = dx_loc**(ndim-1)

  select case (nbor_ilevel(2*idim-1)) ! Gather main characteristics of left neighbour

  case (1)

     if (robin  > zero) then
        C_g(:) = one/robin
     else
        C_g(:) = zero
     endif

     do irad=1,nvar_bicg
        phi_g (irad) = uold(cell_left,ind_bicg(irad))/norm_bicg(irad)
     enddo

     Told = cmp_temp(cell_left)
     rho  = scale_d * max(uold(cell_left,1),smallr)
     do igroup=1,ngrp
        Trold = cal_Teg(uold(cell_left,firstindex_er+igroup)*scale_d*scale_v**2,igroup)
        if(rt_feedback)then
           nu_g(igroup) = rosseland_ana(rho,Told,Told,igroup,in_sink(cell_left)) / scale_kappa
        else
           nu_g(igroup) = rosseland_ana(rho,Told,Trold,igroup,in_sink(cell_left)) / scale_kappa
        end if
        if(nu_g(igroup)*dx_loc .lt. min_optical_depth) nu_g(igroup)=min_optical_depth/dx_loc
     enddo

  case (0)

     do irad=1,nvar_bicg
        phi_g (irad)       = uold(cell_left,ind_bicg(irad))
        C_g   (irad)       = one
     enddo
     do igroup=1,ngrp
        nu_g  (igroup)       = kappaR_bicg(cell_left,igroup)
     enddo

  case (-1)

     C_g(:) = 1.5_dp
     do irad=1,nvar_bicg
        phi_g (irad) = uold(cell_left,ind_bicg(irad))/norm_bicg(irad)
     enddo

     Told = cmp_temp(cell_left)
     rho  = scale_d * max(uold(cell_left,1),smallr)
     do igroup=1,ngrp
        Trold = cal_Teg(uold(cell_left,firstindex_er+igroup)*scale_d*scale_v**2,igroup)
        if(rt_feedback)then
           nu_g  (igroup) = rosseland_ana(rho,Told,Told,igroup,in_sink(cell_left)) / scale_kappa
        else
           nu_g  (igroup) = rosseland_ana(rho,Told,Trold,igroup,in_sink(cell_left)) / scale_kappa
        end if
        if(nu_g(igroup)*2.0d0*dx_loc .lt. min_optical_depth) nu_g(igroup)=min_optical_depth/(2.0d0*dx_loc)
     enddo

  end select

  select case (nbor_ilevel(2*idim)) ! Gather main characteristics of right neighbour

  case (1)

     if (robin  > zero) then
        C_d(:) = one/robin
     else
        C_d(:) = zero
     endif

     do irad=1,nvar_bicg
        phi_d (irad) = uold(cell_right,ind_bicg(irad))/norm_bicg(irad)
     enddo

     Told = cmp_temp(cell_right)
     rho  = scale_d * max(uold(cell_right,1),smallr)
     do igroup=1,ngrp
        Trold = cal_Teg(uold(cell_right,firstindex_er+igroup)*scale_d*scale_v**2,igroup)
        if(rt_feedback)then
           nu_d  (igroup) = rosseland_ana(rho,Told,Told,igroup,in_sink(cell_right)) / scale_kappa
        else
           nu_d  (igroup) = rosseland_ana(rho,Told,Trold,igroup,in_sink(cell_right)) / scale_kappa
        end if
        if(nu_d(igroup)*dx_loc .lt. min_optical_depth) nu_d(igroup)=min_optical_depth/dx_loc
     enddo

  case (0)

     do irad=1,nvar_bicg
        phi_d (irad)  = uold(cell_right,ind_bicg(irad))
        C_d   (irad)  = one
     enddo
     do igroup=1,ngrp
        nu_d  (igroup)  = kappaR_bicg(cell_right,igroup)
     enddo

  case (-1)

     C_d(:) = 1.5_dp
     do irad=1,nvar_bicg
        phi_d (irad) = uold(cell_right,ind_bicg(irad))/norm_bicg(irad)
     enddo

     Told = cmp_temp(cell_right)
     rho  = scale_d * max(uold(cell_right,1),smallr)
     do igroup=1,ngrp
        Trold = cal_Teg(uold(cell_right,firstindex_er+igroup)*scale_d*scale_v**2,igroup)
        if(rt_feedback)then
           nu_d  (igroup) = rosseland_ana(rho,Told,Told,igroup,in_sink(cell_right)) / scale_kappa
        else
           nu_d  (igroup) = rosseland_ana(rho,Told,Trold,igroup,in_sink(cell_right)) / scale_kappa
        end if
        if(nu_d(igroup)*2.0d0*dx_loc .lt. min_optical_depth) nu_d(igroup)=min_optical_depth/(2.0d0*dx_loc)
     enddo

  end select

  do igroup=1,ngrp
     nu_c (igroup) = kappaR_bicg(i,igroup)
     C_g(igroup) = C_g(igroup) * nu_surf(nu_g(igroup),nu_c(igroup), cell_left ,i,dx_loc)
     C_d(igroup) = C_d(igroup) * nu_surf(nu_d(igroup),nu_c(igroup), cell_right,i,dx_loc)

     phi_c(igroup) = uold(i,firstindex_er+igroup)

     if(C_g(igroup) > zero)then
        R = max(1.0e-10_dp,abs (phi_c(igroup)-phi_g(igroup)) /(half*(phi_c(igroup)+phi_g(igroup))))
        R = R / ( C_g(igroup) * dx_loc )

        lambda=lambda_fld(R)
        C_g(igroup) = C_cal*lambda *dt_imp*surf_loc/(dx_loc*C_g(igroup))
     end if

     if(C_d(igroup) > zero)then
        R = max(1.0e-10_dp,abs (phi_c(igroup)-phi_d(igroup)) /(half*(phi_c(igroup)+phi_d(igroup))))
        R = R / ( C_d(igroup) * dx_loc )

        lambda=lambda_fld(R)
        C_d(igroup) = C_cal*lambda *dt_imp*surf_loc/(dx_loc*C_d(igroup))

     end if

     coeff_left (igroup,igroup)=C_g(igroup)
     coeff_right(igroup,igroup)=C_d(igroup)

  enddo

  return

end subroutine compute_coeff_left_right_in_cell
