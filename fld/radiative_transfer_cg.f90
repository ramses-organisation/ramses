subroutine diffusion_cg (ilevel,Nsub)
  use amr_commons
  use hydro_commons
  use cooling_module,ONLY:kB,mH,clight
  use radiation_parameters
  use units_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !===========================================================================
  ! Iterative solver with Conjugate Gradient method and Newton Raphson
  ! to solve A x = b
  !   r1      : stored in unew(i,1)
  !   p1      : stored in unew(i,2)
  !   Diag(A) : stored in unew(i,4)
  !   Ap1     : stored in unew(i,3)
  !  x1(i)    : stored in uold(i,irad)(i.e. new radiative energy at time n+1)
  !  b1(n)    : stored in unew(i,irad)(i.e. radiative energy at time n)
  ! x1(i-1)   : stored in unew(i,7)
  !  Tg(i)    : stored in unew(i,nvar+3)
  !  Tg(n)    : stored in unew(i,5)
  !===========================================================================
  integer,intent(IN)::ilevel,Nsub
  complex*16 :: final_sum
  real(dp)::error,error_ini,epsilon,error_nr,error_nrm1,error_nrm2,error_nrm3
  real(dp)::error_nr_loc,error_nr_all,error_cg_loc,error_cg_all
  real(dp)::alpha_cg,beta_cg,Cv,told,tnew,wdt,rho,dt_exp,wdtB,wdtE,Tr,Trold
  real(dp)::r2_old,r2,pAp,rhs_norm1
  real(dp)::density,planck_ana,rosseland_ana,norm_Er,temp,cal_Teg
  integer :: i,info,ind,iter,iskip,itermax,icpu,igroup,igrp
  integer :: this,iter_nr,nx_loc,nsub_imp,isub,nleaf_all,nleaf_tot
  real(dp)::radiation_source,deriv_radiation_source,rhs,lhs
  real(dp)::dx,dx_loc,surf_loc,vol_loc,scale
  real(dp)::min_ener,min_ener_all,max_ener,max_ener_all
  real(dp),dimension(1:ngrp)::dener
  logical::exist_leaf_cell=.true.
  real(dp)::ambi_heating,ohm_heating,nimhd_heating,bcell2,bx,by,bz,jsquare,jx,jy,jz,etaohmdiss,betaad
  if(verbose)write(*,111)
  if(numbtot(1,ilevel)==0)return

  ! Rescaling factors
  ! Mesh size at level ilevel
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  surf_loc = dx_loc**(ndim-1)
  vol_loc  = dx_loc**ndim

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

  do igroup=1,ngrp
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do icpu=1,ncpu
           do i=1,reception(icpu,ilevel)%ngrid
              rad_flux(reception(icpu,ilevel)%igrid(i)+iskip,igroup)=0.0d0
           end do
           do i=1,reception(icpu,ilevel-1)%ngrid
              rad_flux(reception(icpu,ilevel-1)%igrid(i)+iskip,igroup)=0.0d0
           end do
        end do
     end do
  enddo

  if (nb_ind == 0)then
!!$     print*,'No leaf-cell - myid=',myid
     exist_leaf_cell=.false.
  end if


  nleaf_tot=nb_ind
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(nb_ind,nleaf_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  nleaf_tot=nleaf_all
#endif
     
  if(nleaf_tot .eq. 0)then
     write(*,*)'No leaf cells at level',ilevel,'. Exiting CG'
     deallocate(liste_ind)
     return
  end if

  do i=1,nb_ind
     this = liste_ind(i)
     unew(this,1:firstindex_er+ngrp)=0.0d0
!     urad(this,1:ngrp)=0.0d0
     unew(this,nvar+1)=0.d0
 !    unew(this,nvar+2)=0.d0
     unew(this,nvar+3)=0.0d0
     unew(this,2)=0.0d0
     unew(this,3)=0.0d0
     divu(this)=0.0d0
     enew(this)=0.0d0
     mat_residual_glob(this,1:ngrp,1:ngrp)=0.0d0
     coeff_glob_left(this,1:ngrp,1:ngrp,1:ndim)=0.0d0
     coeff_glob_right(this,1:ngrp,1:ngrp,1:ndim)=0.0d0
  end do

  ! Set constants
  epsilon = epsilon_diff

  !===================================================================
  ! Compute gas temperature stored in unew(i,nvar+3) and in unew(i,5)
  ! Compute Cv_eos stored in unew(i,nvar+1)
  !===================================================================
  call cmp_energy(1)

  !===================================================================
  ! Begin of subcycles....
  !===================================================================
  dt_exp = dtnew(ilevel)
  dt_imp = dtnew(ilevel)

  dener=0.0d0
  do igroup=1,ngrp
     max_ener=0.0d0
     min_ener=1.0d30
     do i=1,nb_ind
        this = liste_ind(i)
        max_ener=max(max_ener, uold(liste_ind(i),firstindex_er+igroup))
        min_ener=min(min_ener, uold(liste_ind(i),firstindex_er+igroup))
     end do

     ! Compute maximum error
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(max_ener,max_ener_all,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
     max_ener=max_ener_all
     call MPI_ALLREDUCE(min_ener,min_ener_all,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,info)
     min_ener=min_ener_all
#endif
     dener(igroup)=max_ener/min_ener
  enddo



  nsub_imp=1
!!$  if(ngrp .gt. 1 .and. maxval(dener) .gt. 1.d4)then
!!$     nsub_imp=int(maxval(dener)**0.25)
!!$     dt_imp=dt_imp/real(nsub_imp)
!!$  endif
  if(ngrp .gt. 1 .and. maxval(dener) .gt. 1.d4)then
     nsub_imp=10
     dt_imp=dt_imp/real(nsub_imp)
  endif
  if(myid==1)then
     do igroup=1,ngrp
        print*,'ilevel',ilevel,'igroup',igroup,'MAXIMUM OF DENER=',dener(igroup),'NSUB_IMP=',nsub_imp
     end do
  end if

  do isub=1,nsub_imp

     ! DO Newton-Raphson
     error_nr=1.0d0
     error_nrm1=1.0d-10
     error_nrm2=1.0d-10
     error_nrm3=1.0d-10
     error_cg_loc=1.0d0

     iter_nr=0

     
!     do while((error_nr .gt. 1.0d-3) .and. (abs((error_nr-error_nrm1)/error_nrm1) .gt. epsilon_diff) .and. (iter_nr .lt. 1000) .and. (abs(error_nr-error_nrm2) .gt. 1.0d-3) .or. iter_nr .lt. ngrp)
        iter_nr=iter_nr+1
        !Initialize NR error
        error_nr_loc=0.0d0

        error_nrm3=error_nrm2
        error_nrm2=error_nrm1
        error_nrm1=error_nr

!!$        do i=1,nb_ind
!!$           this = liste_ind(i)
!!$           do igroup=1,ngrp
!!$              urad(this,igroup)=unew(this,firstindex_er+igroup)
!!$           enddo
!!$        end do

        do igroup=1,ngrp

           !===================================================
           ! Compute thermal coefficient :
           ! Rosseland opacity : stored in divu(indcell(i))
           ! Planck opacity    : stored in enew(indcell(i)) 
           !===================================================
           do i=1,nb_ind
              this = liste_ind(i)
              density = scale_d * max(uold(this,1),smallr)
              temp = unew(this,5)
              ! Compute Rosseland opacity (Compute kappa*rho)
              divu(this)= rosseland_ana(density,temp,temp,igroup) / scale_kappa       
              if(divu(this)*dx_loc .lt. min_optical_depth) divu(this)=min_optical_depth/dx_loc
              ! Store radiative energy
              enew(this)= unew(this,firstindex_er+igroup)
           end do


           ! compute prdivu-nuprdivu ansd store it unew(i,nvar+2) if implicit integration
!           call cmp_Prdivu(ilevel,igroup)

           ! Update boundaries
           call make_virtual_fine_dp(enew(1),ilevel)
           call make_virtual_fine_dp(unew(1,2),ilevel)
           call make_virtual_fine_dp(unew(1,3),ilevel)
           call make_virtual_fine_dp(unew(1,1),ilevel)
           call make_virtual_fine_dp(unew(1,4),ilevel)
           call make_virtual_fine_dp(divu(1),ilevel)
           call make_virtual_fine_dp(unew(1,5),ilevel)
           call make_virtual_fine_dp(unew(1,nvar+3),ilevel)
           call make_virtual_fine_dp(unew(1,nvar+1),ilevel)
!           call make_virtual_fine_dp(unew(1,nvar+2),ilevel)
           do igrp=1,ngrp
              call make_virtual_fine_dp(uold(1,firstindex_er+igrp),ilevel)
              call make_virtual_fine_dp(unew(1,firstindex_er+igrp),ilevel)
!              call make_virtual_fine_dp(urad(1,igrp),ilevel)
           enddo

           call make_boundary_diffusion(ilevel,igroup)
           call cmp_matrix_coeff(ilevel,igroup)
           !===================================================
           ! Compute r1 = b1 - A1x1 and store it into unew(i,1)
           ! Also set p1 = r1 and store it into unew(i,2)
           !===================================================
!ben           call cmp_matrixA (ilevel, 1,igroup)
           call cmp_matrixA2 (ilevel, 1,igroup)

           !        call make_virtual_reverse_dp(unew(1,1),ilevel)
           call make_virtual_fine_dp(unew(1,1),ilevel)

           !        call make_virtual_reverse_dp(unew(1,2),ilevel)
           call make_virtual_fine_dp(unew(1,2),ilevel)

           !===================================
           ! Compute right-hand side norm (max)
           !===================================
           call dot_product(unew(:,1),unew(:,1),rhs_norm1,final_sum)
!!$ben           call dot_product(unew(:,firstindex_er+igroup),unew(:,firstindex_er+igroup),norm_er,final_sum) ! compute radiative energy norm

           !===================================================
           ! Compute Preconditionner M=1/diag(A) and store it in unew(i,4)
           !===================================================
           !ben
           call make_boundary_diffusion(ilevel,igroup)
!ben           call cmp_matrixA (ilevel 3,igroup)
           call cmp_matrixA2 (ilevel, 3,igroup)

           !        call make_virtual_reverse_dp(unew(1,4),ilevel)
           call make_virtual_fine_dp(unew(1,4),ilevel)

           !====================================
           ! MAIN ITERATION LOOP
           !====================================     

           iter=0; itermax=5000


           error_ini=sqrt(rhs_norm1)
!!$ben           norm_er=sqrt(norm_er)!*vol_loc
           !     error_ini=sqrt(real(final_sum))
           error=error_ini
           error_cg_loc=1.0d0

           do while(error/error_ini>epsilon .or. error_cg_loc> epsilon)!error_ini/norm_er .gt. 1.0d-15)
 
              iter=iter+1
              if(iter>itermax)exit
              !============================================
              ! Compute z = Mr and store it into unew(i,3)
              !============================================

              do i=1,nb_ind
                 this = liste_ind(i)
                 unew(this,3) = unew(this,1) * unew(this,4)
              end do

              !====================================
              ! Compute scalar r.z
              !====================================

              call dot_product(unew(:,1),unew(:,3),r2,final_sum)
              r2=r2!real(final_sum)
              !        r2=real(final_sum)

              !====================================
              ! Compute beta factor
              !====================================

              if(iter==1)then
                 beta_cg = 0.0d0
              else
                 beta_cg = r2/r2_old
              end if
              r2_old=r2

              !====================================
              ! Recurrence on p = z + beta*p
              !====================================
              call cX_plus_Y_to_Z (beta_cg,unew(:,2),unew(:,3),unew(:,2))
              ! Update boundaries
              call make_boundary_diffusion(ilevel,igroup)
              call make_virtual_fine_dp(unew(1,2),ilevel)

              !==============================================
              ! Compute q1 = Ap1 and store it into unew(i,3)
              !==============================================
!ben              call cmp_matrixA (ilevel, 2,igroup)
              call cmp_matrixA2 (ilevel, 2,igroup)

              !        call make_virtual_reverse_dp(unew(1,3),ilevel)
              call make_virtual_fine_dp(unew(1,3),ilevel)

              !====================================
              ! Compute p.Ap scalar product
              !====================================

              call dot_product(unew(:,2),unew(:,3),pAp,final_sum)
              pap = pap!real(final_sum) !DDP
              !        pap = real(final_sum) !DDP

              !====================================
              ! Compute alpha factor
              !====================================
              alpha_cg = r2/pAp

              !====================================
              ! Recurrence on x = x + alpha*p
              !====================================
              error_cg_loc=0.0d0
              do i=1,nb_ind
                 this = liste_ind(i)
                 !        unew(liste_ind(i),firstindex_er+igroup) = max(unew(liste_ind(i),firstindex_er+igroup),eray_min/scale_E0)
                 error_cg_loc=max(error_cg_loc, abs((alpha_cg*unew(liste_ind(i),2))/enew(this)))
              end do

              ! Compute maximum error
#ifndef WITHOUTMPI
              call MPI_ALLREDUCE(error_cg_loc,error_cg_all,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
              error_cg_loc=error_cg_all
#endif

              call cX_plus_Y_to_Z (alpha_cg,unew(:,2),unew(:,firstindex_er+igroup),unew(:,firstindex_er+igroup))


              !====================================
              ! Recurrence on r (unew(i,1))
              !====================================

              call cX_plus_Y_to_Z (- alpha_cg ,unew(:,3),unew(:,1),unew(:,1))

              !===================================
              ! Compute right-hand side norm (max)
              !===================================
              call dot_product(unew(:,1),unew(:,1),rhs_norm1,final_sum)

              error=SQRT(rhs_norm1)
              !        error=SQRT(real(final_sum))
              !        error = error_cg_loc

              if(verbose) write(*,112)iter,error,error/error_ini,error_cg_loc,error_ini

           end do
           ! End main iteration loop

           if(iter >= itermax)then
              if(myid==1)write(*,*)'Diffusion fail to converge...'
           end if

           !====================================
           ! Copie des flux
           !====================================
!           call cmp_matrixA (ilevel, 4,igroup)
           call cmp_matrixA2 (ilevel, 4,igroup)

           if(myid==1) write(*,117)ilevel,igroup,iter,error/error_ini,error_cg_loc
!!$           if(myid==1) write(*,*)'igroup: ',igroup,' CG :',iter, 'error L2=',error/error_ini,error_ini,norm_er,'error Linf='error_cg_loc
           niter=niter+iter

           do i=1,nb_ind
              this = liste_ind(i)
              unew(liste_ind(i),firstindex_er+igroup) = max(unew(liste_ind(i),firstindex_er+igroup),eray_min/scale_E0)
              error_nr_loc=max(error_nr_loc, abs((enew(this)-unew(liste_ind(i),firstindex_er+igroup))/enew(this)))
           end do

           ! Compute maximum error
#ifndef WITHOUTMPI
           call MPI_ALLREDUCE(error_nr_loc,error_nr_all,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
           error_nr_loc=error_nr_all
#endif

        end do
        ! End loop over groups

!!$        !====================================
!!$        ! Update gas temperature
!!$        !====================================
!!$        do i=1,nb_ind
!!$
!!$           rho = uold(liste_ind(i),1)
!!$           Cv = unew(liste_ind(i),nvar+1)
!!$           Told= unew(liste_ind(i),5)
!!$
!!$           rhs=0.d0
!!$           lhs=0.d0
!!$           do igrp=1,ngrp
!!$              wdt = C_cal*dt_imp*planck_ana(rho*scale_d,Told,igrp)/scale_kappa
!!$              rhs=rhs-P_cal*wdt*(radiation_source(Told,igrp)/scale_E0-Told*deriv_radiation_source(Told,igrp)/scale_E0 &
!!$                   & -unew(liste_ind(i),firstindex_er+igrp))
!!$
!!$              lhs=lhs+P_cal*wdt*deriv_radiation_source(Told,igrp)/scale_E0
!!$           enddo
!!$
!!$           Tnew = (cv*unew(liste_ind(i),nvar+3)+rhs)/(cv+lhs)
!!$
!!$           error_nr_loc=max(error_nr_loc, abs((Told-Tnew)/Tnew))
!!$
!!$           unew(liste_ind(i),5)=Tnew
!!$
!!$        end do
!!$        ! Compute global norms
!!$#ifndef WITHOUTMPI
!!$        call MPI_ALLREDUCE(error_nr_loc,error_nr_all,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
!!$        error_nr_loc=error_nr_all
!!$#endif

        error_nr=error_nr_loc
!!$     error_nr=1.d-10
        if(ngrp ==1)error_nr=1.0d-10 ! No NR iterations for 1 group
        do igroup=1,ngrp
           call make_boundary_diffusion(ilevel,igroup)
        end do

!!$        if(myid==1)write(*,115)ilevel,iter,error,error/error_ini

        if(myid==1)write(*,116)ilevel,iter_nr,error_nr

     end do
     !End loop over NR iteration

        !====================================
        ! Update gas temperature
        !====================================
        do i=1,nb_ind

           rho = uold(liste_ind(i),1)
           Cv = unew(liste_ind(i),nvar+1)
           Told= unew(liste_ind(i),5)

                 ambi_heating=zero
                 ohm_heating=zero
                 nimhd_heating=zero

#if NIMHD==1
                 if((nmagdiffu .eq. 1 .or. nambipolar .eq.1 .or. nmagdiffu2 .eq. 1 .or. nambipolar2 .eq.1) .and. radiative_nimhdheating_in_cg)then
                    
                    bx=0.5d0*(uold(liste_ind(i),6)+uold(liste_ind(i),nvar+1))
                    by=0.5d0*(uold(liste_ind(i),7)+uold(liste_ind(i),nvar+2))
                    bz=0.5d0*(uold(liste_ind(i),8)+uold(liste_ind(i),nvar+3))
                    bcell2=(bx**2+by**2+bz**2)
                    jx=uold(liste_ind(i),nvar-3)
                    jy=uold(liste_ind(i),nvar-2)
                    jz=uold(liste_ind(i),nvar-1)
                    jsquare=(jx**2+jy**2+jz**2) 
                    ionisrate=default_ionisrate
                    
                    if(nmagdiffu .eq. 1 )ohm_heating=jsquare*etaohmdiss(rho,bcell2,Told,ionisrate)*dt_imp
                    
                    if(nambipolar .eq. 1 )then
                       ambi_heating = (jy*bz-jz*by)**2+(jz*bx-jx*bz)**2+(jx*by-jy*bx)**2
                       ambi_heating = ambi_heating * betaad(rho,bcell2,Told,ionisrate)*dt_imp
                    endif
                    nimhd_heating=nimhd_heating+ohm_heating
                 end if
#endif  

           rhs=0.d0
           lhs=0.d0
           do igrp=1,ngrp
              Trold = cal_Teg(uold(liste_ind(i),firstindex_er+igrp)*scale_E0,igrp)

              wdtB = C_cal*dt_imp*planck_ana(rho*scale_d,Told,Told ,igrp)/scale_kappa
              wdtE = C_cal*dt_imp*planck_ana(rho*scale_d,Told,Told,igrp)/scale_kappa
!!$              wdtE = C_cal*dt_imp*planck_ana(rho*scale_d,Told,Trold,igrp)/scale_kappa

              rhs=rhs-P_cal*wdtB*(radiation_source(Told,igrp)/scale_E0-Told*deriv_radiation_source(Told,igrp)/scale_E0) &
                   & + P_cal*wdtE*unew(liste_ind(i),firstindex_er+igrp)
              
              lhs=lhs+P_cal*wdtB*deriv_radiation_source(Told,igrp)/scale_E0
           enddo

           Tnew = (cv*unew(liste_ind(i),nvar+3)+rhs+nimhd_heating)/(cv+lhs)

           error_nr_loc=max(error_nr_loc, abs((Told-Tnew)/Tnew))

           unew(liste_ind(i),5)=Tnew

        end do

     !update radiative energy and temperature
     do i=1,nb_ind
        do igrp=1,ngrp
           uold(liste_ind(i),firstindex_er+igrp)=unew(liste_ind(i),firstindex_er+igrp)
        end do
        unew(liste_ind(i),nvar+3)=unew(liste_ind(i),5)
     end do


!  end do
  !ENd loop over subcycles

  if(myid == 1)      print*,'niter tot=',niter

  !=============================
  ! Update energy value
  !=============================
  if(static) then
     do igroup=1,ngrp
        do i=1,nb_ind
           uold(liste_ind(i),firstindex_er+igroup) = unew(liste_ind(i),firstindex_er+igroup)*P_cal
        end do
     enddo
  else
     call cmp_energy(2)
  end if

  ! Update boundaries
  do igroup=1,ngrp
     call make_virtual_fine_dp(uold(1,firstindex_er+igroup),ilevel)
     if(ilevel .gt. levelmin)call make_virtual_reverse_dp(rad_flux(1,igroup),ilevel-1)
     !  if(ilevel .gt. levelmin)call make_virtual_fine_dp(rad_flux(1,igroup),ilevel-1)
     call make_virtual_reverse_dp(rad_flux(1,igroup),ilevel)
     !  call make_virtual_fine_dp(rad_flux(1,igroup),ilevel)
  enddo
  call make_virtual_fine_dp(uold(1,5),ilevel)
  call make_virtual_fine_dp(uold(1,nvar),ilevel)
!  call upload_fine(ilevel)
!  call upload_fine(ilevel-1)


111 format('   Entering diffusion_cg')
112 format('   ==> Step=',i5,' Error=',2(1pe10.3,1x),e23.15,es18.5)
115 format('   ==> Level=',i5,' Step=',i5,' Error=',2(1pe10.3,1x))
116 format('   ==> Level=',i5,' Step NR =',i5,' Error=',(1pe10.3,1x))
117 format('   ==> Level=',i5,' igroup=',i3,' Iteration CG=',i5,' Error L2=',(1pe10.3,1x),' Error Linf=',(1pe10.3,1x))

  deallocate(liste_ind)



contains

  !###########################################################
  !###########################################################

  subroutine dot_product(fact1,fact2,pdtvect,local_sum) !                pdtvect = sum(fact1*fact2)
    implicit none
    real(dp),dimension(1:ncoarse+twotondim*ngridmax),intent(IN)::fact1,fact2
    real(dp),intent(OUT)::pdtvect
    complex*16,intent(OUT)::local_sum

    real(dp)::pdtvect_all
    complex*16 ::global_sum
    integer::this

    pdtvect=0.0d0
    local_sum = cmplx(0.0d0,0.0d0)
    global_sum = cmplx(0.0d0,0.0d0)

    do i=1,nb_ind
       this = liste_ind(i)
       !call DDPDD (cmplx(fact1(this)*fact2(this), 0.0,dp), local_sum, 1, itype)
       pdtvect = pdtvect + fact1(this)*fact2(this)
    end do

    ! Compute global norms
#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(pdtvect,pdtvect_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
    pdtvect   = pdtvect_all
    !call MPI_ALLREDUCE(local_sum,global_sum,1,MPI_COMPLEX,MPI_SUMDD,MPI_COMM_WORLD,info)
    !local_sum = global_sum
    local_sum = pdtvect
#endif

  end subroutine dot_product
  
  !###########################################################
  !###########################################################
  subroutine cX_plus_Y_to_Z (cste,vectX,vectY,vectZ)! vectZ = cste*vectX+vectY
    implicit none
    real(dp),dimension(1:ncoarse+twotondim*ngridmax),intent(IN)::vectX,vectY
    real(dp),intent(IN)::cste
    real(dp),dimension(1:ncoarse+twotondim*ngridmax),intent(OUT)::vectZ
    
    
    do i=1,nb_ind
       vectZ(liste_ind(i)) = vectY(liste_ind(i)) + cste*vectX(liste_ind(i)) 
    end do

  end subroutine cX_plus_Y_to_Z


end subroutine diffusion_cg
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmp_matrixA (ilevel,compute,igroup)
  !------------------------------------------------------------------
  ! This routine computes the matrix A to vect_in and create vect_out
  ! compute = 1 : residu           	return B - Ax
  ! compute = 2 : Produit                 return  A.p
  ! compute = 3 : Preconditionner         return diag(A)
  ! compute = 4 : Compute flux in rad_flux
  !------------------------------------------------------------------

  use amr_commons
  use hydro_commons
  use cooling_module,ONLY:kB,mH,clight
  use radiation_parameters
  use units_commons
  use const

  implicit none

  integer,intent(IN)::compute,ilevel,igroup

  integer , dimension(1:nvector,1:2*ndim),save:: nbor_ilevel
  integer , dimension(1:nvector,1:ndim),save::   cell_left , cell_right , big_left, big_right
  integer ,dimension(1:nvector,0:2*ndim),save::  igridn
  integer ,dimension(1:nvector),save ::          ind_cell , ind_grid

  real(dp),dimension(1:nvector),save:: residu,C_g,C_d,nu_g,nu_c,nu_d
  real(dp),dimension(1:nvector),save:: phi_g,phi_c,phi_d,val_g,val_d


  integer :: i,idim,ind,igrid,ngrid,ncache,iskip,igrp,nx_loc
  integer :: supG,sub,supD

  real(dp):: radiation_source,deriv_radiation_source,rhs,lhs,cal_Teg
  real(dp):: dx,dx_loc,surf_loc,vol_loc,scale
  real(dp):: nu_surf,Cv,rho,wdtB,wdtE,Told,Trold,lambda,lambda_fld,R
  real(dp):: cmp_temp,rosseland_ana,planck_ana,Prdivu
  real(dp)::ambi_heating,ohm_heating,nimhd_heating,bcell2,bx,by,bz,jsquare,jx,jy,jz,etaohmdiss,betaad

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel

  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  surf_loc = dx_loc**(ndim-1)
  vol_loc  = dx_loc**ndim

  ! **************************** LOOP OVER CELLS ********************************** !

  residu = 0.0d0

  ! Loop over myid grids by vector sweeps
  ncache = active(ilevel)%ngrid
  do igrid=1,ncache,nvector


     ! Gather nvector grids
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i) = active(ilevel)%igrid(igrid+i-1)
     end do


     do i=1,ngrid
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

        select case (compute)

        case (1)
           ! residu = b - Ax
           do i=1,ngrid
              if(son(ind_cell(i)) == 0 )then
                 rho   = uold (ind_cell(i),1)
                 Cv = unew(ind_cell(i),nvar+1)
                 Prdivu = 0.0d0!unew(ind_cell(i),nvar+2)*dt_imp
                 Told  =  unew(ind_cell(i),5)

                 ambi_heating=zero
                 ohm_heating=zero
                 nimhd_heating=zero

#if NIMHD==1
                 bx=0.5d0*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
                 by=0.5d0*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
                 bz=0.5d0*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))
                 bcell2=(bx**2+by**2+bz**2)
                 jx=uold(ind_cell(i),nvar-3)
                 jy=uold(ind_cell(i),nvar-2)
                 jz=uold(ind_cell(i),nvar-1)
                 jsquare=(jx**2+jy**2+jz**2) 
                 ionisrate=default_ionisrate
                 
                 if((nmagdiffu .eq. 1 .or. nambipolar .eq.1 .or. nmagdiffu2 .eq. 1 .or. nambipolar2 .eq.1) .and. radiative_nimhdheating_in_cg)then
                    
                    if(nmagdiffu .eq. 1 )ohm_heating=jsquare*etaohmdiss(rho,bcell2,Told,ionisrate)*dt_imp
                    
                    if(nambipolar .eq. 1 )then
                       ambi_heating = (jy*bz-jz*by)**2+(jz*bx-jx*bz)**2+(jx*by-jy*bx)**2
                       ambi_heating = ambi_heating * betaad(rho,bcell2,Told,ionisrate)*dt_imp
                    endif
                    nimhd_heating=nimhd_heating+ohm_heating
                 end if
#endif  
                 
                 rhs=0.0d0
                 lhs=0.0d0
                 do igrp=1,ngrp
                    Trold=cal_Teg(uold(ind_cell(i),firstindex_er+igrp)*scale_E0,igrp)
                    wdtB = C_cal*dt_imp*planck_ana(rho*scale_d,Told,Told ,igrp)/scale_kappa
                    wdtE = C_cal*dt_imp*planck_ana(rho*scale_d,Told,Told,igrp)/scale_kappa
!!$                    wdtE = C_cal*dt_imp*planck_ana(rho*scale_d,Told,Trold,igrp)/scale_kappa
                    
                    lhs=lhs+P_cal*wdtB*deriv_radiation_source(Told,igrp)/scale_E0
                    rhs=rhs-P_cal*wdtB*(radiation_source(Told,igrp)/scale_E0-Told*deriv_radiation_source(Told,igrp)/scale_E0)
                 enddo

                 wdtB = C_cal*dt_imp*planck_ana(rho*scale_d,Told,Told ,igroup)/scale_kappa
                 residu(i) = uold(ind_cell(i),firstindex_er+igroup)*vol_loc  &
                      & + vol_loc*wdtB*(radiation_source(Told,igroup)/scale_E0-Told*deriv_radiation_source(Told,igroup)/scale_E0) &
                      & + vol_loc*wdtB*deriv_radiation_source(Told,igroup)/scale_E0*(cv*unew(ind_cell(i),nvar+3)+rhs+nimhd_heating)/(cv+lhs)

                 ! Terms of coupling radiative groups
                 do igrp=1,ngrp
                    if(igrp .ne. igroup) then
                       Trold=cal_Teg(uold(ind_cell(i),firstindex_er+igrp)*scale_E0,igrp)
                       wdtE = C_cal*dt_imp*planck_ana(rho*scale_d,Told,Told,igrp)/scale_kappa
!!$                       wdtE = C_cal*dt_imp*planck_ana(rho*scale_d,Told,Trold,igrp)/scale_kappa
                       residu(i) = residu(i) + vol_loc*wdtB*deriv_radiation_source(Told,igroup)/scale_E0/(cv+lhs) * P_cal*wdtE*unew(ind_cell(i),firstindex_er+igrp)
                    end if
                 enddo

!!$                 residu(i) = residu(i) &
!!$                      &       - (1.0d0+Prdivu+wdt*(1.0d0-deriv_radiation_source(Told,igroup)*P_cal*wdt/scale_E0/(cv+lhs))) *unew(ind_cell(i),firstindex_er+igroup) *vol_loc
                 Trold=cal_Teg(uold(ind_cell(i),firstindex_er+igroup)*scale_E0,igroup)
                 wdtE = C_cal*dt_imp*planck_ana(rho*scale_d,Told,Told,igroup)/scale_kappa
!!$                 wdtE = C_cal*dt_imp*planck_ana(rho*scale_d,Told,Trold,igroup)/scale_kappa
                 residu(i) = residu(i) &
                      &       - (1.0d0+Prdivu+wdtE-wdtB*deriv_radiation_source(Told,igroup)*P_cal*wdtE/scale_E0/(cv+lhs)) *unew(ind_cell(i),firstindex_er+igroup) *vol_loc

                 !compute b
                 residu(i) = residu(i)+(1.0d0-robin)*rad_flux(ind_cell(i),igroup)
              end if
           end do


        case (2)
           ! residu = Ix
           do i=1,ngrid
              if(son(ind_cell(i)) == 0 )then
                 rho  = uold(ind_cell(i),1)
                 Cv = unew(ind_cell(i),nvar+1)
                 Told = unew(ind_cell(i),5)

                 lhs=0.0d0
                 do igrp=1,ngrp
                    wdtB = C_cal*dt_imp*planck_ana(rho*scale_d,Told,Told,igrp)/scale_kappa
                    lhs=lhs+P_cal*wdtB*deriv_radiation_source(Told,igrp)/scale_E0
                 enddo

                 Trold=cal_Teg(uold(ind_cell(i),firstindex_er+igroup)*scale_E0,igroup)
                 wdtB = C_cal*dt_imp*planck_ana(rho*scale_d,Told,Told ,igroup)/scale_kappa
                 wdtE = C_cal*dt_imp*planck_ana(rho*scale_d,Told,Told,igroup)/scale_kappa
!!$                 wdtE = C_cal*dt_imp*planck_ana(rho*scale_d,Told,Trold,igroup)/scale_kappa
                 Prdivu = 0.0d0!unew(ind_cell(i),nvar+2)*dt_imp
                 residu(i) =  (1.0d0+Prdivu+wdtE-wdtB*deriv_radiation_source(Told,igroup)*P_cal*wdtE/scale_E0/(cv+lhs)) *unew(ind_cell(i),2) *vol_loc
              end if
           end do

        case (3)
           ! residu = I
           do i=1,ngrid
              if(son(ind_cell(i)) == 0 )then

                 rho  = uold(ind_cell(i),1)
                 Cv = unew(ind_cell(i),nvar+1)
                 Told = unew(ind_cell(i),5)

                 lhs=0.0d0
                 do igrp=1,ngrp
                    wdtB = C_cal*dt_imp*planck_ana(rho*scale_d,Told,Told,igrp)/scale_kappa
                    lhs=lhs+P_cal*wdtB*deriv_radiation_source(Told,igrp)/scale_E0
                 enddo

                 wdtB = C_cal*dt_imp*planck_ana(rho*scale_d,Told,Told ,igroup)/scale_kappa
                 wdtE = C_cal*dt_imp*planck_ana(rho*scale_d,Told,Told,igroup)/scale_kappa
!!$                 wdtE = C_cal*dt_imp*planck_ana(rho*scale_d,Told,Trold,igroup)/scale_kappa
                 Prdivu = 0.0d0!unew(ind_cell(i),nvar+2)*dt_imp
                 residu(i) =  (1.0d0+Prdivu+wdtE-wdtB*deriv_radiation_source(Told,igroup)*P_cal*wdtE/scale_E0/(cv+lhs)) *vol_loc
              end if
           end do

        case (4)
           ! reinitialize rad_flux for this level
           do i=1,ngrid
              if(son(ind_cell(i)) == 0 )then
                 rad_flux(ind_cell(i),igroup) = 0.0d0
              end if
           end do
        end select
        
        ! Determine the two2ndim and the direction of the grid of neighboors (-1,0,1)
        do idim = 1,ndim
           if (modulo((ind-1)/2**(idim-1),2)==0)then
              supG = (idim-1)*2+1               !direction of left nbor grid
	      supD = 0              		!direction of right nbor grid
              sub = ind + 2**(idim-1)           ! position of nbor in its own grid
	   else
 	      supG = 0              		!direction of left nbor grid
	      supD = (idim-1)*2+2   		!direction of right nbor grid
              sub = ind - 2**(idim-1)           !position of nbor in its own grid
	   end if

           sub = ncoarse + (sub-1)*ngridmax     !nbor indice offset from its own grid

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


	do idim = 1,ndim

           select case (compute)! Getting val_g and val_d
	   case(1)
              do i=1,ngrid
                 if(son(ind_cell(i)) == 0 )then

                    val_g	(i)       = enew(cell_left (i,idim))
                    val_d	(i)       = enew(cell_right(i,idim))

                 end if
              end do

	   case(2)
              do i=1,ngrid
                 if(son(ind_cell(i)) == 0 )then

                    val_g	(i)       = unew(cell_left (i,idim),2)
                    val_d	(i)       = unew(cell_right(i,idim),2)	
                 end if
              end do

	   case(4)

              do i=1,ngrid
                 if(son(ind_cell(i)) == 0 )then


                    val_g	(i)       = unew(cell_left (i,idim),firstindex_er+igroup)
                    val_d	(i)       = unew(cell_right(i,idim),firstindex_er+igroup)

                 end if
              end do

           end select



	   do i=1,ngrid			! Gather main characteristics of left neighboor
              if(son(ind_cell(i)) == 0 )then

                 select case (nbor_ilevel(i,2*idim-1))

                 case (1)
                    phi_g (i)       = max(uold(cell_left(i,idim),firstindex_er+igroup)/P_cal,eray_min/scale_E0)
                    if (compute==2)		  val_g	(i)       = 0.0d0
                    if (compute/=2)		  val_g	(i)       = phi_g(i)!uold(cell_left(i,idim),firstindex_er+igroup) / P_cal
                    nu_g	(i)       = divu(ind_cell(i))
                    if (robin > 0.0d0)		C_g	(i)       = 1.0d0/robin
                    if (robin == 0.0d0)		  C_g	(i)       = 0.0d0

                    Told      	= cmp_temp(cell_left(i,idim))
                    rho      	= scale_d * max(uold(cell_left(i,idim),1),smallr)
                    nu_g  (i)	= rosseland_ana(rho,Told,Told,igroup) / scale_kappa
                    if(nu_g(i)*dx_loc .lt. min_optical_depth) nu_g(i)=min_optical_depth/dx_loc

                 case (0)

                    phi_g (i)       = enew(cell_left(i,idim))
                    val_g	(i)       = val_g (i)
                    nu_g	(i)       = divu(cell_left(i,idim))
                    C_g	(i)       = 1.0d0

                 case (-1)

                    phi_g (i) 	= max(uold(cell_left(i,idim),firstindex_er+igroup)/P_cal,eray_min/scale_E0)
                    Told      	= cmp_temp(cell_left(i,idim))
                    rho      	= scale_d * max(uold(cell_left(i,idim),1),smallr)

                    if (compute==2)		  val_g (i) 	= 0.0d0
                    if (compute/=2) 	  val_g (i)     = phi_g(i)
                    nu_g  (i)	= rosseland_ana(rho,Told,Told,igroup) / scale_kappa
                    if(nu_g(i)*2.0d0*dx_loc .lt. min_optical_depth) nu_g(i)=min_optical_depth/(2.0d0*dx_loc)         
                    C_g	(i) 	= 1.5d0

                 end select
              end if
	   end do




	   do i=1,ngrid				! Gather main characteristics of right neighboor
              if(son(ind_cell(i)) == 0 )then

                 select case (nbor_ilevel(i,2*idim))
                 case (1)
                    phi_d (i)       = max(uold(cell_right(i,idim),firstindex_er+igroup)/P_cal,eray_min/scale_E0)
                    if (compute==2)		  val_d	(i)       = 0.0d0
                    if (compute/=2)		  val_d	(i)       = phi_d(i)!uold(cell_right(i,idim),firstindex_er+igroup) / P_cal
                    nu_d	(i)       = divu(ind_cell(i))
                    if (robin > 0.0d0)		C_d	(i)       = 1.0d0/robin
                    if (robin == 0.0d0)		  C_d	(i)       = 0.0d0

                    Told 		= cmp_temp(cell_right(i,idim))
                    rho 		= scale_d * max(uold(cell_right(i,idim),1),smallr)
                    nu_d  (i) 	= rosseland_ana(rho,Told,Told,igroup) / scale_kappa
                    if(nu_d(i)*1.0d0*dx_loc .lt. min_optical_depth) nu_d(i)=min_optical_depth/(1.0d0*dx_loc)

                 case (0)

                    phi_d (i)       = enew(cell_right(i,idim))
                    val_d	(i)       = val_d (i)
                    nu_d  (i)       = divu(cell_right(i,idim))
                    C_d   (i)       = 1.0d0

                 case (-1)

                    phi_d (i)	= max(uold(cell_right(i,idim),firstindex_er+igroup)/P_cal,eray_min/scale_E0)
                    Told 		= cmp_temp(cell_right(i,idim))
                    rho 		= scale_d * max(uold(cell_right(i,idim),1),smallr)

                    if (compute==2)		  val_d (i) 	= 0.0d0
                    if (compute/=2)		  val_d (i)     = phi_d(i)
                    nu_d  (i) 	= rosseland_ana(rho,Told,Told,igroup) / scale_kappa
                    C_d   (i)	= 1.5d0
                    if(nu_d(i)*2.0d0*dx_loc .lt. min_optical_depth) nu_d(i)=min_optical_depth/(2.0d0*dx_loc)
           
                 end select
              end if
	   end do

	   do i=1,ngrid
              if(son(ind_cell(i)) == 0 )then ! getting nu_surface
                 nu_c (i) = divu(ind_cell(i))

                 C_g(i)           = C_g(i) * nu_surf(nu_g(i),nu_c(i), cell_left(i,idim) ,ind_cell(i),dx_loc)
                 C_d(i)           = C_d(i) * nu_surf(nu_d(i),nu_c(i), cell_right(i,idim),ind_cell(i),dx_loc)

              end if
	   end do

	   do i=1,ngrid
              if(son(ind_cell(i)) == 0 )then
                 phi_c(i) = enew(ind_cell(i))


                 if(C_g(i) > 0.0d0)then

                    R = max(1.0d-10,abs (phi_c(i)-phi_g(i)) /(0.5d0*(phi_c(i)+phi_g(i))))
                    R = R / ( C_g(i) * dx_loc )

                    lambda=lambda_fld(R)
                    C_g(i) = C_cal*lambda *dt_imp*surf_loc/(dx_loc*C_g(i))
                 end if

                 if(C_d(i) > 0.0d0)then
                    R = max(1.0d-10,abs (phi_c(i)-phi_d(i)) /(0.5d0*(phi_c(i)+phi_d(i))))
                    R = R / ( C_d(i) * dx_loc )

                    lambda=lambda_fld(R)
                    C_d(i) = C_cal*lambda *dt_imp*surf_loc/(dx_loc*C_d(i))

                 end if

              end if
	   end do

           if (compute ==4)then		! Computing and saving flux to the coarser ilevel

              do i=1,ngrid
                 if(son(ind_cell(i)) == 0)then

                    if( nbor_ilevel(i,2*idim-1) == -1)then

                       rad_flux(cell_left(i,idim),igroup)  = rad_flux(cell_left(i,idim),igroup)  + &
                            & C_g(i) *( alpha_imp * (unew(ind_cell(i),firstindex_er+igroup) - val_g(i)) + (1.0d0-alpha_imp)*(phi_c(i) - phi_g(i))) 
                    end if

                    if( nbor_ilevel(i,2*idim)   == -1 )then

                       rad_flux(cell_right(i,idim),igroup) = rad_flux(cell_right(i,idim),igroup) + &
                            & C_d(i) *( alpha_imp * (unew(ind_cell(i),firstindex_er+igroup) - val_d(i)) + (1.0d0-alpha_imp)*(phi_c(i) - phi_d(i)))
                    end if

                 end if
              end do
           end if


	   select case (compute)

	   case (1)		 ! compute b-Ax from b-Ix by adding intern flux
              do i=1,ngrid
                 if(son(ind_cell(i)) == 0 )then
                    residu(i) = residu(i) - ((C_g(i)+C_d(i))* enew(ind_cell(i)) - C_g(i)*val_g(i) - C_d(i)*val_d(i))
                 end if
              end do

	   case (2)		! compute Ap from Ip by adding intern flux
              do i=1,ngrid
                 if(son(ind_cell(i)) == 0 )then
                    residu(i) = residu(i) + ((C_g(i)+C_d(i))* unew(ind_cell(i),2) - C_g(i)*val_g(i) - C_d(i)*val_d(i))*alpha_imp
                 end if
              end do

	   case (3)		! compute Diag(A) for preconditionner
              do i=1,ngrid
                 if(son(ind_cell(i)) == 0 )then
                    residu(i) = residu(i) + (C_g(i) + C_d(i))*alpha_imp
                 end if
              end do

	   end select

	end do !ndim


	select case (compute)
    ! get the result out

	case (1)
           do i=1,ngrid
              if(son(ind_cell(i)) == 0 )then
                 unew(ind_cell(i),1) = residu(i)
                 unew(ind_cell(i),2) = residu(i)
              end if
           end do

	case (2)
           do i=1,ngrid
              if(son(ind_cell(i)) == 0 )then
                 unew(ind_cell(i),3) = residu(i)
              end if
           end do

	case (3)
           do i=1,ngrid
              if(son(ind_cell(i)) == 0 )then
                 unew(ind_cell(i),4) = 1.0d0/residu(i)
              end if
           end do

	end select

     end do ! twotodim
  end do	  ! ncache

end subroutine cmp_matrixA

!################################################################
!################################################################
!################################################################ 
!################################################################
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
  use amr_commons
  use hydro_commons
  use cooling_module,ONLY:kB,mH,clight
  use radiation_parameters
  use units_commons
  implicit none
  integer,intent(in) :: Etype ! Etype=1 : beginning ; Etype=2 : end
  integer ::i,idim,this,mvar,igroup
  real(dp)::usquare,Cv,eps,ekin,emag,rho,erad_loc
  real(dp)::Tp_loc

  ! EOS
  real(dp) :: dd,ee,cmp_temp

  do i=1,nb_ind
     this = liste_ind(i)
     rho   = uold(this,1)

     ! Compute total kinetic energy
     usquare=0.0d0
     do idim=1,ndim
        usquare=usquare+(uold(this,idim+1)/uold(this,1))**2
     end do
     ekin  = rho*usquare/2.0d0

     ! Compute total magnetic energy
     emag = 0.0d0
     do mvar=1,3
        emag = emag + 0.125d0*(uold(this,5+mvar)+uold(this,nvar+mvar))**2
     end do

     if(Etype==1)then

        ! Compute total radiative energy
        erad_loc = 0.0d0
        do igroup=1,nener
           erad_loc = erad_loc + uold(this,8+igroup)
        enddo
        
        dd=rho*scale_d

        eps = (uold(this,5)-ekin-emag-erad_loc)
        if(energy_fix)eps = (uold(this,nvar)) !neil : comment this for radiative shock


        Tp_loc = cmp_temp(this)
        
        Cv = eps/Tp_loc
        
        unew(this,nvar+3) = Tp_loc
        unew(this,5)      = unew(this,nvar+3)
        unew(this,nvar+1) = Cv

        do igroup=1,ngrp
           uold(this,firstindex_er+igroup)   = uold(this,firstindex_er+igroup)/P_cal 
           uold(this,firstindex_er+igroup) = max(uold(this,firstindex_er+igroup),eray_min/scale_E0)
           unew(this,firstindex_er+igroup)   = uold(this,firstindex_er+igroup)
        enddo

     elseif(Etype==2)then
        Cv = unew(this,nvar+1) 
        eps = Cv * unew(this,5)

        uold(this,5)    = eps + ekin + emag
        uold(this,nvar) = eps

        do igroup=1,ngrp
           unew(this,firstindex_er+igroup) = max(unew(this,firstindex_er+igroup),eray_min/scale_E0)
           uold(this,firstindex_er+igroup) = unew(this,firstindex_er+igroup)*P_cal
        end do
        do igroup=1,nener
           uold(this,5)    =  uold(this,5) + uold(this,8+igroup)
        enddo

     end if
  end do

end subroutine cmp_energy

!################################################################
!################################################################
!################################################################ 
!################################################################
function cmp_temp(this)
  use amr_commons
  use hydro_commons
  use cooling_module,ONLY:kB,mH,clight
  use radiation_parameters
  use units_commons
  implicit none
  integer,intent(in) ::this
  integer ::idim,mvar,igrp
  real(dp)::usquare,Cv,eps,ekin,emag,rho,erad_loc
  real(dp)::cmp_temp

  ! EOS
  real(dp) :: dd,ee
  integer  :: ht
  
  real(dp)::sum_dust
#if NDUST>0  
  integer::idust
#endif
  rho   = uold(this,1)
!!$  Cv    = rho*kB/(mu_gas*mH*(gamma-1.0d0))/scale_v**2

  ! Compute total kinetic energy
  usquare=0.0d0
  do idim=1,ndim
     usquare=usquare+(uold(this,idim+1)/uold(this,1))**2
  end do
  ekin  = rho*usquare/2.0d0

  ! Compute total magnetic energy
  emag = 0.0d0
  do mvar=1,3
     emag = emag + 0.125d0*(uold(this,5+mvar)+uold(this,nvar+mvar))**2
  end do

  ! Compute total radiative energy
  erad_loc  = 0.0d0
  do igrp=1,nener
     erad_loc = erad_loc + uold(this,8+igrp)
  enddo
  eps = (uold(this,5)-ekin-emag-erad_loc)
  if(energy_fix)eps = (uold(this,nvar)) !neil : comment this for radiative shock
     
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
  real(dp)::nu_surf,nu_harmo,nu_ari,dl,dr

  nu_ari=(Er2+Er1)*half

  nu_harmo = Er2*Er1/nu_ari
  if(nu_harmo < four/three*dx) nu_harmo = four/three*dx  
  nu_harmo=max(Er2*Er1/nu_ari,four/(three*dx))
  nu_surf = nu_ari

  nu_surf=min(nu_harmo,nu_ari)

  return 
end function nu_surf
!################################################################
!################################################################
!################################################################ 
!################################################################

subroutine ind_bound(b_type,inbor,ind_ref)
  implicit none
  integer,intent(IN)::b_type
  integer,intent(OUT)::inbor
  integer,dimension(1:8),intent(OUT)::ind_ref
  integer::D,num,sens,dir

  dir = modulo(b_type,10)
  D = (dir+1)/2
  sens = modulo(dir,2)
  ! Compute direction of reference neighbors
  inbor = dir + 2*sens - 1
  ! Compute index of reference cells
  do num = 1,8
     if (modulo((num-1)/2**(D-1),2)==0)then 
        if(b_type>10)then
           ind_ref(num) = num + (1-sens) * 2**(D-1)
        else
           ind_ref(num) = num + 2**(D-1)
        end if
     else
        if(b_type>10) then
           ind_ref(num) = num - 2**(D-1) + (1-sens) * 2**(D-1)
        else
           ind_ref(num) = num - 2**(D-1) 
        end if
     end if
  end do

end subroutine ind_bound
!################################################################
!################################################################
!################################################################
!################################################################
subroutine cmp_Prdivu(ilevel,igroup)
  !------------------------------------------------------------------
  ! This routine computes the radiative pressure work and store it in unew(nvar+2)
  !------------------------------------------------------------------
  use amr_commons
  use hydro_commons
  use cooling_module,ONLY:kB,mH,clight
  use radiation_parameters
  use units_commons
  implicit none

  integer,intent(IN)::ilevel,igroup

  integer , dimension(1:nvector,1:2*ndim),save:: nbor_ilevel
  integer , dimension(1:nvector,1:ndim),save::   cell_left , cell_right , big_left, big_right
  integer ,dimension(1:nvector,0:2*ndim),save::  igridn
  integer ,dimension(1:nvector),save ::          ind_cell , ind_grid

  real(dp),dimension(1:nvector),save:: residu,nu_c
  real(dp),dimension(1:nvector),save:: phi_c
  real(dp),dimension(1:nvector,1:ndim,1:ndim),save:: vel_g,vel_d
  real(dp),dimension(1:nvector,1:ndim,1:ngrp),save:: Er_g,Er_d
  real(dp),dimension(1:nvector,1:ndim),save:: dx_g,dx_d

  integer :: i,j,k,idim,ind,igrid,ngrid,ncache,iskip,igrp,nx_loc
  integer :: supG,sub,supD

  real(dp):: d_loc,Tp_loc,dx,dx_loc,surf_loc,vol_loc,scale
  real(dp):: rosseland_ana
  real(dp):: kappa_R,gradEr_norm,gradEr_norm2,R,lambda,lambda_fld,chi

  real(dp) ,dimension(1:ndim,1:ngrp)::gradEr
  real(dp) ,dimension(1:ndim,1:ndim)::divu_loc
  real(dp) ,dimension(1:ndim,1:ndim,1:ngrp)::Pg
  real(dp) :: nuPrDivu,nuPr,nuPl,Pr_nu,Pgdivu
  real(dp), dimension(1:5) :: Pr_temp

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel

  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  surf_loc = dx_loc**(ndim-1)
  vol_loc  = dx_loc**ndim

  ! **************************** LOOP OVER CELLS ********************************** !

  residu = 0.0d0

  ! Loop over myid grids by vector sweeps
  ncache = active(ilevel)%ngrid
  do igrid=1,ncache,nvector

     ! Gather nvector grids
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i) = active(ilevel)%igrid(igrid+i-1)
     end do
     do i=1,ngrid
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

           sub = ncoarse + (sub-1)*ngridmax     !nbor indice offset from its own grid

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


	do idim = 1,ndim
	   do i=1,ngrid			! Gather main characteristics of left neighboor
              if(son(ind_cell(i)) == 0 )then

                 do igrp=1,ngrp
                    select case (nbor_ilevel(i,2*idim-1))
                    case (1)
                       Er_g (i,idim,igrp) = max(uold(cell_left(i,idim),firstindex_er+igrp)/P_cal,eray_min/scale_E0)
                       dx_g (i,idim) = dx_loc
                    case (0)
                       Er_g (i,idim,igrp) = unew(cell_left(i,idim),firstindex_er+igrp)
!!$                       Er_g (i,idim,igrp) = urad(cell_left(i,idim),igrp)
                       dx_g (i,idim) = dx_loc
                    case (-1)
                       Er_g (i,idim,igrp) = max(uold(cell_left(i,idim),firstindex_er+igrp)/P_cal,eray_min/scale_E0)
                       dx_g (i,idim) = 1.5d0*dx_loc
                    end select
                    
                    select case (nbor_ilevel(i,2*idim))
                    case (1)
                       Er_d (i,idim,igrp) = max(uold(cell_right(i,idim),firstindex_er+igrp)/P_cal,eray_min/scale_E0)
                       dx_d (i,idim) = dx_loc
                    case (0)
                       Er_d (i,idim,igrp) = unew(cell_right(i,idim),firstindex_er+igrp)
!!$                       Er_d (i,idim,igrp) = urad(cell_right(i,idim),igrp)
                       dx_d (i,idim) = dx_loc
                    case (-1)
                       
                       Er_d (i,idim,igrp) = max(uold(cell_right(i,idim),firstindex_er+igrp)/P_cal,eray_min/scale_E0)
                       dx_g (i,idim) = 1.5d0*dx_loc
                    end select
                 end do
              
                 vel_g(i,idim,1:ndim)=uold(cell_left(i,idim),2:ndim+1)/uold(cell_left(i,idim),1)
                 vel_d(i,idim,1:ndim)=uold(cell_right(i,idim),2:ndim+1)/uold(cell_right(i,idim),1)
              end if
              end do
        end do

        do i=1,ngrid
           if(son(ind_cell(i)) == 0 )then
              phi_c(i) = enew(ind_cell(i))
              nu_c (i) = divu(ind_cell(i))

              !compute divu
              do j=1,ndim
                 do k=1,ndim
                    divu_loc(j,k) = (vel_d(i,j,k)-vel_g(i,j,k))/(dx_g(i,j)+dx_d(i,j))
                 enddo
                 do igrp=1,ngrp
                    gradEr(j,igrp) = (Er_d(i,j,igrp)-Er_g(i,j,igrp))/(dx_g(i,j)+dx_d(i,j))
                 enddo
              enddo

              do igrp=1,ngrp
                 Tp_loc = unew(ind_cell(i),5)
                 d_loc  = uold(ind_cell(i),1)*scale_d

                 gradEr_norm2 = (sum(gradEr(1:ndim,igrp)**2))
                 gradEr_norm  = (gradEr_norm2)**0.5
                 kappa_R=rosseland_ana(d_loc,Tp_loc,Tp_loc,igrp)/scale_kappa
                 R =  gradEr_norm/(unew(ind_cell(i),firstindex_er+igrp)*kappa_R)
!!$                 R =  gradEr_norm/(urad(ind_cell(i),igrp)*kappa_R)
                 lambda = lambda_fld(R)
                 chi = lambda + (lambda*R)**2
                 
                 do j=1,ndim
                    do k=1,ndim
                       Pg(j,k,igrp)=0.0d0
                       if(j .eq. k) Pg(j,k,igrp) = (1.0d0-chi)/2.0d0
                       if(R .gt. 1.d-8)Pg(j,k,igrp) = Pg(j,k,igrp) &
                            & + (3.0d0*chi-1.0d0)/2.0d0*gradEr(j,igrp)*gradEr(k,igrp)/gradEr_norm2
                    enddo
                 enddo
                 Pg(1:ndim,1:ndim,igrp)=Pg(1:ndim,1:ndim,igrp)*unew(ind_cell(i),firstindex_er+igrp)
!!$                 Pg(1:ndim,1:ndim,igrp)=Pg(1:ndim,1:ndim,igrp)*urad(ind_cell(i),igrp)
                 
              end do
                 
              Pgdivu     = 0.0d0
              nuPrDivu   = 0.0d0 ! multig: this is the Doppler shift term
              
              do j=1,ndim
                 do k=1,ndim
                    
                    ! compure Pr:divU term
                    Pgdivu    = Pgdivu    + Pg(j,k,igroup)*divu_loc(j,k)

                    ! fill in Pr_temp array for Doppler shift terms
                    Pr_temp(3) = Pg(j,k,igroup)
                    if(igroup >    1  )then
                       Pr_temp(2) = Pg(j,k,igroup-1)
                    else
                       Pr_temp(2) = Pg(j,k,igroup)
                    endif
                    if(igroup >    2  )then
                       Pr_temp(1) = Pg(j,k,igroup-2)
                    else
                       Pr_temp(1) = Pr_temp(2)
                    endif
                    if(igroup < ngrp  )then
                       Pr_temp(4) = Pg(j,k,igroup+1)
                    else
                       Pr_temp(4) = Pg(j,k,igroup)
                    endif
                    if(igroup < ngrp-1)then
                       Pr_temp(5) = Pg(j,k,igroup+2)
                    else
                       Pr_temp(5) = Pr_temp(4)
                    endif

                    ! compute -[nu P]*Div(u) term
                    if(divu_loc(j,k) > 0.0d0)then

                       if(igroup == ngrp)then
                          nuPr = 0.0d0
                       else
                          nuPr = nu_max_hz(igroup)*Pr_nu(Pr_temp,igroup+1,4,nu_max_hz(igroup))
                       endif
                       if(igroup == 1)then
                          nuPl = 0.0d0
                       else
                          nuPl = nu_min_hz(igroup)*Pr_nu(Pr_temp,igroup,3,nu_min_hz(igroup))
                       endif

                    else

                       if(igroup == ngrp)then
                          nuPr = 0.0d0
                       else
                          nuPr = nu_max_hz(igroup)*Pr_nu(Pr_temp,igroup,3,nu_max_hz(igroup))
                       endif
                       if(igroup == 1)then
                          nuPl = 0.0d0
                       else
                          nuPl = nu_min_hz(igroup)*Pr_nu(Pr_temp,igroup-1,2,nu_min_hz(igroup))
                       endif

                    endif

                    nuPrDivu = nuPrDivu - (nuPr - nuPl)*divu_loc(j,k)

                 enddo
              end do
              unew(ind_cell(i),nvar+2) = 0.0d0!(Pgdivu+nuPrDivu)/urad(ind_cell(i),igroup)
           end if

        end do


     end do ! twotodim
  end do	  ! ncache

end subroutine cmp_Prdivu
!################################################################
!################################################################
!################################################################ 
!################################################################
subroutine cmp_matrix_coeff(ilevel,igroup)
  !------------------------------------------------------------------
  ! This routine computes the matrix A to vect_in and create vect_out
  ! compute = 1 : residu           	return B - Ax
  ! compute = 2 : Produit                 return  A.p
  ! compute = 3 : Preconditionner         return diag(A)
  ! compute = 4 : Compute flux in rad_flux
  !------------------------------------------------------------------
  use amr_commons
  use hydro_commons
  use cooling_module,ONLY:kB,mH,clight
  use radiation_parameters
  use units_commons
  implicit none

  integer,intent(IN)::ilevel,igroup

  integer , dimension(1:nvector,1:2*ndim),save:: nbor_ilevel
  integer , dimension(1:nvector,1:ndim),save::   cell_left , cell_right , big_left, big_right
  integer ,dimension(1:nvector,0:2*ndim),save::  igridn
  integer ,dimension(1:nvector),save ::          ind_cell , ind_grid

  real(dp),dimension(1:nvector),save:: residu,C_g,C_d,nu_g,nu_c,nu_d
  real(dp),dimension(1:nvector),save:: phi_g,phi_c,phi_d,val_g,val_d


  integer :: i,idim,ind,igrid,ngrid,ncache,iskip,igrp,nx_loc
  integer :: supG,sub,supD

  real(dp):: radiation_source,deriv_radiation_source,rhs,lhs
  real(dp):: dx,dx_loc,surf_loc,vol_loc,scale
  real(dp):: nu_surf,Cv,rho,wdt,wdt_igrp,Told,lambda,lambda_fld,R
  real(dp):: cmp_temp,rosseland_ana,planck_ana,Prdivu

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel

  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  surf_loc = dx_loc**(ndim-1)
  vol_loc  = dx_loc**ndim

  ! **************************** LOOP OVER CELLS ********************************** !

  residu = 0.0d0

  ! Loop over myid grids by vector sweeps
  ncache = active(ilevel)%ngrid
  do igrid=1,ncache,nvector


     ! Gather nvector grids
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i) = active(ilevel)%igrid(igrid+i-1)
     end do


     do i=1,ngrid
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

        ! residu = I
        do i=1,ngrid
           if(son(ind_cell(i)) == 0 )then
              
              rho  = uold(ind_cell(i),1)
              Cv = unew(ind_cell(i),nvar+1)
              Told = unew(ind_cell(i),5)
              
              lhs=0.0d0
              do igrp=1,ngrp
                 wdt = C_cal*dt_imp*planck_ana(rho*scale_d,Told,Told,igrp)/scale_kappa
                 lhs=lhs+P_cal*wdt*deriv_radiation_source(Told,igrp)/scale_E0
              enddo
              
              wdt = C_cal*dt_imp*planck_ana(rho*scale_d,Told,Told,igroup)/scale_kappa
              Prdivu = 0.0d0!unew(ind_cell(i),nvar+2)*dt_imp
              mat_residual_glob(ind_cell(i),igroup,igroup) =  (1.0d0+Prdivu+wdt*(1.0d0-deriv_radiation_source(Told,igroup)*P_cal*wdt/scale_E0/(cv+lhs))) *vol_loc
           end if
        end do
        
        ! Determine the two2ndim and the direction of the grid of neighboors (-1,0,1)
        do idim = 1,ndim
           if (modulo((ind-1)/2**(idim-1),2)==0)then
              supG = (idim-1)*2+1               !direction of left nbor grid
	      supD = 0              		!direction of right nbor grid
              sub = ind + 2**(idim-1)           ! position of nbor in its own grid
           else
              supG = 0              		!direction of left nbor grid
	      supD = (idim-1)*2+2   		!direction of right nbor grid
              sub = ind - 2**(idim-1)           !position of nbor in its own grid
           end if
    
           sub = ncoarse + (sub-1)*ngridmax     !nbor indice offset from its own grid

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


        do idim = 1,ndim

!!$           do i=1,ngrid
!!$              if(son(ind_cell(i)) == 0 )then
!!$                 
!!$                 val_g	(i)       = enew(cell_left (i,idim))
!!$                 val_d	(i)       = enew(cell_right(i,idim))
!!$
!!$              end if
!!$           end do
              

	   do i=1,ngrid			! Gather main characteristics of left neighboor
              if(son(ind_cell(i)) == 0 )then

                 select case (nbor_ilevel(i,2*idim-1))

                 case (1)
                    phi_g (i)       = max(uold(cell_left(i,idim),firstindex_er+igroup)/P_cal,eray_min/scale_E0)
!!$                    if (compute==2)		  val_g	(i)       = 0.0d0
!!$                    if (compute/=2)		  val_g	(i)       = phi_g(i)!uold(cell_left(i,idim),firstindex_er+igroup) / P_cal
                    nu_g	(i)       = divu(ind_cell(i))
                    if (robin > 0.0d0)		C_g	(i)       = 1.0d0/robin
                    if (robin == 0.0d0)		  C_g	(i)       = 0.0d0

                    Told      	= cmp_temp(cell_left(i,idim))
                    rho      	= scale_d * max(uold(cell_left(i,idim),1),smallr)
                    nu_g  (i)	= rosseland_ana(rho,Told,Told,igroup) / scale_kappa
                    if(nu_g(i)*dx_loc .lt. min_optical_depth) nu_g(i)=min_optical_depth/dx_loc
                 case (0)

                    phi_g (i)       = enew(cell_left(i,idim))
!!$                    val_g	(i)       = val_g (i)
                    nu_g	(i)       = divu(cell_left(i,idim))
                    C_g	(i)       = 1.0d0

                 case (-1)

                    phi_g (i) 	= max(uold(cell_left(i,idim),firstindex_er+igroup)/P_cal,eray_min/scale_E0)
                    Told      	= cmp_temp(cell_left(i,idim))
                    rho      	= scale_d * max(uold(cell_left(i,idim),1),smallr)

!!$                    if (compute==2)		  val_g (i) 	= 0.0d0
!!$                    if (compute/=2) 	  val_g (i)     = phi_g(i)
                    nu_g  (i)	= rosseland_ana(rho,Told,Told,igroup) / scale_kappa
                    if(nu_g(i)*2.0d0*dx_loc .lt. min_optical_depth) nu_g(i)=min_optical_depth/(2.0d0*dx_loc)         
                    C_g	(i) 	= 1.5d0

                 end select
              end if
	   end do




	   do i=1,ngrid				! Gather main characteristics of right neighboor
              if(son(ind_cell(i)) == 0 )then

                 select case (nbor_ilevel(i,2*idim))
                 case (1)
                    phi_d (i)       = max(uold(cell_right(i,idim),firstindex_er+igroup)/P_cal,eray_min/scale_E0)
!!$                    if (compute==2)		  val_d	(i)       = 0.0d0
!!$                    if (compute/=2)		  val_d	(i)       = phi_d(i)!uold(cell_right(i,idim),firstindex_er+igroup) / P_cal
                    nu_d	(i)       = divu(ind_cell(i))
                    if (robin > 0.0d0)		C_d	(i)       = 1.0d0/robin
                    if (robin == 0.0d0)		  C_d	(i)       = 0.0d0

                    Told 		= cmp_temp(cell_right(i,idim))
                    rho 		= scale_d * max(uold(cell_right(i,idim),1),smallr)
                    nu_d  (i) 	= rosseland_ana(rho,Told,Told,igroup) / scale_kappa
                    if(nu_d(i)*1.0d0*dx_loc .lt. min_optical_depth) nu_d(i)=min_optical_depth/(1.0d0*dx_loc)

                 case (0)

                    phi_d (i)       = enew(cell_right(i,idim))
!!$                    val_d	(i)       = val_d (i)
                    nu_d  (i)       = divu(cell_right(i,idim))
                    C_d   (i)       = 1.0d0

                 case (-1)

                    phi_d (i)	= max(uold(cell_right(i,idim),firstindex_er+igroup)/P_cal,eray_min/scale_E0)
                    Told 		= cmp_temp(cell_right(i,idim))
                    rho 		= scale_d * max(uold(cell_right(i,idim),1),smallr)

!!$                    if (compute==2)		  val_d (i) 	= 0.0d0
!!$                    if (compute/=2)		  val_d (i)     = phi_d(i)
                    nu_d  (i) 	= rosseland_ana(rho,Told,Told,igroup) / scale_kappa
                    C_d   (i)	= 1.5d0
                    if(nu_d(i)*2.0d0*dx_loc .lt. min_optical_depth) nu_d(i)=min_optical_depth/(2.0d0*dx_loc)
           
                 end select
              end if
	   end do

	   do i=1,ngrid
              if(son(ind_cell(i)) == 0 )then ! getting nu_surface
                 nu_c (i) = divu(ind_cell(i))

                 C_g(i)           = C_g(i) * nu_surf(nu_g(i),nu_c(i), cell_left(i,idim) ,ind_cell(i),dx_loc)
                 C_d(i)           = C_d(i) * nu_surf(nu_d(i),nu_c(i), cell_right(i,idim),ind_cell(i),dx_loc)

              end if
	   end do

	   do i=1,ngrid
              if(son(ind_cell(i)) == 0 )then
                 phi_c(i) = enew(ind_cell(i))

                 coeff_glob_left(ind_cell(i),igroup,igroup,idim)=0.0d0
                 if(C_g(i) > 0.0d0)then

                    R = max(1.0d-10,abs (phi_c(i)-phi_g(i)) /(0.5d0*(phi_c(i)+phi_g(i))))
                    R = R / ( C_g(i) * dx_loc )

                    lambda=lambda_fld(R)
                    C_g(i) = C_cal*lambda *dt_imp*surf_loc/(dx_loc*C_g(i))

                    coeff_glob_left(ind_cell(i),igroup,igroup,idim)=C_g(i)

                 end if
                 coeff_glob_right(ind_cell(i),igroup,igroup,idim)=0.0d0
                 if(C_d(i) > 0.0d0)then

                    R = max(1.0d-10,abs (phi_c(i)-phi_d(i)) /(0.5d0*(phi_c(i)+phi_d(i))))
                    R = R / ( C_d(i) * dx_loc )

                    lambda=lambda_fld(R)
                    C_d(i) = C_cal*lambda *dt_imp*surf_loc/(dx_loc*C_d(i))
                    coeff_glob_right(ind_cell(i),igroup,igroup,idim)=C_d(i)
                 end if
              end if
	   end do

       end do !ndim
 

     end do ! twotodim
  end do	  ! ncache
end subroutine cmp_matrix_coeff
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmp_matrixA2 (ilevel,compute,igroup)
  !------------------------------------------------------------------
  ! This routine computes the matrix A to vect_in and create vect_out
  ! compute = 1 : residu           	return B - Ax
  ! compute = 2 : Produit                 return  A.p
  ! compute = 3 : Preconditionner         return diag(A)
  ! compute = 4 : Compute flux in rad_flux
  !------------------------------------------------------------------

  use amr_commons
  use hydro_commons
  use cooling_module,ONLY:kB,mH,clight
  use radiation_parameters
  use units_commons
  implicit none

  integer,intent(IN)::compute,ilevel,igroup

  integer , dimension(1:nvector,1:2*ndim),save:: nbor_ilevel
  integer , dimension(1:nvector,1:ndim),save::   cell_left , cell_right , big_left, big_right
  integer ,dimension(1:nvector,0:2*ndim),save::  igridn
  integer ,dimension(1:nvector),save ::          ind_cell , ind_grid

  real(dp),dimension(1:nvector),save:: residu,C_g,C_d,nu_g,nu_c,nu_d
  real(dp),dimension(1:nvector),save:: phi_g,phi_c,phi_d,val_g,val_d


  integer :: i,idim,ind,igrid,ngrid,ncache,iskip,igrp,nx_loc
  integer :: supG,sub,supD

  real(dp):: radiation_source,deriv_radiation_source,rhs,lhs
  real(dp):: dx,dx_loc,surf_loc,vol_loc,scale
  real(dp):: nu_surf,Cv,rho,wdt,wdt_igrp,Told,lambda,lambda_fld,R
  real(dp):: cmp_temp,rosseland_ana,planck_ana,Prdivu

  ! Mesh size at level ilevel
  dx=0.5D0**ilevel

  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  surf_loc = dx_loc**(ndim-1)
  vol_loc  = dx_loc**ndim

  ! **************************** LOOP OVER CELLS ********************************** !

  residu = 0.0d0

  ! Loop over myid grids by vector sweeps
  ncache = active(ilevel)%ngrid
  do igrid=1,ncache,nvector


     ! Gather nvector grids
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i) = active(ilevel)%igrid(igrid+i-1)
     end do


     do i=1,ngrid
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

        select case (compute)

        case (1)
           ! residu = b - Ax
           do i=1,ngrid
              if(son(ind_cell(i)) == 0 )then
                 rho   = uold (ind_cell(i),1)
                 Cv = unew(ind_cell(i),nvar+1)
                 Prdivu = 0.0d0!unew(ind_cell(i),nvar+2)*dt_imp
                 Told  =  unew(ind_cell(i),5)

                 rhs=0.0d0
                 lhs=0.0d0
                 do igrp=1,ngrp
                    wdt = C_cal*dt_imp*planck_ana(rho*scale_d,Told,Told,igrp)/scale_kappa
                    rhs=rhs-P_cal*wdt*(radiation_source(Told,igrp)/scale_E0-Told*deriv_radiation_source(Told,igrp)/scale_E0)
                    lhs=lhs+P_cal*wdt*deriv_radiation_source(Told,igrp)/scale_E0
                 enddo

                 wdt = C_cal*dt_imp*planck_ana(rho*scale_d,Told,Told,igroup)/scale_kappa
                 residu(i) = uold(ind_cell(i),firstindex_er+igroup)*vol_loc  &
                      & + vol_loc*wdt*(radiation_source(Told,igroup)/scale_E0-Told*deriv_radiation_source(Told,igroup)/scale_E0) &
                      & + vol_loc*wdt*deriv_radiation_source(Told,igroup)/scale_E0*(cv*unew(ind_cell(i),nvar+3)+rhs)/(cv+lhs)

                 ! Terms of coupling radiative groups
                 do igrp=1,ngrp
                    if(igrp .ne. igroup) then
                       wdt_igrp = C_cal*dt_imp*planck_ana(rho*scale_d,Told,Told,igrp)/scale_kappa
                       residu(i) = residu(i) + vol_loc*wdt*deriv_radiation_source(Told,igroup)/scale_E0/(cv+lhs) * P_cal*wdt_igrp*unew(ind_cell(i),firstindex_er+igrp)
!!$                       residu(i) = residu(i) + vol_loc*wdt*deriv_radiation_source(Told,igroup)/scale_E0/(cv+lhs) * P_cal*wdt_igrp*urad(ind_cell(i),igrp)
                    end if
                 enddo

                 residu(i) = residu(i) &
                      &       - (1.0d0+Prdivu+wdt*(1.0d0-deriv_radiation_source(Told,igroup)*P_cal*wdt/scale_E0/(cv+lhs))) *unew(ind_cell(i),firstindex_er+igroup) *vol_loc

                 !compute b
                 residu(i) = residu(i)+(1.0d0-robin)*rad_flux(ind_cell(i),igroup)
              end if
           end do


        case (2)
           ! residu = Ix
           do i=1,ngrid
              if(son(ind_cell(i)) == 0 )then
!!$                 rho  = uold(ind_cell(i),1)
!!$                 Cv = unew(ind_cell(i),nvar+1)
!!$                 Told = unew(ind_cell(i),5)
!!$
!!$                 lhs=0.0d0
!!$                 do igrp=1,ngrp
!!$                    wdt = C_cal*dt_imp*planck_ana(rho*scale_d,Told,igrp)/scale_kappa
!!$                    lhs=lhs+P_cal*wdt*deriv_radiation_source(Told,igrp)/scale_E0
!!$                 enddo
!!$
!!$                 wdt = C_cal*dt_imp*planck_ana(rho*scale_d,Told,igroup)/scale_kappa
!!$                 Prdivu = 0.0d0 !unew(ind_cell(i),nvar+2)*dt_imp
                 residu(i) =  mat_residual_glob(ind_Cell(i),igroup,igroup)*unew(ind_cell(i),2)!(1.0d0+Prdivu+wdt*(1.0d0-deriv_radiation_source(Told,igroup)*P_cal*wdt/scale_E0/(cv+lhs))) *unew(ind_cell(i),2) *vol_loc
              end if
           end do

        case (3)
           ! residu = I
           do i=1,ngrid
              if(son(ind_cell(i)) == 0 )then
!!$
!!$                 rho  = uold(ind_cell(i),1)
!!$                 Cv = unew(ind_cell(i),nvar+1)
!!$                 Told = unew(ind_cell(i),5)
!!$
!!$                 lhs=0.0d0
!!$                 do igrp=1,ngrp
!!$                    wdt = C_cal*dt_imp*planck_ana(rho*scale_d,Told,igrp)/scale_kappa
!!$                    lhs=lhs+P_cal*wdt*deriv_radiation_source(Told,igrp)/scale_E0
!!$                 enddo
!!$
!!$                 wdt = C_cal*dt_imp*planck_ana(rho*scale_d,Told,igroup)/scale_kappa
!!$                 Prdivu = unew(ind_cell(i),nvar+2)*dt_imp
!!$                 residu(i) =  (1.0d0+Prdivu+wdt*(1.0d0-deriv_radiation_source(Told,igroup)*P_cal*wdt/scale_E0/(cv+lhs))) *vol_loc
                 residu(i) =  mat_residual_glob(ind_Cell(i),igroup,igroup)
              end if
           end do

        case (4)
           ! reinitialize rad_flux for this level
           do i=1,ngrid
              if(son(ind_cell(i)) == 0 )then
                 rad_flux(ind_cell(i),igroup) = 0.0d0
              end if
           end do
        end select
        
        ! Determine the two2ndim and the direction of the grid of neighboors (-1,0,1)
        do idim = 1,ndim
           if (modulo((ind-1)/2**(idim-1),2)==0)then
              supG = (idim-1)*2+1               !direction of left nbor grid
	      supD = 0              		!direction of right nbor grid
              sub = ind + 2**(idim-1)           ! position of nbor in its own grid
	   else
 	      supG = 0              		!direction of left nbor grid
	      supD = (idim-1)*2+2   		!direction of right nbor grid
              sub = ind - 2**(idim-1)           !position of nbor in its own grid
	   end if

           sub = ncoarse + (sub-1)*ngridmax     !nbor indice offset from its own grid

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


	do idim = 1,ndim

           select case (compute)! Getting val_g and val_d
	   case(1)
              do i=1,ngrid
                 if(son(ind_cell(i)) == 0 )then

                    val_g	(i)       = enew(cell_left (i,idim))
                    val_d	(i)       = enew(cell_right(i,idim))

                 end if
              end do

	   case(2)
              do i=1,ngrid
                 if(son(ind_cell(i)) == 0 )then

                    val_g	(i)       = unew(cell_left (i,idim),2)
                    val_d	(i)       = unew(cell_right(i,idim),2)	
                 end if
              end do

	   case(4)

              do i=1,ngrid
                 if(son(ind_cell(i)) == 0 )then


                    val_g	(i)       = unew(cell_left (i,idim),firstindex_er+igroup)
                    val_d	(i)       = unew(cell_right(i,idim),firstindex_er+igroup)

                 end if
              end do

           end select

	   do i=1,ngrid			! Gather main characteristics of left neighboor
              if(son(ind_cell(i)) == 0 )then

                 select case (nbor_ilevel(i,2*idim-1))

                 case (1)
                    phi_g (i)       = max(uold(cell_left(i,idim),firstindex_er+igroup)/P_cal,eray_min/scale_E0)
                    if (compute==2)		  val_g	(i)       = 0.0d0
                    if (compute/=2)		  val_g	(i)       = phi_g(i)!uold(cell_left(i,idim),firstindex_er+igroup) / P_cal

                 case (0)

                    phi_g (i)       = enew(cell_left(i,idim))
                    val_g	(i)       = val_g (i)

                 case (-1)

                    phi_g (i) 	= max(uold(cell_left(i,idim),firstindex_er+igroup)/P_cal,eray_min/scale_E0)
                    if (compute==2)		  val_g (i) 	= 0.0d0
                    if (compute/=2) 	  val_g (i)     = phi_g(i)
                 end select
              end if
	   end do




	   do i=1,ngrid				! Gather main characteristics of right neighboor
              if(son(ind_cell(i)) == 0 )then

                 select case (nbor_ilevel(i,2*idim))
                 case (1)
                    phi_d (i)       = max(uold(cell_right(i,idim),firstindex_er+igroup)/P_cal,eray_min/scale_E0)
                    if (compute==2)		  val_d	(i)       = 0.0d0
                    if (compute/=2)		  val_d	(i)       = phi_d(i)!uold(cell_right(i,idim),firstindex_er+igroup) / P_cal
         
                 case (0)

                    phi_d (i)       = enew(cell_right(i,idim))
                    val_d	(i)       = val_d (i)

                 case (-1)

                    phi_d (i)	= max(uold(cell_right(i,idim),firstindex_er+igroup)/P_cal,eray_min/scale_E0)
                    if (compute==2)		  val_d (i) 	= 0.0d0
                    if (compute/=2)		  val_d (i)     = phi_d(i)

                 end select
              end if
	   end do

           if (compute ==4)then		! Computing and saving flux to the coarser ilevel

              do i=1,ngrid
                 if(son(ind_cell(i)) == 0)then

                    phi_c(i) = enew(ind_cell(i))

                    if( nbor_ilevel(i,2*idim-1) == -1)then

                       rad_flux(cell_left(i,idim),igroup)  = rad_flux(cell_left(i,idim),igroup)  + &
                            & coeff_glob_left(ind_cell(i),igroup,igroup,idim)*( alpha_imp * (unew(ind_cell(i),firstindex_er+igroup) - val_g(i)) + (1.0d0-alpha_imp)*(phi_c(i) - phi_g(i))) 
                    end if

                    if( nbor_ilevel(i,2*idim)   == -1 )then

                       rad_flux(cell_right(i,idim),igroup) = rad_flux(cell_right(i,idim),igroup) + &
                            & coeff_glob_right(ind_cell(i),igroup,igroup,idim)*( alpha_imp * (unew(ind_cell(i),firstindex_er+igroup) - val_d(i)) + (1.0d0-alpha_imp)*(phi_c(i) - phi_d(i)))
                    end if

                 end if
              end do
           end if


	   select case (compute)

	   case (1)		 ! compute b-Ax from b-Ix by adding intern flux
              do i=1,ngrid
                 if(son(ind_cell(i)) == 0 )then
                    residu(i) = residu(i) &
                         & - ((coeff_glob_left(ind_cell(i),igroup,igroup,idim)+coeff_glob_right(ind_cell(i),igroup,igroup,idim))* enew(ind_cell(i)) &
                         & - coeff_glob_left(ind_cell(i),igroup,igroup,idim)*val_g(i) &
                         & - coeff_glob_right(ind_cell(i),igroup,igroup,idim)*val_d(i))
                 end if
              end do

	   case (2)		! compute Ap from Ip by adding intern flux
              do i=1,ngrid
                 if(son(ind_cell(i)) == 0 )then
                    residu(i) = residu(i) &
                         & + ((coeff_glob_left(ind_cell(i),igroup,igroup,idim)+coeff_glob_right(ind_cell(i),igroup,igroup,idim))*unew(ind_cell(i),2) &
                         & - coeff_glob_left(ind_cell(i),igroup,igroup,idim)*val_g(i) &
                         & - coeff_glob_right(ind_cell(i),igroup,igroup,idim)*val_d(i))*alpha_imp
                 end if
              end do

	   case (3)		! compute Diag(A) for preconditionner
              do i=1,ngrid
                 if(son(ind_cell(i)) == 0 )then
                    residu(i) = residu(i) &
                         & + (coeff_glob_left(ind_cell(i),igroup,igroup,idim)+coeff_glob_right(ind_cell(i),igroup,igroup,idim))*alpha_imp
                 end if
              end do

	   end select

	end do !ndim


	select case (compute)
    ! get the result out

	case (1)
           do i=1,ngrid
              if(son(ind_cell(i)) == 0 )then
                 unew(ind_cell(i),1) = residu(i)
                 unew(ind_cell(i),2) = residu(i)
              end if
           end do

	case (2)
           do i=1,ngrid
              if(son(ind_cell(i)) == 0 )then
                 unew(ind_cell(i),3) = residu(i)
              end if
           end do

	case (3)
           do i=1,ngrid
              if(son(ind_cell(i)) == 0 )then
                 unew(ind_cell(i),4) = 1.0d0/residu(i)
              end if
           end do

	end select

     end do ! twotodim
  end do	  ! ncache

end subroutine cmp_matrixA2
