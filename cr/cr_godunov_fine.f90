!###############################################################
!###############################################################
subroutine crmom_step(ilevel)
  ! Two-moment cosmic-ray transport driver: subcycled free-streaming P1/M1
  ! advection followed by the implicit scattering/streaming source terms.
  !
  ! Either do one transport step on ilevel, with CR field updates in the
  ! coarser-level neighbours, or, if cr_nsubcycle>1, do many substeps in
  ! ilevel only, using Dirichlet boundary conditions for the level
  ! boundaries.
  use amr_commons
  use cr_parameters
  use cr_hydro_commons
  use constants
  use mpi_mod
  implicit none
  integer,intent(in)::ilevel
  real(dp)::dt_hydro,t_left,dt_cr,t_save
  integer::i_substep,ivar
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  dt_hydro=dtnew(ilevel)                     ! Store hydro timestep length
  t_left=dt_hydro
  ! Shift the time backwards one hydro-dt, in case of CR subcycling. Not
  ! strictly necessary as CR advection has no explicit time dependence.
  t_save=t ; t=t-t_left

  ! Initialise crunew=cruold for the first subcycle (zeroes virtual cells).
  call cr_set_unew(ilevel)

  call get_crmom_courant(dt_cr,ilevel)
  i_substep=0
  do while(t_left>0)                         !                CR sub-cycle
     i_substep=i_substep+1
     ! Temporarily change timestep length to the CR step:
     dtnew(ilevel)=MIN(t_left,dt_cr)
     if(i_substep.gt.cr_nsubcycle)then
        print*,'Doing CR substeps but should not!! ',i_substep,cr_nsubcycle,dt_hydro,t_left
        call clean_stop
     endif
     t=t+dtnew(ilevel)                        ! Shift the time forwards one dt_cr

     if(i_substep>1)call cr_set_unew(ilevel)

     call cr_godunov_fine(ilevel)
     ! Pull in updates from finer level on other cpus -- only needed when
     ! updating coarser level (i.e. not subcycling).
     if(cr_nsubcycle==1)then
        do ivar=1,ncrvar
           call make_virtual_reverse_dp(crunew(1,ivar),ilevel)
        end do
     endif

     call add_cr_source_terms(ilevel)
     call cr_set_uold(ilevel)

     ! Collisional CR cooling (Coulomb/hadronic losses), on cruold.
     if(cr_cooling)call cr_cooling_fine(ilevel)

     do ivar=1,ncrvar
        call make_virtual_fine_dp(cruold(1,ivar),ilevel)
     end do

     if(simple_boundary)call cr_make_boundary_hydro(ilevel)

     t_left=t_left-dtnew(ilevel)
  end do                                     !          End CR subcycle loop
  dtnew(ilevel)=dt_hydro                     ! Restore hydro timestep length
  t=t_save        ! Restore original time (otherwise tiny roundoff error)

  call cr_upload_fine(ilevel)

  if(myid==1 .and. mod(nstep_coarse,ncontrol)==0)then
     call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
     write(*,901)ilevel,i_substep,cr_vmax(ilevel)*scale_v/1e5, &
          cr_vmax(ilevel)*scale_v/c_cgs,dt_cr
  endif

111 format('   Entering crmom_step for level ',I2)
901 format(' Performed level',I3,' CR-step with ',I5,' subcycles, cr_vmax(km/s)=', &
         1pe9.2,' cr_c_fraction=',1pe9.2,',  dt_cr=',1pe9.2)

end subroutine crmom_step

!##########################################################################
!##########################################################################
subroutine add_cr_source_terms(ilevel)
  ! Implicit per-cell cosmic-ray scattering/streaming source terms.
  ! It runs inside crmom_step between cr_godunov_fine
  ! (explicit transport, which writes crunew) and cr_set_uold.
  !
  ! For each cell the implicit solve relaxes (E_cr,F_cr) by the scattering
  ! coefficient sigma=1/Dcr_code(iGrp) via a 4x4 coefficient matrix reduced
  ! with a Schur complement on E_cr; the flux is rotated onto the local B
  ! field with rotatevec and back with invrotatevec.
  !---------------------------------------------------------
  use amr_commons
  use hydro_commons
  use cr_hydro_commons
  use cr_parameters
  use cr_flux_module, only: rotatevec, invrotatevec
  implicit none
  integer::ilevel
  !---------------------------------------------------------
  ! This routine adds the cosmic ray source term .
  !---------------------------------------------------------
  integer::i,ind,iskip,nx_loc
  integer::ncache,igrid,ngrid,idim,id1,ig1,ih1,id2,ig2,ih2
  integer,dimension(1:3,1:2,1:8)::iii,jjj
  real(dp)::scale,dx,dx_loc,B_field(3),dt

  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  integer ,dimension(1:nvector,0:twondim),save::igridn
  integer ,dimension(1:nvector,1:ndim),save::ind_left,ind_right
  real(dp),dimension(1:nvector,1:ndim),save::dx_g,dx_d
  real(dp),dimension(1:nvector,1:ndim,1:ncr_groups),save::gradEcr_loc,gradpcr_loc
  real(dp),dimension(1:nvector,1:3),save::B_field_loc
  real(dp),dimension(1:nvector),save::bdotgradE_loc
  real(dp),dimension(1:nvector,1:3),save::vs_loc
  real(dp),dimension(1:nvector),save::va_loc
  real(dp)::norm,frotx,froty,frotz,bxby,cosp,sinp,cost,sint
  real(dp)::f1,f2,f3
  integer::j,iGrp,icrE
  real(dp),dimension(1:nvector,1:ndim,1:ncr_groups),save::pcrg,pcrd
  real(dp)::coef_11, coef_12, coef_13, coef_14, coef_21, coef_22
  real(dp)::coef_31, coef_33, coef_41, coef_44
  real(dp)::e_coef, new_ec, old_ec, sigma_x, sigma_y, sigma_z, sigma_stream
  real(dp)::rhs1, rhs2, rhs3, rhs4, fred, sqrt3
  real(dp)::v1, v2, v3, vtot1, vtot2, vtot3
  real(dp)::mom_change
  real(dp)::etherm, ekin, emag, err, f_decouple, smallp

  dt=dtnew(ilevel)

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  smallp = smallc**2/gamma

  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5d0**ilevel
  dx_loc=dx*scale

  if(isotropic_pressure)then
     sqrt3=sqrt(3d0)
  else
     sqrt3=1.0d0
  endif

  iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

  ! Loop over myid grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector

     ! Gather nvector grids
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do

     ! Gather neighboring grids
     do i=1,ngrid
        igridn(i,0)=ind_grid(i)
     end do
     do idim=1,ndim
        do i=1,ngrid
           ind_left (i,idim)=nbor(ind_grid(i),2*idim-1)
           ind_right(i,idim)=nbor(ind_grid(i),2*idim  )
           igridn(i,2*idim-1)=son(ind_left (i,idim))
           igridn(i,2*idim  )=son(ind_right(i,idim))
        end do
     end do

     ! Loop over cells
     do ind=1,twotondim

        ! Compute central cell index
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do

        ! Gather all neighboring CR energies
        do idim=1,ndim
           id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
           ih1=ncoarse+(id1-1)*ngridmax
           do i=1,ngrid
           do iGrp=1,ncr_groups
              icrE = iCRu+(ndim+1)*(iGrp-1)  ! starting index of cr variables
              if(igridn(i,ig1)>0)then
                  pcrg(i,idim,iGrp)   = max(cruold(igridn(i,ig1)+ih1,icrE),smallcr)
                  dx_g(i,idim) = dx_loc
              else
                  pcrg(i,idim,iGrp)   = max(cruold(ind_left(i,idim),icrE),smallcr)
                  dx_g(i,idim) = dx_loc*1.5_dp
              end if
           enddo
           enddo
           id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
           ih2=ncoarse+(id2-1)*ngridmax
           do i=1,ngrid
           do iGrp=1,ncr_groups
              icrE = iCRu+(ndim+1)*(iGrp-1)  ! starting index of cr variables
              if(igridn(i,ig2)>0)then
                  pcrd(i,idim,iGrp)  = max(cruold(igridn(i,ig2)+ih2,icrE),smallcr)
                  dx_d(i,idim)=dx_loc
              else
                  pcrd(i,idim,iGrp)  = max(cruold(ind_right(i,idim),icrE),smallcr)
                  dx_d(i,idim)=dx_loc*1.5_dp
              end if
           enddo
           enddo
        end do
        ! End loop over dimensions

        do i=1,ngrid
        do iGrp=1,ncr_groups
           do idim=1,ndim
              gradEcr_loc(i,idim,iGrp) = (pcrd(i,idim,iGrp)-pcrg(i,idim,iGrp)) &
                   &                    / (dx_g(i,idim)     +dx_d(i,idim))
              gradpcr_loc(i,idim,iGrp) = gradEcr_loc(i,idim,iGrp)*(gamma_cr(iGrp)-1.)
           enddo
           vs_loc(i,:)=0.
           va_loc(i)=0.
           do idim=1,3
              B_field_loc(i,idim) = 0.5 * (uold(ind_cell(i), 5+idim) + uold(ind_cell(i), nvar+idim) )
              vs_loc(i,idim) = B_field_loc(i,idim)/sqrt(uold(ind_cell(i),1))
              va_loc(i) = va_loc(i) + vs_loc(i,idim)**2
           enddo
           va_loc(i) = sqrt(va_loc(i))

           if(v_alfven.gt.0.) then
              vs_loc(i,:) = 0.
              vs_loc(i,1) = v_alfven
              va_loc(i) = v_alfven
           endif

           bdotgradE_loc(i)=0.
           ! bdotgrade is needed for eq 3 in Jiang & Oh
           do idim=1,ndim
              bdotgradE_loc(i) = bdotgradE_loc(i) + B_field_loc(i,idim)*gradEcr_loc(i,idim,iGrp)
           enddo
           ! Local streaming velocity, needed for eq 18 in Jiang & Oh
           vs_loc(i,:) = -vs_loc(i,:)  * bdotgradE_loc(i)/max(1d-50,abs(bdotgradE_loc(i))) !! eq 3 in PO17

           ! Source term, eq 18 in Jiang & Oh 2017

            icrE = iCRu+(ndim+1)*(iGrp-1)  ! starting index of cr variables
            norm=0.
            do j=1,3
               ! Reuse the face-averaged B already stored in B_field_loc(i,:)
               ! at the top of this i/iGrp iteration (identical expression).
               B_field(j) = B_field_loc(i,j)
               norm = norm + B_field(j)**2
            end do
            norm = max(sqrt(norm), 1d-30)

            bxby = sqrt(B_field(1)**2+B_field(2)**2)
            if(norm.gt.1e-10) then
               sint = bxby/norm
               cost = B_field(3)/norm
            else
               sint = 1d0
               cost = 0d0
            endif
            if(bxby.gt.1e-10) then
               sinp = B_field(2)/bxby
               cosp = B_field(1)/bxby
            else
               sinp = 0d0
               cosp = 1d0
            endif
            !B_field = B_field/norm

            rhs1 = max(crunew(ind_cell(i),icrE),smallcr)

            ! Flux update ... first rotate flux vector onto B-field,
            ! i.e. describing the flux in the B coordinate system
            f2=0. ; f3=0.
            f1 = crunew(ind_cell(i),icrE+1)
#if NDIM>1
            f2 = crunew(ind_cell(i),icrE+2)
#endif
#if NDIM>2
            f3 = crunew(ind_cell(i),icrE+3)
#endif

            ! Make sure that |F|<=cE/sqrt3 (sqrt3=1 or sqrt(3) for M1 or P1 resp.)
            fred = sqrt(f1**2+f2**2+f3**2)/(cr_vmax(ilevel)*rhs1)*sqrt3
            if(fred .gt. 1d0) then
               f1 = f1/fred ; f2 = f2/fred ; f3 = f3/fred
            endif

            frotx=f1; froty=f2; frotz=f3
            call rotatevec(sint, cost, sinp, cosp, frotx, froty, frotz)
            v1=0. ; v2=0. ; v3=0.
            v1 = uold(ind_cell(i),2)/uold(ind_cell(i),1)
#if NDIM>1
            v2 = uold(ind_cell(i),3)/uold(ind_cell(i),1)
#endif
#if NDIM>2
            v3 = uold(ind_cell(i),4)/uold(ind_cell(i),1)
#endif
            vtot1 = v1 ; vtot2 = v2 ; vtot3 = v3
            if(mom_streaming_heating) then
               vtot1 = vtot1+vs_loc(i,1)
               vtot2 = vtot2+vs_loc(i,2)
               vtot3 = vtot3+vs_loc(i,3)
            endif

            ! Rotate velocity
            call rotatevec(sint, cost, sinp, cosp, v1, v2, v3)
            ! Rotate streaming velocity
            call rotatevec(sint, cost, sinp, cosp, vtot1, vtot2, vtot3)

            rhs2 = frotx
            rhs3 = froty
            rhs4 = frotz
            if(abs(rhs2).lt.1e-10*smallcr) rhs2=0d0
            if(abs(rhs3).lt.1e-10*smallcr) rhs3=0d0
            if(abs(rhs4).lt.1e-10*smallcr) rhs4=0d0

            ! Factor for decoupling CRs from gas at low densities
            f_decouple = MAX(exp(-smallr*cr_smallr_decouple/uold(ind_cell(i),1)),1d-10)

            sigma_stream = max(1./DCRmax_code &
               & ,abs(bdotgradE_loc(i))/3d0 / norm / va_loc(i) / gamma_cr(iGrp) / max(crunew(ind_cell(i),icrE),smallcr))
            if(mom_streaming_diffusion) then
               sigma_x = 1./(Dcr_code(iGrp) + 1./sigma_stream)
            else
               sigma_x = 1./Dcr_code(iGrp)
            endif
            sigma_y = 1./(Dcr_code(iGrp)*Dcr_perp_factor(iGrp))
            sigma_z = 1./(Dcr_code(iGrp)*Dcr_perp_factor(iGrp))

            coef_11 = 1.0 - dt * sigma_x * vtot1 * v1 * gamma_cr(iGrp) *3d0*(gamma_cr(iGrp)-1d0)   &
                          - dt * sigma_y * vtot2 * v2 * gamma_cr(iGrp) *3d0*(gamma_cr(iGrp)-1d0)   &
                          - dt * sigma_z * vtot3 * v3 * gamma_cr(iGrp) *3d0*(gamma_cr(iGrp)-1d0)
            coef_12 = dt * sigma_x * vtot1 *3d0*(gamma_cr(iGrp)-1d0)
            coef_13 = dt * sigma_y * vtot2 *3d0*(gamma_cr(iGrp)-1d0)
            coef_14 = dt * sigma_z * vtot3 *3d0*(gamma_cr(iGrp)-1d0)

            coef_21 = -dt * v1 * sigma_x * gamma_cr(iGrp) * cr_vmax(ilevel)**2
            coef_22 = 1.0 + dt * sigma_x * cr_vmax(ilevel)**2

            coef_31 = -dt * v2 * sigma_y * gamma_cr(iGrp) * cr_vmax(ilevel)**2
            coef_33 = 1.0 + dt * sigma_y * cr_vmax(ilevel)**2

            coef_41 = -dt * v3 * sigma_z * gamma_cr(iGrp) * cr_vmax(ilevel)**2
            coef_44 = 1.0 + dt * sigma_z * cr_vmax(ilevel)**2

            ! newfr1 = (rhs2 - coef21 * newEc)/coef22
            ! newfr2= (rhs3 - coef31 * newEc)/coef33
            ! newfr3 = (rhs4 - coef41 * newEc)/coef44
            ! coef11 - coef21 * coef12 /coef22 - coef13 * coef31 /coef33 - coef41 * coef14 /coef44)* newec
            !    =rhs1 - coef12 *rhs2/coef22 - coef13 * rhs3/coef33 - coef14 * rhs4/coef44

            e_coef = coef_11 - coef_12 * coef_21/coef_22 - coef_13 * coef_31/coef_33 &
                            - coef_14 * coef_41/coef_44
            new_ec = rhs1 - coef_12 * rhs2/coef_22 - coef_13 * rhs3/coef_33 &
                         - coef_14 * rhs4/coef_44
            new_ec = new_ec / e_coef

            old_ec = crunew(ind_cell(i),icrE)
            crunew(ind_cell(i),icrE) = crunew(ind_cell(i),icrE) + (new_ec-old_ec) * f_decouple

            ! Floor the CR energy and update total energy if necessary
            if ( crunew(ind_cell(i),icrE) .lt. smallcr ) crunew(ind_cell(i),icrE) = smallcr
            ! Thermal energy update:
            if(.not. static_gas .and. .not. static) then
               unew(ind_cell(i),5) = unew(ind_cell(i),5) - (crunew(ind_cell(i),icrE) - old_ec)*f_decouple
               unew(ind_cell(i),5) = max(smallp*uold(ind_cell(i),1), unew(ind_cell(i),5))
            endif
            frotx = (rhs2 - coef_21 * new_ec)/coef_22
            froty = (rhs3 - coef_31 * new_ec)/coef_33
            frotz = (rhs4 - coef_41 * new_ec)/coef_44

            ! Make sure that |F|<=cE/sqrt3 (sqrt3=1 or sqrt(3) for M1 or P1 resp.)
            fred = sqrt(frotx**2+froty**2+frotz**2)/(cr_vmax(ilevel)*new_ec)*sqrt3
            if(fred .gt. 1d0) then
               !print*,'maybe a problem with fred'
               frotx = frotx/fred ; froty = froty/fred ; frotz = frotz/fred
            endif

            ! We are missing the perpendicular energy source term which is done here in Athena!!!
            ! Rotate the flux back to the simulation coordinate system
            call invrotatevec(sint, cost, sinp, cosp, frotx, froty, frotz)

            crunew(ind_cell(i),icrE+1) = frotx
            ! Momentum update
            if(.not. static_gas .and. .not. static) then
               if (gradpcr_mom) then
                  mom_change = -gradpcr_loc(i,1,iGrp)*dt
               else
                  mom_change = ( f1 - crunew(ind_cell(i),icrE+1) ) / cr_vmax(ilevel)**2
               endif
               unew(ind_cell(i),2) = unew(ind_cell(i),2) + mom_change*f_decouple
            endif
#if NDIM>1
            crunew(ind_cell(i),icrE+2) = froty
            if(.not. static_gas .and. .not. static) then
               if (gradpcr_mom) then
                  mom_change = -gradpcr_loc(i,2,iGrp)*dt
               else
                  mom_change = ( f2 - crunew(ind_cell(i),icrE+2) ) / cr_vmax(ilevel)**2
               endif
               unew(ind_cell(i),3) = unew(ind_cell(i),3) + mom_change*f_decouple
            endif
#endif
#if NDIM>2
            crunew(ind_cell(i),icrE+3) = frotz
            if(.not. static_gas .and. .not. static) then
               if (gradpcr_mom) then
                  mom_change = -gradpcr_loc(i,3,iGrp)*dt
               else
                  mom_change = ( f3 - crunew(ind_cell(i),icrE+3) ) / cr_vmax(ilevel)**2
               endif
               unew(ind_cell(i),4) = unew(ind_cell(i),4) + mom_change*f_decouple
            endif
#endif

        end do ! ncr_groups

            ekin=0d0
            do idim=1,ndim
               ekin = ekin + 0.5d0*unew(ind_cell(i),idim+1)**2/unew(ind_cell(i),1)
            end do
            err=0d0
            emag=0d0
            do idim=1,ndim
               emag = emag + 0.125d0*(unew(ind_cell(i),idim+5)+unew(ind_cell(i),idim+nvar))**2
            end do

            etherm = unew(ind_cell(i),5) - ekin - emag - err
            if(etherm .lt. smallp*uold(ind_cell(i),1)) then
               unew(ind_cell(i),5) = ekin + emag + err + smallp*uold(ind_cell(i),1)
               if(myid.eq.1) print*,'Noooooo negative energy!!'
               !call clean_stop
            endif

        end do ! ncache

     enddo
     ! End loop over cells
  end do
  ! End loop over grids

111 format('   Entering add_cr_source_terms for level ',i2)

end subroutine add_cr_source_terms

!************************************************************************
SUBROUTINE cr_godunov_fine(ilevel)

! This routine is a wrapper to the grid solver for cosmic-ray transport.
! Small grids (2x2x2) are gathered from level ilevel and sent to the CR
! solver. On entry, CR variables are gathered from array cruold. On exit,
! crunew has been updated.
!------------------------------------------------------------------------
   use amr_commons
   use cr_hydro_commons
   implicit none
   integer::ilevel
   integer::i,igrid,ncache,ngrid
   integer,dimension(1:nvector),save::ind_grid
!------------------------------------------------------------------------
   if(numbtot(1,ilevel)==0)return  ! # of grids at ilevel
   if(verbose)write(*,111)ilevel

   ncache=active(ilevel)%ngrid  ! total # of grids at level ilevel
   do igrid=1,ncache,nvector    ! take steps of nvector grids up to ncache
      ngrid=MIN(nvector,ncache-igrid+1) ! # of grids in each sweep
      do i=1,ngrid              ! collect grid indices for one sweep
         ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
      end do
      call cr_godfine1(ind_grid,ngrid,ilevel)
   end do

   if(verbose)write(*,112)ilevel

111 format('   Entering cr_godunov_fine for level ',i2)
112 format('   Exiting cr_godunov_fine for level ',i2)

END SUBROUTINE cr_godunov_fine

!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE cr_set_unew(ilevel)
   use amr_commons
   use cr_hydro_commons
   use cr_parameters
   implicit none
   integer::ilevel
   !--------------------------------------------------------------------------
   ! This routine sets array crunew to its initial value cruold before
   ! calling the CR advection scheme. crunew is set to zero in virtual
   ! boundaries.
   !--------------------------------------------------------------------------
   integer::i,ivar,ind,icpu,iskip
   real(dp)::fred,sqrt3

   if(numbtot(1,ilevel)==0)return
   if(verbose)write(*,111)ilevel

   if(isotropic_pressure)then
      sqrt3=sqrt(3d0)
   else
      sqrt3=1.0d0
   endif

   if(reduced_CR_flux_correction)then
      do ind=1,twotondim
         iskip=ncoarse+(ind-1)*ngridmax
         do i=1,active(ilevel)%ngrid
            fred=sqrt(sum(cruold(active(ilevel)%igrid(i)+iskip,iCRu+1:iCRu+ndim)**2)) &
                 /(cr_vmax(ilevel)*cruold(active(ilevel)%igrid(i)+iskip,iCRu))*sqrt3
            if(fred>1.0)then
               cruold(active(ilevel)%igrid(i)+iskip,iCRu+1:iCRu+ndim) = &
                    cruold(active(ilevel)%igrid(i)+iskip,iCRu+1:iCRu+ndim)/fred
            endif
         end do
      end do
   endif

   ! Set crunew to cruold for myid cells
   do ind=1,twotondim
      iskip=ncoarse+(ind-1)*ngridmax
      do ivar=1,ncrvar
         do i=1,active(ilevel)%ngrid
            crunew(active(ilevel)%igrid(i)+iskip,ivar)=cruold(active(ilevel)%igrid(i)+iskip,ivar)
         end do
      end do
   end do

   ! Set crunew to 0 for virtual boundary cells
   do icpu=1,ncpu
   do ind=1,twotondim
      iskip=ncoarse+(ind-1)*ngridmax
      do ivar=1,ncrvar
         do i=1,reception(icpu,ilevel)%ngrid
#ifdef LIGHT_MPI_COMM
            crunew(reception(icpu,ilevel)%pcomm%igrid(i)+iskip,ivar)=0.0
#else
            crunew(reception(icpu,ilevel)%igrid(i)+iskip,ivar)=0.0
#endif
         end do
      end do
   end do
   end do

111 format('   Entering cr_set_unew for level ',i2)

END SUBROUTINE cr_set_unew

!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE cr_set_uold(ilevel)
   use amr_commons
   use cr_hydro_commons
   use cr_parameters
   implicit none
   integer::ilevel
   !--------------------------------------------------------------------------
   ! This routine sets array cruold to its new value crunew after the CR
   ! step.
   !
   ! Transport-robustness floor: the per-group CR energy is floored at
   ! smallcr after the update. This makes the explicit transport robust on
   ! its own -- without it the explicit transport of a degenerate state can
   ! drive E_cr slightly negative and trip the Ecr<0 guard in
   ! cmp_cr_flux_tensors. This is pure transport robustness -- no
   ! scattering/streaming/cooling source physics.
   !--------------------------------------------------------------------------
   integer::i,ivar,ind,iskip,iGrp,icrE
   real(dp)::fred,sqrt3

   if(numbtot(1,ilevel)==0)return
   if(verbose)write(*,111)ilevel

   if(isotropic_pressure)then
      sqrt3=sqrt(3d0)
   else
      sqrt3=1.0d0
   endif

   if(reduced_CR_flux_correction)then
      do ind=1,twotondim
         iskip=ncoarse+(ind-1)*ngridmax
         do i=1,active(ilevel)%ngrid
            fred=sqrt(sum(crunew(active(ilevel)%igrid(i)+iskip,iCRu+1:iCRu+ndim)**2)) &
                 /(cr_vmax(ilevel)*crunew(active(ilevel)%igrid(i)+iskip,iCRu))*sqrt3
            if(fred>1.0)then
               crunew(active(ilevel)%igrid(i)+iskip,iCRu+1:iCRu+ndim) = &
                    crunew(active(ilevel)%igrid(i)+iskip,iCRu+1:iCRu+ndim)/fred
            endif
         end do
      end do
   endif

   ! Set cruold to crunew for myid cells
   do ind=1,twotondim
      iskip=ncoarse+(ind-1)*ngridmax
      do ivar=1,ncrvar
         do i=1,active(ilevel)%ngrid
            cruold(active(ilevel)%igrid(i)+iskip,ivar)=crunew(active(ilevel)%igrid(i)+iskip,ivar)
         end do
      end do

      ! Floor each group's CR energy at smallcr (no negative CR densities).
      do iGrp=1,ncr_groups
         icrE=iCRu+(ndim+1)*(iGrp-1)
         do i=1,active(ilevel)%ngrid
            cruold(active(ilevel)%igrid(i)+iskip,icrE)= &
                 max(cruold(active(ilevel)%igrid(i)+iskip,icrE),smallcr)
         end do
      end do
   end do

111 format('   Entering cr_set_uold for level ',i2)

END SUBROUTINE cr_set_uold

!*************************************************************************
SUBROUTINE cr_godfine1(ind_grid, ncache, ilevel)

! This routine first gathers CR variables from neighboring grids to set
! initial conditions in a 6x6x6 grid. It interpolates from coarser level to
! missing grid variables. It then calls the solver that computes fluxes.
! These fluxes are zeroed at coarse-fine boundaries, since contribution from
! finer levels has already been taken into account. Conservative variables
! are updated and stored in array crunew(:), both at the current level and at
! the coarser level if necessary.
!
! in ind_grid: Indexes of grids/octs to solve in
! in ncache:   Length of ind_grid (i.e. number of grids)
! in ilevel:   Level at which the grids are
!-------------------------------------------------------------------------
   use amr_commons
   use hydro_commons
   use cr_hydro_commons
   use cr_parameters
   use cr_flux_module
   use amr_constants, only:i1min,i1max,j1min,j1max,k1min,k1max, &
                        &  i2min,i2max,j2min,j2max,k2min,k2max, &
                        &  i3min,i3max,j3min,j3max,k3min,k3max
   implicit none
   integer::ilevel,ncache
   real(dp)::dt
   integer   ,dimension(1:nvector)::ind_grid
   integer   ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
   integer   ,dimension(1:nvector,0:twondim    ),save::ibuffer_father
   integer   ,dimension(1:nvector,0:twondim    ),save::ind1
   ! Split buffers for the coarser-level neighbour interpolation:
   real(dp),dimension(1:nvector,0:twondim  ,1:nvar+3),save::u1_gas
   real(dp),dimension(1:nvector,1:twotondim,1:nvar+3),save::u2_gas
   real(dp),dimension(1:nvector,0:twondim  ,1:ncrvar),save::u1_cr
   real(dp),dimension(1:nvector,1:twotondim,1:ncrvar),save::u2_cr

   ! Split 6x6x6 stencils:
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3),save::uin_gas
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ncrvar),save::uin_cr
   real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:ncrvar,1:ndim),save::flux
   logical,dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::ok
   integer,dimension(1:nvector),save::igrid_nbor,ind_cell,ind_buffer
   integer,dimension(1:nvector),save::ind_exist,ind_nexist

   integer::i,j,idim,ind_son,ind_father,iskip,nbuffer
   integer::i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,k3,nx_loc,nb_noneigh,nexist
   integer::ivar,iGrp
   real(dp)::dx,scale,oneontwotondim
   real(dp)::fred,sqrt3
!------------------------------------------------------------------------
   oneontwotondim = 1d0/dble(twotondim) ! 1/8 in 3D
   ! Mesh spacing
   nx_loc=icoarse_max-icoarse_min+1    ! =1
   scale=boxlen/dble(nx_loc)           ! length per coarse oct (=boxlen)
   dx=0.5D0**ilevel*scale              ! length per oct/grid at ilevel
   dt=dtnew(ilevel)
   if(isotropic_pressure)then
      sqrt3=sqrt(3d0)
   else
      sqrt3=1.0d0
   endif

   !------------------------------------------
   ! Gather 3^ndim neighboring father cells
   ! of grid (ilevel-1). ind_cell are indexes
   ! in uold/cruold.
   !------------------------------------------
   do i=1,ncache
      ind_cell(i)=father(ind_grid(i))
   end do
   call get3cubefather(ind_cell,nbors_father_cells,ncache,ilevel)

   !---------------------------
   ! Gather 6x6x6 cells stencil
   !---------------------------
   ! Loop over 3x3x3 neighboring father cells (ilevel-1)
   do k1=k1min,k1max
   do j1=j1min,j1max
   do i1=i1min,i1max

      ! Check if neighbor has a grid (or is just a cell)
      nbuffer=0
      nexist=0

      ind_father=1+i1+3*j1+9*k1
      do i=1,ncache
         igrid_nbor(i)=son(nbors_father_cells(i,ind_father))
         if(igrid_nbor(i)>0) then
            nexist=nexist+1
            ind_exist(nexist)=i
         else
            nbuffer=nbuffer+1
            ind_nexist(nbuffer)=i
            ind_buffer(nbuffer)=nbors_father_cells(i,ind_father)
         end if
      end do

      if(nbuffer>0)then
         ! For those nongrid cells we interpolate variables from parent cells.
         call getnborfather(ind_buffer,ibuffer_father,nbuffer,ilevel)
         do j=0,twondim
            ! Gas (rho/vel/B) from uold, interpolated with the MHD-aware
            ! interpol_hydro (preserves div B).
            do ivar=1,nvar+3
               do i=1,nbuffer
                  u1_gas(i,j,ivar)=uold(ibuffer_father(i,j),ivar)
               end do
            end do
            ! CR (E,F) from cruold, interpolated as cell-centered scalars.
            do ivar=1,ncrvar
               do i=1,nbuffer
                  u1_cr(i,j,ivar)=cruold(ibuffer_father(i,j),ivar)
               end do
            end do
            do i=1,nbuffer
               ind1(i,j)=son(ibuffer_father(i,j))
            end do
         end do
         call interpol_hydro(u1_gas,ind1,u2_gas,nbuffer)
         call cr_interpol_hydro(u1_cr,u2_cr,nbuffer)
         if(reduced_CR_flux_correction)then
            do j=1,twotondim
               do i=1,nbuffer
                  ! F<e*c/sqrt(3) for P1 closure
                  fred=sqrt(sum(u2_cr(i,j,iCRu+1:iCRu+ndim)**2)) &
                       /(cr_vmax(ilevel)*u2_cr(i,j,iCRu))*sqrt3
                  if(fred.gt.1.0)then
                     u2_cr(i,j,iCRu+1:iCRu+ndim)=u2_cr(i,j,iCRu+1:iCRu+ndim)/fred
                  endif
               end do
            end do
         endif
      endif

      ! Loop 2x2x2 cells within father cell and add them to the stencils
      do k2=k2min,k2max
      do j2=j2min,j2max
      do i2=i2min,i2max

         ind_son=1+i2+2*j2+4*k2
         iskip=ncoarse+(ind_son-1)*ngridmax
         do i=1,nexist
            ind_cell(i)=iskip+igrid_nbor(ind_exist(i))
         end do

         i3=1; j3=1; k3=1
         if(ndim>0)i3=1+2*(i1-1)+i2 ! From -1 to 4 over outer loop, but
         if(ndim>1)j3=1+2*(j1-1)+j2 ! only over 2x2x2 indexes in inner loop
         if(ndim>2)k3=1+2*(k1-1)+k2

         ! Gather gas variables (rho/vel/B) from uold into uin_gas
         do ivar=1,nvar+3
            do i=1,nexist
               uin_gas(ind_exist(i),i3,j3,k3,ivar)=uold(ind_cell(i),ivar)
            end do
            do i=1,nbuffer
               uin_gas(ind_nexist(i),i3,j3,k3,ivar)=u2_gas(i,ind_son,ivar)
            end do
         end do

         ! Gather CR variables (E,F) from cruold into uin_cr
         do ivar=1,ncrvar
            do i=1,nexist
               uin_cr(ind_exist(i),i3,j3,k3,ivar)=cruold(ind_cell(i),ivar)
            end do
            do i=1,nbuffer
               uin_cr(ind_nexist(i),i3,j3,k3,ivar)=u2_cr(i,ind_son,ivar)
            end do
         end do

         ! Gather refinement flag
         do i=1,nexist
            ok(ind_exist(i),i3,j3,k3)=son(ind_cell(i))>0
         end do
         do i=1,nbuffer
            ok(ind_nexist(i),i3,j3,k3)=.false.
         end do

      end do
      end do
      end do
      ! End loop over cells

   end do
   end do
   end do
   ! End loop over neighboring grids

   !----------------------------------------------------------------------
   ! Compute fluxes of each CR group, using the Eddington tensor
   !----------------------------------------------------------------------
   do iGrp=1,ncr_groups
      call cmp_cr_faces(uin_gas,uin_cr,flux,dx,dt,iGrp,ncache,ilevel)
   end do

   !----------------------------------------------------------------------
   ! Reset flux along direction at refined interface, if not cr-subcycling
   !----------------------------------------------------------------------
   if(cr_nsubcycle==1)then
   do idim=1,ndim
      i0=0; j0=0; k0=0
      if(idim==1)i0=1
      if(idim==2)j0=1
      if(idim==3)k0=1
      do k3=k3min,k3max+k0
      do j3=j3min,j3max+j0
      do i3=i3min,i3max+i0
         do i=1,ncache
            if(ok(i,i3-i0,j3-j0,k3-k0) .or. ok(i,i3,j3,k3))then
               flux(i,i3,j3,k3,:,idim)=0.0d0
            end if
         end do
      end do
      end do
      end do
   end do
   endif

   !--------------------------------------
   ! Conservative update at level ilevel
   !--------------------------------------
   do idim=1,ndim
      i0=0; j0=0; k0=0
      if(idim==1)i0=1
      if(idim==2)j0=1
      if(idim==3)k0=1
      do k2=k2min,k2max      ! all from 0 to 1
      do j2=j2min,j2max      ! => update 2x2x2 cells
      do i2=i2min,i2max      ! i.e. update one grid at level ilevel
         ind_son=1+i2+2*j2+4*k2
         iskip=ncoarse+(ind_son-1)*ngridmax
         do i=1,ncache
            ind_cell(i)=iskip+ind_grid(i)
         end do
         i3=1+i2
         j3=1+j2   ! just because the flux indexes are (1:3), not (0:2)
         k3=1+k2
         ! Update conservative CR variables (new state vector)
         do i=1,ncache
            crunew(ind_cell(i),1:ncrvar)= &
               crunew(ind_cell(i),1:ncrvar)   &
               & +(flux(i,i3   ,j3   ,k3   ,1:ncrvar,idim)  &
               & - flux(i,i3+i0,j3+j0,k3+k0,1:ncrvar,idim))
         end do
      end do
      end do
      end do
   end do

   !--------------------------------------
   ! Conservative update at level ilevel-1
   !--------------------------------------
   if(cr_nsubcycle==1)then
      ! Loop over dimensions
      do idim=1,ndim
         i0=0; j0=0; k0=0
         if(idim==1)i0=1
         if(idim==2)j0=1
         if(idim==3)k0=1

         !----------------------
         ! Left flux at boundary
         !----------------------
         ! Check if grid sits near left boundary
         ! and gather neighbor father cells index
         nb_noneigh=0
         do i=1,ncache
            if(son(nbor(ind_grid(i),2*idim-1))==0)then
               nb_noneigh=nb_noneigh+1
               ind_buffer(nb_noneigh)=nbor(ind_grid(i),2*idim-1)
               ind_cell(nb_noneigh)=i
            end if
         end do
         ! Conservative update of new state variables
         ! Loop over boundary cells
         do k3=k3min,k3max-k0 ! 1 to 1 if dim=3, 1 to 2 otherwise
            do j3=j3min,j3max-j0 ! 1 to 1 if dim=2, 1 to 2 otherwise
               do i3=i3min,i3max-i0 ! 1 to 1 if dim=1, 1 to 2 otherwise
                  do i=1,nb_noneigh
                     crunew(ind_buffer(i),1:ncrvar)= &
                          crunew(ind_buffer(i),1:ncrvar)   &
                          & -flux(ind_cell(i),i3,j3,k3,:,idim)  &
                          & *oneontwotondim
                  end do
               end do
            end do
         end do

         !-----------------------
         ! Right flux at boundary
         !-----------------------
         ! Check if grid sits near right boundary
         ! and gather neighbor father cells index
         nb_noneigh=0
         do i=1,ncache
            if(son(nbor(ind_grid(i),2*idim))==0)then
               nb_noneigh=nb_noneigh+1
               ind_buffer(nb_noneigh)=nbor(ind_grid(i),2*idim)
               ind_cell(nb_noneigh)=i
            end if
         end do
         ! Conservative update of new state variables
         ! Loop over boundary cells
         do k3=k3min+k0,k3max
            do j3=j3min+j0,j3max
               do i3=i3min+i0,i3max
                  do i=1,nb_noneigh
                     crunew(ind_buffer(i),1:ncrvar)= &
                          crunew(ind_buffer(i),1:ncrvar)   &
                          & +flux(ind_cell(i),i3+i0,j3+j0,k3+k0,:,idim)  &
                          & *oneontwotondim
                  end do
               end do
            end do
         end do

      end do
      ! End loop over dimensions
   end if
   ! End if-clause for cr-subcycling

END SUBROUTINE cr_godfine1
