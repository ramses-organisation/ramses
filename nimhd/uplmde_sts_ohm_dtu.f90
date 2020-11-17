!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine diffusion_sts_dtu
   use amr_commons
   use hydro_commons
   use mpi_mod
   implicit none
   integer::ilevel,Nsub,nsubdiff_loc,nsubdiff_all
   integer::icycle,nsubdiff,ivar,i,iskip,ind,info,jcycle,ifirst,isafe
   real(dp)::dx,scale,dx_loc,dtdiffohm,dtdiffamb,norm,nu,dtad,dtohm
   real(dp)::dtdiff,dtdiffold,dtdiffsum,dtdiffcoef,dtdiffcoef_all,dtdiffcoef_loc,sommeinit,dtloc,dtall
   real(dp),dimension(1:3)::skip_loc
   real(dp)::max_diff_loc,max_diff_all,dtdiff_loc,dtdiff_all,dtdiffexpl
   logical::sts
!   if(numbtot(1,ilevel)==0)return 
   ! Rescaling factors
   skip_loc=(/0.0d0,0.0d0,0.0d0/)
   if(ndim>0)skip_loc(1)=dble(icoarse_min)
   if(ndim>1)skip_loc(2)=dble(jcoarse_min)
   if(ndim>2)skip_loc(3)=dble(kcoarse_min)
   scale=dble(icoarse_max-icoarse_min+1)/boxlen

!!$    call set_unew_sts(ilevel,dtloc,dx_loc,0)  ! needed ? probably yes
!!$#ifndef WITHOUTMPI
!!$   call MPI_ALLREDUCE(dtloc     ,dtall      ,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
!!$      &MPI_COMM_WORLD,info)
!!$   dtloc=dtall
!!$#endif
!!$   dtdiffohm=dtmagdiff(ilevel)
!!$   nu=0.0d0
!!$   nsub=1
!!$   do ilevel=levelmin,nlevelmax
!!$      if(numbtot(1,ilevel)==0)cycle 
!!$      dtdiffamb=dtambdiff(ilevel)
!!$      dtdiffohm=dtmagdiff(ilevel)
!!$!   dtdiff=min(dtdiffohm,dtdiffamb)
!!$      dtdiff=1.0d0/dtdiffohm+1.0d0/dtdiffamb
!!$      dtdiff=1.0d0/dtdiff
!!$      dtdiffold=min(dtdiff,dtnew(ilevel))
!!$      dtdiffsum=0.0d0
!!$      !nu=0
!!$      nu = max(1d0/(4d0*nsts(ilevel)**2),nu)
!!$      nsub=max(nsub,nsts(ilevel))
!!$!      if(myid==1)print*,ilevel,dtambdiff(ilevel),dtmagdiff(ilevel),nsts(ilevel),dtnew(ilevel),dtsts(ilevel)
!!$   end do
!!$   call make_virtual_fine_dp(unew(1,1),ilevel)
!!$   call make_virtual_fine_dp(unew(1,5),ilevel)
!!$   call make_virtual_fine_dp(unew(1,6),ilevel)
!!$   call make_virtual_fine_dp(unew(1,7),ilevel)
!!$   call make_virtual_fine_dp(unew(1,8),ilevel)
!!$   call make_virtual_fine_dp(unew(1,nvar+1),ilevel)
!!$   call make_virtual_fine_dp(unew(1,nvar+2),ilevel)
!!$   call make_virtual_fine_dp(unew(1,nvar+3),ilevel)
!!$   call make_virtual_fine_dp(unew(1,nvar),ilevel)

   jcycle=0
!!$   do  jcycle=1,100
   dtdiffcoef=1.0d0
   dtdiffsum=0.0d0
   dtdiffexpl=dtnew(levelmin)
   nsubdiff = 1
   do ilevel=levelmin,nlevelmax
      
      if(numbtot(1,ilevel)==0)cycle 
      dx=0.5D0**ilevel
      dx_loc=dx/scale
      
      call set_unew_sts(ilevel,0,dtohm,dtad)  ! needed ? probably yes
      
      dtdiffohm=dtohm!magdiff(ilevel)
      dtdiffamb=dtad !mbdiff(ilevel)
!!$      dtdiffamb=dtloc!dtambdiff(ilevel)
!!$      dtdiff_loc=min(dtdiffohm,dtdiffamb)
      dtdiff_loc=1.0d0/dtdiffohm+1.0d0/dtdiffamb
      dtdiff_loc=1.0d0/dtdiff_loc
!      print*,dtdiff_loc,dtohm,dtad
      
!      dtdiffold=dtdiff_loc!min(dtdiffohm,dtdiffamb)

#ifndef WITHOUTMPI
      call MPI_ALLREDUCE(dtdiff_loc     ,dtdiff_all      ,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
         &MPI_COMM_WORLD,info)
      dtdiff_loc=dtdiff_all
#endif

!       if(myid==1)print*,dtnew(ilevel),dtdiffsum,dtdiff,nsts(ilevel),dtsts(ilevel)
      nsubdiff_loc=ceiling(sqrt((dtnew(ilevel)-dtdiffsum)/dtdiff_loc))
!!$      if(nsts(ilevel) .gt. nsubdiff)then
!!$!         print*,ilevel,numbtot(1,ilevel),nsts(ilevel),dtnew(ilevel),dtsts(ilevel)
!!$         
!!$         nsubdiff_loc=max(nsubdiff_loc,nsts(ilevel))
!!$         dtdiffcoef_loc = max(dtnew(ilevel)/dtsts(ilevel),dtdiffcoef_loc)
!!$      end if

      if(nsubdiff_loc .gt. nsubdiff)then
         dtdiffexpl = dtdiff_loc
         nsubdiff = nsubdiff_loc
!         dtdtdiffcoef = min(dtnew(ilevel)/dtdiffexpl,1)
      end if

      call make_virtual_fine_dp(unew(1,1),ilevel)
      call make_virtual_fine_dp(unew(1,5),ilevel)
      call make_virtual_fine_dp(unew(1,6),ilevel)
      call make_virtual_fine_dp(unew(1,7),ilevel)
      call make_virtual_fine_dp(unew(1,8),ilevel)
      call make_virtual_fine_dp(unew(1,nvar+1),ilevel)
      call make_virtual_fine_dp(unew(1,nvar+2),ilevel)
      call make_virtual_fine_dp(unew(1,nvar+3),ilevel)
      call make_virtual_fine_dp(unew(1,nvar),ilevel)
      call make_virtual_fine_dp(enew(1),ilevel)
      call make_virtual_fine_dp(divu(1),ilevel)

   enddo!end loop over levels

   dtdiffold=dtdiffexpl
   sommeinit = dtdiffold*nsubdiff/(2.*sqrt(nu_sts))*((1.+sqrt(nu_sts))**(2.*nsubdiff)-(1.-sqrt(nu_sts))**(2.*nsubdiff)) &         ! somme des dtsts, theorique
        & /((1.+sqrt(nu_sts))**(2.*nsubdiff)+(1.-sqrt(nu_sts))**(2.*nsubdiff))
   if  (sommeinit > dtnew(levelmin)-dtdiffsum) then
      dtdiffcoef = (dtnew(levelmin)-dtdiffsum)/sommeinit  ! coefficient pour rescaler les dtsts
   else
      dtdiffcoef = 1.
   end if
!print*,dtnew(levelmin),dtdiffsum,sommeinit
!!$#ifndef WITHOUTMPI
!!$      call MPI_ALLREDUCE(nsubdiff     ,nsubdiff_all     ,1,MPI_DOUBLE_PRECISION,MPI_MAX,&
!!$         &MPI_COMM_WORLD,info)
!!$      nsubdiff = nsubdiff_all
!!$      call MPI_ALLREDUCE(dtdiffcoef_loc     ,dtdiffcoef_all      ,1,MPI_DOUBLE_PRECISION,MPI_MAX,&
!!$         &MPI_COMM_WORLD,info)
!!$      dtdiffcoef=dtdiffcoef_all
!!$#endif
!!$      
!!$      nu = 1d0/(4d0*nsubdiff_loc**2)

      do  icycle=1,nsubdiff
         if (nsubdiff > 0) then
            
!            if(myid==1 .and. (mod(nstep,ncontrol)==0))write(*,*)nu_sts,'dtdiffcoeff',dtdiffcoef,'dtdiffsum',dtdiffsum,'dtdiff',dtdiffold,'dtsts',sommeinit*dtdiffcoef,dtnew(levelmin:nlevelmax)
            
!!$            dtdiff=dtdiffold*((-1d0+nu)*cos((2d0*dble(icycle)-1d0)*acos(-1d0)/(2d0*nsubdiff))+1d0+nu)**(-1)    ! pas de temps du STS
            dtdiff=dtdiffcoef*dtdiffold*((-1d0+nu_sts)*cos((2d0*dble(icycle)-1d0)*acos(-1d0)/(2d0*nsubdiff))+1d0+nu_sts)**(-1)    ! pas de temps du STS
            !!$               if (dtdiffold==dtnew(ilevel)) dtdiff = dtdiffold
!!$               if  (dtdiff+dtdiffsum > dtnew(ilevel))then
         else if (nsubdiff == 0) then
            dtdiff = dtnew(ilevel)
         end if
         !                  write(*,*)'WARNING: leaves STS before last dt_sts',dtdiff,dtdiffsum ,dtnew(ilevel),dtdiff+dtdiffsum - dtnew(ilevel)
         !   stop
         !!$               end if
!!$               if  (nsubdiff == 1) dtdiff = dtnew(ilevel)-dtdiffsum
         dtdiffsum=dtdiffsum+dtdiff
         if(myid==1 .and. (mod(nstep,ncontrol)==0))write(*,*)'subcycling through STS_dtu',icycle,'ilevel',ilevel,'nsubdiff :',nsubdiff,&
              'dtmagdiff :',dtmagdiff(ilevel),'dtambdiff :',dtambdiff(ilevel),&
              'dtnew :',dtnew(ilevel),'dtdiff :',dtdiff,'sum :',dtdiffsum, 'jcycle :',jcycle

         do ilevel=levelmin,nlevelmax
            dx=0.5D0**ilevel
            dx_loc=dx/scale
            
            call set_unew_sts(ilevel,1,dtohm,dtad)   ! needed ? probably yes
         end do

         do ilevel=nlevelmax,levelmin,-1
            call diffusion_fine_sts_dtu(ilevel,dtdiff,icycle)
#ifndef WITHOUTMPI
            call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
            
!            call make_virtual_reverse_dp(unew(1,5),ilevel)
            call make_virtual_reverse_dp(unew(1,6),ilevel)
            call make_virtual_reverse_dp(unew(1,7),ilevel)
            call make_virtual_reverse_dp(unew(1,8),ilevel)
            call make_virtual_reverse_dp(unew(1,nvar+1),ilevel)
            call make_virtual_reverse_dp(unew(1,nvar+2),ilevel)
            call make_virtual_reverse_dp(unew(1,nvar+3),ilevel)
!!$         call make_virtual_reverse_dp(unew(1,nvar),ilevel)

         call set_uold_sts(ilevel,0,dtdiff)
!         call upload_fine(ilevel)!

!!$         call make_virtual_fine_dp(uold(1,1),ilevel)
!         call make_virtual_fine_dp(uold(1,5),ilevel)
         call make_virtual_fine_dp(uold(1,6),ilevel)
         call make_virtual_fine_dp(uold(1,7),ilevel)
         call make_virtual_fine_dp(uold(1,8),ilevel)
         call make_virtual_fine_dp(uold(1,nvar+1),ilevel)
         call make_virtual_fine_dp(uold(1,nvar+2),ilevel)
         call make_virtual_fine_dp(uold(1,nvar+3),ilevel)
!!$         call make_virtual_fine_dp(uold(1,nvar),ilevel)

         if  (nsubdiff == 1) exit 
      end do
      end do

!!$   end if

!!$   if  (jcycle == 100) then
!!$      print *, "PROBLEM IN STS"
!!$      stop
!!$   end if

#ifndef WITHOUTMPI
   call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
   do ilevel=levelmin,nlevelmax
 
   call set_uold_sts(ilevel,1,dtdiffsum)
   call upload_fine(ilevel)
   call make_virtual_fine_dp(uold(1,5),ilevel)
   call make_virtual_fine_dp(uold(1,6),ilevel)
   call make_virtual_fine_dp(uold(1,7),ilevel)
   call make_virtual_fine_dp(uold(1,8),ilevel)
   call make_virtual_fine_dp(uold(1,nvar+1),ilevel)
   call make_virtual_fine_dp(uold(1,nvar+2),ilevel)
   call make_virtual_fine_dp(uold(1,nvar+3),ilevel)
   call make_virtual_fine_dp(uold(1,nvar),ilevel)
end do
end subroutine diffusion_sts_dtu
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine diffusion_fine_sts_dtu(ilevel,dtdiff,icycle)!,first_cmp_current_sts)
   use amr_commons
   use hydro_commons
   implicit none
   integer::ilevel,icycle
   real(dp)::dtdiff

   integer::i,ivar,igrid,ncache,ngrid
   integer,dimension(1:nvector),save::ind_grid

   if(numbtot(1,ilevel)==0)return
   if(verbose)write(*,111)ilevel

   ncache=active(ilevel)%ngrid
   do igrid=1,ncache,nvector
      ngrid=MIN(nvector,ncache-igrid+1)
      do i=1,ngrid
         ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
      end do
      call diffine1_sts_dtu(ind_grid,ngrid,dtdiff,ilevel,icycle)
   end do

   111 format('   Entering diffusion_fine_sts for level ',i2)

 end subroutine diffusion_fine_sts_dtu
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine diffine1_sts_dtu(ind_grid,ncache,dtdiff,ilevel,icycle)
   use amr_commons
   use hydro_commons
   use hydro_parameters, only:iu1,iu2,ju1,ju2,ku1,ku2
   implicit none
   integer::ilevel,ncache,icycle
   real(dp)::dtdiff
   integer,dimension(1:ncache)::ind_grid
   !-------------------------------------------------------------------
   ! This routine gathers first MHD variables from neighboring grids
   ! to set initial conditions in a 6x6x6 grid. It then computes
   ! the current at cell edges. Finally, currents are corrected from finer level
   ! and boundary currents are stored in buffer regions. Updated 
   ! conservative variables are stored in array unew(:).
   !-------------------------------------------------------------------
   integer ,dimension(1:nvector,1:threetondim     ),save::nbors_father_cells
   integer ,dimension(1:nvector,1:twotondim       ),save::nbors_father_grids
   integer ,dimension(1:nvector,0:twondim         ),save::ibuffer_father
   integer ,dimension(1:nvector,0:twondim)         ,save::ind1
   real(dp),dimension(1:nvector,0:twondim  ,1:nvar+3),save::u1
   real(dp),dimension(1:nvector,1:twotondim,1:nvar+3),save::u2
   real(dp),dimension(1:nvector,0:twondim  ,1:ndim),save::g1=0.0d0
   real(dp),dimension(1:nvector,1:twotondim,1:ndim),save::g2=0.0d0

   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3),save::uloc
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim),save::gloc=0.0d0
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3),save::jcell=0d0 

   real(dp),dimension(1:nvector,1:twotondim,1:ndim),save::v2
   real(dp),dimension(1:nvector,1:ndim),save::vv,xx

   logical ,dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::ok

   REAL(dp),DIMENSION(1:nvector,1:3,1:3,1:3),save::emfx
   REAL(dp),DIMENSION(1:nvector,1:3,1:3,1:3),save::emfy
   REAL(dp),DIMENSION(1:nvector,1:3,1:3,1:3),save::emfz

   real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:3),save::flux
   real(dp),dimension(1:nvector),save :: dB

   integer,dimension(1:nvector),save::igrid_nbor,ind_cell,ind_buffer,igrid
   integer,dimension(1:nvector),save::ind_exist,ind_nexist
   logical,dimension(1:nvector),save::exist_nbor
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::jcentersquare,jxbsquare

   real(dp),dimension(1:3)::skip_loc
   integer::i,j,ivar,idim,ind_son,ind_father,iskip,nbuffer,ibuffer
   integer::ind,ix,iy,iz,nx_loc,nb_noneigh,nexist
   integer::i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,k3
   integer::i1min,i1max,j1min,j1max,k1min,k1max
   integer::i2min,i2max,j2min,j2max,k2min,k2max
   integer::i3min,i3max,j3min,j3max,k3min,k3max
   integer::ind_father1,ind_father2,ind_father3
   integer::ind_buffer1,ind_buffer2,ind_buffer3
   integer::interpol_type_old,ivar1,ivar2,ivar3,ivar4,ivar5,ivar6
   real(dp)::dflux_x,dflux_y,dflux_z
   real(dp)::dx_loc,scale,oneontwotondim
   real(dp)::emag,emagold
   real(dp)::dflux,weight,A,B,C
   integer ::mvar

   integer::neul=5

   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

   flux=0.0d0
   emfx=0d0
   emfy=0d0
   emfz=0d0

   oneontwotondim = 1d0/dble(twotondim)

   ! Mesh spacing in that level
   nx_loc=icoarse_max-icoarse_min+1
   scale=boxlen/dble(nx_loc)
   dx_loc=0.5D0**ilevel*scale

   ! Integer constants
   i1min=0; i1max=0; i2min=0; i2max=0; i3min=1; i3max=1
   j1min=0; j1max=0; j2min=0; j2max=0; j3min=1; j3max=1
   k1min=0; k1max=0; k2min=0; k2max=0; k3min=1; k3max=1
   if(ndim>0)then
      i1max=2; i2max=1; i3max=2
   end if
   if(ndim>1)then
      j1max=2; j2max=1; j3max=2
   end if
   if(ndim>2)then
      k1max=2; k2max=1; k3max=2
   end if

   ivar1=6; ivar2=7; ivar3=8
   ivar4=nvar+1; ivar5=nvar+2; ivar6=nvar+3

   ! Gather 3^ndim neighboring father cells
   do i=1,ncache
      ind_cell(i)=father(ind_grid(i))
   end do
   call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ncache,ilevel)

   !---------------------------
   ! Gather 6x6x6 cells stencil
   !---------------------------
   ! Loop over 3x3x3 neighboring father cells
   do k1=k1min,k1max
      do j1=j1min,j1max
         do i1=i1min,i1max

            ! Check if neighboring grid exists
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

            ! If not, interpolate variables from parent cells
            if(nbuffer>0)then
               call getnborfather(ind_buffer,ibuffer_father,nbuffer,ilevel)
               do j=0,twondim
                  do ivar=1,nvar+3
                     do i=1,nbuffer
                        u1(i,j,ivar)=uold(ibuffer_father(i,j),ivar)
                        if(ivar .eq. 2)then
                           u1(i,j,ivar)=enew(ibuffer_father(i,j))
                        end if
                        if(ivar .eq. 3)then
                           u1(i,j,ivar)=divu(ibuffer_father(i,j))
                        end if
                     end do
                  end do
                  do i=1,nbuffer
                     ind1(i,j)=son(ibuffer_father(i,j))
                  end do
               end do
               call interpol_hydro_sts(u1,ind1,u2,nbuffer)
            endif

            ! Loop over 2x2x2 cells
            do k2=k2min,k2max
               do j2=j2min,j2max
                  do i2=i2min,i2max
                     ind_son=1+i2+2*j2+4*k2
                     iskip=ncoarse+(ind_son-1)*ngridmax
                     do i=1,nexist
                        ind_cell(i)=iskip+igrid_nbor(ind_exist(i))
                     end do

                     i3=1; j3=1; k3=1
                     if(ndim>0)i3=1+2*(i1-1)+i2
                     if(ndim>1)j3=1+2*(j1-1)+j2
                     if(ndim>2)k3=1+2*(k1-1)+k2


                     ! Gather hydro variables
                     do ivar=1,nvar+3
                        do i=1,nexist
                           uloc(ind_exist(i),i3,j3,k3,ivar)=uold(ind_cell(i),ivar)
                           if(ivar==2)uloc(ind_exist(i),i3,j3,k3,ivar)=enew(ind_cell(i))
                           if(ivar==3)uloc(ind_exist(i),i3,j3,k3,ivar)=divu(ind_cell(i))
                        end do
                        do i=1,nbuffer
                           uloc(ind_nexist(i),i3,j3,k3,ivar)=u2(i,ind_son,ivar)
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

   !----------------
   ! Compute current
   !----------------
   call cmp_current_sts_dtu(uloc,emfx,emfy,emfz,flux,ncache,dx_loc,dx_loc,dx_loc,dtdiff,jxbsquare,jcentersquare,ok,jcell)

   !-------------------------------------------------
   ! Reset current along direction x at refined edges
   !-------------------------------------------------
   do k3=1,3
      do j3=1,3
         do i3=1,2
            do i=1,ncache
               if(ok(i,i3,j3  ,k3  ) .or. ok(i,i3,j3  ,k3-1) .or.  &
                  & ok(i,i3,j3-1,k3  ) .or. ok(i,i3,j3-1,k3-1))then
               emfx(i,i3,j3,k3)=0.0d0
            end if
         end do
      end do
   end do
end do
!-------------------------------------------------
! Reset current along direction y at refined edges
!-------------------------------------------------
do k3=1,3
   do j3=1,2
      do i3=1,3
         do i=1,ncache
            if(ok(i,i3  ,j3,k3  ) .or. ok(i,i3  ,j3,k3-1) .or.  &
               & ok(i,i3-1,j3,k3  ) .or. ok(i,i3-1,j3,k3-1))then
            emfy(i,i3,j3,k3)=0.0d0
         end if
      end do
   end do
end do
  end do
  !-------------------------------------------------
  ! Reset current along direction z at refined edges
  !-------------------------------------------------
  do k3=1,2
     do j3=1,3
        do i3=1,3
           do i=1,ncache
              if(ok(i,i3  ,j3  ,k3) .or. ok(i,i3  ,j3-1,k3) .or.  &
                 & ok(i,i3-1,j3  ,k3) .or. ok(i,i3-1,j3-1,k3))then
              emfz(i,i3,j3,k3)=0.0d0
           end if
        end do
     end do
  end do
  end do

  !------------------------------------
  ! Conservative update at level ilevel
  !------------------------------------
  do k3=1,2
     do j3=1,2
        do i3=1,2
           ind_son=i3+2*(j3-1)+4*(k3-1)
           iskip=ncoarse+(ind_son-1)*ngridmax
           do i=1,ncache
              ind_cell(i)=iskip+ind_grid(i)
                end do
                ! Update Bx using constraint transport
                do i=1,ncache
   !              if(son(ind_cell(i))==0)then
              
                 !JB
              !                if (.not. barotrop) then
              ! Remove emag to etot
                    !!$                    unew(ind_cell(i),5)=unew(ind_cell(i),5)-0.125d0*&
                    !!$                         &((unew(ind_cell(i),6)+unew(ind_cell(i),nvar+1))**2&
                    !!$                         &+(unew(ind_cell(i),7)+unew(ind_cell(i),nvar+2))**2&
                    !!$                         &+(unew(ind_cell(i),8)+unew(ind_cell(i),nvar+3))**2)
!                 end if

              ! Update Bx using constraint transport
              dflux_x=( emfy(i,i3,j3,k3)-emfy(i,i3,j3,k3+1) ) &
                   &    -( emfz(i,i3,j3,k3)-emfz(i,i3,j3+1,k3) )
              unew(ind_cell(i),6)=unew(ind_cell(i),6)+dflux_x*dtdiff/dx_loc
              dflux_x=( emfy(i,i3+1,j3,k3)-emfy(i,i3+1,j3,k3+1) ) &
                   &    -( emfz(i,i3+1,j3,k3)-emfz(i,i3+1,j3+1,k3) )  
              unew(ind_cell(i),nvar+1)=unew(ind_cell(i),nvar+1)+dflux_x*dtdiff/dx_loc
           end do
        end do
     end do
     end do
     
  do k3=1,2
     do j3=1,2
        do i3=1,2
           ind_son=i3+2*(j3-1)+4*(k3-1)
           iskip=ncoarse+(ind_son-1)*ngridmax
           do i=1,ncache
              ind_cell(i)=iskip+ind_grid(i)
           end do
           ! Update By using constraint transport
           do i=1,ncache
                 ! Update By using constraint transport
                 dflux_y=( emfz(i,i3,j3,k3)-emfz(i,i3+1,j3,k3) ) &
                    &    -( emfx(i,i3,j3,k3)-emfx(i,i3,j3,k3+1) )
                 unew(ind_cell(i),7)=unew(ind_cell(i),7)+dflux_y*dtdiff/dx_loc
                 dflux_y=( emfz(i,i3,j3+1,k3)-emfz(i,i3+1,j3+1,k3) ) &
                    &    -( emfx(i,i3,j3+1,k3)-emfx(i,i3,j3+1,k3+1) )
                 unew(ind_cell(i),nvar+2)=unew(ind_cell(i),nvar+2)+dflux_y*dtdiff/dx_loc
              end do
           end do
        end do
     end do

  do k3=1,2
     do j3=1,2
        do i3=1,2
           ind_son=i3+2*(j3-1)+4*(k3-1)
           iskip=ncoarse+(ind_son-1)*ngridmax
           do i=1,ncache
              ind_cell(i)=iskip+ind_grid(i)
           end do
           do i=1,ncache
              ! Update Bz using constraint transport
                 dflux_z=( emfx(i,i3,j3,k3)-emfx(i,i3,j3+1,k3) ) &
                    &    -( emfy(i,i3,j3,k3)-emfy(i,i3+1,j3,k3) )
                 unew(ind_cell(i),8)=unew(ind_cell(i),8)+dflux_z*dtdiff/dx_loc
                 dflux_z=( emfx(i,i3,j3,k3+1)-emfx(i,i3,j3+1,k3+1) ) &
                    &    -( emfy(i,i3,j3,k3+1)-emfy(i,i3+1,j3,k3+1) )
                 unew(ind_cell(i),nvar+3)=unew(ind_cell(i),nvar+3)+dflux_z*dtdiff/dx_loc
              end do
           end do
        end do
     end do
!if(radiative_nimhdheating)then
  do k3=1,2
     do j3=1,2
        do i3=1,2
           ind_son=i3+2*(j3-1)+4*(k3-1)
           iskip=ncoarse+(ind_son-1)*ngridmax
           do i=1,ncache
              ind_cell(i)=iskip+ind_grid(i)
           end do
           do i=1,ncache
!!$                 if (fld) then
                    ! update jcenter
                    do idim=1,3
                       unew(ind_cell(i),nvar-3+idim)=jcell(i,i3   ,j3   ,k3   ,idim)
                    end do

!!$                 else
!!$                    unew(ind_cell(i),nvar)=unew(ind_cell(i),nvar)+&
!!$                       &jxbsquare(i,i3,j3,k3) +jcentersquare(i,i3,j3,k3)
!!$                 end if
!print*,jcentersquare(i,i3,j3,k3),'toto'
                    !JB
                    !!$                    unew(ind_cell(i),5)=unew(ind_cell(i),5)+0.125d0*&
                    !!$                      &((unew(ind_cell(i),6)+unew(ind_cell(i),nvar+1))**2&
                    !!$                      &+(unew(ind_cell(i),7)+unew(ind_cell(i),nvar+2))**2&
                    !!$                      &+(unew(ind_cell(i),8)+unew(ind_cell(i),nvar+3))**2)


                    !!$               unew(ind_cell(i),5)=unew(ind_cell(i),5)-&
                    !!$                    &((flux(i,i3+1,j3,k3,1)-flux(i,i3,j3,k3,1))+&
                    !!$                    &(flux(i,i3,j3+1,k3,2)-flux(i,i3,j3,k3,2))+&
                    !!$                    &(flux(i,i3,j3,k3+1,3)-flux(i,i3,j3,k3,3)))*dtdiff/dx_loc
                    !            end if
                    !         
                    !            if (.not. barotrop) then


                    !!$                    unew(ind_cell(i),5) = unew(ind_cell(i),5) + jxbsquare(i,i3,j3,k3) + jcentersquare(i,i3,j3,k3)

                    !BEN

!!$                    unew(ind_cell(i),nvar)=unew(ind_cell(i),nvar)+&
!!$                       &jxbsquare(i,i3,j3,k3) +jcentersquare(i,i3,j3,k3)

                    !!$

                    !                    unew(ind_cell(i),5)=unew(ind_cell(i),5)-&
                    !                         &((flux(i,i3+1,j3,k3,1)-flux(i,i3,j3,k3,1))+&
                    !                         &(flux(i,i3,j3+1,k3,2)-flux(i,i3,j3,k3,2))+&
                    !                         &(flux(i,i3,j3,k3+1,3)-flux(i,i3,j3,k3,3)))*dtdiff/dx_loc
                    !                 end if
                    !!$          end do
!                 endif
!              end if
              end do!
              !$         end if
           end do
        end do
     end do
!  end if

  if(ilevel>levelmin)then

     !-----------------------------------------------------------
     ! Conservative update at level ilevel-1 for the Euler system
     !-----------------------------------------------------------
     !!$  ! Loop over dimensions
     !!$     do idim=1,ndim
     !!$        i0=0; j0=0; k0=0
     !!$        if(idim==1)i0=1
     !!$        if(idim==2)j0=1
     !!$        if(idim==3)k0=1
     !!$        
     !!$        !----------------------
     !!$        ! Left flux at boundary
     !!$        !----------------------     
     !!$        ! Check if grids sits near left boundary
     !!$        ! and gather neighbor father cells index
     !!$        nb_noneigh=0
     !!$        do i=1,ncache
     !!$           if (son(nbor(ind_grid(i),2*idim-1))==0) then
     !!$              nb_noneigh = nb_noneigh + 1
     !!$              ind_buffer(nb_noneigh) = nbor(ind_grid(i),2*idim-1)
     !!$              ind_cell(nb_noneigh) = i
     !!$           end if
     !!$        end do
     !!$        ! Conservative update of etot
     !!$        ! Loop over boundary cells
     !!$        do k3=k3min,k3max-k0
     !!$           do j3=j3min,j3max-j0
     !!$              do i3=i3min,i3max-i0
     !!$                 do i=1,nb_noneigh
     !!$                    unew(ind_buffer(i),5)=unew(ind_buffer(i),5) &
     !!$                         & -flux(ind_cell(i),i3,j3,k3,idim)*oneontwotondim
     !!$                 end do
     !!$              end do
     !!$           end do
     !!$        end do
     !!$     
     !!$        !-----------------------
     !!$        ! Right flux at boundary
     !!$        !-----------------------     
     !!$        ! Check if grids sits near right boundary
     !!$        ! and gather neighbor father cells index
     !!$        nb_noneigh=0
     !!$        do i=1,ncache
     !!$           if (son(nbor(ind_grid(i),2*idim))==0) then
     !!$              nb_noneigh = nb_noneigh + 1
     !!$              ind_buffer(nb_noneigh) = nbor(ind_grid(i),2*idim)
     !!$              ind_cell(nb_noneigh) = i
     !!$           end if
     !!$        end do
     !!$        ! Conservative update of new state variables
     !!$        ! Loop over boundary cells
     !!$        do k3=k3min+k0,k3max
     !!$           do j3=j3min+j0,j3max
     !!$              do i3=i3min+i0,i3max
     !!$                 do i=1,nb_noneigh
     !!$                    unew(ind_buffer(i),5)=unew(ind_buffer(i),5) &
     !!$                         & +flux(ind_cell(i),i3+i0,j3+j0,k3+k0,idim)*oneontwotondim
     !!$              end do
     !!$           end do
     !!$        end do
     !!$     end do
     !!$
     !!$  end do
     !!$  ! End loop over dimensions

     !--------------------------------------
     ! Conservative update at level ilevel-1
     !--------------------------------------
     i1=1; j1=1; k1=1

     !--------------------------------------
     ! Deal with 4 EMFx edges
     !--------------------------------------

     ! Update coarse By and Bz using fine EMFx on Y=0 and Z=0 grid edge
     ind_father1=1+(i1  )+3*(j1  )+9*(k1-1)
     ind_father2=1+(i1  )+3*(j1-1)+9*(k1-1)
     ind_father3=1+(i1  )+3*(j1-1)+9*(k1  )
     do i=1,ncache
        ind_buffer1=nbors_father_cells(i,ind_father1)
        ind_buffer2=nbors_father_cells(i,ind_father2)
        ind_buffer3=nbors_father_cells(i,ind_father3)
        weight=1.0d0
        if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
        if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5d0
        dflux=(emfx(i,1,1,1)+emfx(i,2,1,1))*0.25d0*weight*dtdiff/dx_loc
        unew(ind_buffer1,7)=unew(ind_buffer1,7)+dflux
        unew(ind_buffer2,nvar+2)=unew(ind_buffer2,nvar+2)+dflux
        unew(ind_buffer2,nvar+3)=unew(ind_buffer2,nvar+3)-dflux
        unew(ind_buffer3,8)=unew(ind_buffer3,8)-dflux
        if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
           unew(ind_buffer1,3+nvar)=unew(ind_buffer1,3+nvar)+dflux*0.5d0
           unew(ind_buffer3,2+nvar)=unew(ind_buffer3,2+nvar)-dflux*0.5d0
        endif
     end do

     ! Update coarse By and Bz using fine EMFx on Y=0 and Z=1 grid edge
     ind_father1=1+(i1  )+3*(j1-1)+9*(k1  )
     ind_father2=1+(i1  )+3*(j1-1)+9*(k1+1)
     ind_father3=1+(i1  )+3*(j1  )+9*(k1+1)
     do i=1,ncache
        ind_buffer1=nbors_father_cells(i,ind_father1)
        ind_buffer2=nbors_father_cells(i,ind_father2)
        ind_buffer3=nbors_father_cells(i,ind_father3)
        weight=1.0d0
        if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
        if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5d0
        dflux=(emfx(i,1,1,3)+emfx(i,2,1,3))*0.25d0*weight*dtdiff/dx_loc
        unew(ind_buffer1,nvar+3)=unew(ind_buffer1,nvar+3)-dflux
        unew(ind_buffer2,8)=unew(ind_buffer2,8)-dflux
        unew(ind_buffer2,nvar+2)=unew(ind_buffer2,nvar+2)-dflux
        unew(ind_buffer3,7)=unew(ind_buffer3,7)-dflux
        if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
           unew(ind_buffer1,2+nvar)=unew(ind_buffer1,2+nvar)+dflux*0.5d0
           unew(ind_buffer3,8)=unew(ind_buffer3,8)+dflux*0.5d0
        endif
     end do

     ! Update coarse By and Bz using fine EMFx on Y=1 and Z=1 grid edge
     ind_father1=1+(i1  )+3*(j1  )+9*(k1+1)
     ind_father2=1+(i1  )+3*(j1+1)+9*(k1+1)
     ind_father3=1+(i1  )+3*(j1+1)+9*(k1  )
     do i=1,ncache
        ind_buffer1=nbors_father_cells(i,ind_father1)
        ind_buffer2=nbors_father_cells(i,ind_father2)
        ind_buffer3=nbors_father_cells(i,ind_father3)
        weight=1.0d0
        if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
        if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5d0
        dflux=(emfx(i,1,3,3)+emfx(i,2,3,3))*0.25d0*weight*dtdiff/dx_loc
        unew(ind_buffer1,nvar+2)=unew(ind_buffer1,nvar+2)-dflux
        unew(ind_buffer2,7)=unew(ind_buffer2,7)-dflux
        unew(ind_buffer2,8)=unew(ind_buffer2,8)+dflux
        unew(ind_buffer3,nvar+3)=unew(ind_buffer3,nvar+3)+dflux
        if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
           unew(ind_buffer3,7)=unew(ind_buffer3,7)+dflux*0.5d0
           unew(ind_buffer1,8)=unew(ind_buffer1,8)-dflux*0.5d0
        endif
     end do

     ! Update coarse By and Bz using fine EMFx on Y=1 and Z=0 grid edge
     ind_father1=1+(i1  )+3*(j1+1)+9*(k1  )
     ind_father2=1+(i1  )+3*(j1+1)+9*(k1-1)
     ind_father3=1+(i1  )+3*(j1  )+9*(k1-1)
     do i=1,ncache
        ind_buffer1=nbors_father_cells(i,ind_father1)
        ind_buffer2=nbors_father_cells(i,ind_father2)
        ind_buffer3=nbors_father_cells(i,ind_father3)
        weight=1.0d0
        if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
        if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5d0
        dflux=(emfx(i,1,3,1)+emfx(i,2,3,1))*0.25d0*weight*dtdiff/dx_loc
        unew(ind_buffer1,8)=unew(ind_buffer1,8)+dflux
        unew(ind_buffer2,nvar+3)=unew(ind_buffer2,nvar+3)+dflux
        unew(ind_buffer2,7)=unew(ind_buffer2,7)+dflux
        unew(ind_buffer3,nvar+2)=unew(ind_buffer3,nvar+2)+dflux
        if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
           unew(ind_buffer3,3+nvar)=unew(ind_buffer3,3+nvar)-dflux*0.5d0
           unew(ind_buffer1,7)=unew(ind_buffer1,7)-dflux*0.5d0
        endif
     end do

     !--------------------------------------
     ! Deal with 4 EMFy edges
     !--------------------------------------

     ! Update coarse Bx and Bz using fine EMFy on X=0 and Z=0 grid edge
     ind_father1=1+(i1  )+3*(j1  )+9*(k1-1)
     ind_father2=1+(i1-1)+3*(j1  )+9*(k1-1)
     ind_father3=1+(i1-1)+3*(j1  )+9*(k1  )
     do i=1,ncache
        ind_buffer1=nbors_father_cells(i,ind_father1)
        ind_buffer2=nbors_father_cells(i,ind_father2)
        ind_buffer3=nbors_father_cells(i,ind_father3)
        weight=1.0d0
        if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
        if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5d0
        dflux=(emfy(i,1,1,1)+emfy(i,1,2,1))*0.25d0*weight*dtdiff/dx_loc
        unew(ind_buffer1,6)=unew(ind_buffer1,6)-dflux
        unew(ind_buffer2,nvar+1)=unew(ind_buffer2,nvar+1)-dflux
        unew(ind_buffer2,nvar+3)=unew(ind_buffer2,nvar+3)+dflux
        unew(ind_buffer3,8)=unew(ind_buffer3,8)+dflux
        if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
           unew(ind_buffer3,1+nvar)=unew(ind_buffer3,1+nvar)+dflux*0.5d0
           unew(ind_buffer1,3+nvar)=unew(ind_buffer1,3+nvar)-dflux*0.5d0
        endif
     end do

     ! Update coarse Bx and Bz using fine EMFy on X=0 and Z=1 grid edge
     ind_father1=1+(i1-1)+3*(j1  )+9*(k1  )
     ind_father2=1+(i1-1)+3*(j1  )+9*(k1+1)
     ind_father3=1+(i1  )+3*(j1  )+9*(k1+1)
     do i=1,ncache
        ind_buffer1=nbors_father_cells(i,ind_father1)
        ind_buffer2=nbors_father_cells(i,ind_father2)
        ind_buffer3=nbors_father_cells(i,ind_father3)
        weight=1.0d0
        if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
        if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5d0
        dflux=(emfy(i,1,1,3)+emfy(i,1,2,3))*0.25d0*weight*dtdiff/dx_loc
        unew(ind_buffer1,nvar+3)=unew(ind_buffer1,nvar+3)+dflux
        unew(ind_buffer2,8)=unew(ind_buffer2,8)+dflux
        unew(ind_buffer2,nvar+1)=unew(ind_buffer2,nvar+1)+dflux
        unew(ind_buffer3,6)=unew(ind_buffer3,6)+dflux
        if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
           unew(ind_buffer3,8)=unew(ind_buffer3,8)-dflux*0.5d0
           unew(ind_buffer1,1+nvar)=unew(ind_buffer1,1+nvar)-dflux*0.5d0
        endif
     end do

     ! Update coarse Bx and Bz using fine EMFy on X=1 and Z=1 grid edge
     ind_father1=1+(i1  )+3*(j1  )+9*(k1+1)
     ind_father2=1+(i1+1)+3*(j1  )+9*(k1+1)
     ind_father3=1+(i1+1)+3*(j1  )+9*(k1  )
     do i=1,ncache
        ind_buffer1=nbors_father_cells(i,ind_father1)
        ind_buffer2=nbors_father_cells(i,ind_father2)
        ind_buffer3=nbors_father_cells(i,ind_father3)
        weight=1.0d0
        if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
        if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5d0
        dflux=(emfy(i,3,1,3)+emfy(i,3,2,3))*0.25d0*weight*dtdiff/dx_loc
        unew(ind_buffer1,nvar+1)=unew(ind_buffer1,nvar+1)+dflux
        unew(ind_buffer2,6)=unew(ind_buffer2,6)+dflux
        unew(ind_buffer2,8)=unew(ind_buffer2,8)-dflux
        unew(ind_buffer3,nvar+3)=unew(ind_buffer3,nvar+3)-dflux
        if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
           unew(ind_buffer3,6)=unew(ind_buffer3,6)-dflux*0.5d0
           unew(ind_buffer1,8)=unew(ind_buffer1,8)+dflux*0.5d0
        endif
     end do

     ! Update coarse Bx and Bz using fine EMFx on X=1 and Z=0 grid edge
     ind_father1=1+(i1+1)+3*(j1  )+9*(k1  )
     ind_father2=1+(i1+1)+3*(j1  )+9*(k1-1)
     ind_father3=1+(i1  )+3*(j1  )+9*(k1-1)
     do i=1,ncache
        ind_buffer1=nbors_father_cells(i,ind_father1)
        ind_buffer2=nbors_father_cells(i,ind_father2)
        ind_buffer3=nbors_father_cells(i,ind_father3)
        weight=1.0d0
        if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
        if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5d0
        dflux=(emfy(i,3,1,1)+emfy(i,3,2,1))*0.25d0*weight*dtdiff/dx_loc
        unew(ind_buffer1,8)=unew(ind_buffer1,8)-dflux
        unew(ind_buffer2,nvar+3)=unew(ind_buffer2,nvar+3)-dflux
        unew(ind_buffer2,6)=unew(ind_buffer2,6)-dflux
        unew(ind_buffer3,nvar+1)=unew(ind_buffer3,nvar+1)-dflux
        if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
           unew(ind_buffer3,3+nvar)=unew(ind_buffer3,3+nvar)+dflux*0.5d0
           unew(ind_buffer1,6)=unew(ind_buffer1,6)+dflux*0.5d0  
        end if
     end do

     !--------------------------------------
     ! Deal with 4 EMFz edges
     !--------------------------------------

     ! Update coarse Bx and By using fine EMFz on X=0 and Y=0 grid edge
     ind_father1=1+(i1  )+3*(j1-1)+9*(k1  )
     ind_father2=1+(i1-1)+3*(j1-1)+9*(k1  )
     ind_father3=1+(i1-1)+3*(j1  )+9*(k1  )
     do i=1,ncache
        ind_buffer1=nbors_father_cells(i,ind_father1)
        ind_buffer2=nbors_father_cells(i,ind_father2)
        ind_buffer3=nbors_father_cells(i,ind_father3)
        weight=1.0d0
        if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
        if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5d0
        dflux=(emfz(i,1,1,1)+emfz(i,1,1,2))*0.25d0*weight*dtdiff/dx_loc
        unew(ind_buffer1,6)=unew(ind_buffer1,6)+dflux
        unew(ind_buffer2,nvar+1)=unew(ind_buffer2,nvar+1)+dflux
        unew(ind_buffer2,nvar+2)=unew(ind_buffer2,nvar+2)-dflux
        unew(ind_buffer3,7)=unew(ind_buffer3,7)-dflux
        if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
           unew(ind_buffer3,1+nvar)=unew(ind_buffer3,1+nvar)-dflux*0.5d0
           unew(ind_buffer1,2+nvar)=unew(ind_buffer1,2+nvar)+dflux*0.5d0
        endif
     end do

     ! Update coarse Bx and By using fine EMFz on X=0 and Y=1 grid edge
     ind_father1=1+(i1-1)+3*(j1  )+9*(k1  )
     ind_father2=1+(i1-1)+3*(j1+1)+9*(k1  )
     ind_father3=1+(i1  )+3*(j1+1)+9*(k1  )
     do i=1,ncache
        ind_buffer1=nbors_father_cells(i,ind_father1)
        ind_buffer2=nbors_father_cells(i,ind_father2)
        ind_buffer3=nbors_father_cells(i,ind_father3)
        weight=1.0d0
        if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
        if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5d0
        dflux=(emfz(i,1,3,1)+emfz(i,1,3,2))*0.25d0*weight*dtdiff/dx_loc
        unew(ind_buffer1,nvar+2)=unew(ind_buffer1,nvar+2)-dflux
        unew(ind_buffer2,7)=unew(ind_buffer2,7)-dflux
        unew(ind_buffer2,nvar+1)=unew(ind_buffer2,nvar+1)-dflux
        unew(ind_buffer3,6)=unew(ind_buffer3,6)-dflux
        if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
           unew(ind_buffer3,7)=unew(ind_buffer3,7)+dflux*0.5d0
           unew(ind_buffer1,1+nvar)=unew(ind_buffer1,1+nvar)+dflux*0.5d0
        endif
     end do

     ! Update coarse Bx and By using fine EMFz on X=1 and Y=1 grid edge
     ind_father1=1+(i1  )+3*(j1+1)+9*(k1  )
     ind_father2=1+(i1+1)+3*(j1+1)+9*(k1  )
     ind_father3=1+(i1+1)+3*(j1  )+9*(k1  )
     do i=1,ncache
        ind_buffer1=nbors_father_cells(i,ind_father1)
        ind_buffer2=nbors_father_cells(i,ind_father2)
        ind_buffer3=nbors_father_cells(i,ind_father3)
        weight=1.0d0
        if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
        if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5d0
        dflux=(emfz(i,3,3,1)+emfz(i,3,3,2))*0.25d0*weight*dtdiff/dx_loc
        unew(ind_buffer1,nvar+1)=unew(ind_buffer1,nvar+1)-dflux
        unew(ind_buffer2,6)=unew(ind_buffer2,6)-dflux
        unew(ind_buffer2,7)=unew(ind_buffer2,7)+dflux
        unew(ind_buffer3,nvar+2)=unew(ind_buffer3,nvar+2)+dflux
        if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
           unew(ind_buffer3,6)=unew(ind_buffer3,6)+dflux*0.5d0
           unew(ind_buffer1,7)=unew(ind_buffer1,7)-dflux*0.5d0
        endif
     end do

     ! Update coarse Bx and By using fine EMFz on X=1 and Y=0 grid edge
     ind_father1=1+(i1+1)+3*(j1  )+9*(k1  )
     ind_father2=1+(i1+1)+3*(j1-1)+9*(k1  )
     ind_father3=1+(i1  )+3*(j1-1)+9*(k1  )
     do i=1,ncache
        ind_buffer1=nbors_father_cells(i,ind_father1)
        ind_buffer2=nbors_father_cells(i,ind_father2)
        ind_buffer3=nbors_father_cells(i,ind_father3)
        weight=1.0d0
        if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
        if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5d0
        dflux=(emfz(i,3,1,1)+emfz(i,3,1,2))*0.25d0*weight*dtdiff/dx_loc
        unew(ind_buffer1,7)=unew(ind_buffer1,7)+dflux
        unew(ind_buffer2,nvar+2)=unew(ind_buffer2,nvar+2)+dflux
        unew(ind_buffer2,6)=unew(ind_buffer2,6)+dflux
        unew(ind_buffer3,nvar+1)=unew(ind_buffer3,nvar+1)+dflux
        if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
           unew(ind_buffer3,2+nvar)=unew(ind_buffer3,2+nvar)-dflux*0.5d0
           unew(ind_buffer1,6)=unew(ind_buffer1,6)-dflux*0.5d0
        endif
     end do

  endif

end subroutine diffine1_sts_dtu
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmp_current_sts_dtu(u,Ex_arete,Ey_arete,Ez_arete,fluxni, &
      &ngrid,dx,dy,dz,dt,jxbsquare,jcentersquare,ok,jcell)
   use amr_parameters
   use const
   use hydro_parameters
   implicit none
   real(dp) :: dBx_arete_dy,dBx_arete_dz
   real(dp) :: dBy_arete_dx,dBy_arete_dz
   real(dp) :: dBz_arete_dx,dBz_arete_dy
   integer  :: ic,i,j,k,im1,jm1,km1
   integer  :: ii,jj,kk

   real(dp) :: etaohmdiss, etaodx, etaody, etaodz
   real(dp) :: pressurex, pressurey, pressurez
   real(dp) :: rhox, rhoy, rhoz

   integer :: ngrid,ivar
   real(dp)::dx,dy,dz,dt
   real(dp),dimension(1:nvector,1:3,1:3,1:3) :: Ex_arete
   real(dp),dimension(1:nvector,1:3,1:3,1:3) :: Ey_arete
   real(dp),dimension(1:nvector,1:3,1:3,1:3) :: Ez_arete

   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::u
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),save::q
   logical ,dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::ok

   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3),save::bemfx,bemfy,bemfz
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3),save::jemfx,jemfy,jemfz
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3),save::florentzx,florentzy,florentzz
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3),save::fluxmd,fluxh,fluxad
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3),save::emfambdiff,fluxambdiff
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3),save::emfohmdiss,fluxohm
   real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:3)::fluxni

   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::jcentersquare,jxbsquare
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3,1:3),save::bmagij

   integer::ilo,ihi,jlo,jhi,klo,khi

  ! Output courant vector in the cell
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:3)::jcell
jcell=0.0d0

   emfohmdiss=0.0d0
   fluxohm=0.0d0
   emfambdiff=0.0d0
   fluxambdiff=0.0d0

   bmagij=0d0
   emfambdiff=0d0
   fluxambdiff=0d0
   emfohmdiss=0d0
   fluxohm=0d0
   jcentersquare=0d0
   jxbsquare=0d0
   fluxmd=0d0
   fluxh=0d0
   fluxad=0d0

   fluxni=0.0d0
   Ex_arete=0.0d0 ; Ey_arete=0.0d0 ; Ez_arete=0.0d0 

   bemfx=0d0
   bemfy=0d0
   bemfz=0d0
   jemfx=0d0
   jemfy=0d0
   jemfz=0d0
   florentzx=0d0
   florentzy=0d0
   florentzz=0d0

   q=0.0d0

   call ctoprim_sts(u,q,ngrid)
   if((nambipolar2.eq.1).or.(nmagdiffu2.eq.1)) then
      call computejb2(u,q,ngrid,dx,dy,dz,dt,bemfx,bemfy,bemfz,jemfx,jemfy,jemfz,bmagij,florentzx,florentzy,florentzz,fluxmd,fluxh,fluxad,jcell)
      !     call computejb(u,q,ngrid,dx,dy,dz,dt,bemfx,bemfy,bemfz,jemfx,jemfy,jemfz,bmagij,florentzx,florentzy,florentzz,fluxmd,fluxh,fluxad)
   end if
   if (nmagdiffu2.eq.1) then
      call computdifmag(u,q,ngrid,dx,dy,dz,dt,bemfx,bemfy,bemfz,jemfx,jemfy,jemfz,bmagij,fluxmd,emfohmdiss,fluxohm,jcentersquare)
   endif
   if (nambipolar2.eq.1)  then
      call computambip(u,q,ngrid,dx,dy,dz,dt,bemfx,bemfy,bemfz,florentzx,florentzy,florentzz,fluxad,bmagij,emfambdiff,fluxambdiff,jxbsquare)
   end if

   !ben
   ilo=MIN(1,iu1+2); ihi=MAX(1,iu2-2)
   jlo=MIN(1,ju1+2); jhi=MAX(1,ju2-2)
   klo=MIN(1,ku1+2); khi=MAX(1,ku2-2)

   ! Aretes paralleles a l'axe des x
   do k=kf1,kf2
      do j=jf1,jf2
         do i=ilo,ihi
            do ic=1,ngrid
               Ex_arete(ic,i,j,k)=emfohmdiss(ic,i,j,k,1)+emfambdiff(ic,i,j,k,1)
               fluxni(ic,i,j,k,1)=fluxohm(ic,i,j,k,1)+fluxambdiff(ic,i,j,k,1)
               fluxni(ic,3,j,k,1)=fluxohm(ic,3,j,k,1)+fluxambdiff(ic,3,j,k,1)
            enddo
         enddo
      enddo
   enddo

#if NDIM>1
   ! Aretes paralleles a l'axe des y
   do k=kf1,kf2
      do j=jlo,jhi
         do i=if1,if2
            do ic=1,ngrid
               Ey_arete(ic,i,j,k)=emfohmdiss(ic,i,j,k,2)+emfambdiff(ic,i,j,k,2)
               fluxni(ic,i,j,k,2)=fluxohm(ic,i,j,k,2)+fluxambdiff(ic,i,j,k,2)
               fluxni(ic,i,3,k,2)=fluxohm(ic,i,3,k,2)+fluxambdiff(ic,i,3,k,2)
            enddo
         enddo
      enddo
   enddo

   ! Aretes paralleles a l'axe des z
   do k=klo,khi
      do j=jf1,jf2
         do i=if1,if2
            do ic=1,ngrid
               Ez_arete(ic,i,j,k)=emfohmdiss(ic,i,j,k,3)+emfambdiff(ic,i,j,k,3)
               fluxni(ic,i,j,k,3)=fluxohm(ic,i,j,k,3)+fluxambdiff(ic,i,j,k,3)
               fluxni(ic,i,j,3,3)=fluxohm(ic,i,j,3,3)+fluxambdiff(ic,i,j,3,3)
            enddo
         enddo
      enddo
   enddo
#endif

end subroutine cmp_current_sts_dtu
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine ctoprim_sts(uin,q,ngrid)
   use amr_parameters
   use hydro_parameters
   use const
   implicit none

   integer ::ngrid
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::uin
   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q  

   integer ::i, j, k, l, n, idim
   real(dp)::eint, smalle, smallp, etot
   real(dp),dimension(1:nvector),save::eken,emag

   ! EOS
   real(dp)  :: pp_eos

   real(dp):: small_er=1.0d-30       ! minimum rad energy inside frequency group in code units

   smalle = smallc**2/gamma/(gamma-one)
   smallp = smallr*smallc**2/gamma

   ! Convert to primitive variable
   do k = ku1, ku2
      do j = ju1, ju2
         do i = iu1, iu2

            ! Compute density
            do l = 1, ngrid
               q(l,i,j,k,1) = max(uin(l,i,j,k,1),smallr)
            end do

            ! Debug
            if(debug)then
               do l = 1, ngrid
                  if(uin(l,i,j,k,1).le.smallr)then
                     write(*,*)'negative density'
                     write(*,*)uin(l,i,j,k,1)
                     stop
                  end if
               end do
            end if

            ! Compute velocities
            do l = 1, ngrid
               q(l,i,j,k,2) = uin(l,i,j,k,2)/uin(l,i,j,k,1)
               q(l,i,j,k,3) = uin(l,i,j,k,3)/uin(l,i,j,k,1)
               q(l,i,j,k,4) = uin(l,i,j,k,4)/uin(l,i,j,k,1)
            end do

            ! Compute cell centered magnetic field
            DO l = 1, ngrid
               q(l,i,j,k,6) = (uin(l,i,j,k,6)+uin(l,i,j,k,nvar+1))*half
               q(l,i,j,k,7) = (uin(l,i,j,k,7)+uin(l,i,j,k,nvar+2))*half
               q(l,i,j,k,8) = (uin(l,i,j,k,8)+uin(l,i,j,k,nvar+3))*half
            END DO

#if NENER>0
            ! Compute radiative energy
            do l = 1, ngrid
               erad_loc(l)=0.0d0
               do n = 1, nener
                  q(l,i,j,k,8+n) = max(uin(l,i,j,k,8+n),small_er)
                  erad_loc(l) = erad_loc(l) + q(l,i,j,k,8+n)
               end do
            end do
#endif

            ! Compute specific kinetic energy and magnetic energy
            do l = 1, ngrid
               eken(l) = half*(q(l,i,j,k,2)**2+q(l,i,j,k,3)**2+q(l,i,j,k,4)**2)
               emag(l) = half*(q(l,i,j,k,6)**2+q(l,i,j,k,7)**2+q(l,i,j,k,8)**2)
            end do

            ! Compute thermal pressure through EOS
            do l = 1, ngrid
               etot = uin(l,i,j,k,5) - emag(l)
#if NENER>0
               etot = etot - erad_loc(l)
#endif
               eint = etot/uin(l,i,j,k,1)-eken(l)
!!$               if(eos)then
!!$                  !                  eint = uin(l,i,j,k,nvar)!eint*q(l,i,j,k,1)   ! volumic 
                  eint = eint*q(l,i,j,k,1)   ! volumic 

                  call pressure_eos(q(l,i,j,k,1),eint,pp_eos)
                  q(l,i,j,k,5)=MAX(pp_eos,smallp)
!!$               else
!!!!BENOIT
!                  q(l,i,j,k,5)=MAX((gamma-one)*q(l,i,j,k,1)*eint,smallp)
!!$                  !                 q(l,i,j,k,5)=MAX((gamma-one)*eint,smallp)
!!$               endif
            end do

         end do
      end do
   end do
   ! Passive scalar (and extinction and internal energy and rad fluxes in M1) !!!!!!
#if NVAR>8+NENER
  do n = 9+nener, nvar
     do k = ku1, ku2
        do j = ju1, ju2
           do i = iu1, iu2
              do l = 1, ngrid
                 q(l,i,j,k,n) = uin(l,i,j,k,n)/q(l,i,j,k,1)
              end do
           end do
        end do
     end do
  end do
#endif

 end subroutine ctoprim_sts
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine interpol_hydro_sts(u1,ind1,u2,nn)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::nn,n
  real(dp),dimension(1:nvector,0:twondim  ,1:nvar+3)::u1
  integer ,dimension(1:nvector,0:twondim)           ::ind1
  real(dp),dimension(1:nvector,1:twotondim,1:nvar+3)::u2
  !----------------------------------------------------------
  ! This routine performs a prolongation without interpolation
  ! operation for buffer cells 
  ! The interpolated variables are:
  ! interpol_var=0: rho, rho u and E
  ! interpol_var=1: rho, rho u and rho epsilon
  ! The interpolation method is:
  ! interpol_type=0 straight injection
  !----------------------------------------------------------
  integer::i,j,ivar,idim,ind,ix,iy,iz,neul=5

  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,0:twondim),save::a
  real(dp),dimension(1:nvector,1:ndim),save::w
  real(dp),dimension(1:nvector),save::ekin,emag,erad_loc
  real(dp),dimension(1:nvector,0:twondim  ,1:6),save::B1
  real(dp),dimension(1:nvector,1:twotondim,1:6),save::B2

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)
  end do

  ! If necessary, convert father total energy into internal energy
  if(interpol_var==1)then
     do j=0,twondim
        ekin(1:nn)=0.0d0
        do idim=1,3
           do i=1,nn
              ekin(i)=ekin(i)+0.5d0*u1(i,j,idim+1)**2/u1(i,j,1)
           end do
        end do
        emag(1:nn)=0.0d0
        do idim=1,3
           do i=1,nn
              emag(i)=emag(i)+0.125d0*(u1(i,j,idim+neul)+u1(i,j,idim+nvar))**2
           end do
        end do
        erad_loc(1:nn)=0.0d0
#if NENER>0
        do n=1,nener
           do i=1,nn
              erad_loc(i)=erad_loc(i) + u1(i,j,8+n)
           end do
        end do
#endif
        do i=1,nn
           u1(i,j,neul)=u1(i,j,neul)-ekin(i)-emag(i)-erad_loc(i)
        end do
     end do
  end if

  !------------------------------------------------
  ! Loop over cell-centered interpolation variables
  !------------------------------------------------
  do ivar=1,nvar
  if(ivar<=neul.or.ivar>neul+ndim)then

     ! Load father variable
     do j=0,twondim
        do i=1,nn 
           a(i,j)=u1(i,j,ivar)
        end do
     end do

     ! Reset gradient
     w(1:nn,1:ndim)=0.0D0

!BEN
!!$     ! Compute gradient with chosen limiter
!!$     if(interpol_type==1)call compute_limiter_minmod(a,w,nn)
!!$     if(interpol_type==2)call compute_limiter_central(a,w,nn)
!!$     if(interpol_type==3)call compute_central(a,w,nn)

     ! Interpolate over children cells
     do ind=1,twotondim
        u2(1:nn,ind,ivar)=a(1:nn,0)
        do idim=1,ndim
           do i=1,nn
              u2(i,ind,ivar)=u2(i,ind,ivar)+w(i,idim)*xc(ind,idim)
           end do
        end do
     end do

  end if
  end do
  ! End loop over cell-centered variables

  ! Update cell centered magnetic field also in redundant array
#if NDIM<2
  do ind=1,twotondim
     do i=1,nn
        u2(i,ind,2+nvar)=u2(i,ind,2+neul)
     end do
  end do
#endif
#if NDIM<3
  do ind=1,twotondim
     do i=1,nn
        u2(i,ind,3+nvar)=u2(i,ind,3+neul)
     end do
  end do
#endif

  !------------------------------------------------
  ! Loop over face-centered interpolation variables
  !------------------------------------------------
  do j=0,twondim
     do i=1,nn
        B1(i,j,1)=u1(i,j,neul+1)
        B1(i,j,2)=u1(i,j,neul+2)
        B1(i,j,3)=u1(i,j,neul+3)
        B1(i,j,4)=u1(i,j,nvar+1)
        B1(i,j,5)=u1(i,j,nvar+2)
        B1(i,j,6)=u1(i,j,nvar+3)
     end do
  end do
  call interpol_mag_sts(B1,ind1,B2,nn)
  do ind=1,twotondim
     do i=1,nn
        u2(i,ind,neul+1)=B2(i,ind,1)
        u2(i,ind,nvar+1)=B2(i,ind,4)
#if NDIM>1        
        u2(i,ind,neul+2)=B2(i,ind,2)
        u2(i,ind,nvar+2)=B2(i,ind,5)
#endif
#if NDIM>2
        u2(i,ind,neul+3)=B2(i,ind,3)
        u2(i,ind,nvar+3)=B2(i,ind,6)
#endif
     end do
  end do

  ! If necessary, convert children internal energy into total energy
  if(interpol_var==1)then
     do ind=1,twotondim
        ekin(1:nn)=0.0d0
        do idim=1,3
           do i=1,nn
              ekin(i)=ekin(i)+0.5d0*u2(i,ind,idim+1)**2/u2(i,ind,1)
           end do
        end do
        emag(1:nn)=0.0d0
        do idim=1,3
           do i=1,nn
              emag(i)=emag(i)+0.125d0*(u2(i,ind,idim+neul)+u2(i,ind,idim+nvar))**2
           end do
        end do
        erad_loc(1:nn)=0.0d0
#if NENER>0 
       do n=1,nener
           do i=1,nn
              erad_loc(i)=erad_loc(i) + u2(i,ind,8+n)
           end do
        end do
#endif
        do i=1,nn
           u2(i,ind,neul)=u2(i,ind,neul)+ekin(i)+emag(i)+erad_loc(i)
        end do
     end do
  end if

end subroutine interpol_hydro_sts
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine interpol_mag_sts(B1,ind1,B2,nn)
  use amr_commons
  implicit none
  integer::nn
  real(dp),dimension(1:nvector,0:twondim  ,1:6)::B1
  integer ,dimension(1:nvector,0:twondim)      ::ind1
  real(dp),dimension(1:nvector,1:twotondim,1:6)::B2
  !----------------------------------------------------------
  ! This routine performs a prolongation (interpolation)
  ! operation for newly refined cells or buffer cells.
  ! The interpolated variables are Bx, By and Bz.
  ! Divergence free is garanteed.
  ! The scheme is the one invented by Toth and Balsara.
  ! interpol_mag_type=0: straight injection
  ! interpol_mag_type=1: linear interpolation with MinMod slope
  ! interpol_mag_type=2: linear interpolation with Monotonized Central slope
  ! interpol_mag_type=3: linear interpolation without limiters
  !----------------------------------------------------------
  integer::i,j,k,ind,l,idim,imax,jmax,kmax
  real(dp),dimension(1:nvector,-1:1,0:1,0:1),save::u
  real(dp),dimension(1:nvector,0:1,-1:1,0:1),save::v
  real(dp),dimension(1:nvector,0:1,0:1,-1:1),save::w

  imax=1; jmax=0; kmax=0
#if NDIM>1
  jmax=1
#endif
#if NDIM>2
  kmax=1
#endif

  ! Compute interpolated fine B over coarse side faces
  call interpol_faces_sts(B1,u,v,w,nn)
  
  ! Get fine B from refined faces, if any
  call copy_from_refined_faces(B1,ind1,u,v,w,nn)
 
  ! Compute interpolated fine B inside coarse cell.
  call cmp_central_faces(u,v,w,nn)

  ! Scatter results
  do i=0,imax
  do j=0,jmax
  do k=0,kmax
     ind=1+i+2*j+4*k
     do l=1,nn
        B2(l,ind,1)=u(l,i-1,j,k)
        B2(l,ind,2)=v(l,i,j-1,k)
        B2(l,ind,3)=w(l,i,j,k-1)
        B2(l,ind,4)=u(l,i,j,k)
        B2(l,ind,5)=v(l,i,j,k)
        B2(l,ind,6)=w(l,i,j,k)
     end do
  end do
  end do
  end do

end subroutine interpol_mag_sts
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine interpol_faces_sts(b1,u,v,w,nn)
  use amr_commons
  use hydro_commons, ONLY: interpol_mag_type
  implicit none
  integer::nn
  real(dp),dimension(1:nvector,0:twondim,1:6)::b1
  real(dp),dimension(1:nvector,-1:1,0:1,0:1)::u
  real(dp),dimension(1:nvector,0:1,-1:1,0:1)::v
  real(dp),dimension(1:nvector,0:1,0:1,-1:1)::w

  ! TVD interpolation from coarse faces
  integer::i,j,k,l,imax,jmax,kmax
  real(dp),dimension(1:nvector,0:4),save::b
  real(dp),dimension(1:nvector,1:2),save::s

  imax=1; jmax=0; kmax=0
#if NDIM>1
  jmax=1
#endif
#if NDIM>2
  kmax=1
#endif

  ! Left face along direction x (interpolate Bx)
  do l=1,nn
     b(l,0)=b1(l,0,1)
  end do
#if NDIM>1     
  do l=1,nn
     b(l,1)=b1(l,3,1)
     b(l,2)=b1(l,4,1)
  end do
#endif
#if NDIM>2
  do l=1,nn
     b(l,3)=b1(l,5,1)
     b(l,4)=b1(l,6,1)
  end do
#endif

  s(1:nn,1:2)=0.0
!BEN
!!$#if NDIM==2
!!$  if(interpol_mag_type>0)call compute_1d_tvd(b,s,nn)
!!$#endif
!!$#if NDIM==3
!!$  if(interpol_mag_type>0)call compute_2d_tvd(b,s,nn)
!!$#endif
  do j=0,jmax
  do k=0,kmax
     do l=1,nn
        u(l,-1,j,k)=b(l,0)+0.5*s(l,1)*(dble(j)-0.5)+0.5*s(l,2)*(dble(k)-0.5)
     end do
  end do
  end do

  ! Right face along direction x (interpolate Bx)
  do l=1,nn
     b(l,0)=b1(l,0,4)
  end do
#if NDIM>1     
  do l=1,nn
     b(l,1)=b1(l,3,4)
     b(l,2)=b1(l,4,4)
  end do
#endif
#if NDIM>2
  do l=1,nn
     b(l,3)=b1(l,5,4)
     b(l,4)=b1(l,6,4)
  end do
#endif

  s(1:nn,1:2)=0.0
!BEN
!!$#if NDIM==2
!!$  if(interpol_mag_type>0)call compute_1d_tvd(b,s,nn)
!!$#endif
!!$#if NDIM==3
!!$  if(interpol_mag_type>0)call compute_2d_tvd(b,s,nn)
!!$#endif
  do j=0,jmax
  do k=0,kmax
     do l=1,nn
        u(l,+1,j,k)=b(l,0)+0.5*s(l,1)*(dble(j)-0.5)+0.5*s(l,2)*(dble(k)-0.5)
     end do
  end do
  end do

#if NDIM>1
  ! Left face along direction y (interpolate By)
  do l=1,nn
     b(l,0)=b1(l,0,2)
  end do
  do l=1,nn
     b(l,1)=b1(l,1,2)
     b(l,2)=b1(l,2,2)
  end do
#if NDIM>2     
  do l=1,nn
     b(l,3)=b1(l,5,2)
     b(l,4)=b1(l,6,2)
  end do
#endif

  s(1:nn,1:2)=0.0
!BEN
!!$#if NDIM==2
!!$  if(interpol_mag_type>0)call compute_1d_tvd(b,s,nn)
!!$#endif
!!$#if NDIM==3
!!$  if(interpol_mag_type>0)call compute_2d_tvd(b,s,nn)
!!$#endif
  do i=0,imax
  do k=0,kmax
     do l=1,nn
        v(l,i,-1,k)=b(l,0)+0.5*s(l,1)*(dble(i)-0.5)+0.5*s(l,2)*(dble(k)-0.5)
     end do
  end do
  end do

  ! Right face along direction y (interpolate By)
  do l=1,nn
     b(l,0)=b1(l,0,5)
  end do
  do l=1,nn
     b(l,1)=b1(l,1,5)
     b(l,2)=b1(l,2,5)
  end do
#if NDIM>2     
  do l=1,nn
     b(l,3)=b1(l,5,5)
     b(l,4)=b1(l,6,5)
  end do
#endif

  s(1:nn,1:2)=0.0
!BEN
!!$#if NDIM==2
!!$  if(interpol_mag_type>0)call compute_1d_tvd(b,s,nn)
!!$#endif
!!$#if NDIM==3
!!$  if(interpol_mag_type>0)call compute_2d_tvd(b,s,nn)
!!$#endif
  do i=0,imax
  do k=0,kmax
     do l=1,nn
        v(l,i,+1,k)=b(l,0)+0.5*s(l,1)*(dble(i)-0.5)+0.5*s(l,2)*(dble(k)-0.5)
     end do
  end do
  end do
#endif

#if NDIM>2
  ! Left face along direction z (interpolate Bz)
  do l=1,nn
     b(l,0)=b1(l,0,3)
     b(l,1)=b1(l,1,3)
     b(l,2)=b1(l,2,3)
     b(l,3)=b1(l,3,3)
     b(l,4)=b1(l,4,3)
  end do

  s(1:nn,1:2)=0.0
!!$BEN  if(interpol_mag_type>0)call compute_2d_tvd(b,s,nn)
  do i=0,1
     do j=0,1
        do l=1,nn
           w(l,i,j,-1)=b(l,0)+0.5*s(l,1)*(dble(i)-0.5)+0.5*s(l,2)*(dble(j)-0.5)
        end do
     end do
  end do

  ! Right face along direction z (interpolate Bz)
  do l=1,nn
     b(l,0)=b1(l,0,6)
     b(l,1)=b1(l,1,6)
     b(l,2)=b1(l,2,6)
     b(l,3)=b1(l,3,6)
     b(l,4)=b1(l,4,6)
  end do

  s(1:nn,1:2)=0.0
!!$BEN  if(interpol_mag_type>0)call compute_2d_tvd(b,s,nn)
  do i=0,1
     do j=0,1
        do l=1,nn
           w(l,i,j,+1)=b(l,0)+0.5*s(l,1)*(dble(i)-0.5)+0.5*s(l,2)*(dble(j)-0.5)
        end do
     end do
  end do
#endif

end subroutine interpol_faces_sts
!###########################################################
!###########################################################
!###########################################################
!###########################################################