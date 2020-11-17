!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine diffusion_sts(ilevel,nsub)
   use amr_commons
   use hydro_commons
   use mpi_mod
   implicit none
   integer,intent(IN)::ilevel,Nsub
   integer::icycle,nsubdiff,ivar,i,iskip,ind,info,jcycle,ifirst,isafe
   real(dp)::dx,scale,dx_loc,dtdiffohm,dtdiffamb,norm,nu,dtohm,dtad
   real(dp)::dtdiff,dtdiffold,dtdiffsum,dtdiffcoef,sommeinit,dtloc,dtall
   real(dp),dimension(1:3)::skip_loc
   real(dp)::max_diff_loc,max_diff_all
   logical::sts
   if(numbtot(1,ilevel)==0)return 
   ! Rescaling factors
   dx=0.5D0**ilevel
   skip_loc=(/0.0d0,0.0d0,0.0d0/)
   if(ndim>0)skip_loc(1)=dble(icoarse_min)
   if(ndim>1)skip_loc(2)=dble(jcoarse_min)
   if(ndim>2)skip_loc(3)=dble(kcoarse_min)
   scale=dble(icoarse_max-icoarse_min+1)/boxlen
   dx_loc=dx/scale

   call set_unew_sts(ilevel,0,dtohm,dtad)  ! needed ? probably yes
#ifndef WITHOUTMPI
   call MPI_ALLREDUCE(dtloc     ,dtall      ,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
      &MPI_COMM_WORLD,info)
   dtloc=dtall
#endif
   dtdiffohm=dtmagdiff(ilevel)
   dtdiffamb=dtambdiff(ilevel)
!!$   dtdiff=min(dtdiffohm,dtdiffamb)
   dtdiff=1.0d0/dtdiffohm+1.0d0/dtdiffamb
   dtdiff=1.0d0/dtdiff
   dtdiffold=min(dtdiff,dtnew(ilevel))
   dtdiffsum=0.0d0
   !nu=0
   nu = nu_sts!1d0/(4d0*nsts(ilevel)**2)

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

   call set_unew_sts(ilevel,0,dtohm,dtad)  ! needed ? probably yes
#ifndef WITHOUTMPI
   call MPI_ALLREDUCE(dtloc     ,dtall      ,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
        &MPI_COMM_WORLD,info)
   dtloc=dtall
#endif

      dtdiffohm=dtmagdiff(ilevel)
      dtdiffamb=dtambdiff(ilevel)
      !!$      dtdiffamb=dtloc!dtambdiff(ilevel)
!!$      dtdiff=min(dtdiffohm,dtdiffamb)
      dtdiff=1.0d0/dtdiffohm+1.0d0/dtdiffamb
      dtdiff=1.0d0/dtdiff

      dtdiffold=dtdiff!min(dtdiffohm,dtdiffamb)

      call make_virtual_fine_dp(unew(1,1),ilevel)
      call make_virtual_fine_dp(unew(1,5),ilevel)
      call make_virtual_fine_dp(unew(1,6),ilevel)
      call make_virtual_fine_dp(unew(1,7),ilevel)
      call make_virtual_fine_dp(unew(1,8),ilevel)
      call make_virtual_fine_dp(unew(1,nvar+1),ilevel)
      call make_virtual_fine_dp(unew(1,nvar+2),ilevel)
      call make_virtual_fine_dp(unew(1,nvar+3),ilevel)
      call make_virtual_fine_dp(unew(1,nvar),ilevel)
      !      call make_virtual_fine_dp(enew(1),ilevel)

      nsts(ilevel)=min(100,floor(sqrt(dtnew(ilevel)/dtdiffold))+1)
      nsubdiff=max(1,nsts(ilevel))
      dtsts(ilevel) = dtdiffold*nsts(ilevel)/(2.*sqrt(nu_sts))*((1.+sqrt(nu_sts))**(2.*nsts(ilevel))-(1.-sqrt(nu_sts))**(2.*nsts(ilevel))) &         ! somme des dtsts, theorique
             & /((1.+sqrt(nu_sts))**(2.*nsts(ilevel))+(1.-sqrt(nu_sts))**(2.*nsts(ilevel)))
      dtdiffcoef= dtnew(ilevel)/dtsts(ilevel)

      do  icycle=1,nsubdiff
         if (nsts(ilevel) > 0) then
            if(myid==1 .and. (mod(nstep,ncontrol)==0))write(*,*)nu,'dtdiffcoeff',dtdiffcoef,'sommeinit',sommeinit,'dtdiffsum',dtdiffsum,'dtdiffold',dtdiffold
!!$            dtdiff=dtdiffold*((-1d0+nu)*cos((2d0*dble(icycle)-1d0)*acos(-1d0)/(2d0*nsubdiff))+1d0+nu)**(-1)    ! pas de temps du STS
            dtdiff=dtdiffcoef*dtdiffold*((-1d0+nu)*cos((2d0*dble(icycle)-1d0)*acos(-1d0)/(2d0*nsubdiff))+1d0+nu)**(-1)    ! pas de temps du STS
         else if (nsts(ilevel) == 0) then
            dtdiff = dtnew(ilevel)
         end if

         dtdiffsum=dtdiffsum+dtdiff
         if(myid==1 .and. (mod(nstep,ncontrol)==0))write(*,*)'subcycling through STS',icycle,'ilevel',ilevel,'nsubdiff :',nsubdiff,&
            'dtmagdiff :',dtmagdiff(ilevel),'dtambdiff :',dtambdiff(ilevel),&
            'dtnew :',dtnew(ilevel),'dtdiff :',dtdiff,'sum :',dtdiffsum, 'jcycle :',jcycle

         call set_unew_sts(ilevel,1,dtohm,dtad)  ! needed ? probably yes

         call diffusion_fine_sts(ilevel,dtdiff,icycle)
#ifndef WITHOUTMPI
         call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif

!         call make_virtual_reverse_dp(unew(1,5),ilevel)
         call make_virtual_reverse_dp(unew(1,6),ilevel)
         call make_virtual_reverse_dp(unew(1,7),ilevel)
         call make_virtual_reverse_dp(unew(1,8),ilevel)
         call make_virtual_reverse_dp(unew(1,nvar+1),ilevel)
         call make_virtual_reverse_dp(unew(1,nvar+2),ilevel)
         call make_virtual_reverse_dp(unew(1,nvar+3),ilevel)
!         call make_virtual_reverse_dp(unew(1,nvar),ilevel)

         call set_uold_sts(ilevel,0,dtdiff)
         !            call upload_fine(ilevel)!

!         call make_virtual_fine_dp(uold(1,1),ilevel)
!         call make_virtual_fine_dp(uold(1,5),ilevel)
         call make_virtual_fine_dp(uold(1,6),ilevel)
         call make_virtual_fine_dp(uold(1,7),ilevel)
         call make_virtual_fine_dp(uold(1,8),ilevel)
         call make_virtual_fine_dp(uold(1,nvar+1),ilevel)
         call make_virtual_fine_dp(uold(1,nvar+2),ilevel)
         call make_virtual_fine_dp(uold(1,nvar+3),ilevel)
!         call make_virtual_fine_dp(uold(1,nvar),ilevel)

         if  (nsubdiff == 1) exit 

      end do

#ifndef WITHOUTMPI
   call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif

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

end subroutine diffusion_sts
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine diffusion_fine_sts(ilevel,dtdiff,icycle)!,first_cmp_current_sts)
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
      call diffine1_sts(ind_grid,ngrid,dtdiff,ilevel,icycle)
   end do

   111 format('   Entering diffusion_fine_sts for level ',i2)

end subroutine diffusion_fine_sts
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine diffine1_sts(ind_grid,ncache,dtdiff,ilevel,icycle)
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
   real(dp)::emag,emagold,ekin,u,v,w,d
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
                     end do
                  end do
                  do i=1,nbuffer
                     ind1(i,j)=son(ibuffer_father(i,j))
                  end do
               end do
               call interpol_hydro(u1,ind1,u2,nbuffer)
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

                     if(icycle==1)then
                        do i=1,nexist
                           A=0.5d0*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
                           B=0.5d0*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
                           C=0.5d0*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))
                           enew(ind_cell(i))=A**2+B**2+C**2

                           d=max(uold(ind_cell(i),1),smallr)
                           u=uold(ind_cell(i),2)/d
                           v=uold(ind_cell(i),3)/d
                           w=uold(ind_cell(i),4)/d
                           ekin=0.5d0*d*(u**2+v**2+w**2)
#if NENER>0
                           do ivar=1,nener
                              ekin=ekin+uold(ind_cell(i),8+ivar)
                           end do
#endif
                           emag=0.5d0*(A**2+B**2+C**2)
                           divu(ind_cell(i))=uold(ind_cell(i),5)-ekin-emag
                        end do

                     end if

                     ! Gather hydro variables
                     do ivar=1,nvar+3
                        do i=1,nexist
                           uloc(ind_exist(i),i3,j3,k3,ivar)=uold(ind_cell(i),ivar)
                           if(ivar==2)uloc(ind_exist(i),i3,j3,k3,ivar)=enew(ind_cell(i))
                           if(ivar==3)uloc(ind_exist(i),i3,j3,k3,ivar)=divu(ind_cell(i))
                           end do
                        do i=1,nbuffer
                           uloc(ind_nexist(i),i3,j3,k3,ivar)=u2(i,ind_son,ivar)
                           if(ivar==2)then
                              uloc(ind_nexist(i),i3,j3,k3,ivar)=  &
                                 &  (0.5d0*(u2(i,ind_son,6)+u2(i,ind_son,nvar+1)))**2 &
                                 & +(0.5d0*(u2(i,ind_son,7)+u2(i,ind_son,nvar+2)))**2 &
                                 & +(0.5d0*(u2(i,ind_son,8)+u2(i,ind_son,nvar+3)))**2
                           end if
                           if(ivar==3)then
                              uloc(ind_nexist(i),i3,j3,k3,ivar)=u2(i,ind_son,nvar)
                           end if
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
   call cmp_current_sts(uloc,emfx,emfy,emfz,flux,ncache,dx_loc,dx_loc,dx_loc,dtdiff,jxbsquare,jcentersquare,ok,jcell)

   !------------------------------------------------
   ! Reset flux along direction at refined interface    
   !------------------------------------------------
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
                     flux(i,i3,j3,k3,idim)=0.0d0
                  end if
               end do
            end do
         end do
      end do
   end do

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
              if(son(ind_cell(i))==0)then

                 ! Update Bx using constraint transport
                 dflux_x=( emfy(i,i3,j3,k3)-emfy(i,i3,j3,k3+1) ) &
                    &    -( emfz(i,i3,j3,k3)-emfz(i,i3,j3+1,k3) )
                 unew(ind_cell(i),ivar1)=unew(ind_cell(i),ivar1)+dflux_x*dtdiff/dx_loc
                 dflux_x=( emfy(i,i3+1,j3,k3)-emfy(i,i3+1,j3,k3+1) ) &
                    &    -( emfz(i,i3+1,j3,k3)-emfz(i,i3+1,j3+1,k3) )  
                 unew(ind_cell(i),ivar4)=unew(ind_cell(i),ivar4)+dflux_x*dtdiff/dx_loc

                 ! Update By using constraint transport
                 dflux_y=( emfz(i,i3,j3,k3)-emfz(i,i3+1,j3,k3) ) &
                    &    -( emfx(i,i3,j3,k3)-emfx(i,i3,j3,k3+1) )
                 unew(ind_cell(i),ivar2)=unew(ind_cell(i),ivar2)+dflux_y*dtdiff/dx_loc
                 dflux_y=( emfz(i,i3,j3+1,k3)-emfz(i,i3+1,j3+1,k3) ) &
                    &    -( emfx(i,i3,j3+1,k3)-emfx(i,i3,j3+1,k3+1) )
                 unew(ind_cell(i),ivar5)=unew(ind_cell(i),ivar5)+dflux_y*dtdiff/dx_loc

                 ! Update Bz using constraint transport
                 dflux_z=( emfx(i,i3,j3,k3)-emfx(i,i3,j3+1,k3) ) &
                    &    -( emfy(i,i3,j3,k3)-emfy(i,i3+1,j3,k3) )
                 unew(ind_cell(i),ivar3)=unew(ind_cell(i),ivar3)+dflux_z*dtdiff/dx_loc
                 dflux_z=( emfx(i,i3,j3,k3+1)-emfx(i,i3,j3+1,k3+1) ) &
                    &    -( emfy(i,i3,j3,k3+1)-emfy(i,i3+1,j3,k3+1) )
                 unew(ind_cell(i),ivar6)=unew(ind_cell(i),ivar6)+dflux_z*dtdiff/dx_loc

!!$                 if ((.not. barotrop).and. (.not.radiative_nimhdheating)) then
!!$                    unew(ind_cell(i),nvar)=unew(ind_cell(i),nvar)+&
!!$                       &jxbsquare(i,i3,j3,k3) +jcentersquare(i,i3,j3,k3)
!!$                 endif
!!$
!!$                 if(radiative_nimhdheating)then
                    do idim=1,3
                       unew(ind_cell(i),nvar-3+idim)=jcell(i,i3   ,j3   ,k3   ,idim)
                    end do
!!$                 end if
              end if
           end do!
        end do
     end do
  end do
  
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
        weight=1.0
        if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
        if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
        dflux=(emfx(i,1,1,1)+emfx(i,2,1,1))*0.25*weight*dtdiff/dx_loc
        unew(ind_buffer1,ivar2)=unew(ind_buffer1,ivar2)+dflux
        unew(ind_buffer2,ivar5)=unew(ind_buffer2,ivar5)+dflux
        unew(ind_buffer2,ivar6)=unew(ind_buffer2,ivar6)-dflux
        unew(ind_buffer3,ivar3)=unew(ind_buffer3,ivar3)-dflux
        if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
           unew(ind_buffer1,3+nvar)=unew(ind_buffer1,3+nvar)+dflux*0.5
           unew(ind_buffer3,2+nvar)=unew(ind_buffer3,2+nvar)-dflux*0.5
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
        weight=1.0
        if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
        if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
        dflux=(emfx(i,1,1,3)+emfx(i,2,1,3))*0.25*weight*dtdiff/dx_loc
        unew(ind_buffer1,ivar6)=unew(ind_buffer1,ivar6)-dflux
        unew(ind_buffer2,ivar3)=unew(ind_buffer2,ivar3)-dflux
        unew(ind_buffer2,ivar5)=unew(ind_buffer2,ivar5)-dflux
        unew(ind_buffer3,ivar2)=unew(ind_buffer3,ivar2)-dflux
        if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
           unew(ind_buffer1,2+nvar)=unew(ind_buffer1,2+nvar)+dflux*0.5
           unew(ind_buffer3,3+neul)=unew(ind_buffer3,3+neul)+dflux*0.5
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
        weight=1.0
        if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
        if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
        dflux=(emfx(i,1,3,3)+emfx(i,2,3,3))*0.25*weight*dtdiff/dx_loc
        unew(ind_buffer1,ivar5)=unew(ind_buffer1,ivar5)-dflux
        unew(ind_buffer2,ivar2)=unew(ind_buffer2,ivar2)-dflux
        unew(ind_buffer2,ivar3)=unew(ind_buffer2,ivar3)+dflux
        unew(ind_buffer3,ivar6)=unew(ind_buffer3,ivar6)+dflux
        if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
           unew(ind_buffer3,2+neul)=unew(ind_buffer3,2+neul)+dflux*0.5
           unew(ind_buffer1,3+neul)=unew(ind_buffer1,3+neul)-dflux*0.5
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
        weight=1.0
        if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
        if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
        dflux=(emfx(i,1,3,1)+emfx(i,2,3,1))*0.25*weight*dtdiff/dx_loc
        unew(ind_buffer1,ivar3)=unew(ind_buffer1,ivar3)+dflux
        unew(ind_buffer2,ivar6)=unew(ind_buffer2,ivar6)+dflux
        unew(ind_buffer2,ivar2)=unew(ind_buffer2,ivar2)+dflux
        unew(ind_buffer3,ivar5)=unew(ind_buffer3,ivar5)+dflux
        if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
           unew(ind_buffer3,3+nvar)=unew(ind_buffer3,3+nvar)-dflux*0.5
           unew(ind_buffer1,2+neul)=unew(ind_buffer1,2+neul)-dflux*0.5
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
        weight=1.0
        if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
        if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
        dflux=(emfy(i,1,1,1)+emfy(i,1,2,1))*0.25*weight*dtdiff/dx_loc
        unew(ind_buffer1,ivar1)=unew(ind_buffer1,ivar1)-dflux
        unew(ind_buffer2,ivar4)=unew(ind_buffer2,ivar4)-dflux
        unew(ind_buffer2,ivar6)=unew(ind_buffer2,ivar6)+dflux
        unew(ind_buffer3,ivar3)=unew(ind_buffer3,ivar3)+dflux
        if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
           unew(ind_buffer3,1+nvar)=unew(ind_buffer3,1+nvar)+dflux*0.5
           unew(ind_buffer1,3+nvar)=unew(ind_buffer1,3+nvar)-dflux*0.5
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
        weight=1.0
        if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
        if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
        dflux=(emfy(i,1,1,3)+emfy(i,1,2,3))*0.25*weight*dtdiff/dx_loc
        unew(ind_buffer1,ivar6)=unew(ind_buffer1,ivar6)+dflux
        unew(ind_buffer2,ivar3)=unew(ind_buffer2,ivar3)+dflux
        unew(ind_buffer2,ivar4)=unew(ind_buffer2,ivar4)+dflux
        unew(ind_buffer3,ivar1)=unew(ind_buffer3,ivar1)+dflux
        if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
           unew(ind_buffer3,3+neul)=unew(ind_buffer3,3+neul)-dflux*0.5
           unew(ind_buffer1,1+nvar)=unew(ind_buffer1,1+nvar)-dflux*0.5
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
        weight=1.0
        if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
        if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
        dflux=(emfy(i,3,1,3)+emfy(i,3,2,3))*0.25*weight*dtdiff/dx_loc
        unew(ind_buffer1,ivar4)=unew(ind_buffer1,ivar4)+dflux
        unew(ind_buffer2,ivar1)=unew(ind_buffer2,ivar1)+dflux
        unew(ind_buffer2,ivar3)=unew(ind_buffer2,ivar3)-dflux
        unew(ind_buffer3,ivar6)=unew(ind_buffer3,ivar6)-dflux
        if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
           unew(ind_buffer3,1+neul)=unew(ind_buffer3,1+neul)-dflux*0.5
           unew(ind_buffer1,3+neul)=unew(ind_buffer1,3+neul)+dflux*0.5
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
        weight=1.0
        if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
        if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
        dflux=(emfy(i,3,1,1)+emfy(i,3,2,1))*0.25*weight*dtdiff/dx_loc
        unew(ind_buffer1,ivar3)=unew(ind_buffer1,ivar3)-dflux
        unew(ind_buffer2,ivar6)=unew(ind_buffer2,ivar6)-dflux
        unew(ind_buffer2,ivar1)=unew(ind_buffer2,ivar1)-dflux
        unew(ind_buffer3,ivar4)=unew(ind_buffer3,ivar4)-dflux
        if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
           unew(ind_buffer3,3+nvar)=unew(ind_buffer3,3+nvar)+dflux*0.5
           unew(ind_buffer1,1+neul)=unew(ind_buffer1,1+neul)+dflux*0.5  
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
        weight=1.0
        if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
        if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
        dflux=(emfz(i,1,1,1)+emfz(i,1,1,2))*0.25*weight*dtdiff/dx_loc
        unew(ind_buffer1,ivar1)=unew(ind_buffer1,ivar1)+dflux
        unew(ind_buffer2,ivar4)=unew(ind_buffer2,ivar4)+dflux
        unew(ind_buffer2,ivar5)=unew(ind_buffer2,ivar5)-dflux
        unew(ind_buffer3,ivar2)=unew(ind_buffer3,ivar2)-dflux
        if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
           unew(ind_buffer3,1+nvar)=unew(ind_buffer3,1+nvar)-dflux*0.5
           unew(ind_buffer1,2+nvar)=unew(ind_buffer1,2+nvar)+dflux*0.5
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
        weight=1.0
        if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
        if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
        dflux=(emfz(i,1,3,1)+emfz(i,1,3,2))*0.25*weight*dtdiff/dx_loc
        unew(ind_buffer1,ivar5)=unew(ind_buffer1,ivar5)-dflux
        unew(ind_buffer2,ivar2)=unew(ind_buffer2,ivar2)-dflux
        unew(ind_buffer2,ivar4)=unew(ind_buffer2,ivar4)-dflux
        unew(ind_buffer3,ivar1)=unew(ind_buffer3,ivar1)-dflux
        if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
           unew(ind_buffer3,2+neul)=unew(ind_buffer3,2+neul)+dflux*0.5
           unew(ind_buffer1,1+nvar)=unew(ind_buffer1,1+nvar)+dflux*0.5
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
        weight=1.0
        if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
        if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
        dflux=(emfz(i,3,3,1)+emfz(i,3,3,2))*0.25*weight*dtdiff/dx_loc
        unew(ind_buffer1,ivar4)=unew(ind_buffer1,ivar4)-dflux
        unew(ind_buffer2,ivar1)=unew(ind_buffer2,ivar1)-dflux
        unew(ind_buffer2,ivar2)=unew(ind_buffer2,ivar2)+dflux
        unew(ind_buffer3,ivar5)=unew(ind_buffer3,ivar5)+dflux
        if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
           unew(ind_buffer3,1+neul)=unew(ind_buffer3,1+neul)+dflux*0.5
           unew(ind_buffer1,2+neul)=unew(ind_buffer1,2+neul)-dflux*0.5
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
        weight=1.0
        if(son(ind_buffer1)>0.and.son(ind_buffer3)>0) cycle
        if(son(ind_buffer1)>0.or.son(ind_buffer2)>0.or.son(ind_buffer3)>0)weight=0.5
        dflux=(emfz(i,3,1,1)+emfz(i,3,1,2))*0.25*weight*dtdiff/dx_loc
        unew(ind_buffer1,ivar2)=unew(ind_buffer1,ivar2)+dflux
        unew(ind_buffer2,ivar5)=unew(ind_buffer2,ivar5)+dflux
        unew(ind_buffer2,ivar1)=unew(ind_buffer2,ivar1)+dflux
        unew(ind_buffer3,ivar4)=unew(ind_buffer3,ivar4)+dflux
        if(son(ind_buffer1)==0.and.son(ind_buffer2)==0.and.son(ind_buffer3)==0) then
           unew(ind_buffer3,2+nvar)=unew(ind_buffer3,2+nvar)-dflux*0.5
           unew(ind_buffer1,1+neul)=unew(ind_buffer1,1+neul)-dflux*0.5
        endif
     end do

  endif

end subroutine diffine1_sts
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmp_current_sts(u,Ex_arete,Ey_arete,Ez_arete,fluxni, &
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

end subroutine cmp_current_sts
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine set_unew_sts(ilevel,iupdate,dtohm,dtad)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine sets array unew to its initial value uold before calling
  ! the hydro scheme. unew is set to zero in virtual boundaries.
  !--------------------------------------------------------------------------
  integer::i,j,ivar,ind,icpu,iskip,iupdate
  real(dp)::d,u,v,w,e,A,B,C,e_r,Cv,dtohm,dtad,dx,dtohmb,dtadb
  real(dp)::xx,tcell,B2,betaad,eps,scale,etaohmdiss
  real(dp),dimension(1:3)::skip_loc

   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

   ! Rescaling factors
   skip_loc=(/0.0d0,0.0d0,0.0d0/)
   if(ndim>0)skip_loc(1)=dble(icoarse_min)
   if(ndim>1)skip_loc(2)=dble(jcoarse_min)
   if(ndim>2)skip_loc(3)=dble(kcoarse_min)
   scale=dble(icoarse_max-icoarse_min+1)/boxlen

  dtohm=1d36
  dtad=1d36
  dx=0.5D0**ilevel/scale

  ! Set unew to uold for myid cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do ivar=1,nvar+3
        do i=1,active(ilevel)%ngrid
!           if(iupdate==0)enew(active(ilevel)%igrid(i)+iskip)=0.0d0
!           if(son(active(ilevel)%igrid(i)+iskip)==0)then
              unew(active(ilevel)%igrid(i)+iskip,ivar) = uold(active(ilevel)%igrid(i)+iskip,ivar)
!           endif
        end do
     end do
     do i=1,active(ilevel)%ngrid
!        if(iupdate==0)enew(active(ilevel)%igrid(i)+iskip)=0.0d0
!        if(son(active(ilevel)%igrid(i)+iskip)==0)then
           d=uold(active(ilevel)%igrid(i)+iskip,1)
           u=uold(active(ilevel)%igrid(i)+iskip,2)/d
           v=uold(active(ilevel)%igrid(i)+iskip,3)/d
           w=uold(active(ilevel)%igrid(i)+iskip,4)/d
           A=0.5d0*(uold(active(ilevel)%igrid(i)+iskip,6)+uold(active(ilevel)%igrid(i)+iskip,nvar+1))
           B=0.5d0*(uold(active(ilevel)%igrid(i)+iskip,7)+uold(active(ilevel)%igrid(i)+iskip,nvar+2))
           C=0.5d0*(uold(active(ilevel)%igrid(i)+iskip,8)+uold(active(ilevel)%igrid(i)+iskip,nvar+3))
           e_r=0.0d0
#if NENER>0
           do j=1,nener
              e_r=e_r+uold(active(ilevel)%igrid(i)+iskip,8+j)
           end do
#endif
           e=uold(active(ilevel)%igrid(i)+iskip,5)-0.5d0*d*(u**2+v**2+w**2)-0.5d0*(A**2+B**2+C**2)-e_r

           if(iupdate==0)then
              enew(active(ilevel)%igrid(i)+iskip) = (A**2+B**2+C**2)

              eps   = e 
              divu(active(ilevel)%igrid(i)+iskip) = eps !uold(active(ilevel)%igrid(i)+iskip,nvar)

              ! Compute gas temperature in cgs
              call temperature_eos(d,eps,Tcell)

              B2 = A**2+B**2+C**2
              !compute nimhd timesteps
              ! Ohmic dissipation
              if (nmagdiffu2.eq.1) then
                 xx=etaohmdiss(d,B2,tcell)
                 if(xx.gt.0d0) then
                    dtohmb=coefohm*dx*dx/xx
                 else
                    dtohmb=1d35
                 endif
                 dtohm=min(dtohmb,dtohm)
              end if
                            
              ! ambipolar diffusion
              if (nambipolar2.eq.1)then
                 dtad=1d36
                 xx=B2*betaad(d,B2,tcell) 
                 if (xx.gt.0d0) then
                    !! WARNING RHOAD mandatory because rho(k) is not density cf lines above
                    dtadb=coefad*dx*dx/xx
                 else
                    dtadb=1d36
                 endif
                 dtad=min(dtadb,dtad)
              endif  ! diffusion ambipolaire
              
           end if
           
        end do
     end do
   
  ! Set unew to 0 for virtual boundary cells
  do icpu=1,ncpu
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do ivar=1,nvar+3
        do i=1,reception(icpu,ilevel)%ngrid
              unew(reception(icpu,ilevel)%igrid(i)+iskip,ivar)=0.0d0
        end do
     end do
!!$        do i=1,reception(icpu,ilevel)%ngrid
!!$              enew(reception(icpu,ilevel)%igrid(i)+iskip)=0.0d0
!!$              divu(reception(icpu,ilevel)%igrid(i)+iskip)=0.0d0
!!$        end do
  end do
  end do

111 format('   Entering set_unew for level ',i2)

end subroutine set_unew_sts
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine set_uold_sts(ilevel,iend,dtloc)
  use amr_commons
  use hydro_commons
  use poisson_commons
  use pm_commons
  use mpi_mod
  use constants
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine sets array uold to its new value unew after the
  ! hydro step.
  !--------------------------------------------------------------------------
  integer::i,j,k,ivar,ind,iskip,nx_loc,info,iend,jdim,icpu
  real(dp)::scale,d,u,v,w,A,B,C,d_old
  real(dp)::e_mag,e_kin,e_cons,e_prim,e_trunc,div,dx,fact,e_r

  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  integer ,dimension(1:nvector,0:twondim),save::igridn
  integer ,dimension(1:nvector,1:ndim),save::ind_left,ind_right
  real(dp),dimension(1:nvector,1:ndim,1:ndim),save::velg,veld
  real(dp),dimension(1:nvector,1:ndim,1:ndim),save::Bg,Bd
  real(dp),dimension(1:nvector,1:ndim)::dx_g,dx_d
  real(dp)::Pgdivu,u_square
  real(dp)::d_loc,Tp_loc

  integer::ncache,igrid,ngrid,idim,id1,ig1,ih1,id2,ig2,ih2,igroup
  integer  ,dimension(1:3,1:2,1:8)::iii,jjj
  real(dp)::dx_loc,surf_loc,vol_loc,usquare,emag,erad_loc,ekin,eps,cv,pp_eos
  real(dp)::kappa_R,gradEr_norm,gradEr_norm2,R,lambda,lambda_fld,chi,PgmErdivu,gradEru
  real(dp) ,dimension(1:3)::skip_loc
  real(dp) ,dimension(1:ndim,1:ndim)::divu_loc
  real(dp) ,dimension(1:ndim       )::u_loc
  real(dp) :: nuPrDivu,nuPr,nuPl,Pr_nu
  real(dp), dimension(1:5) :: Pr_temp

  !  EOS
  real(dp) :: dd,ee
  
  !SINK
  real(dp),dimension(1:twotondim,1:3)::xc  
 
  real(dp)::ekinold,emagold,erold,etotold,tcell,xy,b2,dtloc,dtlocb,betaad
  real(dp)::norm_jcenter,norm_jxb,eta_ad,eta_ohm,Bx,By,Bz
  real(dp)::Ohm_heating,etaohmdiss,jsquare
  real(dp)::ambi_heating,nimhd_heating,jx,jy,jz,bcell2

   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5d0**ilevel

  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale ! Warning: scale factor already done in dx
  vol_loc=dx_loc**ndim
  surf_loc=dx_loc**(ndim-1)

  iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

  ! Set uold to unew for myid cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,active(ilevel)%ngrid
        do ivar=1,nvar+3
              uold(active(ilevel)%igrid(i)+iskip,ivar) = unew(active(ilevel)%igrid(i)+iskip,ivar)
           end do
           d=uold(active(ilevel)%igrid(i)+iskip,1)
           u=uold(active(ilevel)%igrid(i)+iskip,2)/d
           v=uold(active(ilevel)%igrid(i)+iskip,3)/d
           w=uold(active(ilevel)%igrid(i)+iskip,4)/d
           A=0.5d0*(uold(active(ilevel)%igrid(i)+iskip,6)+uold(active(ilevel)%igrid(i)+iskip,nvar+1))
           B=0.5d0*(uold(active(ilevel)%igrid(i)+iskip,7)+uold(active(ilevel)%igrid(i)+iskip,nvar+2))
           C=0.5d0*(uold(active(ilevel)%igrid(i)+iskip,8)+uold(active(ilevel)%igrid(i)+iskip,nvar+3))
           e_kin=0.5d0*d*(u**2+v**2+w**2)
           e_mag=0.5d0*(A**2+B**2+C**2)
           e_r=0.0D0
#if NENER>0
           do j=1,nener
              e_r=e_r+uold(active(ilevel)%igrid(i)+iskip,8+j)
           end do      
#endif     
           e_cons=divu(active(ilevel)%igrid(i)+iskip)!uold(active(ilevel)%igrid(i)+iskip,5)-e_kin-e_mag-e_r

           if(iend==1)then
              nimhd_heating=0.0d0
              ambi_heating=0.0d0
              ohm_heating=0.0d0
               
              bcell2=(A**2+B**2+C**2)
              jx=uold(active(ilevel)%igrid(i)+iskip,nvar-3)
              jy=uold(active(ilevel)%igrid(i)+iskip,nvar-2)
              jz=uold(active(ilevel)%igrid(i)+iskip,nvar-1)
              jsquare=(jx**2+jy**2+jz**2)
              if(ntestDADM.eq.1)then
                 tcell=1.0d0
              else 
                 call temperature_eos(d,uold(active(ilevel)%igrid(i)+iskip,nvar),tcell)
              end if

              if(nmagdiffu2 .eq. 1 )ohm_heating=jsquare*etaohmdiss(d,bcell2,tcell)*dtnew(ilevel)!*vol_loc

              if(nambipolar2 .eq. 1 )then
                 ambi_heating = (jy*C-jz*B)**2+(jz*A-jx*C)**2+(jx*B-jy*A)**2
                 ambi_heating = ambi_heating * betaad(d,bcell2,tcell)*dtnew(ilevel)
              endif
              nimhd_heating=ambi_heating+ohm_heating
              uold(active(ilevel)%igrid(i)+iskip,5)=e_cons+e_kin+e_mag+e_r+nimhd_heating
              uold(active(ilevel)%igrid(i)+iskip,nvar)= e_cons+nimhd_heating
           end if

           ! Compute gas temperature in cgs
           dd=d *scale_d
           ee=e_cons *scale_d*scale_v**2
           e_prim=e_prim *scale_d*scale_v**2
           
           ! Compute temperature for perfect gas to prevent crash in the interpolation routine of the EOS
           Cv= dd*kB/(mu_gas*mH*(gamma-1.0d0))!(cgs)
           Tp_loc =  ee/(Cv)
           e_prim = e_prim/CV
           
           if(e_cons < 0.0)then
              write(*,*) 'uold(nvar) < 0',uold(active(ilevel)%igrid(i)+iskip,nvar)
              write(*,*)'new', uold(active(ilevel)%igrid(i)+iskip,5)-e_kin-e_mag-e_r
              !            write(*,*)'old', etotold-ekinold-emagold-erold,iend
           endif
     end do
  end do
  
111 format('   Entering set_uold_sts for level ',i2)

end subroutine set_uold_sts
!###########################################################
!###########################################################
!###########################################################
!###########################################################