#if NDIM==3
subroutine compute_clump_properties(xx)
  use amr_commons
  use hydro_commons, ONLY:uold
  use clfind_commons
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  real(dp),dimension(1:ncoarse+ngridmax*twotondim)::xx
  !----------------------------------------------------------------------------
  ! this subroutine performs a loop over all cells above the threshold and
  ! collects the  relevant information. After some MPI communications,
  ! all necessary peak-patch properties are computed
  !----------------------------------------------------------------------------
  integer::ipart,grid,peak_nr,ilevel,global_peak_id,ipeak,plevel
  real(dp)::zero=0
  !variables needed temporarily store cell properties
  real(dp)::d=0, vol=0
  ! variables related to the size of a cell on a given level
  real(dp)::dx,dx_loc,scale,vol_loc
  real(dp),dimension(1:nlevelmax)::volume
  real(dp),dimension(1:3)::skip_loc,xcell
  real(dp),dimension(1:twotondim,1:3)::xc
  integer::nx_loc,ind,ix,iy,iz,idim
  logical,dimension(1:ndim)::period
  logical::periodic

#ifndef WITHOUTMPI
  integer::i
  real(dp)::tot_mass_tot
#endif

  period(1)=(nx==1)
  period(2)=(ny==1)
  period(3)=(nz==1)

  periodic=period(1)
  periodic=periodic.or.period(2)
  periodic=periodic.or.period(3)

  !peak-patch related arrays before sharing information with other cpus

  min_dens=huge(zero)
  n_cells=0; n_cells_halo=0
  halo_mass=0d0; clump_mass=0d0; clump_vol=0d0
  center_of_mass=0d0; clump_velocity=0d0
  peak_pos=0d0

  if(verbose)write(*,*)'Entering compute clump properties'
  !------------------------------------------
  ! compute volume of a cell in a given level
  !------------------------------------------
  do ilevel=1,nlevelmax
     ! Mesh spacing in that level
     dx=0.5D0**ilevel
     nx_loc=(icoarse_max-icoarse_min+1)
     scale=boxlen/dble(nx_loc)
     dx_loc=dx*scale
     vol_loc=dx_loc**ndim
     volume(ilevel)=vol_loc
  end do

  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)

  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     xc(ind,1)=(dble(ix)-0.5D0)
     xc(ind,2)=(dble(iy)-0.5D0)
     xc(ind,3)=(dble(iz)-0.5D0)
  end do

  !--------------------------------------------------------------------------
  ! loop over all cells above the threshold
  !--------------------------------------------------------------------------
  do ipart=1,ntest
     global_peak_id=flag2(icellp(ipart))
     if (global_peak_id /=0 ) then
        call get_local_peak_id(global_peak_id,peak_nr)

        ! Cell coordinates
        ind=(icellp(ipart)-ncoarse-1)/ngridmax+1 ! cell position
        grid=icellp(ipart)-ncoarse-(ind-1)*ngridmax ! grid index
        dx=0.5D0**levp(ipart)
        xcell(1:ndim)=(xg(grid,1:ndim)+xc(ind,1:ndim)*dx-skip_loc(1:ndim))*scale

        ! gas density
        if(ivar_clump==0 .or. ivar_clump==-1)then
           d=xx(icellp(ipart))
        else
           if(hydro)then
              d=uold(icellp(ipart),1)
           endif
        endif

        ! Cell volume
        vol=volume(levp(ipart))

        ! Number of leaf cells per clump
        n_cells(peak_nr)=n_cells(peak_nr)+1

        ! Min density
        min_dens(peak_nr)=min(d,min_dens(peak_nr))

        ! Clump mass
        clump_mass(peak_nr)=clump_mass(peak_nr)+vol*d

        ! Clump volume
        clump_vol(peak_nr)=clump_vol(peak_nr)+vol

        ! Clump center of mass location
        center_of_mass(peak_nr,1:3)=center_of_mass(peak_nr,1:3)+vol*d*xcell(1:3)

        ! Clump center of mass velocity
        if (hydro)clump_velocity(peak_nr,1:3)=clump_velocity(peak_nr,1:3)+vol*uold(icellp(ipart),2:4)

     end if
  end do

  !--------------------------------------------------------------------------
  ! Loop over local peaks and identify true peak positions
  !--------------------------------------------------------------------------
  do ipeak=1,npeaks
     ! Peak cell coordinates
     ind=(peak_cell(ipeak)-ncoarse-1)/ngridmax+1    ! cell position
     grid=peak_cell(ipeak)-ncoarse-(ind-1)*ngridmax ! grid index
     plevel=peak_cell_level(ipeak)
     dx=0.5D0**plevel
     xcell(1:ndim)=(xg(grid,1:ndim)+xc(ind,1:ndim)*dx-skip_loc(1:ndim))*scale
     call true_max(xcell(1),xcell(2),xcell(3),plevel)
     peak_pos(ipeak,1:3)=xcell(1:3)
  end do
  call build_peak_communicator

#ifndef WITHOUTMPI
  ! Collect results from all MPI domains
  call virtual_peak_int(n_cells,'sum')
  call virtual_peak_dp(min_dens,'min')
  call virtual_peak_dp(clump_mass,'sum')
  call virtual_peak_dp(clump_vol,'sum')
  do i=1,ndim
     call virtual_peak_dp(center_of_mass(1,i),'sum')
     call virtual_peak_dp(clump_velocity(1,i),'sum')
  end do
#endif

  ! Compute specific quantities
  do ipeak=1,npeaks
     if (relevance(ipeak)>0..and.n_cells(ipeak)>0)then
        center_of_mass(ipeak,1:3)=center_of_mass(ipeak,1:3)/clump_mass(ipeak)
        clump_velocity(ipeak,1:3)=clump_velocity(ipeak,1:3)/clump_mass(ipeak)
     end if
  end do

#ifndef WITHOUTMPI
  ! Scatter results to all MPI domains
  do i=1,ndim
     call boundary_peak_dp(peak_pos(1,i))
     call boundary_peak_dp(center_of_mass(1,i))
     call boundary_peak_dp(clump_velocity(1,i))
  end do
#endif

  ! Initialize halo mass to clump mass
  halo_mass=clump_mass

  ! Calculate total mass above threshold
  tot_mass=sum(clump_mass(1:npeaks))

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(tot_mass,tot_mass_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  tot_mass=tot_mass_tot
#endif

  ! Compute further properties of the clumps
  do ipeak=1,npeaks
     if (relevance(ipeak)>0..and.n_cells(ipeak)>0)then
        av_dens(ipeak)=clump_mass(ipeak)/clump_vol(ipeak)
     end if
  end do

#ifndef WITHOUTMPI
  ! Scatter results to all MPI domains
  call boundary_peak_dp(av_dens)
#endif

  ! For periodic boxes, recompute center of mass relative to peak position
  if(periodic)then
     center_of_mass=0d0;
     do ipart=1,ntest
        global_peak_id=flag2(icellp(ipart))
        if (global_peak_id /=0 ) then
           call get_local_peak_id(global_peak_id,peak_nr)

           ! Cell coordinates
           ind=(icellp(ipart)-ncoarse-1)/ngridmax+1 ! cell position
           grid=icellp(ipart)-ncoarse-(ind-1)*ngridmax ! grid index
           dx=0.5D0**levp(ipart)
           xcell(1:ndim)=(xg(grid,1:ndim)+xc(ind,1:ndim)*dx-skip_loc(1:ndim))*scale

           do idim=1,ndim
              if (period(idim) .and. (xcell(idim)-peak_pos(peak_nr,idim))>boxlen*0.5)xcell(idim)=xcell(idim)-boxlen
              if (period(idim) .and. (xcell(idim)-peak_pos(peak_nr,idim))<boxlen*(-0.5))xcell(idim)=xcell(idim)+boxlen
           end do

           ! gas density
           if(ivar_clump==0 .or. ivar_clump==-1)then
              d=xx(icellp(ipart))
           else
              if(hydro)then
                 d=uold(icellp(ipart),1)
              endif
           endif

           ! Cell volume
           vol=volume(levp(ipart))

           ! Clump center of mass location
           center_of_mass(peak_nr,1:3)=center_of_mass(peak_nr,1:3)+vol*d*xcell(1:3)

        end if
     end do
     call build_peak_communicator

#ifndef WITHOUTMPI
     ! Collect results from all MPI domains
     do i=1,ndim
        call virtual_peak_dp(center_of_mass(1,i),'sum')
     end do
#endif

     ! Compute specific quantity
     do ipeak=1,npeaks
        if (relevance(ipeak)>0..and.n_cells(ipeak)>0)then
           center_of_mass(ipeak,1:3)=center_of_mass(ipeak,1:3)/clump_mass(ipeak)
        end if
     end do

#ifndef WITHOUTMPI
     ! Scatter results to all MPI domains
     do i=1,ndim
        call boundary_peak_dp(center_of_mass(1,i))
     end do
#endif

  end if
end subroutine compute_clump_properties
!################################################################
!################################################################
!################################################################
!################################################################
subroutine write_clump_properties(to_file)
  use amr_commons
  use pm_commons,ONLY:mp
  use hydro_commons,ONLY:mass_sph
  use clfind_commons
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer,parameter::tag=1101
  integer::dummy_io,info,info2
#endif
  logical::to_file
  !---------------------------------------------------------------------------
  ! this routine writes the clump properties to screen and to file
  !---------------------------------------------------------------------------

  integer::i,j,jj,ilun,ilun2,n_rel,n_rel_tot,nx_loc
  real(dp)::rel_mass,rel_mass_tot,scale,particle_mass=0
  character(LEN=80)::fileloc,filedir
  character(LEN=5)::nchar,ncharcpu
  real(dp),dimension(1:npeaks)::peakd
  integer,dimension(1:npeaks)::ind_sort

#ifndef WITHOUTMPI
  real(dp)::particle_mass_tot
#endif

  if (.not. to_file)return

  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  if(ivar_clump==0 .or. ivar_clump==-1)then
     particle_mass=MINVAL(mp, MASK=(mp > 0))
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(particle_mass,particle_mass_tot,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,info)
     particle_mass=particle_mass_tot
#endif
  else
     if(hydro)then
        particle_mass=mass_sph
     endif
  endif

  ! sort clumps by peak density in ascending order
  do i=1,npeaks
     peakd(i)=max_dens(i)
     ind_sort(i)=i
  end do
  call quick_sort_dp(peakd,ind_sort,npeaks)

  if(to_file)then
     ilun=20
     ilun2=22
  else
     ilun=6
     ilun2=6
  end if

  ! print results in descending order to screen/file
  rel_mass=0
  n_rel=0

  if (to_file .eqv. .true.) then
     ! first create directories
     call title(ifout,nchar)
     filedir='output_'//TRIM(nchar)
     call create_output_dirs(filedir)
     ! Wait for the token
#ifndef WITHOUTMPI
     if(IOGROUPSIZE>0) then
        if (mod(myid-1,IOGROUPSIZE)/=0) then
           call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
                & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
        end if
     endif
#endif

     if(IOGROUPSIZEREP>0)then
        call title(((myid-1)/IOGROUPSIZEREP)+1,ncharcpu)
        fileloc=TRIM(filedir)//'/group_'//TRIM(ncharcpu)//'/clump_'//TRIM(nchar)//'.txt'
     else
        fileloc=TRIM(filedir)//'/clump_'//TRIM(nchar)//'.txt'
     endif
     call title(myid,nchar)
     fileloc=TRIM(fileloc)//TRIM(nchar)
     open(unit=ilun,file=fileloc,form='formatted')

     if(saddle_threshold>0)then
        call title(ifout,nchar)
        if(IOGROUPSIZEREP>0)then
           call title(((myid-1)/IOGROUPSIZEREP)+1,ncharcpu)
           fileloc=TRIM(filedir)//'/group_'//TRIM(ncharcpu)//'/halo_'//TRIM(nchar)//'.txt'
        else
           fileloc=TRIM(filedir)//'/halo_'//TRIM(nchar)//'.txt'
        endif
        call title(myid,nchar)
        fileloc=TRIM(fileloc)//TRIM(nchar)
        open(unit=ilun2,file=fileloc,form='formatted')
     endif
  end if

  if (to_file .or. myid==1)then
     write(ilun,'(144A)')'   index  halo   lev   parent      ncell    peak_x             peak_y             peak_z     '//&
          '        rho-               rho+               rho_av             mass_cl            relevance   '
     if(saddle_threshold>0)then
        write(ilun2,'(135A)')'     index      ncell    peak_x             peak_y             peak_z     '//&
             '        rho+               mass      '
     endif
  end if

  if (particlebased_clump_output) then ! write particle based data
    do j=npeaks,1,-1
      jj=ind_sort(j)
      if (relevance(jj) > relevance_threshold .and. clmp_mass_pb(jj) > mass_threshold*particle_mass)then
        write(ilun,'(I8,X,I2,X,I10,X,I10,8(X,1PE18.9E2))')&
               jj+ipeak_start(myid)&
               ,lev_peak(jj)&
               ,new_peak(jj)&
               ,n_cells(jj)&
#ifdef UNBINDINGCOM
               ,clmp_com_pb(jj,1)&
               ,clmp_com_pb(jj,2)&
               ,clmp_com_pb(jj,3)&
#else
               ,peak_pos(jj,1)&
               ,peak_pos(jj,2)&
               ,peak_pos(jj,3)&
#endif
               ,min_dens(jj)&
               ,max_dens(jj)&
               ,clmp_mass_pb(jj)/clump_vol(jj)&
               ,clmp_mass_pb(jj)&
               ,relevance(jj)
         rel_mass=rel_mass+clmp_mass_exclusive(jj)
         n_rel=n_rel+1
      end if

      if(saddle_threshold>0)then
        if(ind_halo(jj).EQ.jj+ipeak_start(myid).AND.clmp_mass_pb(jj) > mass_threshold*particle_mass)then
           write(ilun2,'(I10,X,I10,5(X,1PE18.9E2))')&
                  jj+ipeak_start(myid)&
                  ,n_cells_halo(jj)&
#ifdef UNBINDINGCOM
                  ,clmp_com_pb(jj,1)&
                  ,clmp_com_pb(jj,2)&
                  ,clmp_com_pb(jj,3)&
#else
                  ,peak_pos(jj,1)&
                  ,peak_pos(jj,2)&
                  ,peak_pos(jj,3)&
#endif
                  ,max_dens(jj)&
                  ,clmp_mass_pb(jj)
        endif
      endif
    end do

  else ! write cell based data

    do j=npeaks,1,-1
       jj=ind_sort(j)
       if (relevance(jj) > relevance_threshold .and. halo_mass(jj) > mass_threshold*particle_mass)then
          write(ilun,'(I8,X,I8,1X,I2,X,I10,X,I10,8(X,1PE18.9E2))')&
               jj+ipeak_start(myid)&
               ,ind_halo(jj)&
               ,lev_peak(jj)&
               ,new_peak(jj)&
               ,n_cells(jj)&
               ,peak_pos(jj,1)&
               ,peak_pos(jj,2)&
               ,peak_pos(jj,3)&
               ,min_dens(jj)&
               ,max_dens(jj)&
               ,clump_mass(jj)/clump_vol(jj)&
               ,clump_mass(jj)&
               ,relevance(jj)
          rel_mass=rel_mass+clump_mass(jj)
          n_rel=n_rel+1
       end if
       if(saddle_threshold>0)then
          if(ind_halo(jj).EQ.jj+ipeak_start(myid).AND.halo_mass(jj) > mass_threshold*particle_mass)then
             write(ilun2,'(I10,X,I10,5(X,1PE18.9E2))')&
                  jj+ipeak_start(myid)&
                  ,n_cells_halo(jj)&
                  ,peak_pos(jj,1)&
                  ,peak_pos(jj,2)&
                  ,peak_pos(jj,3)&
                  ,max_dens(jj)&
                  ,halo_mass(jj)
          endif
       endif
    end do
  end if

  if (to_file)then
     close(ilun)
     if(saddle_threshold>0)then
        close(ilun2)
     endif
  end if

     ! Send the token
#ifndef WITHOUTMPI
     if(IOGROUPSIZE>0) then
        if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
           dummy_io=1
           call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
                & MPI_COMM_WORLD,info2)
        end if
     endif
#endif

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(n_rel,n_rel_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  n_rel=n_rel_tot
  call MPI_ALLREDUCE(rel_mass,rel_mass_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  rel_mass=rel_mass_tot
#else
  n_rel_tot = n_rel
  rel_mass_tot = rel_mass
#endif
  if(myid==1)then
     if(clinfo)write(*,'(A,1PE12.5)')' Total mass [code units] above threshold =',tot_mass
     if(clinfo)write(*,'(A,I10,A,1PE12.5)')' Total mass [code units] in',n_rel_tot,' listed clumps =',rel_mass_tot
  endif

end subroutine write_clump_properties
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine merge_clumps(action)
  use amr_commons
  use clfind_commons
  use mpi_mod
  implicit none
  character(len=9)::action
#ifndef WITHOUTMPI
  integer::info
#endif
  !---------------------------------------------------------------------------
  ! This routine merges the irrelevant clumps
  ! -clumps are sorted by ascending max density
  ! -irrelevent clumps are merged to most relevant neighbor
  !---------------------------------------------------------------------------

  integer::j,i,merge_to,ipart
  integer::current,nmove,ipeak,jpeak,iter
  integer::nsurvive,nzero,idepth
  integer::ilev,global_peak_id
  real(dp)::value_iij,zero=0,relevance_peak
  integer,dimension(1:npeaks_max)::alive,ind_sort
  real(dp),dimension(1:npeaks_max)::peakd
  logical::do_merge=.false.

#ifndef WITHOUTMPI
  integer::mergelevel_max_global
  integer::nmove_all,nsurvive_all,nzero_all
#endif

  if (verbose)then
     if(action.EQ.'relevance')then
        write(*,*)'Now merging irrelevant clumps'
     endif
     if(action.EQ.'saddleden')then
        write(*,*)'Now merging clumps into halos'
     endif
  endif

  ! Initialize new_peak array to global peak id
  ! All peaks are alive at the start
  do i=1,npeaks
     new_peak(i)=ipeak_start(myid)+i
     if(action.EQ.'relevance')then
        alive(i)=1
     endif
     if(action.EQ.'saddleden')then
        if(relevance(i)>relevance_threshold)then
           alive(i)=1
        else
           alive(i)=0
        endif
     endif
  end do

  ! Sort peaks by maximum peak density in ascending order
  do i=1,npeaks
     peakd(i)=max_dens(i)
     ind_sort(i)=i
  end do
  call quick_sort_dp(peakd,ind_sort,npeaks)

  ! Loop over peak levels
  nzero=npeaks_tot
  idepth=0
  do while(nzero>0)

     ! Compute maximum saddle density for each clump
     do i=1,hfree-1
        call get_max(i,sparse_saddle_dens)
     end do

#ifndef WITHOUTMPI
     ! Create new local duplicated peaks and update communicator
     call virtual_saddle_max
     call build_peak_communicator

     ! Set up bounday values
     call boundary_peak_dp(sparse_saddle_dens%maxval)
     call boundary_peak_int(sparse_saddle_dens%maxloc)
     call boundary_peak_dp(max_dens)
     call boundary_peak_int(new_peak)
     call boundary_peak_int(alive)
#endif

     ! Merge peaks
     nmove=npeaks_tot
     iter=0
     do while(nmove>0)
        nmove=0
        do i=npeaks,1,-1
           ipeak=ind_sort(i)
           merge_to=new_peak(ipeak)
           if(alive(ipeak)>0)then
              if(action.EQ.'relevance')then
                 if(sparse_saddle_dens%maxval(ipeak)>0)then
                    relevance_peak=max_dens(ipeak)/sparse_saddle_dens%maxval(ipeak)
                 else
                    relevance_peak=max_dens(ipeak)/density_threshold
                 end if
                 do_merge=(relevance_peak<relevance_threshold)
              endif
              if(action.EQ.'saddleden')then
                 do_merge=(sparse_saddle_dens%maxval(ipeak)>saddle_threshold)
              endif
              if(do_merge)then
                 if(sparse_saddle_dens%maxloc(ipeak)>0)then
                    call get_local_peak_id(sparse_saddle_dens%maxloc(ipeak),jpeak)
                    if(max_dens(jpeak)>max_dens(ipeak))then
                       merge_to=new_peak(jpeak)
                    else if(max_dens(jpeak)==max_dens(ipeak))then
                       merge_to=MIN(new_peak(ipeak),new_peak(jpeak))
                    endif
                 endif
              endif
           endif
           if(new_peak(ipeak).NE.merge_to)then
              nmove=nmove+1
              new_peak(ipeak)=merge_to
           endif
        end do
        ! Update boundary conditions for new_peak array
        call boundary_peak_int(new_peak)
        iter=iter+1
#ifndef WITHOUTMPI
        call MPI_ALLREDUCE(nmove,nmove_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
        nmove=nmove_all
#endif
        if(verbose)write(*,*)'niter=',iter,'nmove=',nmove
     end do

     ! Transfer matrix elements of merged peaks to surviving peaks
     ! Create new local duplicated peaks and update communicator
     do ipeak=1,hfree-1
        if(alive(ipeak)>0)then
           merge_to=new_peak(ipeak)
           if(ipeak.LE.npeaks)then
              global_peak_id=ipeak_start(myid)+ipeak
           else
              global_peak_id=gkey(ipeak)
           endif
           if(merge_to.NE.global_peak_id)then
              call get_local_peak_id(merge_to,jpeak)
              current=sparse_saddle_dens%first(ipeak) ! first element of line ipeak
              do while(current>0) ! walk the line
                 j=sparse_saddle_dens%col(current)
                 value_iij=sparse_saddle_dens%val(current) ! value of the matrix
                 ! Copy the value of density only if larger
                 if(value_iij>get_value(jpeak,j,sparse_saddle_dens))then
                    call set_value(jpeak,j,value_iij,sparse_saddle_dens)
                    call set_value(j,jpeak,value_iij,sparse_saddle_dens)
                 end if
                 current=sparse_saddle_dens%next(current)
              end do
              call set_value(jpeak,jpeak,zero,sparse_saddle_dens)
           end if
        endif
     end do
     call build_peak_communicator

     ! Set alive to zero for newly merged peaks
     nzero=0
     nsurvive=0
     do ipeak=1,npeaks
        if(alive(ipeak)>0)then
           merge_to=new_peak(ipeak)
           if(merge_to.NE.(ipeak_start(myid)+ipeak))then
              alive(ipeak)=0
              lev_peak(ipeak)=idepth
              nzero=nzero+1
           else
              nsurvive=nsurvive+1
           end if
        endif
     end do
     call boundary_peak_int(alive)

     ! Remove all matrix elements corresponding to merged peaks
     do ipeak=1,hfree-1
        current=sparse_saddle_dens%first(ipeak) ! first element of line ipeak
        do while(current>0) ! walk the line
           j=sparse_saddle_dens%col(current)
           if(alive(ipeak)==0 .or. alive(j)==0)then
              call set_value(ipeak,j,zero,sparse_saddle_dens)
           endif
           current=sparse_saddle_dens%next(current)
        end do
     end do

#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(nzero,nzero_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
     nzero=nzero_all
     call MPI_ALLREDUCE(nsurvive,nsurvive_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
     nsurvive=nsurvive_all
#endif
     if(verbose)write(*,*)'level=',idepth,'nmove=',nzero,'survived=',nsurvive
     idepth=idepth+1

  end do
  ! End loop over peak levels

  mergelevel_max=idepth-2 ! last level has no more clumps, also idepth=idepth+1 still happens on last level.
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(mergelevel_max,mergelevel_max_global,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
  mergelevel_max=mergelevel_max_global
#endif

  ! Compute maximum saddle density for each surviving clump
  ! Create new local duplicated peaks and update communicator
  do i=1,hfree-1
     call get_max(i,sparse_saddle_dens)
  end do
#ifndef WITHOUTMPI
  call virtual_saddle_max
  call build_peak_communicator
#endif

  ! Set up bounday values
  call boundary_peak_dp(sparse_saddle_dens%maxval)
  call boundary_peak_int(sparse_saddle_dens%maxloc)
  call boundary_peak_dp(max_dens)
  call boundary_peak_int(new_peak)
  call boundary_peak_int(alive)

  ! Count surviving peaks
  nsurvive=0
  do ipeak=1,npeaks
     if(alive(ipeak)>0)then
        nsurvive=nsurvive+1
     endif
  end do
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(nsurvive,nsurvive_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  nsurvive=nsurvive_all
#endif
  if(myid==1)then
     if(action.EQ.'relevance')then
        write(*,*)'Found',nsurvive,' relevant peaks'
     endif
     if(action.EQ.'saddleden')then
        write(*,*)'Found',nsurvive,' halos'
     endif
  endif

  if(action.EQ.'relevance')then

     ! Compute relevance
     do ipeak=1,npeaks
        if(alive(ipeak)>0)then
           if (sparse_saddle_dens%maxval(ipeak)>0)then
              relevance_peak=max_dens(ipeak)/sparse_saddle_dens%maxval(ipeak)
           else
              relevance_peak=max_dens(ipeak)/density_threshold
           end if
           relevance(ipeak)=relevance_peak
        else
           relevance(ipeak)=0
        endif
     end do

     ! Merge all peaks to deepest level
     do ilev=idepth-2,0,-1
        do ipeak=1,npeaks
           if(lev_peak(ipeak)==ilev)then
              merge_to=new_peak(ipeak)
              call get_local_peak_id(merge_to,jpeak)
              new_peak(ipeak)=new_peak(jpeak)
           endif
        end do
        call build_peak_communicator
        call boundary_peak_int(new_peak)
     end do

     ! Update flag2 field
     do ipart=1,ntest
        if (flag2(icellp(ipart))>0)then
           call get_local_peak_id(flag2(icellp(ipart)),ipeak)
           merge_to=new_peak(ipeak)
           call get_local_peak_id(merge_to,jpeak)
           flag2(icellp(ipart))=merge_to
        end if
     end do
     call build_peak_communicator

  endif

  if(action.EQ.'saddleden')then

     ! Compute peak index for the halo
     do ipeak=1,npeaks
        ind_halo(ipeak)=new_peak(ipeak)
     end do
     call boundary_peak_int(ind_halo)
     do ilev=idepth-2,0,-1
        do ipeak=1,npeaks
           if(lev_peak(ipeak)==ilev)then
              merge_to=ind_halo(ipeak)
              call get_local_peak_id(merge_to,jpeak)
              ind_halo(ipeak)=ind_halo(jpeak)
           endif
        end do
        call build_peak_communicator
        call boundary_peak_int(ind_halo)
     end do

     ! Compute halo masses
     halo_mass=0
     n_cells_halo=0
     do ipeak=1,npeaks
        merge_to=ind_halo(ipeak)
        call get_local_peak_id(merge_to,jpeak)
        halo_mass(jpeak)=halo_mass(jpeak)+clump_mass(ipeak)
        n_cells_halo(jpeak)=n_cells_halo(jpeak)+n_cells(ipeak)
     end do
     call build_peak_communicator
     call virtual_peak_dp(halo_mass,'sum')
     call boundary_peak_dp(halo_mass)
     call virtual_peak_int(n_cells_halo,'sum')
     call boundary_peak_int(n_cells_halo)
     ! Assign back halo mass to peak
     do ipeak=1,npeaks
        merge_to=ind_halo(ipeak)
        call get_local_peak_id(merge_to,jpeak)
        halo_mass(ipeak)=halo_mass(jpeak)
     end do

  endif

end subroutine merge_clumps
!################################################################
!################################################################
!################################################################
!################################################################
subroutine get_max(i,mat)
  use amr_commons,ONLY:myid
  use sparse_matrix
  use clfind_commons,ONLY: npeaks,ipeak_start,gkey
  type(sparse_mat)::mat
  integer::i
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! get maximum in i-th line by walking the linked list
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer::current,icol

  mat%maxval(i)=0
  mat%maxloc(i)=0

  ! walk the line...
  current=mat%first(i)
  do  while( current /= 0 )
     if(mat%maxval(i)<mat%val(current))then
        mat%maxval(i)=mat%val(current)
        icol=mat%col(current)
        if(icol<=npeaks)then
           mat%maxloc(i)=ipeak_start(myid)+icol
        else
           mat%maxloc(i)=gkey(icol)
        endif
     end if
     current=mat%next(current)
  end do
end subroutine get_max
!################################################################
!################################################################
!################################################################
!################################################################
subroutine allocate_peak_patch_arrays
  use amr_commons, ONLY:ndim,dp
  use clfind_commons
  use sparse_matrix
  implicit none

  integer::bit_length,ncode,ipart,peak_nr
  integer::ipeak

  !-------------------------------
  ! Allocate peak-patch_properties
  !-------------------------------
  allocate(n_cells(1:npeaks_max))
  allocate(n_cells_halo(1:npeaks_max))
  allocate(lev_peak(1:npeaks_max))
  allocate(new_peak(npeaks_max))
  allocate(clump_size(1:npeaks_max,1:ndim))
  allocate(peak_pos(1:npeaks_max,1:ndim))
  allocate(center_of_mass(1:npeaks_max,1:ndim))
  allocate(min_dens(1:npeaks_max))
  allocate(av_dens(1:npeaks_max))
  allocate(ind_halo(1:npeaks_max))
  allocate(halo_mass(1:npeaks_max))
  allocate(clump_mass(1:npeaks_max))
  allocate(clump_vol(1:npeaks_max))
  allocate(relevance(1:npeaks_max))
  call sparse_initialize(npeaks_max,sparse_saddle_dens)

  ! These arrays are not used by the clump finder
  allocate(clump_velocity(1:npeaks_max,1:ndim))
  allocate(clump_mass4(npeaks_max))
  allocate(clump_star4(npeaks_max))
  allocate(kinetic_support(npeaks_max))
  allocate(thermal_support(npeaks_max))
  allocate(magnetic_support(npeaks_max))
  allocate(Psurf(npeaks_max))
  allocate(MagPsurf(npeaks_max))
  allocate(MagTsurf(npeaks_max))
  allocate(grav_term(npeaks_max))
  allocate(rad_term(npeaks_max))
  allocate(contracting(npeaks_max))
  allocate(Icl(npeaks_max))
  allocate(Icl_d(npeaks_max))
  allocate(Icl_dd(npeaks_max))
  allocate(Icl_d_3by3(npeaks_max,1:3,1:3))
  allocate(Icl_3by3(npeaks_max,1:3,1:3))

  !--------------------
  ! Allocate hash table
  !--------------------
  ncode=npeaks_max-npeaks
  do bit_length=1,32
     ncode=ncode/2
     if(ncode<=1) exit
  end do
  nhash=prime(bit_length+1)
  allocate(hkey(1:nhash))
  hfree=npeaks+1
  hcollision=0
  allocate(gkey(npeaks+1:npeaks_max))
  allocate(nkey(npeaks+1:npeaks_max))
  hkey=0; gkey=0; nkey=0


  !------------------------------------------------
  ! Initialize the hash table with interior patches
  !------------------------------------------------
  do ipart=1,ntest
     peak_nr=flag2(icellp(ipart)) ! global peak id
     if (peak_nr>0)then
        call get_local_peak_id(peak_nr,ipeak)
     end if
  end do

  !------------------------------------------------
  ! Initialize all peak based arrays for clump finder
  !------------------------------------------------
  lev_peak=0; new_peak=0; ind_halo=0; relevance=1

  !------------------------------------------------
  ! Allocate peak communicator arrays
  !------------------------------------------------
  allocate(peak_send_cnt(1:ncpu),peak_send_oft(1:ncpu))
  allocate(peak_recv_cnt(1:ncpu),peak_recv_oft(1:ncpu))

end subroutine allocate_peak_patch_arrays
!################################################################
!################################################################
!################################################################
!################################################################
subroutine deallocate_all
  use clfind_commons
  use amr_commons, only:smbh
  use sparse_matrix
  implicit none

  deallocate(n_cells)
  deallocate(n_cells_halo)
  deallocate(lev_peak)
  deallocate(new_peak)
  deallocate(clump_size)
  deallocate(peak_pos)
  deallocate(center_of_mass)
  deallocate(min_dens)
  deallocate(av_dens)
  deallocate(max_dens)
  deallocate(peak_cell, peak_cell_level)
  deallocate(ind_halo)
  deallocate(halo_mass)
  deallocate(clump_mass)
  deallocate(clump_vol)
  deallocate(relevance)
  call sparse_kill(sparse_saddle_dens)

  deallocate(clump_mass4)
  deallocate(clump_star4)
  deallocate(clump_velocity)
  deallocate(grav_term,rad_term)
  deallocate(thermal_support,kinetic_support,magnetic_support)
  deallocate(Psurf,MagPsurf,MagTsurf)
  deallocate(contracting)
  deallocate(Icl_dd,Icl_d,Icl,Icl_d_3by3,Icl_3by3)

  deallocate(hkey,gkey,nkey)

  deallocate(peak_send_cnt,peak_send_oft)
  deallocate(peak_recv_cnt,peak_recv_oft)

end subroutine deallocate_all
!################################################################
!################################################################
!################################################################
!################################################################
subroutine get_local_peak_id(global_peak_id,local_peak_id)
  use amr_commons
  use clfind_commons
  implicit none
  integer::global_peak_id,local_peak_id

  integer::ihash,ikey,jkey

  if(    global_peak_id> ipeak_start(myid) .and. &
       & global_peak_id<=ipeak_start(myid)+npeaks_per_cpu(myid))then
     local_peak_id=global_peak_id-ipeak_start(myid)
  else
     ihash=MOD(global_peak_id,nhash)+1 ! compute the simple prime hash key
     if(hkey(ihash)==0)then ! hash table is empty
        hkey(ihash)=hfree
        local_peak_id=hfree
        gkey(hfree)=global_peak_id
        hfree=hfree+1
        if(hfree.eq.npeaks_max)then
           write(*,*)'Too many peaks'
           write(*,*)'Increase npeaks_max'
           stop
        endif
     else
        ikey=hkey(ihash) ! collision in the hash table
        do while(ikey>0)
           jkey=ikey
           if(gkey(ikey)==global_peak_id)exit
           ikey=nkey(ikey)
        end do
        if(ikey==0)then ! peak doesn't already exist
           nkey(jkey)=hfree
           local_peak_id=hfree
           gkey(hfree)=global_peak_id
           hfree=hfree+1
           hcollision=hcollision+1
           if(hfree.eq.npeaks_max)then
              write(*,*)'Too many peaks'
              write(*,*)'Increase npeaks_max'
              stop
           endif
        else            ! peak already exists
           local_peak_id=ikey
        end if
     end if
  end if

end subroutine get_local_peak_id
!################################################################
!################################################################
!################################################################
!################################################################
subroutine get_local_peak_cpu(local_peak_id,peak_cpu)
  use amr_commons
  use clfind_commons
  implicit none
  integer::local_peak_id,peak_cpu
  integer::icpu,global_peak_id

  ! get the mpi-domain a peak belongs to from its LOCAL id

  if(local_peak_id <= npeaks) then
     peak_cpu=myid
  else
     global_peak_id=gkey(local_peak_id)
     peak_cpu=ncpu
     do icpu=1,ncpu
        if(    global_peak_id> ipeak_start(icpu) .and. &
             & global_peak_id<=ipeak_start(icpu)+npeaks_per_cpu(icpu))then
           peak_cpu=icpu
        endif
     end do
  end if

end subroutine get_local_peak_cpu
!################################################################
!################################################################
!################################################################
!################################################################
subroutine build_peak_communicator
  use amr_commons
  use clfind_commons
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info,ipeak,icpu
  integer,dimension(1:ncpu)::ipeak_alltoall

  peak_send_cnt=0; peak_recv_cnt=0
  do ipeak=npeaks+1,hfree-1
     call get_local_peak_cpu(ipeak,icpu)
     peak_send_cnt(icpu)=peak_send_cnt(icpu)+1
  end do
  call MPI_ALLTOALL(peak_send_cnt,1,MPI_INTEGER,peak_recv_cnt,1,MPI_INTEGER,MPI_COMM_WORLD,info)

  peak_send_oft=0; peak_send_tot=0
  peak_recv_oft=0; peak_recv_tot=0
  do icpu=1,ncpu
     peak_send_tot=peak_send_tot+peak_send_cnt(icpu)
     peak_recv_tot=peak_recv_tot+peak_recv_cnt(icpu)
     if(icpu<ncpu)then
        peak_send_oft(icpu+1)=peak_send_oft(icpu)+peak_send_cnt(icpu)
        peak_recv_oft(icpu+1)=peak_recv_oft(icpu)+peak_recv_cnt(icpu)
     endif
  end do

  if(allocated(peak_send_buf))then
     deallocate(peak_send_buf,peak_recv_buf)
  endif
  allocate(peak_send_buf(1:peak_send_tot))
  allocate(peak_recv_buf(1:peak_recv_tot))
  ipeak_alltoall=0
  do ipeak=npeaks+1,hfree-1
     call get_local_peak_cpu(ipeak,icpu)
     ipeak_alltoall(icpu)=ipeak_alltoall(icpu)+1
     peak_send_buf(peak_send_oft(icpu)+ipeak_alltoall(icpu))=gkey(ipeak)
  end do
  call MPI_ALLTOALLV(peak_send_buf,peak_send_cnt,peak_send_oft,MPI_INTEGER, &
       &             peak_recv_buf,peak_recv_cnt,peak_recv_oft,MPI_INTEGER,MPI_COMM_WORLD,info)
#endif
end subroutine build_peak_communicator
!################################################################
!################################################################
!################################################################
!################################################################
subroutine virtual_peak_int(xx,action)
  use amr_commons
  use clfind_commons
  use mpi_mod
  implicit none
  integer,dimension(1:npeaks_max)::xx
  character(len=3)::action
#ifndef WITHOUTMPI
  integer,allocatable,dimension(:)::int_peak_send_buf,int_peak_recv_buf
  integer::ipeak,icpu,info,j
  integer,dimension(1:ncpu)::ipeak_alltoall

  allocate(int_peak_send_buf(1:peak_send_tot))
  allocate(int_peak_recv_buf(1:peak_recv_tot))
  ipeak_alltoall=0
  do ipeak=npeaks+1,hfree-1
     call get_local_peak_cpu(ipeak,icpu)
     ipeak_alltoall(icpu)=ipeak_alltoall(icpu)+1
     int_peak_send_buf(peak_send_oft(icpu)+ipeak_alltoall(icpu))=xx(ipeak)
  end do
  call MPI_ALLTOALLV(int_peak_send_buf,peak_send_cnt,peak_send_oft,MPI_INTEGER, &
       &             int_peak_recv_buf,peak_recv_cnt,peak_recv_oft,MPI_INTEGER,MPI_COMM_WORLD,info)
  select case (action)
  case('sum')
     do j=1,peak_recv_tot
        ipeak=peak_recv_buf(j)-ipeak_start(myid)
        xx(ipeak)=xx(ipeak)+int_peak_recv_buf(j)
     end do
  case('min')
     do j=1,peak_recv_tot
        ipeak=peak_recv_buf(j)-ipeak_start(myid)
        xx(ipeak)=MIN(xx(ipeak),int_peak_recv_buf(j))
     end do
  case('max')
     do j=1,peak_recv_tot
        ipeak=peak_recv_buf(j)-ipeak_start(myid)
        xx(ipeak)=MAX(xx(ipeak),int_peak_recv_buf(j))
     end do
  end select
  deallocate(int_peak_send_buf,int_peak_recv_buf)
#endif
end subroutine virtual_peak_int
!################################################################
!################################################################
!################################################################
!################################################################
subroutine virtual_peak_dp(xx,action)
  use amr_commons
  use clfind_commons
  use mpi_mod
  implicit none
  real(dp),dimension(1:npeaks_max)::xx
  character(len=3)::action
#ifndef WITHOUTMPI
  real(kind=8),allocatable,dimension(:)::dp_peak_send_buf,dp_peak_recv_buf
  integer::ipeak,icpu,info,j
  integer,dimension(1:ncpu)::ipeak_alltoall

  allocate(dp_peak_send_buf(1:peak_send_tot))
  allocate(dp_peak_recv_buf(1:peak_recv_tot))
  ipeak_alltoall=0
  do ipeak=npeaks+1,hfree-1
     call get_local_peak_cpu(ipeak,icpu)
     ipeak_alltoall(icpu)=ipeak_alltoall(icpu)+1
     dp_peak_send_buf(peak_send_oft(icpu)+ipeak_alltoall(icpu))=xx(ipeak)
  end do
  call MPI_ALLTOALLV(dp_peak_send_buf,peak_send_cnt,peak_send_oft,MPI_DOUBLE_PRECISION, &
       &             dp_peak_recv_buf,peak_recv_cnt,peak_recv_oft,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,info)
  select case (action)
  case('sum')
     do j=1,peak_recv_tot
        ipeak=peak_recv_buf(j)-ipeak_start(myid)
        xx(ipeak)=xx(ipeak)+dp_peak_recv_buf(j)
     end do
  case('min')
     do j=1,peak_recv_tot
        ipeak=peak_recv_buf(j)-ipeak_start(myid)
        xx(ipeak)=MIN(xx(ipeak),dp_peak_recv_buf(j))
     end do
  case('max')
     do j=1,peak_recv_tot
        ipeak=peak_recv_buf(j)-ipeak_start(myid)
        xx(ipeak)=MAX(xx(ipeak),dp_peak_recv_buf(j))
     end do
  end select
  deallocate(dp_peak_send_buf,dp_peak_recv_buf)
#endif
end subroutine virtual_peak_dp
!################################################################
!################################################################
!################################################################
!################################################################
subroutine virtual_saddle_max
  use amr_commons
  use clfind_commons
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info,icpu
  real(kind=8),allocatable,dimension(:)::dp_peak_send_buf,dp_peak_recv_buf
  integer,allocatable,dimension(:)::int_peak_send_buf,int_peak_recv_buf
  integer::ipeak,jpeak,j
  integer,dimension(1:ncpu)::ipeak_alltoall

  allocate(int_peak_send_buf(1:peak_send_tot))
  allocate(int_peak_recv_buf(1:peak_recv_tot))
  allocate(dp_peak_send_buf(1:peak_send_tot))
  allocate(dp_peak_recv_buf(1:peak_recv_tot))
  ipeak_alltoall=0
  do ipeak=npeaks+1,hfree-1
     call get_local_peak_cpu(ipeak,icpu)
     ipeak_alltoall(icpu)=ipeak_alltoall(icpu)+1
     dp_peak_send_buf(peak_send_oft(icpu)+ipeak_alltoall(icpu))=sparse_saddle_dens%maxval(ipeak)
     int_peak_send_buf(peak_send_oft(icpu)+ipeak_alltoall(icpu))=sparse_saddle_dens%maxloc(ipeak)
  end do
  call MPI_ALLTOALLV(dp_peak_send_buf,peak_send_cnt,peak_send_oft,MPI_DOUBLE_PRECISION, &
       &             dp_peak_recv_buf,peak_recv_cnt,peak_recv_oft,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,info)
  call MPI_ALLTOALLV(int_peak_send_buf,peak_send_cnt,peak_send_oft,MPI_INTEGER, &
       &             int_peak_recv_buf,peak_recv_cnt,peak_recv_oft,MPI_INTEGER,MPI_COMM_WORLD,info)
  do j=1,peak_recv_tot
     ipeak=peak_recv_buf(j)-ipeak_start(myid)
     if(sparse_saddle_dens%maxval(ipeak)<dp_peak_recv_buf(j))then
        sparse_saddle_dens%maxval(ipeak)=dp_peak_recv_buf(j)
        sparse_saddle_dens%maxloc(ipeak)=int_peak_recv_buf(j)
        call get_local_peak_id(int_peak_recv_buf(j),jpeak)
     endif
  end do
  deallocate(dp_peak_send_buf,dp_peak_recv_buf)
  deallocate(int_peak_send_buf,int_peak_recv_buf)
#endif
end subroutine virtual_saddle_max
!################################################################
!################################################################
!################################################################
!################################################################
subroutine boundary_peak_int(xx)
  use amr_commons
  use clfind_commons
  use mpi_mod
  implicit none
  integer,dimension(1:npeaks_max)::xx
#ifndef WITHOUTMPI
  integer,allocatable,dimension(:)::int_peak_send_buf,int_peak_recv_buf
  integer::ipeak,icpu,info,j
  integer,dimension(1:ncpu)::ipeak_alltoall

  allocate(int_peak_send_buf(1:peak_send_tot))
  allocate(int_peak_recv_buf(1:peak_recv_tot))
  do j=1,peak_recv_tot
     ipeak=peak_recv_buf(j)-ipeak_start(myid)
     int_peak_recv_buf(j)=xx(ipeak)
  end do
  call MPI_ALLTOALLV(int_peak_recv_buf,peak_recv_cnt,peak_recv_oft,MPI_INTEGER, &
       &             int_peak_send_buf,peak_send_cnt,peak_send_oft,MPI_INTEGER,MPI_COMM_WORLD,info)
  ipeak_alltoall=0
  do ipeak=npeaks+1,hfree-1
     call get_local_peak_cpu(ipeak,icpu)
     ipeak_alltoall(icpu)=ipeak_alltoall(icpu)+1
     xx(ipeak)=int_peak_send_buf(peak_send_oft(icpu)+ipeak_alltoall(icpu))
  end do
  deallocate(int_peak_send_buf,int_peak_recv_buf)
#endif
end subroutine boundary_peak_int
!################################################################
!################################################################
!################################################################
!################################################################
subroutine boundary_peak_dp(xx)
  use amr_commons
  use clfind_commons
  use mpi_mod
  implicit none
  real(dp),dimension(1:npeaks_max)::xx
#ifndef WITHOUTMPI
  real(kind=8),allocatable,dimension(:)::dp_peak_send_buf,dp_peak_recv_buf
  integer::ipeak,icpu,info,j
  integer,dimension(1:ncpu)::ipeak_alltoall

  allocate(dp_peak_send_buf(1:peak_send_tot))
  allocate(dp_peak_recv_buf(1:peak_recv_tot))
  do j=1,peak_recv_tot
     ipeak=peak_recv_buf(j)-ipeak_start(myid)
     dp_peak_recv_buf(j)=xx(ipeak)
  end do
  call MPI_ALLTOALLV(dp_peak_recv_buf,peak_recv_cnt,peak_recv_oft,MPI_DOUBLE_PRECISION, &
       &             dp_peak_send_buf,peak_send_cnt,peak_send_oft,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,info)
  ipeak_alltoall=0
  do ipeak=npeaks+1,hfree-1
     call get_local_peak_cpu(ipeak,icpu)
     ipeak_alltoall(icpu)=ipeak_alltoall(icpu)+1
     xx(ipeak)=dp_peak_send_buf(peak_send_oft(icpu)+ipeak_alltoall(icpu))
  end do
  deallocate(dp_peak_send_buf,dp_peak_recv_buf)
#endif
end subroutine boundary_peak_dp
!################################################################
!################################################################
!################################################################
!################################################################
subroutine write_clump_map
  use amr_commons
  use clfind_commons
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::dummy_io,info2
  integer,parameter::tag=1102
#endif
  !---------------------------------------------------------------------------
  ! This routine writes a csv-file of cell center coordinates and clump number
  ! for each cell which is in a clump. Makes only sense to be called when the
  ! clump finder is called at output-writing and not for sink-formation.
  !---------------------------------------------------------------------------

  integer::ind,grid,ix,iy,iz,ipart,nx_loc,peak_nr
  real(dp)::scale,dx
  real(dp),dimension(1:3)::xcell,skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  character(LEN=5)::myidstring,nchar,ncharcpu

  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)

  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     xc(ind,1)=(dble(ix)-0.5D0)
     xc(ind,2)=(dble(iy)-0.5D0)
     xc(ind,3)=(dble(iz)-0.5D0)
  end do

  !prepare file output for peak map
  ! Wait for the token
#ifndef WITHOUTMPI
     if(IOGROUPSIZE>0) then
        if (mod(myid-1,IOGROUPSIZE)/=0) then
           call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
                & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
        end if
     endif
#endif

  call title(ifout,nchar)
  call title(myid,myidstring)
  if(IOGROUPSIZEREP>0)then
     call title(((myid-1)/IOGROUPSIZEREP)+1,ncharcpu)
     open(unit=20,file=TRIM('output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/clump_map.csv'//myidstring),form='formatted')
  else
     open(unit=20,file=TRIM('output_'//TRIM(nchar)//'/clump_map.csv'//myidstring),form='formatted')
  endif
  !loop parts
  do ipart=1,ntest
     peak_nr=flag2(icellp(ipart))
     if (peak_nr /=0 ) then
        ! Cell coordinates
        ind=(icellp(ipart)-ncoarse-1)/ngridmax+1 ! cell position
        grid=icellp(ipart)-ncoarse-(ind-1)*ngridmax ! grid index
        dx=0.5D0**levp(ipart)
        xcell(1:ndim)=(xg(grid,1:ndim)+xc(ind,1:ndim)*dx-skip_loc(1:ndim))*scale
        !peak_map
        write(20,'(1PE18.9E2,A,1PE18.9E2,A,1PE18.9E2A,I4,A,I8)')xcell(1),',',xcell(2),',',xcell(3),',',levp(ipart),',',peak_nr

     end if
  end do
  close(20)

  ! Send the token
#ifndef WITHOUTMPI
     if(IOGROUPSIZE>0) then
        if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
           dummy_io=1
           call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
                & MPI_COMM_WORLD,info2)
        end if
     endif
#endif


end subroutine write_clump_map
!################################################################
!################################################################
!################################################################
subroutine analyze_peak_memory
  use amr_commons
  use clfind_commons
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  integer::i,j
  integer,dimension(1:ncpu)::npeak_all,npeak_tot
  integer,dimension(1:ncpu)::hfree_all,hfree_tot
  integer,dimension(1:ncpu)::sparse_all,sparse_tot
  integer,dimension(1:ncpu)::coll_all,coll_tot

  npeak_all=0
  npeak_all(myid)=npeaks
  coll_all=0
  coll_all(myid)=hcollision
  hfree_all=0
  hfree_all(myid)=hfree-npeaks
  sparse_all=0
  sparse_all(myid)=sparse_saddle_dens%used
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(npeak_all,npeak_tot,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(coll_all,coll_tot,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(hfree_all,hfree_tot,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(sparse_all,sparse_tot,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#else
  npeak_tot=npeak_all
  coll_tot=coll_all
  hfree_tot=hfree_all
  sparse_tot=sparse_all
#endif
  if(myid==1)then
     write(*,*)'peaks per cpu'
     do i=0,ncpu-1,10
        write(*,'(256(I8,1X))')(npeak_tot(j),j=i+1,min(i+10,ncpu))
     end do
     write(*,*)'ghost peaks per cpu'
     do i=0,ncpu-1,10
        write(*,'(256(I8,1X))')(hfree_tot(j),j=i+1,min(i+10,ncpu))
     end do
     write(*,*)'hash table collisions'
     do i=0,ncpu-1,10
        write(*,'(256(I8,1X))')(coll_tot(j),j=i+1,min(i+10,ncpu))
     end do
     write(*,*)'sparse matrix used'
     do i=0,ncpu-1,10
        write(*,'(256(I8,1X))')(sparse_tot(j),j=i+1,min(i+10,ncpu))
     end do
  end if
end subroutine analyze_peak_memory
!################################################################
!################################################################
!################################################################
!################################################################
#endif
