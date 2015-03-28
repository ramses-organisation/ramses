subroutine compute_clump_properties(xx)
  use amr_commons
  use hydro_commons, ONLY:uold
  use clfind_commons
  use poisson_commons, ONLY:phi,f
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  real(dp),dimension(1:ncoarse+ngridmax*twotondim)::xx
  !----------------------------------------------------------------------------
  ! this subroutine performs a loop over all cells above the threshold and 
  ! collects the  relevant information. After some MPI communications,
  ! all necessary peak-patch properties are computed
  !----------------------------------------------------------------------------
  integer::ipart,grid,info,i,j,peak_nr,ilevel,global_peak_id,ipeak
  real(dp)::zero=0.
  !variables needed temporarily store cell properties
  real(dp)::d,vol
  real(dp),dimension(1:3)::vd
  ! variables related to the size of a cell on a given level
  real(dp)::dx,dx_loc,scale,vol_loc,tot_mass_tot
  real(dp),dimension(1:nlevelmax)::volume
  real(dp),dimension(1:3)::skip_loc,xcell
  real(dp),dimension(1:twotondim,1:3)::xc
  integer::nx_loc,ind,ix,iy,iz,idim
  logical,dimension(1:ndim)::period
  logical::periodic

  period(1)=(nx==1)
#if NDIM>1
  if(ndim>1)period(2)=(ny==1)
#endif
#if NDIM>2
  if(ndim>2)period(3)=(nz==1)
#endif

  periodic=period(1)
#if NDIM>1
  if(ndim>1)periodic=periodic.or.period(2)
#endif
#if NDIM>2
  if(ndim>2)periodic=periodic.or.period(3)
#endif
  !peak-patch related arrays before sharing information with other cpus

  min_dens=huge(zero); max_dens=0.d0; av_dens=0d0
  n_cells=0; n_cells_halo=0 
  halo_mass=0d0; clump_mass=0.d0; clump_vol=0.d0
  center_of_mass=0.d0; clump_velocity=0.d0; clump_force=0.d0 
  peak_pos=0.

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
        if(ivar_clump==0)then
           d=xx(icellp(ipart))
        endif
        if(hydro)then
           d=uold(icellp(ipart),1)
        endif

        ! Cell volume
        vol=volume(levp(ipart))

        ! Number of leaf cells per clump
        n_cells(peak_nr)=n_cells(peak_nr)+1
        
        ! Min density
        min_dens(peak_nr)=min(d,min_dens(peak_nr))

        ! Max density and peak location
        if(d>=max_dens(peak_nr))then
           max_dens(peak_nr)=d
           av_dens(peak_nr)=d
           peak_pos(peak_nr,1:ndim)=xcell(1:ndim)
        end if

        ! Clump mass
        clump_mass(peak_nr)=clump_mass(peak_nr)+vol*d

        ! Clump volume
        clump_vol(peak_nr)=clump_vol(peak_nr)+vol
        
        ! center of mass location
        center_of_mass(peak_nr,1:3)=center_of_mass(peak_nr,1:3)+vol*d*xcell(1:3)

        ! center of mass velocity
        if (hydro)clump_velocity(peak_nr,1:3)=clump_velocity(peak_nr,1:3)+vol*uold(icellp(ipart),2:4)

        ! average grav acceleration
        clump_force(peak_nr,1:3)=clump_force(peak_nr,1:3)+vol*d*f(icellp(ipart),1:3)

     end if
  end do

  call build_peak_communicator
  ! MPI communication to collect the results from the different cpus
#ifndef WITHOUTMPI     
  call virtual_peak_int(n_cells,'sum')
  call virtual_peak_dp(min_dens,'min')
  call virtual_peak_dp(max_dens,'max')
  call virtual_peak_dp(clump_mass,'sum')
  call virtual_peak_dp(clump_vol,'sum')
  call boundary_peak_dp(max_dens)
  ! Reset position for false local peaks
  do ipeak=1,hfree-1
     if(av_dens(ipeak)<max_dens(ipeak))then
        peak_pos(ipeak,1:ndim)=0.d0 
     endif
  end do
  do i=1,ndim
     call virtual_peak_dp(peak_pos(1,i),'sum')
  end do

  do i=1,ndim
     call virtual_peak_dp(center_of_mass(1,i),'sum')
     call virtual_peak_dp(clump_velocity(1,i),'sum')
     call virtual_peak_dp(clump_force(1,i),'sum')
  end do
  do ipeak=1,npeaks
     if (relevance(ipeak)>0.)then
        center_of_mass(ipeak,1:3)=center_of_mass(ipeak,1:3)/clump_mass(ipeak)
        clump_velocity(ipeak,1:3)=clump_velocity(ipeak,1:3)/clump_mass(ipeak)
        clump_force(ipeak,1:3)=clump_force(ipeak,1:3)/clump_mass(ipeak)
     end if
  end do
  do i=1,ndim
     call boundary_peak_dp(center_of_mass(1,i))
     call boundary_peak_dp(clump_velocity(1,i))
     call boundary_peak_dp(clump_force(1,i))
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
     if (relevance(ipeak)>0.)then
        av_dens(ipeak)=clump_mass(ipeak)/clump_vol(ipeak)
     end if
  end do

  !for periodic boxes the center of mass can be meaningless at that stage
  ! -> recompute center of mass relative to peak position
  
  if(periodic)then
     
     center_of_mass=0.d0;
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
           if(ivar_clump==0)then
              d=xx(icellp(ipart))
           endif
           if(hydro)then
              d=uold(icellp(ipart),1)
           endif

           ! Cell volume
           vol=volume(levp(ipart))

           ! center of mass location
           center_of_mass(peak_nr,1:3)=center_of_mass(peak_nr,1:3)+vol*d*xcell(1:3)

        end if
     end do

     call build_peak_communicator
     ! MPI communication to collect the results from the different cpus
#ifndef WITHOUTMPI     
     do i=1,ndim
        call virtual_peak_dp(center_of_mass(1,i),'sum')
     end do
     do ipeak=1,npeaks
        if (relevance(ipeak)>0.)then
           center_of_mass(ipeak,1:3)=center_of_mass(ipeak,1:3)/clump_mass(ipeak)
        end if
     end do
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
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  logical::to_file

  !---------------------------------------------------------------------------
  ! this routine writes the clump properties to screen and to file
  !---------------------------------------------------------------------------

  integer::i,j,jj,ilun,ilun2,n_rel,n_rel_tot,info,nx_loc
  real(dp)::rel_mass,rel_mass_tot,scale,particle_mass,particle_mass_tot
  character(LEN=80)::fileloc
  character(LEN=5)::nchar
  real(dp),dimension(1:npeaks)::peakd
  integer,dimension(1:npeaks)::ind_sort

  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  if(hydro)then
     particle_mass=mass_sph
  else
     particle_mass=MINVAL(mp, MASK=(mp.GT.0.))
#ifndef WITHOUTMPI  
     call MPI_ALLREDUCE(particle_mass,particle_mass_tot,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,info)
     particle_mass=particle_mass_tot  
#endif
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
  rel_mass=0.
  n_rel=0
  if (to_file .eqv. .true.) then
     call title(ifout-1,nchar)
     fileloc=TRIM('output_'//TRIM(nchar)//'/clump_'//TRIM(nchar)//'.txt')
     call title(myid,nchar)
     fileloc=TRIM(fileloc)//TRIM(nchar)
     open(unit=ilun,file=fileloc,form='formatted')
     if(saddle_threshold>0)then
        call title(ifout-1,nchar)
        fileloc=TRIM('output_'//TRIM(nchar)//'/halo_'//TRIM(nchar)//'.txt')
        call title(myid,nchar)
        fileloc=TRIM(fileloc)//TRIM(nchar)
        open(unit=ilun2,file=fileloc,form='formatted')
     endif
  end if
  
  if (to_file .or. myid==1)then
     write(ilun,'(135A)')'   index  lev   parent      ncell    peak_x             peak_y             peak_z     '//&
          '        rho-               rho+               rho_av             mass_cl            relevance   '
     if(saddle_threshold>0)then
        write(ilun2,'(135A)')'     index      ncell    peak_x             peak_y             peak_z     '//&
             '        rho+               mass      '
     endif
  end if
  
  do j=npeaks,1,-1
     jj=ind_sort(j)
     if (relevance(jj) > relevance_threshold .and. halo_mass(jj) > mass_threshold*particle_mass)then           
        write(ilun,'(I8,X,I2,X,I10,X,I10,8(X,1PE18.9E2))')&
             jj+ipeak_start(myid)&
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
#ifndef WITHOUTMPI  
  call MPI_ALLREDUCE(n_rel,n_rel_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  n_rel=n_rel_tot  
  call MPI_ALLREDUCE(rel_mass,rel_mass_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  rel_mass=rel_mass_tot  
#endif
  if(myid==1)then
     if(clinfo)write(*,'(A,1PE12.5)')' Total mass above threshold =',tot_mass
     if(clinfo)write(*,'(A,I10,A,1PE12.5)')' Total mass in',n_rel,' listed clumps =',rel_mass
  endif
  if (to_file)then
     close(ilun)
     if(saddle_threshold>0)then
        close(ilun2)
     endif
  end if


end subroutine write_clump_properties
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine merge_clumps(action)
  use amr_commons
  use clfind_commons
  implicit none
  character(len=9)::action
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif  
  !---------------------------------------------------------------------------
  ! This routine merges the irrelevant clumps 
  ! -clumps are sorted by ascending max density
  ! -irrelevent clumps are merged to most relevant neighbor
  !---------------------------------------------------------------------------

  integer::info,j,i,ii,merge_count,final_peak,merge_to,ipart,saddle_max_host
  integer::peak,next_peak,current,isearch,nmove,nmove_all,ipeak,jpeak,iter
  integer::nsurvive,nsurvive_all,nzero,nzero_all,idepth
  integer::jmerge,ilev,global_peak_id
  real(dp)::value_iij,zero=0.,relevance_peak,dens_max
  real(dp)::mass_threshold_sm,mass_threshold_uu
  integer,dimension(1:npeaks_max)::alive,ind_sort
  real(dp),dimension(1:npeaks_max)::peakd
  real(dp),dimension(1:2)::max_loc_input,max_loc_output
  logical::do_merge

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
        if(myid==1.and.clinfo)write(*,*)'niter=',iter,'nmove=',nmove
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
           if(alive(ipeak)==0.OR.alive(j)==0)then
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
     if(myid==1.and.clinfo)write(*,*)'level=',idepth,'nmove=',nzero,'survived=',nsurvive
     idepth=idepth+1
     
  end do
  ! End loop over peak levels
  
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
           relevance(ipeak)=0.
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
     halo_mass=0.0
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
  
  mat%maxval(i)=0.
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

  real(dp)::zero=0.
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
  allocate(saddle_max(1:npeaks_max))
  allocate(relevance(1:npeaks_max))
  call sparse_initialize(npeaks_max,npeaks_max,sparse_saddle_dens)

  ! These arrays are not used by the clump finder
  allocate(clump_velocity(1:npeaks_max,1:ndim))
  allocate(clump_force(1:npeaks_max,1:ndim))
  allocate(e_kin_int(npeaks_max))
  allocate(e_thermal(npeaks_max))
  allocate(Psurf(npeaks_max))
  allocate(clump_check(npeaks_max))
  allocate(grav_term(npeaks_max))
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

  !---------------------------------
  ! Initialize all peak based arrays
  !---------------------------------
  n_cells=0; n_cells_halo=0; lev_peak=0; new_peak=0; ind_halo=0
  saddle_max=0.; relevance=1.; clump_size=0.
  min_dens=huge(zero)
  av_dens=0.; halo_mass=0.
  clump_mass=0.; clump_vol=0.; peak_pos=0.; center_of_mass=0.

  clump_force=0.
  clump_velocity=0.
  e_kin_int=0.
  e_thermal=0.
  grav_term=0.d0
  clump_check=-1.
  contracting=.false.
  Icl=0.; Icl_d=0.; Icl_dd=0.; Icl_d_3by3=0.; Icl_3by3=0.

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
  deallocate(ind_halo)
  deallocate(halo_mass)
  deallocate(clump_mass)
  deallocate(clump_vol)
  deallocate(saddle_max)
  deallocate(relevance)
  call sparse_kill(sparse_saddle_dens)

  deallocate(clump_force)
  deallocate(clump_velocity)
  deallocate(e_kin_int)
  deallocate(grav_term)
  deallocate(e_thermal)
  deallocate(Psurf)
  deallocate(clump_check)
  deallocate(contracting)
  deallocate(Icl_dd,Icl_d,Icl,Icl_d_3by3,Icl_3by3)

  deallocate(hkey,gkey,nkey)


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
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  integer::ipeak,icpu,info,j
  integer,dimension(1:ncpu,1:ncpu)::npeak_alltoall
  integer,dimension(1:ncpu,1:ncpu)::npeak_alltoall_tot
  integer,dimension(1:ncpu)::ipeak_alltoall
#ifndef WITHOUTMPI
  npeak_alltoall=0
  do ipeak=npeaks+1,hfree-1
     call get_local_peak_cpu(ipeak,icpu)
     npeak_alltoall(myid,icpu)=npeak_alltoall(myid,icpu)+1
  end do
  call MPI_ALLREDUCE(npeak_alltoall,npeak_alltoall_tot,ncpu*ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  npeak_alltoall=npeak_alltoall_tot
!!$  if(myid==1)then
!!$     write(*,'(129(I3,1X))')0,(j,j=1,ncpu)
!!$     do icpu=1,ncpu
!!$        write(*,'(129(I3,1X))')icpu,(npeak_alltoall(icpu,j),j=1,ncpu)
!!$     end do
!!$  end if
  if(.not. allocated(peak_send_cnt))then
     allocate(peak_send_cnt(1:ncpu),peak_send_oft(1:ncpu))
     allocate(peak_recv_cnt(1:ncpu),peak_recv_oft(1:ncpu))
  endif
  peak_send_cnt=0; peak_send_oft=0; peak_send_tot=0
  peak_recv_cnt=0; peak_recv_oft=0; peak_recv_tot=0
  do icpu=1,ncpu
     peak_send_cnt(icpu)=npeak_alltoall(myid,icpu)
     peak_recv_cnt(icpu)=npeak_alltoall(icpu,myid)
     peak_send_tot=peak_send_tot+peak_send_cnt(icpu)
     peak_recv_tot=peak_recv_tot+peak_recv_cnt(icpu)
     if(icpu<ncpu)then
        peak_send_oft(icpu+1)=peak_send_oft(icpu)+npeak_alltoall(myid,icpu)
        peak_recv_oft(icpu+1)=peak_recv_oft(icpu)+npeak_alltoall(icpu,myid)
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
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif  
  integer,dimension(1:npeaks_max)::xx
  character(len=3)::action

  integer,allocatable,dimension(:)::int_peak_send_buf,int_peak_recv_buf
  integer::ipeak,icpu,info,j
  integer,dimension(1:ncpu)::ipeak_alltoall
#ifndef WITHOUTMPI
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
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  real(dp),dimension(1:npeaks_max)::xx
  character(len=3)::action
  
  real(kind=8),allocatable,dimension(:)::dp_peak_send_buf,dp_peak_recv_buf
  integer::ipeak,icpu,info,j
  integer,dimension(1:ncpu)::ipeak_alltoall

#ifndef WITHOUTMPI
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
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  
  real(kind=8),allocatable,dimension(:)::dp_peak_send_buf,dp_peak_recv_buf
  integer,allocatable,dimension(:)::int_peak_send_buf,int_peak_recv_buf
  integer::ipeak,icpu,info,j,jpeak
  integer,dimension(1:ncpu)::ipeak_alltoall

#ifndef WITHOUTMPI
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
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif  
  integer,dimension(1:npeaks_max)::xx

  integer,allocatable,dimension(:)::int_peak_send_buf,int_peak_recv_buf
  integer::ipeak,icpu,info,j
  integer,dimension(1:ncpu)::ipeak_alltoall
#ifndef WITHOUTMPI
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
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  real(dp),dimension(1:npeaks_max)::xx
  
  real(kind=8),allocatable,dimension(:)::dp_peak_send_buf,dp_peak_recv_buf
  integer::ipeak,icpu,info,j
  integer,dimension(1:ncpu)::ipeak_alltoall
#ifndef WITHOUTMPI
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
  implicit none

  !---------------------------------------------------------------------------
  ! This routine writes a csv-file of cell center coordinates and clump number
  ! for each cell which is in a clump. Makes only sense to be called when the 
  ! clump finder is called at output-writing and not for sink-formation.
  !---------------------------------------------------------------------------

  integer::ind,grid,ix,iy,iz,ipart,nx_loc,peak_nr
  real(dp)::scale,dx
  real(dp),dimension(1:3)::xcell,skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  character(LEN=5)::myidstring,nchar 

  nx_loc=(icoarse_max-icoarse_min+1)
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
  call title(ifout-1,nchar)
  call title(myid,myidstring)
  open(unit=20,file=TRIM('output_'//TRIM(nchar)//'/clump_map.csv'//myidstring),form='formatted')

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
        write(20,'(F11.8,A,F11.8,A,F11.8,A,I8)')xcell(1),',',xcell(2),',',xcell(3),',',peak_nr
     end if
  end do
  close(20)
end subroutine write_clump_map
!################################################################
!################################################################
!################################################################
subroutine analyze_peak_memory
  use amr_commons
  use clfind_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::icpu,info,i,j
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
