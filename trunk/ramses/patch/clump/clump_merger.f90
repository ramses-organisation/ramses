subroutine compute_clump_properties(xx,ntest)
  use amr_commons
  use hydro_commons, ONLY:uold
  use clfind_commons
  use poisson_commons, ONLY:phi,f
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ntest
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
  integer::nx_loc,ind,ix,iy,iz
  !peak-patch related arrays before sharing information with other cpus

  min_dens=huge(zero); max_dens=0.d0; av_dens=0d0
  n_cells=0; phi_min=huge(zero); second_moments=0.
  halo_mass=0d0; clump_mass=0.d0; clump_vol=0.d0; clump_momentum=0.d0
  center_of_mass=0.d0; clump_force=0.d0 
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
     ! peak number after merge
     global_peak_id=flag2(icellp(ipart)) 
     call get_local_peak_id(global_peak_id,peak_nr)

     if (peak_nr /=0 ) then

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

     end if
  end do
  
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
  av_dens(1:npeaks)=clump_mass(1:npeaks)/(clump_vol(1:npeaks)+tiny(0.d0))

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

  integer::i,j,jj,ilun,n_rel,n_rel_tot,info,nx_loc
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

  If(to_file)then
     ilun=20
  else 
     ilun=6
  end if

  ! print results in descending order to screen/file
  rel_mass=0.
  n_rel=0
  if (to_file .eqv. .true.) then
     call title(ifout-1,nchar)
     fileloc=TRIM('output_'//TRIM(nchar)//'/clump_'//TRIM(nchar)//'.txt')
     call title(myid,nchar)
     fileloc=TRIM(fileloc)//TRIM(nchar)
     open(unit=20,file=fileloc,form='formatted')
  end if
  
  write(ilun,'(135A)')'   index  lev   parent      ncell    peak_x             peak_y             peak_z     '//&
       '        rho-               rho+               rho_av             mass_cl            relevance   '

  do j=npeaks,1,-1
     jj=ind_sort(j)
     if (relevance(jj) > relevance_threshold .and. halo_mass(jj) > mass_threshold*particle_mass)then           
        write(ilun,'(I8,X,I2,X,I10,X,I10,10(X,1PE18.9E2),5(X,1PE11.2E2))')&
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
     close(20)
  end if

#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif

end subroutine write_clump_properties
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine merge_clumps(ntest,action)
  use amr_commons
  use clfind_commons
  implicit none
  integer::ntest
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
     ! Create new local duplicated peaks and update communicator
     call virtual_saddle_max
     call build_peak_communicator

     ! Set up bounday values
     call boundary_peak_dp(sparse_saddle_dens%maxval)
     call boundary_peak_int(sparse_saddle_dens%maxloc)
     call boundary_peak_dp(max_dens)
     call boundary_peak_int(new_peak)
          
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
        if(myid==1)write(*,*)'niter=',iter,'nmove=',nmove
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
     if(myid==1)write(*,*)'level=',idepth,'nmove=',nzero,'survived=',nsurvive
     idepth=idepth+1
     
  end do
  ! End loop over peak levels
  
  ! Compute maximum saddle density for each surviving clump
  ! Create new local duplicated peaks and update communicator
  do i=1,hfree-1
     call get_max(i,sparse_saddle_dens)
  end do
  call virtual_saddle_max
  call build_peak_communicator

  ! Set up bounday values
  call boundary_peak_dp(sparse_saddle_dens%maxval)
  call boundary_peak_int(sparse_saddle_dens%maxloc)
  call boundary_peak_dp(max_dens)
  call boundary_peak_int(new_peak)
  
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
        call get_local_peak_id(flag2(icellp(ipart)),ipeak)
        merge_to=new_peak(ipeak)
        call get_local_peak_id(merge_to,jpeak)
        flag2(icellp(ipart))=merge_to
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
     do ipeak=1,npeaks
        merge_to=ind_halo(ipeak)
        call get_local_peak_id(merge_to,jpeak)
        halo_mass(jpeak)=halo_mass(jpeak)+clump_mass(ipeak)
     end do
     call build_peak_communicator
     call virtual_peak_dp(halo_mass,'sum')
     call boundary_peak_dp(halo_mass)
     ! Assign halo mass to peak
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
  ! get maximum in line by walking the linked list
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
subroutine allocate_peak_patch_arrays(ntest)
  use amr_commons, ONLY:ndim,dp
  use clfind_commons
  use sparse_matrix
  implicit none
  integer::ntest

  real(dp)::zero=0.
  integer::bit_length,ncode,ipart,peak_nr
  integer::ipeak

  !-------------------------------
  ! Allocate peak-patch_properties
  !-------------------------------
  allocate(n_cells(1:npeaks_max))
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
  allocate(clump_momentum(1:npeaks_max,1:ndim))
  allocate(clump_force(1:npeaks_max,1:ndim))
  allocate(second_moments(1:npeaks_max,1:ndim,1:ndim)) 
  allocate(clump_mass4(1:npeaks_max))
  allocate(e_kin_int(npeaks_max))
  allocate(e_bind(npeaks_max))
  allocate(e_thermal(npeaks_max))
  allocate(phi_min(npeaks_max))
  allocate(minmatch(npeaks_max))
  allocate(phi_ref(npeaks_max))
  allocate(Psurf(npeaks_max))
  allocate(v_therm(npeaks_max))
  allocate(m4(npeaks_max))
  allocate(bulk_momentum(1:npeaks_max,1:ndim))
  allocate(e_kin_iso(npeaks_max))
  allocate(e_bind_iso(npeaks_max))
  allocate(e_therm_iso(npeaks_max))
  allocate(peak_check(npeaks_max))
  allocate(isodens_check(npeaks_max))
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
     call get_local_peak_id(peak_nr,ipeak)
  end do

  !---------------------------------
  ! Initialize all peak based arrays
  !---------------------------------
  n_cells=0; lev_peak=0; new_peak=0; ind_halo=0
  saddle_max=0.; relevance=1.; clump_size=0.
  min_dens=huge(zero)
  av_dens=0.; halo_mass=0.
  clump_mass=0.; clump_vol=0.; peak_pos=0.; center_of_mass=0.

  clump_force=0.
  second_moments=0.
  clump_momentum=0.
  clump_mass4=0.
  e_kin_int=0.
  e_bind=0.
  e_thermal=0.
  phi_min=0.
  minmatch=1
  phi_ref=huge(zero)
  Psurf=0.
  grav_term=0.d0
  isodens_check=-1.;  clump_check=-1.; peak_check=-1.; 
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
  deallocate(second_moments)
  deallocate(clump_momentum)
  deallocate(e_kin_int)
  deallocate(e_bind,grav_term)
  deallocate(e_thermal)
  deallocate(phi_min)
  deallocate(minmatch)
  deallocate(phi_ref)
  deallocate(Psurf)
  deallocate(v_therm)
  deallocate(m4,bulk_momentum)
  deallocate(e_kin_iso,e_bind_iso,e_therm_iso)
  deallocate(peak_check,isodens_check,clump_check)
  deallocate(contracting)
  deallocate(Icl_dd,Icl_d,Icl,Icl_d_3by3,Icl_3by3)
  deallocate(clump_mass4)

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
  call MPI_ALLREDUCE(npeak_all,npeak_tot,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(coll_all,coll_tot,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(hfree_all,hfree_tot,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(sparse_all,sparse_tot,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
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

end subroutine boundary_peak_dp
!################################################################
!################################################################ 
!################################################################
!################################################################     
subroutine get_phi_ref(ntest)
  use amr_commons
  use hydro_commons
  use pm_commons
  use clfind_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ntest

  !---------------------------------------------------------------
  ! This subroutine finds the minimum potential on the saddle 
  ! surface of the peak patch by looping over all "test-particles"
  ! This loop is also used to compute the surface pressure term.
  !---------------------------------------------------------------

  integer::info   
  integer::ipart,ip,ilevel,next_level,dummyzero
  integer,dimension(1:nvector)::ind_cell

  Psurf=0.
  ! loop 'testparts', pass the information of nvector parts to neighborsearch 
  ip=0
  do ipart=1,ntest
     ip=ip+1
     ilevel=levp(testp_sort(ipart)) ! level
     if (verbose.and.ilevel/=nlevelmax)print*,'not all particles in max level',ilevel
     next_level=0
     if(ipart<ntest)next_level=levp(testp_sort(ipart+1)) !level of next particle
     ind_cell(ip)=icellp(testp_sort(ipart))
     if(ip==nvector .or. next_level /= ilevel)then
        call neighborsearch(ind_cell,ip,0,dummyzero,ilevel,5)
        call surface_int(ind_cell,ip,ilevel)
        ip=0
     endif
  end do
  if (ip>0)then 
     call neighborsearch(ind_cell,ip,0,dummyzero,ilevel,5)
     call surface_int(ind_cell,ip,ilevel)
  endif
   
#ifndef WITHOUTMPI     
!!$  call MPI_ALLREDUCE(phi_ref,phi_ref_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,info)
!!$  call MPI_ALLREDUCE(Psurf,Psurf_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI     
!!$  phi_ref_tot=phi_ref
!!$  Psurf_tot=Psurf
#endif

end subroutine get_phi_ref
!################################################################
!################################################################ 
!################################################################
!################################################################     
subroutine trim_clumps(ntest)
  use amr_commons
  use clfind_commons
  use pm_commons, only:ir_cloud
  implicit none
  integer::ntest

  !---------------------------------------------------------------------------
  ! this routine trims the clumps down to the intersection of the clump with 
  ! the accretion zone of the sink. Cells that are too far away from the peak
  ! are removed from the clump by setting flag2 to 0.
  !---------------------------------------------------------------------------

  integer::ipart,nx_loc,ind
  real(dp)::dx,scale,dx_loc,r2
  integer ::ix,iy,iz,grid,peak_nr

  real(dp),dimension(1:3)::skip_loc,xcell
  real(dp),dimension(1:twotondim,1:3)::xc

  ! Mesh spacing in max level
  dx=0.5D0**nlevelmax
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     xc(ind,1)=(dble(ix)-0.5D0)
     xc(ind,2)=(dble(iy)-0.5D0)
     xc(ind,3)=(dble(iz)-0.5D0)
  end do

  !update flag 2
  do ipart=1,ntest
     peak_nr=flag2(icellp(ipart))
     if (peak_nr /=0 ) then
        ind=(icellp(ipart)-ncoarse-1)/ngridmax+1 ! cell position
        grid=icellp(ipart)-ncoarse-(ind-1)*ngridmax ! grid index
        dx=0.5D0**levp(ipart)
        xcell(1:ndim)=(xg(grid,1:ndim)+xc(ind,1:ndim)*dx-skip_loc(1:ndim))*scale
        r2=(peak_pos(peak_nr,1)-xcell(1))**2&
             +(peak_pos(peak_nr,2)-xcell(2))**2&
             +(peak_pos(peak_nr,3)-xcell(3))**2.
        if (r2 > (ir_cloud*dx_loc)**2.)then        
           !remove cell from clump
           flag2(icellp(ipart))=0
        end if
     end if
  end do

end subroutine trim_clumps
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine jacobi(A,x,err2)
  use amr_commons, only:myid,dp
  implicit none
  real(dp)::err2
  real(dp),dimension(3,3)::A,x

  !---------------------------------------------------------------------------
  ! Compute eigenvalues and eigenvectors using the jacobi-Method 
  ! as for example described in Numerical Recipes. 
  ! Returns eigenvalues as diagonal elements of A
  !---------------------------------------------------------------------------

  integer::n
  integer::i,j,k
  real(dp)::b2, bar
  real(dp)::beta, coeff, c, s, cs, sc
  
  n=3
  ! x is identity matrix initially
  x = 0.0
  do i=1,n
     x(i,i) = 1.0
  end do

  ! sum all squared off-diagonal elements 
  b2 = 0.0
  do i=1,n
     do j=1,n
        if (i.ne.j) b2 = b2 + A(i,j)**2
     end do
  end do

  if (b2 <= err2) then
     if (myid==1)write(*,*), 'returning. maybe err2 too small? ',err2
     return
  endif

  ! average for off-diagonal elements /2
  bar = 0.5*b2/9.

  do while (b2 > err2)
     do i=1,n-1
        do j=i+1,n
           if (A(j,i)**2 <= bar) cycle  ! do not touch small elements
           b2 = b2 - 2.0*A(j,i)**2
           bar = 0.5*b2/9.
           ! calculate coefficient c and s for Givens matrix
           beta = (A(j,j)-A(i,i))/(2.0*A(j,i))
           coeff = 0.5*beta*(1.0+beta**2)**(-0.5)
           s = (max(0.5+coeff,0.0))**0.5
           c = (max(0.5-coeff,0.0))**0.5
           ! update rows i and j
           do k=1,n
              cs =  c*A(i,k)+s*A(j,k)
              sc = -s*A(i,k)+c*A(j,k)
              A(i,k) = cs
              A(j,k) = sc
           end do
           ! find new matrix A_{k+1} 
           do k=1,n
              cs =  c*A(k,i)+s*A(k,j)
              sc = -s*A(k,i)+c*A(k,j)
              A(k,i) = cs
              A(k,j) = sc
              cs =  c*x(k,i)+s*x(k,j)
              sc = -s*x(k,i)+c*x(k,j)
              x(k,i) = cs
              x(k,j) = sc
           end do
        end do
     end do
  end do
end subroutine jacobi
!################################################################
!################################################################
!################################################################
!################################################################
subroutine write_clump_map(ntest)
  use amr_commons
  use clfind_commons
  implicit none
  integer::ntest

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
!################################################################
subroutine compute_clump_properties_round2(xx,ntest,all_bound)
  use amr_commons
  use hydro_commons, ONLY:uold,gamma
  use poisson_commons, ONLY:phi,f
  use clfind_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  logical::all_bound
  integer::ntest
  real(dp),dimension(1:ncoarse+ngridmax*twotondim)::xx
  !----------------------------------------------------------------------------
  ! This subroutine performs another loop over all particles and collects 
  ! more information like binding energies, etc, that can not be created by
  ! just summing up cell properties.
  !----------------------------------------------------------------------------

  integer::ipart,ilevel,info,i,peak_nr,j,ii,jj
  integer::grid,nx_loc,ix,iy,iz,ind
  real(dp)::d,vol,M,ekk,phi_rel,de,c_sound,d0,v_bulk2,p
  real(dp)::t_larson1,cont_speed=0.
  real(dp)::dx,dx_loc,scale,vol_loc,abs_err,A1=0.,A2=0.,A3=0.
  real(dp),dimension(1:nlevelmax)::volume
  real(dp),dimension(1:3)::vd,xcell,xpeak,v_cl,rrel,vrel,frel,skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:3,1:3)::eigenv,a
  
  !  first, get minimum potential on saddle surface
!  call get_phi_ref(ntest)
  
  !initialize arrays
  e_kin_int=0.d0; clump_size=0.d0; e_bind=0.d0; e_thermal=0.d0
  clump_mass4=0.d0
  v_therm=0.; bulk_momentum=0.; m4=0.
  e_kin_iso=0.; e_bind_iso=0.; e_therm_iso=0.
  grav_term=0.d0; Icl_d=0.d0; Icl=0.; Icl_dd=0.
  Icl_3by3=0.;  Icl_d_3by3=0.
  contracting=.false.

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

  !---------------------------------------------------------------------------
  ! loop over all test particles to collect information from the cells
  !---------------------------------------------------------------------------
  do ipart=1,ntest     
     ! peak number after merge
     peak_nr=flag2(icellp(ipart)) 

     if (peak_nr /=0 ) then
        
        ! Cell coordinates
        ind=(icellp(ipart)-ncoarse-1)/ngridmax+1 ! cell position
        grid=icellp(ipart)-ncoarse-(ind-1)*ngridmax ! grid index
        dx=0.5D0**levp(ipart)
        xcell(1:ndim)=(xg(grid,1:ndim)+xc(ind,1:ndim)*dx-skip_loc(1:ndim))*scale

        ! gas density and energydensity
        if(ivar_clump==0)then
           d=xx(icellp(ipart))
           de=0.0
           do i=1,ndim
              vd(i)=0.0
              xpeak(i)=peak_pos(peak_nr,i)
           end do
        endif
        if(hydro)then
           d=uold(icellp(ipart),1)
           de=uold(icellp(ipart),ndim+2) 
           do i=1,ndim
              vd(i)=uold(icellp(ipart),i+1)
              xpeak(i)=peak_pos(peak_nr,i)
           end do
        endif

        ! Cell volume
        vol=volume(levp(ipart))                  


        M=clump_mass(peak_nr)
        v_cl(1:ndim)=clump_momentum(peak_nr,1:ndim)/M
        
        !properties of the cell relative to center of mass
        rrel=xcell(1:3)-center_of_mass(peak_nr,1:3)
        vrel=vd(1:3)/d-v_cl(1:3)
        frel=f(icellp(ipart),1:3)-clump_force(peak_nr,1:3)

        do i=1,ndim
           ! size relative to center of mass
           clump_size(peak_nr,i)=clump_size(peak_nr,i)+rrel(i)**2 * vol

           ! internal kinetic energy
           e_kin_int(peak_nr)=e_kin_int(peak_nr)+vrel(i)**2*d*vol*0.5
        end do

        ! potential energy using the acutal phi W= 0.5*int phi_rel*rho
        phi_rel=(phi(icellp(ipart))-phi_ref(peak_nr))*scale
        e_bind(peak_nr)=e_bind(peak_nr)-phi_rel*d*vol*5.d-1
                
        ! thermal energy
        ekk=0.
        do i=1,3 
           ekk=ekk+0.5*vd(i)**2/d                          
        end do
        p=(de-ekk)*(gamma-1)
        e_thermal(peak_nr)=e_thermal(peak_nr)+1.5*vol*p

        ! sound speed
        c_sound=(de-ekk)/d*gamma/(gamma-1)

        !Mass weighted thermal Velocity
        v_therm(peak_nr)=v_therm(peak_nr)+c_sound*d*vol/M        
                
        !properties for regions close to peak (4 cells away)
        if (((xpeak(1)-xcell(1))**2.+(xpeak(2)-xcell(2))**2.+(xpeak(3)-xcell(3))**2.) .LE. 16.*volume(nlevelmax)**(2./3.))then
           do i=1,3
              bulk_momentum(peak_nr,i)=bulk_momentum(peak_nr,i)+(vd(i)/d-v_cl(i))*vol*(d-d0)
           end do
           m4(peak_nr)=m4(peak_nr)+(d-d0)*vol
           clump_mass4(peak_nr)=clump_mass4(peak_nr)+d*vol           
        end if

        !properties for region enclosed by isopotential surface 
        if (phi_rel<0.)then
           do i=1,3
              !not strictly correct since v_cl is av. vel of WHOLE clump
              e_kin_iso(peak_nr)=e_kin_iso(peak_nr)+(vd(i)/d-v_cl(i))**2*d*vol*0.5
           end do
           e_bind_iso(peak_nr)=e_bind_iso(peak_nr)-phi_rel*d*vol*0.5
           e_therm_iso(peak_nr)=e_therm_iso(peak_nr)+1.5*p*vol
        endif

        !terms for virial theorem analysis
        do i=1,3
           grav_term(peak_nr) = grav_term(peak_nr) + frel(i) * rrel(i) * vol*d
           Icl_d(peak_nr)     = Icl_d(peak_nr)     + vrel(i) * rrel(i) * vol*d
           Icl(peak_nr)       = Icl(peak_nr)       + rrel(i) * rrel(i) * vol*d
           do j=1,3
              Icl_d_3by3(peak_nr,i,j)=  Icl_d_3by3(peak_nr,i,j)   + ( vrel(j) * rrel(i)  +  vrel(i) * rrel(j) )   * vol*d
              Icl_3by3(peak_nr,i,j)  =  Icl_3by3(peak_nr,i,j)     +   rrel(j) * rrel(i)                           * vol*d
           end do
        end do
     end if
  end do
  !---------------------------------------------------------------------------
  ! a lot of MPI communication to collect the results from the different cpu's
  !---------------------------------------------------------------------------
#ifndef WITHOUTMPI     
!!$  call MPI_ALLREDUCE(e_kin_int,e_kin_int_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
!!$  call MPI_ALLREDUCE(e_bind,e_bind_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
!!$  call MPI_ALLREDUCE(clump_size,clump_size_tot,3*npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
!!$  call MPI_ALLREDUCE(m4,m4_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
!!$  call MPI_ALLREDUCE(bulk_momentum,bulk_momentum_tot,3*npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
!!$  call MPI_ALLREDUCE(e_thermal,e_thermal_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
!!$  call MPI_ALLREDUCE(v_therm,v_therm_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
!!$  call MPI_ALLREDUCE(e_therm_iso,e_therm_iso_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
!!$  call MPI_ALLREDUCE(e_kin_iso,e_kin_iso_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
!!$  call MPI_ALLREDUCE(e_bind_iso,e_bind_iso_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
!!$  call MPI_ALLREDUCE(grav_term,grav_term_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
!!$  call MPI_ALLREDUCE(Icl_d,Icl_d_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
!!$  call MPI_ALLREDUCE(Icl,Icl_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
!!$  call MPI_ALLREDUCE(Icl_d_3by3,Icl_d_3by3_tot,9*npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
!!$  call MPI_ALLREDUCE(Icl_3by3,Icl_3by3_tot,9*npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
!!$  call MPI_ALLREDUCE(clump_mass4,clump_mass_tot4,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI     
!!$  e_kin_int_tot=e_kin_int
!!$  e_bind_tot=e_bind
!!$  clump_size_tot=clump_size
!!$  m4_tot=m4
!!$  bulk_momentum_tot=bulk_momentum
!!$  e_thermal_tot=e_thermal
!!$  v_therm_tot=v_therm
!!$  e_bind_iso_tot=e_bind_iso
!!$  e_kin_iso_tot=e_kin_iso
!!$  e_therm_iso_tot=e_therm_iso
!!$  grav_term_tot=grav_term
!!$  Icl_d_tot=Icl_d
!!$  Icl_tot=Icl
!!$  Icl_d_3by3_tot=Icl_d_3by3
!!$  Icl_3by3_tot=Icl_3by3
!!$  clump_mass_tot4=clump_mass4
#endif

  if(myid==1)write(*,*)'Loop over cells done.'

#ifdef TOTO
  !second time derivative of I
  Icl_dd(1:npeaks)=2.*(grav_term(1:npeaks)-Psurf(1:npeaks)+2*e_kin_int(1:npeaks)+2*e_thermal(1:npeaks))

  all_bound=.true.

  do j=npeaks,1,-1
     if (relevance(j)>0.)then

        !compute eigenvalues and eigenvectors of Icl_d_3by3
        a=Icl_3by3(j,1:3,1:3)
        abs_err=1.d-8*Icl(j)**2+1.d-40
!        call jacobi(a,eigenv,abs_err)
        A1=a(1,1); A2=a(2,2); A3=a(3,3)

        !compute the contractions along the eigenvectors of Icl
        contractions(j,1:3)=0.
        do ii=1,3
           do jj=1,3
              contractions(j,1)=contractions(j,1)+Icl_d_3by3(j,ii,jj)*eigenv(1,ii)*eigenv(1,jj)
              contractions(j,2)=contractions(j,2)+Icl_d_3by3(j,ii,jj)*eigenv(2,ii)*eigenv(2,jj)
              contractions(j,3)=contractions(j,3)+Icl_d_3by3(j,ii,jj)*eigenv(3,ii)*eigenv(3,jj)
           end do
        end do

        !Check wether clump is contracting fast enough along all axis
        if (Icl(j)>0)then 
           contracting(j)=.true.
           contracting(j)=contracting(j) .and. contractions(j,1)/(A1+tiny(0.d0)) < cont_speed 
           contracting(j)=contracting(j) .and. contractions(j,2)/(A2+tiny(0.d0)) < cont_speed 
           contracting(j)=contracting(j) .and. contractions(j,3)/(A3+tiny(0.d0)) < cont_speed 
        end if
        
        !compute peak check for smbh sink formation
        v_rms=2.*e_kin_int(j)/clump_mass(j)
        v_bulk2=(bulk_momentum(j,1)**2+bulk_momentum(j,2)**2&
             +bulk_momentum(j,3)**2)/(m4(j)**2+tiny(0.d0))     
        peak_check(j)=scale*(phi_ref(j)-phi_min(j))/((v_therm(j)**2+v_rms+v_bulk2)*0.5+tiny(0.d0))

        !compute other checks (currently not needed for sink formation)
        isodens_check(j)=scale*e_bind_iso(j)/(tiny(0.d0)+2*e_kin_iso(j)+2*e_therm_iso(j))
        clump_check(j)=(-1.*grav_term(j)+Psurf(j))/(tiny(0.d0)+2*e_kin_int(j)+2*e_thermal(j))
        
        !update the all_bound property
        all_bound=all_bound.and.(isodens_check(j)>1.)

     endif
  end do

  !write to the log file some information that could be of interest for debugging etc.
  if(myid==1 .and. clinfo .and. .not. smbh .and. sink)then 
     write(*,'(135A)')'==========================================================================================='
     write(*,'(135A)')'Cl_N     t1[y]      t2[y]      t3[y] |I_d|/I_dd[y] tidal_Fg   Psurf      e_kin      e_therm'
     write(*,'(135A)')'==========================================================================================='
     do j=npeaks,1,-1
        if (relevance(j)>0.)then
           write(*,'(I4,2X,8(E8.2E2,3X))'),j&
                ,A1/(contractions(j,1)+tiny(0.d0))*cty,A2/(contractions(j,2)+tiny(0.d0))*cty,A3/(contractions(j,3)+tiny(0.d0))*cty&
                ,abs(Icl_d(j))/Icl_dd(j)*cty&
                ,grav_term(j),-1.*Psurf(j)&
                ,e_kin_int(j),e_thermal(j)
        end if
     end do
     write(*,'(135A)')'==========================================================================================='
  end if
     

#endif

end subroutine compute_clump_properties_round2
