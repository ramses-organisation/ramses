subroutine compute_clump_properties(ntest)
  use amr_commons
  use hydro_commons, ONLY:uold
  use clfind_commons
  use poisson_commons, ONLY:phi
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
   
  integer::ntest

  !----------------------------------------------------------------------------
  ! this subroutine performs one loop over all particles and collects the 
  ! relevant information from the cells where the particles sit. after a lot
  ! of mpi-communication, all necessary peak-patch properties can be computed
  !----------------------------------------------------------------------------

  ! variables used for the loop over all particles
  integer::ipart,ilevel,grid
  integer::info,i,j,peak_nr

  !variables needed temporarily store cell properties
  real(kind=8)::d,vol
  real(kind=8),dimension(1:3)::vd


  ! variables related to the size of a cell on a given level
  real(kind=8)::dx,dx_loc,scale,vol_loc
  real(kind=8),dimension(1:nlevelmax)::volume
  real(dp),dimension(1:3)::skip_loc,xcell
  real(dp),dimension(1:twotondim,1:3)::xc
  integer::nx_loc,ind
  integer ::ix,iy,iz

  !peak-patch related arrays before sharing information with other cpus
  real(kind=8),dimension(1:npeaks_tot)::max_dens,phi_min
  real(kind=8),dimension(1:npeaks_tot)::min_dens,clump_mass,clump_vol
  real(kind=8),dimension(1:npeaks_tot,1:3)::center_of_mass,clump_momentum,peak_pos
  integer,dimension(1:npeaks_tot)::n_cells

  min_dens=huge(0.d0);  max_dens=0.d0; n_cells=0; phi_min=huge(0.d0)
  clump_mass=0.d0; clump_vol=0.d0; clump_momentum=0.d0; center_of_mass=0.d0

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
        
        ! Cell density and energy- and momentum density
        d=uold(icellp(ipart),1)
        do i=1,ndim
           vd(i)=uold(icellp(ipart),i+1)
        end do

        ! Cell volume
        vol=volume(levp(ipart))

        ! number of leaf cells per clump
        n_cells(peak_nr)=n_cells(peak_nr)+1
        
        ! find min density and potential
        min_dens(peak_nr)=min(d,min_dens(peak_nr))
        phi_min(peak_nr)=min(phi(icellp(ipart)),phi_min(peak_nr))


        ! find max density and peak location
        if(d>=max_dens(peak_nr))then
           max_dens(peak_nr)=d
           peak_pos(peak_nr,1:ndim)=xcell(1:ndim)
        end if

        ! find clump mass
        clump_mass(peak_nr)=clump_mass(peak_nr)+vol*d
        
        ! center of mass velocity
        do i=1,3
           clump_momentum(peak_nr,i)=clump_momentum(peak_nr,i)+vd(i)*vol
        end do
        
        ! clump size (maybe take center of mass instead of peak as reference point)
        do i=1,ndim
           center_of_mass(peak_nr,i)=center_of_mass(peak_nr,i)+xcell(i)*vol*d
           
           ! compute second order moments
           do j=1,ndim
              second_moments(peak_nr,i,j)=second_moments(peak_nr,i,j)+xcell(i)*xcell(j)*vol*d
           end do
        end do
        
        ! clump volume
        clump_vol(peak_nr)=clump_vol(peak_nr)+vol
        
     end if
  end do

  !---------------------------------------------------------------------------
  ! a lot of MPI communication to collect the results from the different cpu's
  !---------------------------------------------------------------------------
#ifndef WITHOUTMPI     
  call MPI_ALLREDUCE(n_cells,n_cells_tot,npeaks_tot,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(min_dens,min_dens_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(phi_min,phi_min_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(max_dens,max_dens_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(clump_mass,clump_mass_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(clump_vol,clump_vol_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(clump_momentum,clump_momentum_tot,3*npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(center_of_mass,center_of_mass_tot,3*npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(second_moments,second_moments_tot,9*npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI     
  n_cells_tot=n_cells
  min_dens_tot=min_dens
  phi_min_tot=phi_min
  max_dens_tot=max_dens
  clump_mass_tot=clump_mass
  clump_vol_tot=clump_vol
  clump_momentum_tot=clump_momentum
  center_of_mass_tot=center_of_mass
  second_moments_tot=second_moments
#endif
  do j=1,npeaks_tot
     do i=1,ndim
        if (max_dens(j)<max_dens_tot(j))peak_pos(j,i)=0.
     end do
  end do
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(peak_pos,peak_pos_tot,3*npeaks_tot,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  peak_pos_tot=peak_pos
#endif
  !---------------------------------------------------------------------------
  ! compute further information related to the clumps
  !---------------------------------------------------------------------------
  
  !calculate total mass above threshold
  tot_mass=sum(clump_mass_tot)
  av_dens_tot(1:npeaks_tot)=clump_mass_tot(1:npeaks_tot)/clump_vol_tot(1:npeaks_tot)
  do i=1,ndim
     center_of_mass_tot(1:npeaks_tot,i)=center_of_mass_tot(1:npeaks_tot,i)/clump_mass_tot(1:npeaks_tot)
  end do

end subroutine compute_clump_properties
!################################################################
!################################################################
!################################################################
!################################################################
subroutine compute_clump_properties_round2(ntest,map,all_bound)
  use amr_commons
  use hydro_commons, ONLY:uold,gamma
  use poisson_commons, ONLY:phi
  use clfind_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  logical::map,all_bound
  integer::ntest

  !----------------------------------------------------------------------------
  ! this subroutine performs another loop over all particles and collects the 
  ! information more information like the velocity of a cell realtive to the 
  ! center of mass
  ! this routine is also used to write a peak map of the merged peak patches 
  ! to disk
  !----------------------------------------------------------------------------

  integer::ipart,ilevel,info,i,peak_nr,j

  !variables needed temporarily store cell properties
  real(dp)::d,vol,M,ekk,phi_rel,de,c_sound,d0,v_bulk2
  real(dp),dimension(1:3)::vd,xcell,xpeak,v_cl

  ! variables to be used with vector-sweeps
  integer::grid

  ! variables related to the size of a cell on a given level
  real(kind=8)::dx,dx_loc,scale,vol_loc
  real(kind=8),dimension(1:nlevelmax)::volume
  integer::nx_loc,ix,iy,iz,ind
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc

  !peak-patch related arrays before sharing information with other cpus
  real(kind=8),dimension(1:npeaks_tot)::e_kin_int,e_bind,e_thermal,e_kin_int4,e_bind4,e_thermal4,v_therm,v_rms,m4
  real(kind=8),dimension(1:npeaks_tot)::clump_mass4,E_kin_iso,E_bind_iso,E_therm_iso
  real(kind=8),dimension(1:npeaks_tot,1:3)::clump_size,bulk_momentum
  real(dp)::Ggrav=1.d0
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v

  !strings for file output           
  character(LEN=5)::myidstring,nchar 

  !  first, get minimum potential on saddle surface
  call get_phi_ref(ntest)
  call get_phi_ref2(ntest)
  
  !initialize arrays
  e_kin_int=0.d0; clump_size=0.d0; e_bind=0.d0; e_thermal=0.d0; e_bind4=0.d0; e_thermal4=0.d0; e_kin_int4=0.d0
  clump_mass4=0.d0
  v_therm=0.; v_rms=0.; bulk_momentum=0.; m4=0.
  E_kin_iso=0.; E_bind_iso=0.; E_therm_iso=0.
  
  ! Conversion factor from user units to cgs units                                             
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  d0 = density_threshold/scale_nH;
  if(cosmo)d0=d0/aexp**3

  !prepare file output for peak map
  if (map)then
     call title(ifout-1,nchar)
     call title(myid,myidstring)
     open(unit=20,file=TRIM('output_'//TRIM(nchar)//'/clump_map.csv'//myidstring),form='formatted')
  end if

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

        !peak_map
        if (map)write(20,'(F11.8,A,F11.8,A,F11.8,A,I8)')xcell(1),',',xcell(2),',',xcell(3),',',peak_nr
        
        ! Cell density and energy- and momentum density
        d=uold(icellp(ipart),1)
        de=uold(icellp(ipart),ndim+2)
        do i=1,ndim
           vd(i)=uold(icellp(ipart),i+1)
           xpeak(i)=peak_pos_tot(peak_nr,i)
        end do
        
        vol=volume(levp(ipart))                  

        phi_rel=(phi(icellp(ipart))-phi_min_tot(peak_nr))*scale

        M=clump_mass_tot(peak_nr)
        v_cl(1:ndim)=clump_momentum_tot(peak_nr,1:ndim)/M
       
        ! kinetic energy
        ekk=0.
        do i=1,3 
           ekk=ekk+0.5*vd(i)**2/d                          
        end do
        
        ! internal kinetic energy
        do i=1,3
           e_kin_int(peak_nr)=e_kin_int(peak_nr)+ &
                (vd(i)/d-clump_momentum_tot(peak_nr,i)/M)**2*d*vol*0.5
        end do
        
        ! potential energy using the acutal phi W= 0.5*int phi*rho
        e_bind(peak_nr)=e_bind(peak_nr)-phi_rel*d*vol*5.d-1
        
        ! size relative to center of mass
        do i=1,ndim
           clump_size(peak_nr,i)=clump_size(peak_nr,i)+(xcell(i)-center_of_mass_tot(peak_nr,i))**2.d0*vol
        end do
        
        ! thermal energy
        e_thermal(peak_nr)=e_thermal(peak_nr)+1.5*(de-ekk)*vol*(gamma-1)

        ! sound speed
        c_sound=(de-ekk)/d*gamma/(gamma-1)

        !Mass weighted thermal Velocity
        v_therm(peak_nr)=v_therm(peak_nr)+c_sound*d*vol/M        
                
        !MS Velocity
        do i=1,ndim
           v_rms(peak_nr)=v_rms(peak_nr)+(vd(i)/d-v_cl(i))**2*vol*d/M
        end do
        

        !for smaller region if cell is close enough (4 cells away)
        if (((xpeak(1)-xcell(1))**2.+(xpeak(2)-xcell(2))**2.+(xpeak(3)-xcell(3))**2.) .LE. 16.*volume(nlevelmax)**(2./3.))then
           do i=1,3
              e_kin_int4(peak_nr)=e_kin_int4(peak_nr)+(vd(i)/d-clump_momentum_tot(peak_nr,i)/M)**2*d*vol*0.5
              bulk_momentum(peak_nr,i)=bulk_momentum(peak_nr,i)+(vd(i)/d-v_cl(i))*vol*(d-d0)
           end do

           m4(peak_nr)=m4(peak_nr)+(d-d0)*vol
           
           e_bind4(peak_nr)=e_bind4(peak_nr)-phi_rel*d*vol*0.5
           e_thermal4(peak_nr)=e_thermal4(peak_nr)+1.5*(de-ekk)*vol*(gamma-1)

           clump_mass4(peak_nr)=clump_mass4(peak_nr)+d*vol
        end if

        !repeat for region enclosed by isopotential surface 
        if (phi_rel<0.)then
           do i=1,3
              E_kin_iso(peak_nr)=E_kin_iso(peak_nr)+(vd(i)/d-clump_momentum_tot(peak_nr,i)/M)**2*d*vol*0.5
           end do
           E_bind_iso(peak_nr)=E_bind_iso(peak_nr)-phi_rel*d*vol*0.5
           E_therm_iso(peak_nr)=E_therm_iso(peak_nr)+1.5*(de-ekk)*vol*(gamma-1)
        endif

     end if
  end do
  if (map)close(20)
  !---------------------------------------------------------------------------
  ! a lot of MPI communication to collect the results from the different cpu's
  !---------------------------------------------------------------------------
#ifndef WITHOUTMPI     
  call MPI_ALLREDUCE(e_kin_int,e_kin_int_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(e_kin_int4,e_kin_int_tot4,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(e_bind,e_bind_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(e_bind4,e_bind_tot4,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(clump_size,clump_size_tot,3*npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(m4,m4_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(bulk_momentum,bulk_momentum_tot,3*npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(e_thermal,e_thermal_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(e_thermal4,e_thermal_tot4,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(v_therm,v_therm_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(v_rms,v_rms_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(E_therm_iso,E_therm_iso_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(E_kin_iso,E_kin_iso_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(E_bind_iso,E_bind_iso_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI     
  e_kin_int_tot=e_kin_int
  e_kin_int_tot4=e_kin_int4
  e_bind_tot=e_bind
  e_bind_tot4=e_bind4
  clump_size_tot=clump_size
  m4_tot=m4
  bulk_momentum_tot=bulk_momentum
  e_thermal_tot=e_thermal
  e_thermal_tot4=e_thermal4
  v_therm_tot=v_therm
  v_rms_tot=v_rms
  E_bind_iso_tot=E_bind_iso
  E_kin_iso_tot=E_kin_iso
  E_therm_iso_tot=E_therm_iso
#endif
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(clump_mass4,clump_mass_tot4,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  clump_mass_tot4=clump_mass4
#endif

  all_bound=.true.
  do j=npeaks_tot,1,-1
     if (relevance_tot(j)>0.)then
        !compute all the checks
        v_bulk2=(bulk_momentum_tot(j,1)**2+bulk_momentum_tot(j,2)**2&
             +bulk_momentum_tot(j,3)**2)/(m4_tot(j)**2+tiny(0.d0))     
        peak_check(j)=scale*(phi_ref_tot(j)-phi_min_tot(j))/((v_therm_tot(j)**2+v_rms_tot(j)+v_bulk2)*0.5+tiny(0.d0))
        ball4_check(j)=scale*e_bind_tot4(j)/(tiny(0.d0)+2*e_thermal_tot4(j)+2*e_kin_int_tot4(j))
        isodens_check(j)=scale*E_bind_iso_tot(j)/(tiny(0.d0)+2*E_kin_iso_tot(j)+2*E_therm_iso_tot(j))
        clump_check(j)=(scale*e_bind_tot(j)+Psurf_tot(j))/(tiny(0.d0)+2*e_kin_int_tot(j)+2*e_thermal_tot(j))     
        all_bound=all_bound.and.(isodens_check(j)>1.)
     endif
  end do

end subroutine compute_clump_properties_round2
!################################################################
!################################################################
!################################################################
!################################################################
subroutine write_clump_properties(to_file)
  use amr_commons
  use clfind_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  logical::to_file

  !---------------------------------------------------------------------------
  ! this routine writes the clump properties to screen and to file
  !---------------------------------------------------------------------------

  integer::j,jj,ilun,n_rel,info,nx_loc
  real(kind=8)::rel_mass,scale
  character(LEN=5)::nchar



  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  
  !sort clumps by peak density in ascending order
  call heapsort_index(max_dens_tot,sort_index,npeaks_tot)

  if(to_file)then
     ilun=20
  else 
     ilun=6
  end if

  !print results in descending order to screen/file
  if(myid==1)then
     
     rel_mass=0.
     n_rel=0
     if (to_file .eqv. .true.) then
        call title(ifout-1,nchar)
        open(unit=20,file=TRIM('output_'//TRIM(nchar)//'/clump_info.txt'),form='formatted')
        open(unit=21,file=TRIM('output_'//TRIM(nchar)//'/clump_masses.txt'),form='formatted')
     end if
     if(smbh)then
     if(verbose)write(ilun,'(135A)')'Cl_N #leaf-cells  peak_x [uu] peak_y [uu] peak_z [uu] size_x [cm] size_y [cm] size_z [cc] |v|_CM [u.u.] rho- [H/cc] rho+ [H/cc] rho_av [H/cc] M_cl [M_sol] V_cl [AU^3] rel.  peak_check   ball4_check   isodens_check   clump_check '
     do j=npeaks_tot,1,-1
        jj=sort_index(j)
        if (relevance_tot(jj) > 0)then
           if(verbose)write(ilun,'(I6,X,I10,17(1X,1PE14.7))')jj&          
                ,n_cells_tot(jj)&
                ,peak_pos_tot(jj,1),peak_pos_tot(jj,2),peak_pos_tot(jj,3)&
                ,(5.*clump_size_tot(jj,1)/clump_vol_tot(jj))**0.5*scale_l &
                ,(5.*clump_size_tot(jj,2)/clump_vol_tot(jj))**0.5*scale_l &
                ,(5.*clump_size_tot(jj,3)/clump_vol_tot(jj))**0.5*scale_l &
                ,(clump_momentum_tot(jj,1)**2+clump_momentum_tot(jj,2)**2+ &
                clump_momentum_tot(jj,3)**2)**0.5/clump_mass_tot(jj)*scale_l/scale_t&
                ,min_dens_tot(jj)*scale_nH,max_dens_tot(jj)*scale_nH&
                ,clump_mass_tot(jj)/clump_vol_tot(jj)*scale_nH&
                ,clump_mass_tot(jj)*scale_d*dble(scale_l)**3/1.98892d33&
                ,clump_vol_tot(jj)*(scale_l)**3&
                ,relevance_tot(jj)&
                ,peak_check(jj)&
                ,ball4_check(jj)&
                ,isodens_check(jj)&
                ,clump_check(jj)



           rel_mass=rel_mass+clump_mass_tot(jj)*scale_d*dble(scale_l)**3/1.98892d33
           n_rel=n_rel+1
        end if
     end do

  else
     if(verbose)write(ilun,'(135A)')'Cl_N #leaf-cells  peak_x [uu] peak_y [uu] peak_z [uu] size_x [AU] size_y [AU]'//&
          ' size_z [AU]  |v|_CM [u.u.]  rho- [H/cc]  rho+ [H/cc]  rho_av [H/cc] M_cl [M_sol] V_cl [AU^3]   rel.  peak_check   ball4_check   isodens_check   clump_check phi_ref_tot phi_ref2_tot'
     do j=npeaks_tot,1,-1
        jj=sort_index(j)

        
        if (relevance_tot(jj) > 0)then
           if(verbose)write(ilun,'(I6,X,I10,3(X,F11.5),3(X,F11.5),X,F13.5,3(X,E12.3E2),2(X,E11.2E2),X,E11.2E2,6(2X,E11.2E2))')&
                jj&
                ,n_cells_tot(jj)&
                ,peak_pos_tot(jj,1)&
                ,peak_pos_tot(jj,2)&
                ,peak_pos_tot(jj,3)&
                ,(5.*clump_size_tot(jj,1)/clump_vol_tot(jj))**0.5*(scale_l/1.496d13)&
                ,(5.*clump_size_tot(jj,2)/clump_vol_tot(jj))**0.5*(scale_l/1.496d13)&
                ,(5.*clump_size_tot(jj,3)/clump_vol_tot(jj))**0.5*scale_l/1.496d13&
                ,(clump_momentum_tot(jj,1)**2+clump_momentum_tot(jj,2)**2+clump_momentum_tot(jj,3)**2)**0.5/clump_mass_tot(jj)&
                ,min_dens_tot(jj)*scale_nH&
                ,max_dens_tot(jj)*scale_nH&
                ,clump_mass_tot(jj)/clump_vol_tot(jj)*scale_nH&
                ,clump_mass_tot(jj)*scale_d*dble(scale_l)**3/1.98892d33&
                ,clump_vol_tot(jj)*(scale_l/1.496d13)**3&
                ,relevance_tot(jj)&
                ,peak_check(jj)&
                ,ball4_check(jj)&
                ,isodens_check(jj)&
                ,clump_check(jj)&
                ,phi_ref_tot(jj)&
                ,phi_ref2_tot(jj)
                
           rel_mass=rel_mass+clump_mass_tot(jj)*scale_d*dble(scale_l)**3/1.98892d33
           n_rel=n_rel+1
        end if
     end do
     end if
     if(to_file)then
        write(21,*)n_rel
        do j=npeaks_tot,1,-1
           jj=sort_index(j)
           if (relevance_tot(jj)>0)write(21,*)clump_mass_tot(jj)*scale_d*dble(scale_l)**3/1.98892d33
        end do
     else
        write(ilun,'(A,1PE12.5)')'total mass above threshold =',tot_mass*scale_d*dble(scale_l)**3/1.98892d33
        write(ilun,'(A,I6,A,1PE12.5)')'total mass in',n_rel,' listed clumps =',rel_mass
     endif
     if (to_file)then
        close(20)
        close(21)
     end if
  end if

end subroutine write_clump_properties
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine saddlepoint_search(ntest)
  use amr_commons
  use hydro_commons, ONLY: uold
  use clfind_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  integer::ntest

  !---------------------------------------------------------------------------
  ! subroutine which creates a npeaks_tot**2 sized array of saddlepoint densities
  ! by looping over all testparticles and passing them to neighborcheck
  ! with case 4, which means that saddlecheck will be called for each neighboring
  ! leaf cell. There it is checked, wether the two cells (original cell and 
  ! neighboring cell) are connected by a new densest saddle.
  !---------------------------------------------------------------------------

  integer::ipart,ip,ilevel,next_level,ntr
  integer::i,j,info,dummyint,icpu
  integer,dimension(1:nvector)::ind_cell
  integer::tag,last_imat_free_ext,ii,jj
  integer,allocatable,dimension(:)::tempint
  real(kind=8),allocatable,dimension(:)::temp,temp_tot
  integer,dimension(MPI_STATUS_SIZE)::stat

  ! saddle point array for 1 cpu
  allocate(saddle_dens(1:nmatmax))
  saddle_dens=0.
  last_imat_free=0
  last_imat_free_ext=0

  ! loop 'testparts', pass the information of nvector parts to neighborsearch 
  ip=0
  do ipart=1,ntest
     ip=ip+1
     ilevel=levp(testp_sort(ipart)) ! level
     next_level=0
     if(ipart<ntest)next_level=levp(testp_sort(ipart+1)) !level of next particle
     ind_cell(ip)=icellp(testp_sort(ipart))
     if (flag2(ind_cell(ip))==0)print*,'alert neighborsearch',ind_cell(ip),myid,ilevel
     if(ip==nvector .or. next_level /= ilevel)then
        call neighborsearch(ind_cell,ip,dummyint,ilevel,4)
        ip=0
     endif
  end do
  if (ip>0)call neighborsearch(ind_cell,ip,dummyint,ilevel,4)

!!! OLD VERSION
!  ! share the results among MPI domains  
!#ifndef WITHOUTMPI
!  ! OLDEST VERSION
!  ! call MPI_ALLREDUCE(saddle_dens,saddle_dens_tot,(npeaks_tot**2),MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
!  allocate(temp(1:npeaks_tot))
!  allocate(temp_tot(1:npeaks_tot))
!  do i=1,npeaks_tot
!     temp(1:npeaks_tot)=saddle_dens(1:npeaks_tot,i)
!     temp_tot=0.d0
!     call MPI_ALLREDUCE(temp,temp_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info) 
!     saddle_dens_tot(1:npeaks_tot,i)=temp_tot(1:npeaks_tot)
!  end do
!  deallocate(temp)
!  deallocate(temp_tot)
!#endif
!#ifdef WITHOUTMPI
!  saddle_dens_tot=saddle_dens
!#endif
!!!!!!!!!!!!!!!!!

  !!!if(last_imat_free>0)then
  !!!   do i=1,last_imat_free
  !!!      write(*,*)myid,icurrent(i),imat_prev(i),imat_next(i),saddle_dens(i)
  !!!   end do
  !!!end if

  ! fantastic tree loop to communicate the saddle point linked list
#ifndef WITHOUTMPI
  if(myid==1)write(*,*)'MPI communication - saddle density matrix'
  if(myid==1)call merge_lists(last_imat_free,last_imat_free)
  do icpu=2,ncpu
     tag=icpu
     if(myid==1)call MPI_RECV(last_imat_free_ext,1,MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,stat,info)
     if(myid==icpu)call MPI_SEND(last_imat_free,1,MPI_INTEGER,0,tag,MPI_COMM_WORLD,info)
     if(last_imat_free+last_imat_free_ext>nmatmax)write(*,*)'ERROR: Increase nmatmax!'
     allocate(tempint(1:nmatmax))
     tempint=0
     if(myid==1)call MPI_RECV(tempint,nmatmax,MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,stat,info)
     if(myid==icpu)call MPI_SEND(imat_next,nmatmax,MPI_INTEGER,0,tag,MPI_COMM_WORLD,info)
     if(myid==1.and.last_imat_free_ext>0)imat_next(last_imat_free+1:last_imat_free+last_imat_free_ext)=tempint(1:last_imat_free_ext)
     tempint=0
     if(myid==1)call MPI_RECV(tempint,nmatmax,MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,stat,info)
     if(myid==icpu)call MPI_SEND(imat_prev,nmatmax,MPI_INTEGER,0,tag,MPI_COMM_WORLD,info)
     if(myid==1.and.last_imat_free_ext>0)imat_prev(last_imat_free+1:last_imat_free+last_imat_free_ext)=tempint(1:last_imat_free_ext)
     tempint=0
     if(myid==1)call MPI_RECV(tempint,nmatmax,MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,stat,info)
     if(myid==icpu)call MPI_SEND(icurrent,nmatmax,MPI_INTEGER,0,tag,MPI_COMM_WORLD,info)
     if(myid==1.and.last_imat_free_ext>0)icurrent(last_imat_free+1:last_imat_free+last_imat_free_ext)=tempint(1:last_imat_free_ext)
     deallocate(tempint)
     allocate(temp(1:nmatmax))
     if(myid==1)call MPI_RECV(temp,nmatmax,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,stat,info)
     if(myid==icpu)call MPI_SEND(saddle_dens,nmatmax,MPI_DOUBLE_PRECISION,0,tag,MPI_COMM_WORLD,info)
     if(myid==1.and.last_imat_free_ext>0)saddle_dens(last_imat_free+1:last_imat_free+last_imat_free_ext)=temp(1:last_imat_free_ext)
     deallocate(temp)
     ! Merging of the two linked lists
     if(myid==1.and.last_imat_free_ext>0)call merge_lists(last_imat_free+last_imat_free_ext,last_imat_free)
  end do
  call MPI_BARRIER(MPI_COMM_WORLD,info)     
  ! Broadcast the final linked list to all cores
  call MPI_BCAST (last_imat_free,1,MPI_INTEGER,0,MPI_COMM_WORLD,info)
  call MPI_BCAST (icurrent,nmatmax,MPI_INTEGER,0,MPI_COMM_WORLD,info)
  call MPI_BCAST (imat_next,nmatmax,MPI_INTEGER,0,MPI_COMM_WORLD,info)
  call MPI_BCAST (imat_prev,nmatmax,MPI_INTEGER,0,MPI_COMM_WORLD,info)
  call MPI_BCAST (saddle_dens,nmatmax,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
  if(myid==1)write(*,*)'Saddle density matrix communication done'   
  saddle_dens_tot=saddle_dens   
#endif
#ifdef WITHOUTMPI
  saddle_dens_tot=saddle_dens
#endif

! VERSION 1 - To be tested
!#ifndef WITHOUTMPI   
!  if(myid==1)write(*,*)'MPI Communication'
!  if(ncpu/2*2==ncpu)then ! even number of cores
!     ntr=1 ! level of the tree
!     do while (int(ncpu)/int(2**ntr)>0) ! loop on the tree level 
!        do icpu=1,ncpu-1,2**ntr ! loop on the cores at the same level
!           tag=icpu
!           !if(myid==1.and.myid==icpu)write(*,*)'Core', icpu+2**(ntr-1), 'sends to Core', icpu
!           if(myid==icpu)call MPI_RECV(last_imat_free_ext,1,MPI_INTEGER,(icpu+2**(ntr-1))-1,tag,MPI_COMM_WORLD,stat,info)
!           if(myid==(icpu+2**(ntr-1)))call MPI_SEND(last_imat_free,1,MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,info)
!           if(last_imat_free+last_imat_free_ext>nmatmax)write(*,*)'ERROR: Increase nmatmax!'
!           allocate(tempint(1:nmatmax))
!           tempint=0
!           if(myid==icpu)call MPI_RECV(tempint,nmatmax,MPI_INTEGER,(icpu+2**(ntr-1))-1,tag,MPI_COMM_WORLD,stat,info)
!           if(myid==(icpu+2**(ntr-1)))call MPI_SEND(imat_next,nmatmax,MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,info)
!           imat_next(last_imat_free+1:last_imat_free+last_imat_free_ext)=tempint(1:last_imat_free_ext)
!           tempint=0
!           if(myid==icpu)call MPI_RECV(tempint,nmatmax,MPI_INTEGER,(icpu+2**(ntr-1))-1,tag,MPI_COMM_WORLD,stat,info)
!           if(myid==(icpu+2**(ntr-1)))call MPI_SEND(imat_prev,nmatmax,MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,info)
!           imat_prev(last_imat_free+1:last_imat_free+last_imat_free_ext)=tempint(1:last_imat_free_ext)
!           tempint=0
!           if(myid==icpu)call MPI_RECV(tempint,nmatmax,MPI_INTEGER,(icpu+2**(ntr-1))-1,tag,MPI_COMM_WORLD,stat,info)
!           if(myid==(icpu+2**(ntr-1)))call MPI_SEND(icurrent,nmatmax,MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,info)
!           icurrent(last_imat_free+1:last_imat_free+last_imat_free_ext)=tempint(1:last_imat_free_ext)
!           deallocate(tempint)
!           allocate(temp(1:nmatmax))
!           if(myid==icpu)call MPI_RECV(temp,nmatmax,MPI_DOUBLE_PRECISION,(icpu+2**(ntr-1))-1,tag,MPI_COMM_WORLD,stat,info)
!           if(myid==(icpu+2**(ntr-1)))call MPI_SEND(saddle_dens,nmatmax,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,info)
!           saddle_dens(last_imat_free+1:last_imat_free+last_imat_free_ext)=temp(1:last_imat_free_ext)
!           deallocate(temp)
!           ! Merging of the two linked lists
!           if(myid==icpu)call merge_lists(last_imat_free+last_imat_free_ext,last_imat_free)
!        end do
!        call MPI_BARRIER(MPI_COMM_WORLD,info)
!        ntr=ntr+1
!     end do
!     ! Broadcast the final linked list to all cores
!     call MPI_BCAST (last_imat_free,1,MPI_INTEGER,0,MPI_COMM_WORLD,info)
!     call MPI_BCAST (icurrent,nmatmax,MPI_INTEGER,0,MPI_COMM_WORLD,info)
!     call MPI_BCAST (imat_next,nmatmax,MPI_INTEGER,0,MPI_COMM_WORLD,info)
!     call MPI_BCAST (imat_prev,nmatmax,MPI_INTEGER,0,MPI_COMM_WORLD,info)
!     call MPI_BCAST (saddle_dens,nmatmax,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
!  else ! odd numer of cores
!     ntr=1 ! level of the tree
!     do while (int(ncpu-1)/int(2**ntr)>0) ! loop on the tree level  
!        do i=1,ncpu-2,2**ntr ! loop on the cores at the same level
!           ! point to point communication - the external linked list is appended to the local one 
!           tag=icpu
!           if(myid==icpu)call MPI_RECV(last_imat_free_ext,1,MPI_INTEGER,(icpu+2**(ntr-1))-1,tag,MPI_COMM_WORLD,stat,info)
!           if(myid==(icpu+2**(ntr-1)))call MPI_SEND(last_imat_free,1,MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,info)
!           if(last_imat_free+last_imat_free_ext>nmatmax)write(*,*)'ERROR: Increase nmatmax!'
!           allocate(tempint(1:nmatmax))
!           tempint=0
!           if(myid==icpu)call MPI_RECV(tempint,nmatmax,MPI_INTEGER,(icpu+2**(ntr-1))-1,tag,MPI_COMM_WORLD,stat,info)
!           if(myid==(icpu+2**(ntr-1)))call MPI_SEND(imat_next,nmatmax,MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,info)
!           imat_next(last_imat_free+1:last_imat_free+last_imat_free_ext)=tempint(1:last_imat_free_ext)
!           tempint=0
!           if(myid==icpu)call MPI_RECV(tempint,nmatmax,MPI_INTEGER,(icpu+2**(ntr-1))-1,tag,MPI_COMM_WORLD,stat,info)
!           if(myid==(icpu+2**(ntr-1)))call MPI_SEND(imat_prev,nmatmax,MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,info)
!           imat_prev(last_imat_free+1:last_imat_free+last_imat_free_ext)=tempint(1:last_imat_free_ext)
!           tempint=0
!           if(myid==icpu)call MPI_RECV(tempint,nmatmax,MPI_INTEGER,(icpu+2**(ntr-1))-1,tag,MPI_COMM_WORLD,stat,info)
!           if(myid==(icpu+2**(ntr-1)))call MPI_SEND(icurrent,nmatmax,MPI_INTEGER,icpu-1,tag,MPI_COMM_WORLD,info)
!           icurrent(last_imat_free+1:last_imat_free+last_imat_free_ext)=tempint(1:last_imat_free_ext)
!           deallocate(tempint)
!           allocate(temp(1:nmatmax))
!           if(myid==icpu)call MPI_RECV(temp,nmatmax,MPI_DOUBLE_PRECISION,(icpu+2**(ntr-1))-1,tag,MPI_COMM_WORLD,stat,info)
!           if(myid==(icpu+2**(ntr-1)))call MPI_SEND(saddle_dens,nmatmax,MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,info)
!           saddle_dens(last_imat_free+1:last_imat_free+last_imat_free_ext)=temp(1:last_imat_free_ext)
!           deallocate(temp)
!           ! Merging of the two linked lists 
!           if(myid==icpu)call merge_lists(last_imat_free+last_imat_free_ext,last_imat_free)
!        end do
!        call MPI_BARRIER(MPI_COMM_WORLD,info)
!        ntr=ntr+1
!     end do
!     ! final point to point communication with the last core
!     tag=1
!     if(myid==1)call MPI_RECV(last_imat_free_ext,1,MPI_INTEGER,ncpu-1,tag,MPI_COMM_WORLD,stat,info)
!     if(myid==ncpu)call MPI_SEND(last_imat_free,1,MPI_INTEGER,0,tag,MPI_COMM_WORLD,info)
!     if(last_imat_free+last_imat_free_ext>nmatmax)write(*,*)'ERROR: Increase nmatmax!'
!     allocate(tempint(1:nmatmax))
!     tempint=0
!     if(myid==1)call MPI_RECV(tempint,nmatmax,MPI_INTEGER,ncpu-1,tag,MPI_COMM_WORLD,stat,info)
!     if(myid==ncpu)call MPI_SEND(imat_next,nmatmax,MPI_INTEGER,0,tag,MPI_COMM_WORLD,info)
!     imat_next(last_imat_free+1:last_imat_free+last_imat_free_ext)=tempint(1:last_imat_free_ext)
!     tempint=0
!     if(myid==1)call MPI_RECV(tempint,nmatmax,MPI_INTEGER,ncpu-1,tag,MPI_COMM_WORLD,stat,info)
!     if(myid==ncpu)call MPI_SEND(imat_prev,nmatmax,MPI_INTEGER,0,tag,MPI_COMM_WORLD,info)
!     imat_prev(last_imat_free+1:last_imat_free+last_imat_free_ext)=tempint(1:last_imat_free_ext)
!     tempint=0
!     if(myid==1)call MPI_RECV(tempint,nmatmax,MPI_INTEGER,ncpu-1,tag,MPI_COMM_WORLD,stat,info)
!     if(myid==ncpu)call MPI_SEND(icurrent,nmatmax,MPI_INTEGER,0,tag,MPI_COMM_WORLD,info)
!     icurrent(last_imat_free+1:last_imat_free+last_imat_free_ext)=tempint(1:last_imat_free_ext)
!     deallocate(tempint)
!     allocate(temp(1:nmatmax))
!     if(myid==1)call MPI_RECV(temp,nmatmax,MPI_DOUBLE_PRECISION,ncpu-1,tag,MPI_COMM_WORLD,stat,info)
!     if(myid==ncpu)call MPI_SEND(saddle_dens,nmatmax,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
!     saddle_dens(last_imat_free+1:last_imat_free+last_imat_free_ext)=temp(1:last_imat_free_ext)
!     deallocate(temp)
!     ! Merging of the two linked lists 
!     if(myid==icpu)call merge_lists(last_imat_free+last_imat_free_ext,last_imat_free)
!     call MPI_BARRIER(MPI_COMM_WORLD,info)
!     ! Broadcast the final linked list to all cores                                                                                 
!     call MPI_BCAST (last_imat_free,1,MPI_INTEGER,0,MPI_COMM_WORLD,info)
!     call MPI_BCAST (icurrent,nmatmax,MPI_INTEGER,0,MPI_COMM_WORLD,info)
!     call MPI_BCAST (imat_next,nmatmax,MPI_INTEGER,0,MPI_COMM_WORLD,info)
!     call MPI_BCAST (imat_prev,nmatmax,MPI_INTEGER,0,MPI_COMM_WORLD,info)
!     call MPI_BCAST (saddle_dens,nmatmax,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
!  end if
!  if(myid==1)write(*,*)'Saddle density matrix communication done'
!  saddle_dens_tot=saddle_dens
!#endif
!#ifdef WITHOUTMPI    
!  saddle_dens_tot=saddle_dens
!#endif

  ! OLD VERSION  
  !! check symmetry (could be done only in case of verbose == true)
  !do i=1,npeaks_tot
  !   do j=1,i
  !      if(saddle_dens_tot(i,j)/=saddle_dens_tot(j,i).and.myid==1)then 
  !         write(*,*),'Alert! asymmetric saddle point array!',i,j,saddle_dens_tot(i,j),saddle_dens_tot(j,i)
  !      endif
  !   end do
  !end do

  ! compute saddle_max value and relevance
  ! OLD VERSION
  !saddle_max_tot=maxval(saddle_dens_tot,dim=1)

  saddle_max_tot=0.d0
  do i=1,last_imat_free
     jj=(icurrent(i)-1)/npeaks_tot+1
     ii=icurrent(i)-(jj-1)*npeaks_tot
     saddle_max_tot(ii)=max(saddle_dens_tot(i),saddle_max_tot(ii))
  end do

  do i=1,npeaks_tot
     if (saddle_max_tot(i)>0.)then
        relevance_tot(i)=max_dens_tot(i)/saddle_max_tot(i)
     else
        relevance_tot(i)=max_dens_tot(i)/min_dens_tot(i)
     end if
  end do
  
  !from here only saddle_dens_tot is used
  deallocate(saddle_dens)
  
end subroutine saddlepoint_search
!######################################################################### 
!######################################################################### 
!######################################################################### 
!######################################################################### 
subroutine merge_lists(last_imat_tot,last_imat_merged)
  use amr_commons
  use hydro_commons
  use clfind_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  integer::last_imat_tot,last_imat_merged
  !--------------------------------------------------------------------------- 
  ! Subroutine that merges two linked lists containing the non-zero elements
  ! of the saddlepoint matrix. The two linked lists are saved into a unique
  ! array whose elements will be sorted. Multiple entries in the final linked 
  ! list are finally removed.
  !---------------------------------------------------------------------------

  integer::i,j
  integer,allocatable,dimension(:)::ord,tempint
  real(kind=8),allocatable,dimension(:)::temp,ticurrent

  ! Sorting array by increasing icurrent
  allocate(ord(1:last_imat_tot))
  allocate(ticurrent(1:last_imat_tot))
  do i=1,last_imat_tot
     ord(i)=i
     ticurrent(i)=dble(icurrent(i))
  end do
  call quick_sort(ticurrent(1),ord(1),last_imat_tot)
  !call heapsort_index(ticurrent,ord,last_imat_tot)
  allocate(tempint(1:last_imat_tot))
  do i=1,last_imat_tot
     tempint(i)=int(ticurrent(i))
  end do
  icurrent(1:last_imat_tot)=tempint(1:last_imat_tot)
  deallocate(tempint,ticurrent)
  allocate(temp(1:last_imat_tot))
  do i=1,last_imat_tot
     temp(i)=saddle_dens(ord(i))
  end do
  saddle_dens(1:last_imat_tot)=temp(1:last_imat_tot)
  deallocate(temp)
  deallocate(ord)

  ! Build the non compact linked list 
  do i=1,last_imat_tot
     if (i==1)then
        imat_prev(i)=0
        imat_next(i)=i+1
     else
        if(icurrent(i)==icurrent(i-1))then
           saddle_dens(i-1)=max(saddle_dens(i),saddle_dens(i-1))
           saddle_dens(i)=saddle_dens(i-1)
           imat_next(i)=0
           imat_prev(i)=0
           imat_next(i-1)=i+1
           imat_prev(i+1)=i-1
           icurrent(i)=0
           if(i+1>last_imat_tot)then
              imat_next(i-1)=0
              imat_prev(i+1)=0
           end if
        else
           imat_prev(i)=i-1
           imat_next(i)=i+1
           imat_next(i-1)=i
           imat_prev(i+1)=i
           if(i+1>last_imat_tot)then
              imat_next(i)=0
           end if
        end if
     end if
  end do

  ! Make a compact form of the linked list
  allocate(tempint(1:last_imat_tot))
  tempint=0
  j=0
  do i=1,last_imat_tot
     if(icurrent(i).ne.0)then
        j=j+1
        tempint(j)=icurrent(i)
     end if
  end do
  icurrent=0
  icurrent(1:j)=tempint(1:j)
  deallocate(tempint)
  allocate(temp(1:last_imat_tot))
  temp=0
  j=0
  do i=1,last_imat_tot
     if(icurrent(i).ne.0)then
        j=j+1
        temp(j)=saddle_dens(i)
     end if
  end do
  saddle_dens=0.d0
  saddle_dens(1:j)=temp(1:j)
  deallocate(temp)

  imat_prev=0
  imat_next=0
  last_imat_merged=j
  do i=1,last_imat_merged
     if(i==1)then
        imat_prev(i)=0
        imat_next(i)=i+1
     else 
        if(i==last_imat_merged) then
           imat_prev(i)=i-1
           imat_next(i)=0
        else
           imat_prev(i)=i-1
           imat_next(i)=i+1
        end if
     end if
  end do

end subroutine merge_lists
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine merge_clumps(ntest)
  use amr_commons
  use clfind_commons
  implicit none
  integer::ntest
  !---------------------------------------------------------------------------
  ! This routine merges clumps according to different criteria
  !---------------------------------------------------------------------------
  integer::imat,jloc,iloc,j,i,ii,merge_count,final_peak,merge_to,ipart
  real(kind=8)::max_val
  integer,dimension(1:npeaks_tot)::old_peak
  integer::peak,next_peak,index1,index2
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,d0
  real(dp),allocatable,dimension(:)::temp_saddle

  if (verbose)write(*,*)'Now merging clumps'

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  d0 = density_threshold/scale_nH;
  if(cosmo)d0=d0/aexp**3

  ! Sort clumps by peak density in ascending order
  call heapsort_index(max_dens_tot,sort_index,npeaks_tot)

  if (.not.smbh) then
     do i=1,npeaks_tot
        ii=sort_index(i)
        new_peak(ii)=ii
        
        ! If the relevance is below the threshold -> merge
        if (relevance_tot(ii)<relevance_threshold.and.relevance_tot(ii)>=1.) then
           
           ! Go through the ii-th line in the saddle point array to find the neighbor to merge to
           ! OLD VERSION
           !merge_to=0; max_val=0.
           !do j=1,npeaks_tot
           !   if (saddle_dens_tot(ii,j)>max_val)then
           !      merge_to=j
           !      max_val=saddle_dens_tot(ii,j)
           !   end if
           !end do
           merge_to=0; max_val=0.
           do imat=1,last_imat_free
              jloc=(icurrent(imat)-1)/npeaks_tot+1
              iloc=icurrent(imat)-(jloc-1)*npeaks_tot
              if(saddle_dens_tot(imat)>max_val.and.ii==iloc)then
                 merge_to=jloc
                 max_val=saddle_dens_tot(imat)
              end if
           end do
           
           ! Store new peak index
           new_peak(ii)=merge_to
           if(verbose .and. myid==1)then
              if(merge_to>0)then
                 write(*,*)'clump ',ii,'merged to ',merge_to
              endif
           endif
           
           ! Update clump properties
           if (merge_to>0)then
              do j=1,ndim
                 clump_momentum_tot(merge_to,j)=&
                      clump_momentum_tot(merge_to,j)+clump_momentum_tot(ii,j)
                 clump_momentum_tot(ii,j)=0.
                 center_of_mass_tot(merge_to,j)=&
                      (clump_mass_tot(merge_to)*center_of_mass_tot(merge_to,j) &
                      +clump_mass_tot(ii)*center_of_mass_tot(ii,j)) &
                      /(clump_mass_tot(merge_to)+clump_mass_tot(ii))
                 center_of_mass_tot(ii,j)=0.
              end do
              n_cells_tot(merge_to)=n_cells_tot(merge_to)+n_cells_tot(ii)
              clump_vol_tot(merge_to)=clump_vol_tot(ii)+clump_vol_tot(merge_to)
              max_dens_tot(merge_to)=max(max_dens_tot(merge_to),max_dens_tot(ii))
              min_dens_tot(merge_to)=min(min_dens_tot(merge_to),min_dens_tot(ii))
              clump_mass_tot(merge_to)=clump_mass_tot(merge_to)+clump_mass_tot(ii)
           end if
           n_cells_tot(ii)=0
           clump_vol_tot(ii)=0.
           max_dens_tot(ii)=0.
           min_dens_tot(ii)=0.
           clump_mass_tot(ii)=0.
           
           ! Update saddle point array
           ! OLD VERSION
           !do j=1,npeaks_tot
           !   if (merge_to>0)then
           !      if(saddle_dens_tot(ii,j)>saddle_dens_tot(merge_to,j))then
           !         saddle_dens_tot(merge_to,j)=saddle_dens_tot(ii,j)
           !         saddle_dens_tot(j,merge_to)=saddle_dens_tot(ii,j)
           !      end if
           !      saddle_dens_tot(merge_to,merge_to)=0.
           !   end if
           !   saddle_dens_tot(ii,j)=0.        
           !   saddle_dens_tot(j,ii)=0.
           !end do
           allocate(temp_saddle(1:npeaks_tot))
           temp_saddle=0.d0
           do imat=1,last_imat_free
              if(merge_to>0.and.icurrent(imat)>0)then
                 jloc=(icurrent(imat)-1)/npeaks_tot+1
                 iloc=icurrent(imat)-(jloc-1)*npeaks_tot
                 if(iloc==ii)then
                    temp_saddle(jloc)=saddle_dens_tot(imat)
                    icurrent(imat)=merge_to+(jloc-1)*npeaks_tot
                 end if
                 if(jloc==ii)then
                    temp_saddle(iloc)=saddle_dens_tot(imat)
                    icurrent(imat)=iloc+(merge_to-1)*npeaks_tot
                 end if
              end if
           end do
           do imat=1,last_imat_free
              if (merge_to>0.and.icurrent(imat)>0)then
                 jloc=(icurrent(imat)-1)/npeaks_tot+1
                 iloc=icurrent(imat)-(jloc-1)*npeaks_tot
                 if(iloc==merge_to)then
                    temp_saddle(jloc)=max(temp_saddle(jloc),saddle_dens_tot(imat))
                 end if
              end if
           end do
           do imat=1,last_imat_free
              if(merge_to>0.and.icurrent(imat)>0)then
                 jloc=(icurrent(imat)-1)/npeaks_tot+1
                 iloc=icurrent(imat)-(jloc-1)*npeaks_tot
                 if(iloc==merge_to.or.jloc==merge_to.and.temp_saddle(jloc)>0)then
                    saddle_dens_tot(imat)=temp_saddle(jloc)
                 else
                    if(iloc==merge_to.or.jloc==merge_to.and.temp_saddle(jloc)==0.d0)then
                       icurrent(imat)=npeaks_tot**2
                       saddle_dens_tot(imat)=0.d0
                    end if
                 end if
              end if
           end do
           deallocate(temp_saddle)
           do imat=1,last_imat_free
              if(icurrent(imat)>0)then
                 jloc=(icurrent(imat)-1)/npeaks_tot+1
                 iloc=icurrent(imat)-(jloc-1)*npeaks_tot
                 if(iloc==merge_to.and.jloc==merge_to)saddle_dens_tot(imat)=0.
                 if(iloc==ii.or.jloc==ii)saddle_dens_tot(imat)=0.
              end if
           end do
           
           ! Update saddle_max value
           ! OLD VERSION
           !if (merge_to>0)then
           !   saddle_max_tot(merge_to)=0
           !   do j=1,npeaks_tot
           !      if (saddle_dens_tot(merge_to,j)>saddle_max_tot(merge_to))then
           !         saddle_max_tot(merge_to)=saddle_dens_tot(merge_to,j)
           !      end if
           !   end do
           !end if
           if (merge_to>0)then
              saddle_max_tot(merge_to)=0
              do imat=1,last_imat_free
                 jloc=(icurrent(imat)-1)/npeaks_tot+1
                 iloc=icurrent(imat)-(jloc-1)*npeaks_tot
                 if (iloc==merge_to.and.saddle_dens_tot(imat)>saddle_max_tot(merge_to))then
                    saddle_max_tot(merge_to)=saddle_dens_tot(imat)
                 end if
              end do
           end if
           
           ! Update relevance of clumps
           if (merge_to>0)then
              if (saddle_max_tot(merge_to)>1.d-40)then
                 relevance_tot(merge_to)=max_dens_tot(merge_to)/saddle_max_tot(merge_to)
              else 
                 relevance_tot(merge_to)=max_dens_tot(merge_to)/min_dens_tot(merge_to)
              end if
           end if
           relevance_tot(ii)=0.
        end if
     end do
  else
     ! If SMBH merge all peaks
     do i=1,npeaks_tot
        ii=sort_index(i)
        new_peak(ii)=ii

        ! Go through the ii-th line in the saddle point array to find the neighbor to merge to
        ! OLD VERSION
        !merge_to=0; max_val=0.
        !do j=1,npeaks_tot
        !   if (saddle_dens_tot(ii,j)>max_val)then
        !      merge_to=j
        !      max_val=saddle_dens_tot(ii,j)
        !   end if
        !end do
        merge_to=0; max_val=0.
        do imat=1,last_imat_free
           if(icurrent(imat)>0)then
              jloc=(icurrent(imat)-1)/npeaks_tot+1
              iloc=icurrent(imat)-(jloc-1)*npeaks_tot
              if(saddle_dens_tot(imat)>max_val.and.ii==iloc)then
                 merge_to=jloc
                 max_val=saddle_dens_tot(imat)
              end if
           end if
        end do

        ! Store new peak index
        if(merge_to>0)new_peak(ii)=merge_to
        if(verbose .and. myid==1)then
           if(merge_to>0)then
              write(*,*)'clump ',ii,'merged to ',merge_to
           endif
        endif

        ! Update clump properties
        if (merge_to>0)then
           do j=1,ndim
              clump_momentum_tot(merge_to,j)=&
                   clump_momentum_tot(merge_to,j)+clump_momentum_tot(ii,j)
              clump_momentum_tot(ii,j)=0.
              center_of_mass_tot(merge_to,j)=&
                   (clump_mass_tot(merge_to)*center_of_mass_tot(merge_to,j) &
                   +clump_mass_tot(ii)*center_of_mass_tot(ii,j)) &
                   /(clump_mass_tot(merge_to)+clump_mass_tot(ii))
              center_of_mass_tot(ii,j)=0.
           end do
           n_cells_tot(merge_to)=n_cells_tot(merge_to)+n_cells_tot(ii)
           clump_vol_tot(merge_to)=clump_vol_tot(ii)+clump_vol_tot(merge_to)
           max_dens_tot(merge_to)=max(max_dens_tot(merge_to),max_dens_tot(ii))
           min_dens_tot(merge_to)=min(min_dens_tot(merge_to),min_dens_tot(ii))
           clump_mass_tot(merge_to)=clump_mass_tot(merge_to)+clump_mass_tot(ii)
           n_cells_tot(ii)=0
           clump_vol_tot(ii)=0.
           max_dens_tot(ii)=0.
           min_dens_tot(ii)=0.
           clump_mass_tot(ii)=0.
        end if

        ! Update saddle point array
        ! OLD VERSION
        !do j=1,npeaks_tot
        !   if (merge_to>0)then
        !      if(saddle_dens_tot(ii,j)>saddle_dens_tot(merge_to,j))then
        !         saddle_dens_tot(merge_to,j)=saddle_dens_tot(ii,j)
        !         saddle_dens_tot(j,merge_to)=saddle_dens_tot(ii,j)
        !      end if
        !      saddle_dens_tot(merge_to,merge_to)=0.
        !      saddle_dens_tot(ii,j)=0.
        !      saddle_dens_tot(j,ii)=0.
        !   end if
        !end do
        allocate(temp_saddle(1:npeaks_tot))
        temp_saddle=0.d0
        do imat=1,last_imat_free
           if(merge_to>0.and.icurrent(imat)>0)then
              jloc=(icurrent(imat)-1)/npeaks_tot+1
              iloc=icurrent(imat)-(jloc-1)*npeaks_tot
              if(iloc==ii)then
                 temp_saddle(jloc)=saddle_dens_tot(imat)
                 icurrent(imat)=merge_to+(jloc-1)*npeaks_tot
              end if
              if(jloc==ii)then
                 temp_saddle(iloc)=saddle_dens_tot(imat)
                 icurrent(imat)=iloc+(merge_to-1)*npeaks_tot
              end if
           end if
        end do
        do imat=1,last_imat_free
           if (merge_to>0.and.icurrent(imat)>0)then
              jloc=(icurrent(imat)-1)/npeaks_tot+1
              iloc=icurrent(imat)-(jloc-1)*npeaks_tot
              if(iloc==merge_to)then
                 temp_saddle(jloc)=max(temp_saddle(jloc),saddle_dens_tot(imat))
              end if
           end if
        end do
        do imat=1,last_imat_free
           if(merge_to>0.and.icurrent(imat)>0)then
              jloc=(icurrent(imat)-1)/npeaks_tot+1
              iloc=icurrent(imat)-(jloc-1)*npeaks_tot
              if(iloc==merge_to.or.jloc==merge_to.and.temp_saddle(jloc)>0)then
                 saddle_dens_tot(imat)=temp_saddle(jloc)
              else 
                 if(iloc==merge_to.or.jloc==merge_to.and.temp_saddle(jloc)==0.d0)then
                    icurrent(imat)=npeaks_tot**2
                    saddle_dens_tot(imat)=0.d0
                 end if
              end if
           end if
        end do
        deallocate(temp_saddle)
        do imat=1,last_imat_free
           if(merge_to>0.and.icurrent(imat)>0)then
              jloc=(icurrent(imat)-1)/npeaks_tot+1
              iloc=icurrent(imat)-(jloc-1)*npeaks_tot
              if(iloc==merge_to.and.jloc==merge_to)saddle_dens_tot(imat)=0.
              if(iloc==ii.or.jloc==ii)saddle_dens_tot(imat)=0.              
           end if
        end do

        ! Update saddle_max value
        ! OLD VERSION
        !if (merge_to>0)then
        !   saddle_max_tot(merge_to)=0
        !   do j=1,npeaks_tot
        !      if (saddle_dens_tot(merge_to,j)>saddle_max_tot(merge_to))then
        !         saddle_max_tot(merge_to)=saddle_dens_tot(merge_to,j)
        !      end if
        !   end do
        !end if
        if (merge_to>0)then
           saddle_max_tot(merge_to)=0
           do imat=1,last_imat_free
              if(icurrent(imat)>0)then
                 jloc=(icurrent(imat)-1)/npeaks_tot+1
                 iloc=icurrent(imat)-(jloc-1)*npeaks_tot
                 if (iloc==merge_to.and.saddle_dens_tot(imat)>saddle_max_tot(merge_to))then
                    saddle_max_tot(merge_to)=saddle_dens_tot(imat)
                 end if
              end if
           end do
        end if
        
        ! Update relevance of clumps
        if (merge_to>0)then
           if (saddle_max_tot(merge_to)>1.d-40)then
              relevance_tot(merge_to)=max_dens_tot(merge_to)/saddle_max_tot(merge_to)
           else
              relevance_tot(merge_to)=max_dens_tot(merge_to)/min_dens_tot(merge_to)
           end if
           relevance_tot(ii)=0.
        end if
     end do
  end if

  if (verbose)write(*,*)'Done merging clumps 0'

! Change new_peak so that it points to the end point of the merging 
! history and not only to the clump it has been merged to in first place 
  do i=1,npeaks_tot
     peak=i
     merge_count=0
     next_peak=new_peak(peak)
     do while(peak.NE.next_peak.AND.next_peak>0)
        !construct old_peak to walk the tree from trunk to leafs
        old_peak(next_peak)=peak
        !go to next peak
        merge_count=merge_count+1
        peak=next_peak
        next_peak=new_peak(peak)
     end do
     final_peak=next_peak
     !now we walk the other way, so next_peak will always be the previous one
     next_peak=peak
     do j=merge_count,1,-1
        next_peak=old_peak(next_peak)
        new_peak(next_peak)=final_peak
     end do
  end do

  if (verbose)write(*,*)'Done merging clumps I'

! Remove peaks using HOP-based criterion
  if (smbh)then
     do i=1,npeaks_tot
        ii=sort_index(i)
        if( max_dens_tot(new_peak(ii)) < 3.0*d0 )then
           new_peak(ii)=0
           n_cells_tot(ii)=0
           clump_vol_tot(ii)=0.
           max_dens_tot(ii)=0.
           min_dens_tot(ii)=0.
           clump_mass_tot(ii)=0.
           ! Update saddle point array
           ! OLD VERSION
           !do j=1,npeaks_tot
           !   saddle_dens_tot(ii,j)=0.
           !   saddle_dens_tot(j,ii)=0.
           !end do
           do imat=1,last_imat_free
              jloc=(icurrent(imat)-1)/npeaks_tot+1
              iloc=icurrent(imat)-(jloc-1)*npeaks_tot
              if(iloc==ii.or.jloc==ii)saddle_dens_tot(imat)=0.
           end do
           relevance_tot(ii)=0.        
        end if
     end do
  endif

  !update flag 2
  do ipart=1,ntest
     if (flag2(icellp(ipart))>0)flag2(icellp(ipart))=new_peak(flag2(icellp(ipart)))
  end do

  if (verbose)write(*,*)'Done merging clumps II'  
end subroutine merge_clumps
!################################################################                 
!################################################################ 
!################################################################                 
!################################################################     
!subroutine update_saddle(iold,merge_to,saddle_new)
!  use amr_commons
!  use clfind_commons
!  implicit none
!
!  integer::iold,merge_to,iloc,jloc,imat
!  real(kind=8),dimension(1:npeaks_tot)::saddle_new
  !---------------------------------------------------------------------------  
  ! This routine updates the saddle point
  !--------------------------------------------------------------------------- 
  
  
!end subroutine update_saddle
!################################################################
!################################################################
!################################################################
!################################################################
subroutine allocate_peak_patch_arrays
  use amr_commons, ONLY:ndim
  use clfind_commons
  implicit none

  ! Allocate peak-patch_properties
  allocate(n_cells_tot(1:npeaks_tot))
  allocate(clump_size_tot(1:npeaks_tot,1:ndim))
  allocate(peak_pos_tot(1:npeaks_tot,1:ndim))
  allocate(center_of_mass_tot(1:npeaks_tot,1:ndim))
  allocate(second_moments(1:npeaks_tot,1:ndim,1:ndim)) 
  allocate(second_moments_tot(1:npeaks_tot,1:ndim,1:ndim))
  allocate(min_dens_tot(1:npeaks_tot))
  allocate(av_dens_tot(1:npeaks_tot))
  allocate(max_dens_tot(1:npeaks_tot))
  allocate(clump_mass_tot(1:npeaks_tot))
  allocate(clump_mass_tot4(1:npeaks_tot))
  allocate(clump_vol_tot(1:npeaks_tot))
  allocate(saddle_max_tot(1:npeaks_tot))
  allocate(relevance_tot(1:npeaks_tot))
  allocate(sort_index(1:npeaks_tot))
  allocate(saddle_dens_tot(1:nmatmax))
  allocate(icurrent(1:nmatmax))
  allocate(imat_next(1:nmatmax))
  allocate(imat_prev(1:nmatmax))
  allocate(clump_momentum_tot(1:npeaks_tot,1:ndim))
  allocate(e_kin_int_tot(npeaks_tot))
  allocate(e_kin_int_tot4(npeaks_tot))!all the "4" variables are for the 4cell ball properties...
  allocate(e_bind_tot(npeaks_tot))
  allocate(e_bind_tot4(npeaks_tot))
  allocate(e_thermal_tot(npeaks_tot))
  allocate(e_thermal_tot4(npeaks_tot)) 
  allocate(phi_min_tot(npeaks_tot))
  allocate(minmatch_tot(npeaks_tot))
  allocate(new_peak(npeaks_tot))
  allocate(phi_ref(npeaks_tot))
  allocate(phi_ref_tot(npeaks_tot))
  allocate(Psurf(npeaks_tot))
  allocate(Psurf_tot(npeaks_tot))
  allocate(v_therm_tot(npeaks_tot))
  allocate(v_rms_tot(npeaks_tot))
  allocate(m4_tot(npeaks_tot))
  allocate(bulk_momentum_tot(1:npeaks_tot,1:ndim))
  allocate(E_kin_iso_tot(npeaks_tot))
  allocate(E_bind_iso_tot(npeaks_tot))
  allocate(E_therm_iso_tot(npeaks_tot))
  allocate(peak_check(npeaks_tot))
  allocate(ball4_check(npeaks_tot))
  allocate(isodens_check(npeaks_tot))
  allocate(clump_check(npeaks_tot))
  allocate(phi_ref2_tot(npeaks_tot))


  !initialize all peak based arrays
  n_cells_tot=0
  saddle_max_tot=0.
  relevance_tot=1.
  clump_size_tot=0.
  min_dens_tot=huge(0.d0)
  max_dens_tot=0.
  av_dens_tot=0.
  clump_mass_tot=0.
  clump_mass_tot4=0.
  clump_vol_tot=0.
  peak_pos_tot=0.
  center_of_mass_tot=0.
  second_moments=0.; second_moments_tot=0.
  saddle_dens_tot=0.
  icurrent=0
  imat_next=0
  imat_prev=0
  clump_momentum_tot=0.
  e_kin_int_tot=0.; e_kin_int_tot4=0.
  e_bind_tot=0.; e_bind_tot4=0.
  e_thermal_tot=0.; e_thermal_tot4=0.
  phi_min_tot=0.
  minmatch_tot=1
  new_peak=0
  phi_ref=huge(0.d0)
  Psurf=0.
  

end subroutine allocate_peak_patch_arrays
!################################################################                 
!################################################################ 
!################################################################                 
!################################################################     
subroutine deallocate_all
  use clfind_commons
  implicit none

  deallocate(n_cells_tot)
  deallocate(clump_size_tot)
  deallocate(peak_pos_tot)
  deallocate(center_of_mass_tot)
  deallocate(second_moments)
  deallocate(second_moments_tot)
  deallocate(min_dens_tot)
  deallocate(av_dens_tot)
  deallocate(max_dens_tot)
  deallocate(clump_mass_tot)
  deallocate(clump_vol_tot)
  deallocate(saddle_max_tot)
  deallocate(relevance_tot)
  deallocate(sort_index)
  deallocate(saddle_dens_tot)
  deallocate(icurrent)
  deallocate(imat_prev)
  deallocate(imat_next)
  deallocate(clump_momentum_tot)
  deallocate(e_kin_int_tot,e_kin_int_tot4)
  deallocate(e_bind_tot,e_bind_tot4)
  deallocate(e_thermal_tot,e_thermal_tot4)
  deallocate(phi_min_tot)
  deallocate(minmatch_tot)
  deallocate(new_peak)
  deallocate(phi_ref,phi_ref_tot,phi_ref2_tot)
  deallocate(Psurf,Psurf_tot)
  deallocate(v_therm_tot)
  deallocate(v_rms_tot)
  deallocate(m4_tot,bulk_momentum_tot)
  deallocate(E_kin_iso_tot,E_bind_iso_tot,E_therm_iso_tot)
  deallocate(peak_check,ball4_check,isodens_check,clump_check)

end subroutine deallocate_all
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
  ! surface of the peak patch
  !---------------------------------------------------------------
  integer::info
   
  integer::ipart,ip,ilevel,next_level
  integer,dimension(1:nvector)::ind_cell



  ! loop 'testparts', pass the information of nvector parts to neighborsearch 
  ip=0
  do ipart=1,ntest
     ip=ip+1
     ilevel=levp(testp_sort(ipart)) ! level
     next_level=0
     if(ipart<ntest)next_level=levp(testp_sort(ipart+1)) !level of next particle
     ind_cell(ip)=icellp(testp_sort(ipart))
     if(ip==nvector .or. next_level /= ilevel)then
        call neighborsearch(ind_cell,ip,0,ilevel,5)
        call surface_int(ind_cell,ip,ilevel)
        ip=0
     endif
  end do
  if (ip>0)then 
     call neighborsearch(ind_cell,ip,0,ilevel,5)
     call surface_int(ind_cell,ip,ilevel)
  endif
   
#ifndef WITHOUTMPI     
  call MPI_ALLREDUCE(phi_ref,phi_ref_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(Psurf,Psurf_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI     
  phi_ref_tot=phi_ref
  Psurf_tot=Psurf
#endif

end subroutine get_phi_ref









subroutine get_phi_ref2
  use amr_commons
  use hydro_commons
  use pm_commons
  use clfind_commons
  use poisson_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !---------------------------------------------------------------
  ! This subroutine checks wheter a density maximum corresponds to
  ! a potential minimum -> pot_min(peak_nr)=ok in this case.
  ! Furthermore, the minimum potential is stored for every peak.
  !---------------------------------------------------------------
  integer::k1,j1,i1,jj,info,N,i,itest
  integer,dimension(1:nvector)::cell_index,cell_levl,ind_cell,lev_cell,cc
  real(dp),dimension(1:nvector,1:3)::pos,xtest,xtest_cpu,xtest_ind
  real(dp)::dx,dx_loc,scale,phim,x,y,z,r
  integer::nx_loc

  real(dp),dimension(1:npeaks_tot)::phi_ref2
  integer,dimension(1:npeaks_tot)::minmatch

  phi_ref2=0.d0; minmatch=1

  dx=0.5D0**nlevelmax
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale/aexp  

  N=10
  do jj=1,npeaks_tot
     if (relevance_tot(jj)>0)then
        pos(1,1:3)=peak_pos_tot(jj,1:3)
        
        itest=0
        !construct sample points on 4 cell radius ball surface                                                                                                                                                                         
        do k1=-N,N
           do j1=-N,N
              do i1=-N,N
                 x=dble(k1)/N*4.
                 y=dble(j1)/N*4.
                 z=dble(i1)/N*4.
                 r=(x**2.+y**2.+z**2.)**0.5
                 !reject all points with radius > 4dx
                 if (r<=4. .and. r>0)then
                    itest=itest+1
                    ! project the other points to the sphere surface
                    x=x/r*4.
                    y=y/r*4.
                    z=z/r*4.
                    xtest_cpu(itest,1)=pos(1,1)+x*dx_loc
                    xtest_cpu(itest,2)=pos(1,2)+y*dx_loc
                    xtest_cpu(itest,3)=pos(1,3)+z*dx_loc
                    if (itest==nvector)then
                       call cmp_cpumap(xtest_cpu,cc,itest)
                       do i=1,nvector
                          if(cc(i) == myid) then
                             xtest_ind(1,1:3)=xtest_cpu(i,1:3)
                             call get_cell_index(cell_index,cell_levl,xtest_ind,nlevelmax,1)
                             phi_ref2(jj)=min(phi_ref2(jj),phi(cell_index(1)))
                          endif
                       end do
                       itest=0
                    endif
                 endif
              end do
           end do
        end do

        call cmp_cpumap(xtest_cpu,cc,itest)
        do i=1,itest
           if(cc(i) == myid) then
              xtest_ind(1,1:3)=xtest_cpu(i,1:3)
              call get_cell_index(cell_index,cell_levl,xtest_ind,nlevelmax,1)
              phi_ref2(jj)=min(phi_ref2(jj),phi(cell_index(1)))
           endif
        end do

     end if ! end if relevance > than..
  end do ! end loop over clumps


#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(phi_ref2,phi_ref2_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  phi_ref2_tot=phi_ref2
#endif
end subroutine get_phi_ref2



