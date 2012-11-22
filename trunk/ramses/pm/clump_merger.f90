subroutine compute_clump_properties(ntest)
  use amr_commons
  use hydro_commons, ONLY:uold
  use clfind_commons
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
  real(kind=8),dimension(1:npeaks_tot)::max_dens
  real(kind=8),dimension(1:npeaks_tot)::min_dens,clump_mass,clump_vol
  real(kind=8),dimension(1:npeaks_tot,1:3)::center_of_mass,clump_momentum,peak_pos
  integer,dimension(1:npeaks_tot)::n_cells

  min_dens=huge(0.d0);  max_dens=0.d0; n_cells=0
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
        
        ! find min density
        if(d<=min_dens(peak_nr))min_dens(peak_nr)=d
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
#endif
#ifdef WITHOUTMPI     
  n_cells_tot=n_cells
#endif
#ifndef WITHOUTMPI     
  call MPI_ALLREDUCE(min_dens,min_dens_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI     
  min_dens_tot=min_dens
#endif
#ifndef WITHOUTMPI     
  call MPI_ALLREDUCE(max_dens,max_dens_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI     
  max_dens_tot=max_dens
#endif
#ifndef WITHOUTMPI     
  call MPI_ALLREDUCE(clump_mass,clump_mass_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI     
  clump_mass_tot=clump_mass
#endif
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(clump_vol,clump_vol_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  clump_vol_tot=clump_vol
#endif
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(clump_momentum,clump_momentum_tot,3*npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  clump_momentum_tot=clump_momentum
#endif
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(center_of_mass,center_of_mass_tot,3*npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  center_of_mass_tot=center_of_mass
#endif
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(second_moments,second_moments_tot,9*npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
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
subroutine compute_clump_properties_round2(ntest,map)
  use amr_commons
  use hydro_commons, ONLY:uold,gamma
  use poisson_commons, ONLY:phi
  use clfind_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  logical::map
  integer::ntest

  !----------------------------------------------------------------------------                
  ! this subroutine performs another loop over all particles and collects the 
  ! information more information like the velocity of a cell realtive to the 
  ! center of mass
  ! this routine is also used to write a peak map of the merged peak patches 
  ! to disk
  !----------------------------------------------------------------------------

  integer::ipart,ilevel,info,i,peak_nr

  !variables needed temporarily store cell properties
  real(dp)::d,vol,M,ekk,phi_rel,de
  real(dp),dimension(1:3)::vd,xcell,xpeak

  ! variables to be used with vector-sweeps
  integer::grid

  ! variables related to the size of a cell on a given level
  real(kind=8)::dx,dx_loc,scale,vol_loc
  real(kind=8),dimension(1:nlevelmax)::volume
  integer::nx_loc,ix,iy,iz,ind
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc

  !peak-patch related arrays before sharing information with other cpus
  real(kind=8),dimension(1:npeaks_tot)::e_kin_int,e_bind,e_thermal,e_kin_int4,e_bind4,e_thermal4
  real(kind=8),dimension(1:npeaks_tot,1:3)::clump_size

  !strings for file output
  character(LEN=5)::myidstring,nchar
  

  !  first, check wether clump density max and potential minimum match
  call clump_phi
  
  !initialize arrays
  e_kin_int=0.d0; clump_size=0.d0; e_bind=0.d0; e_thermal=0.d0; e_bind4=0.d0; e_thermal4=0.d0; e_kin_int4=0.d0
  
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
        phi_rel=(phi(icellp(ipart))-phi_min_tot(peak_nr))
        M=clump_mass_tot(peak_nr)
       
        ! Radius and kinetic energy
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
!        e_bind(peak_nr)=e_bind(peak_nr)-phi_rel*d*vol*5.d-1
        ! add gravitational self energy of every cell
!        e_bind(peak_nr)=e_bind(peak_nr)+9.d-1*d**2.d0*vol**(5.d0/3.d0)
        
        
        ! size relative to center of mass 
        do i=1,ndim
           clump_size(peak_nr,i)=clump_size(peak_nr,i)+(xcell(i)-center_of_mass_tot(peak_nr,i))**2.d0*vol
        end do
        
        ! thermal energy
        e_thermal(peak_nr)=e_thermal(peak_nr)+1.5*(de-ekk)*vol*(gamma-1)
        
        ! repeat same for smaller region if cell is close enough (4 cells away)
        if (((xpeak(1)-xcell(1))**2.+(xpeak(2)-xcell(2))**2.+(xpeak(3)-xcell(3))**2.) .LE. 16.*volume(nlevelmax)**(2./3.))then
           do i=1,3
              e_kin_int4(peak_nr)=e_kin_int4(peak_nr)+(vd(i)/d-clump_momentum_tot(peak_nr,i)/M)**2*d*vol*0.5
           end do
           
           e_bind4(peak_nr)=e_bind4(peak_nr)-phi_rel*d*vol*0.5
 !          e_bind4(peak_nr)=e_bind4(peak_nr)+0.9*d**2.*vol**(5./3.)
           e_thermal4(peak_nr)=e_thermal4(peak_nr)+1.5*(de-ekk)*vol*(gamma-1)
        end if
     end if
  end do
  if (map)close(20)
  !---------------------------------------------------------------------------
  ! a lot of MPI communication to collect the results from the different cpu's
  !---------------------------------------------------------------------------
#ifndef WITHOUTMPI     
  call MPI_ALLREDUCE(e_kin_int,e_kin_int_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(e_kin_int4,e_kin_int_tot4,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI     
  e_kin_int_tot=e_kin_int
  e_kin_int_tot4=e_kin_int4
#endif
#ifndef WITHOUTMPI     
  call MPI_ALLREDUCE(e_bind,e_bind_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(e_bind4,e_bind_tot4,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI     
  e_bind_tot=e_bind
  e_bind_tot4=e_bind4
#endif
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(clump_size,clump_size_tot,3*npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  clump_size_tot=clump_size
#endif
#ifndef WITHOUTMPI     
  call MPI_ALLREDUCE(e_thermal,e_thermal_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(e_thermal4,e_thermal_tot4,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI     
  e_thermal_tot=e_thermal
  e_thermal_tot4=e_thermal4
#endif

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

  integer::j,jj,ilun,n_rel
  real(kind=8)::rel_mass

  character(LEN=5)::nchar

  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  !sort clumps by peak density in ascending order
  call heapsort_index(max_dens_tot,sort_index,npeaks_tot)

  if(to_file)then
     ilun=20
  else 
     ilun=6
  end if

  !print results in descending order to screen/file
  if(myid==1) then 
     rel_mass=0.
     n_rel=0
     if (to_file .eqv. .true.) then
        call title(ifout-1,nchar)
        open(unit=20,file=TRIM('output_'//TRIM(nchar)//'/clump_info.txt'),form='formatted')
        open(unit=21,file=TRIM('output_'//TRIM(nchar)//'/clump_masses.txt'),form='formatted')
     end if
     if(smbh)then
     write(ilun,'(135A)')'Cl_N #leaf-cells  peak_x [uu] peak_y [uu] peak_z [uu] size_x [cm] size_y [cm] size_z [AU] |v|_CM [u.u.] rho- [H/cc] rho+ [H/cc] rho_av [H/cc] M_cl [M_sol] V_cl [AU^3] rel. V/(U+Q) V4/(U4+Q4) m_match'
     do j=npeaks_tot,1,-1
        jj=sort_index(j)
        if (relevance_tot(jj) > 0)then
           write(ilun,'(I6,X,I10,15(1X,1PE14.7),1X,I1)')jj&          
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
                ,e_bind_tot(jj)/(e_thermal_tot(jj)+e_kin_int_tot(jj)+tiny(0.d0))&
                ,e_bind_tot4(jj)/(e_thermal_tot4(jj)+e_kin_int_tot4(jj)+tiny(0.d0))&
                ,minmatch_tot(jj)

           rel_mass=rel_mass+clump_mass_tot(jj)*scale_d*dble(scale_l)**3/1.98892d33
           n_rel=n_rel+1
        end if
     end do

     else
     write(ilun,*)'Cl_N #leaf-cells  peak_x [uu] peak_y [uu] peak_z [uu] size_x [AU] size_y [AU]'//&
          ' size_z [AU] |v|_CM [u.u.] rho- [H/cc] rho+ [H/cc] rho_av [H/cc] M_cl [M_sol] V_cl [AU^3] rel. V/(U+Q) V4/(U4+Q4) m_match'
     do j=npeaks_tot,1,-1
        jj=sort_index(j)
        if (relevance_tot(jj) > 0)then
           write(ilun,'(I6,X,I10,3(X,F11.5),3(X,F11.5),X,F13.5,3(XE21.12E2),X,F13.5,XE11.2E2,X,F7.3,1X,F6.3,3X,F6.3,4X,I1)')jj&          
                ,n_cells_tot(jj)&
                ,peak_pos_tot(jj,1),peak_pos_tot(jj,2),peak_pos_tot(jj,3)&
                ,(5.*clump_size_tot(jj,1)/clump_vol_tot(jj))**0.5*(scale_l/1.496d13)&
                ,(5.*clump_size_tot(jj,2)/clump_vol_tot(jj))**0.5*(scale_l/1.496d13)&
                ,(5.*clump_size_tot(jj,3)/clump_vol_tot(jj))**0.5*scale_l/1.496d13&
                ,(clump_momentum_tot(jj,1)**2+clump_momentum_tot(jj,2)**2+ &
                clump_momentum_tot(jj,3)**2)**0.5/clump_mass_tot(jj)&
                ,min_dens_tot(jj)*scale_nH,max_dens_tot(jj)*scale_nH&
                ,clump_mass_tot(jj)/clump_vol_tot(jj)*scale_nH&
                ,clump_mass_tot(jj)*scale_d*dble(scale_l)**3/1.98892d33&
                ,clump_vol_tot(jj)*(scale_l/1.496d13)**3&
                ,relevance_tot(jj)&
                ,e_bind_tot(jj)/(e_thermal_tot(jj)+e_kin_int_tot(jj)+tiny(0.d0))&
                ,e_bind_tot4(jj)/(e_thermal_tot4(jj)+e_kin_int_tot4(jj)+tiny(0.d0))&
                ,minmatch_tot(jj)

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
     end if
     write(ilun,'(A,F9.3)')'total mass above threshold =',tot_mass*scale_d*dble(scale_l)**3/1.98892d33
     write(ilun,'(A,I6,A,F9.3)')'total mass in',n_rel,' listed clumps =',rel_mass
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

  integer::ipart,ip,ilevel,next_level
  integer::i,j,info,dummyint
  integer,dimension(1:nvector)::ind_cell


  ! saddle point array for 1 cpu
  allocate(saddle_dens(1:npeaks_tot,1:npeaks_tot))
  saddle_dens=0.

  ! loop 'testparts', pass the information of nvector parts to neighborsearch 
  ip=0
  do ipart=1,ntest
     ip=ip+1
     ilevel=levp(testp_sort(ipart)) ! level
     next_level=0
     if(ipart<ntest)next_level=levp(testp_sort(ipart+1)) !level of next particle
     ind_cell(ip)=icellp(testp_sort(ipart))
     if(ip==nvector .or. next_level /= ilevel)then
        call neighborsearch(ind_cell,ip,dummyint,ilevel,4)
        ip=0
     endif
  end do
  if (ip>0)call neighborsearch(ind_cell,ip,dummyint,ilevel,4)

  ! share the results among MPI domains  
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(saddle_dens,saddle_dens_tot,(npeaks_tot**2),MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  saddle_dens_tot=saddle_dens
#endif

  ! check symmetry (could be done only in case of verbose == true)
  do i=1,npeaks_tot
     do j=1,i
        if(saddle_dens_tot(i,j)/=saddle_dens_tot(j,i).and.myid==1)then 
           write(*,*),'Alert! asymmetric saddle point array!',i,j,saddle_dens_tot(i,j),saddle_dens_tot(j,i)
        endif
     end do
  end do

  ! compute saddle_max value and relevance
  saddle_max_tot=maxval(saddle_dens_tot,dim=1)
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
subroutine merge_clumps(ntest)
  use amr_commons
  use clfind_commons
  implicit none
  integer::ntest
  !---------------------------------------------------------------------------
  ! 
  !---------------------------------------------------------------------------
  integer::j,i,ii,merge_count,final_peak,merge_to,ipart
  real(kind=8)::max_val
  integer,dimension(1:npeaks_tot)::old_peak
  integer::peak,next_peak

  if (verbose)write(*,*)'Now merging clumps'

  ! Sort clumps by peak density in ascending order
  call heapsort_index(max_dens_tot,sort_index,npeaks_tot)

  do i=1,npeaks_tot
     ii=sort_index(i)
     new_peak(ii)=ii

     ! If the relevance is below the threshold -> merge
     if (relevance_tot(ii)<relevance_threshold) then
        
        ! Go through the ii-th line in the saddle point array to find the neighbor to merge to
        merge_to=0; max_val=0.
        do j=1,npeaks_tot
           if (saddle_dens_tot(ii,j)>max_val)then
              merge_to=j
              max_val=saddle_dens_tot(ii,j)
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
        do j=1,npeaks_tot
           if (merge_to>0)then
              if(saddle_dens_tot(ii,j)>saddle_dens_tot(merge_to,j))then
                 saddle_dens_tot(merge_to,j)=saddle_dens_tot(ii,j)
                 saddle_dens_tot(j,merge_to)=saddle_dens_tot(ii,j)
              end if
              saddle_dens_tot(merge_to,merge_to)=0.
           end if
           saddle_dens_tot(ii,j)=0.        
           saddle_dens_tot(j,ii)=0.
        end do
        
        ! Update saddle_max value
        if (merge_to>0)then
           saddle_max_tot(merge_to)=0
           do j=1,npeaks_tot
              if (saddle_dens_tot(merge_to,j)>saddle_max_tot(merge_to))then
                 saddle_max_tot(merge_to)=saddle_dens_tot(merge_to,j)
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

  !update flag 2
   do ipart=1,ntest
     flag2(icellp(ipart))=new_peak(flag2(icellp(ipart)))
  end do

  if (verbose)write(*,*)'Done merging clumps'  
end subroutine merge_clumps
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
  allocate(clump_vol_tot(1:npeaks_tot))
  allocate(saddle_max_tot(1:npeaks_tot))
  allocate(relevance_tot(1:npeaks_tot))
  allocate(sort_index(1:npeaks_tot))
  allocate(saddle_dens_tot(1:npeaks_tot,1:npeaks_tot))
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

  !initialize all peak based arrays
  n_cells_tot=0
  saddle_max_tot=0.
  relevance_tot=1.
  clump_size_tot=0.
  min_dens_tot=huge(0.d0)
  max_dens_tot=0.
  av_dens_tot=0.
  clump_mass_tot=0.
  clump_vol_tot=0.
  peak_pos_tot=0.
  center_of_mass_tot=0.
  second_moments=0.; second_moments_tot=0.
  saddle_dens_tot=0.
  clump_momentum_tot=0.
  e_kin_int_tot=0.; e_kin_int_tot4=0.
  e_bind_tot=0.; e_bind_tot4=0.
  e_thermal_tot=0.; e_thermal_tot4=0.
  phi_min_tot=0.
  minmatch_tot=1
  new_peak=0

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
  deallocate(clump_momentum_tot)
  deallocate(e_kin_int_tot,e_kin_int_tot4)
  deallocate(e_bind_tot,e_bind_tot4)
  deallocate(e_thermal_tot,e_thermal_tot4)
  deallocate(phi_min_tot)
  deallocate(minmatch_tot)
  deallocate(saddle_dens_tot)
  deallocate(peak_pos_tot)
  deallocate(new_peak)

end subroutine deallocate_all
!################################################################                 
!################################################################ 
!################################################################                 
!################################################################     
subroutine clump_phi
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

  real(dp),dimension(1:npeaks_tot)::phi_min
  integer,dimension(1:npeaks_tot)::minmatch

  phi_min=0.d0; minmatch=1

  dx=0.5D0**nlevelmax 
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale  


  N=10
  do jj=1,npeaks_tot
     if (relevance_tot(jj)>0)then
        pos(1,1:3)=peak_pos_tot(jj,1:3)
        call cmp_cpumap(pos,cc,1)
        if(cc(1) == myid) then
           call get_cell_index(cell_index,cell_levl,pos,nlevelmax,1)
           phim=phi(cell_index(1))

           ! Check for neighbors

           ! one cell offset
           do k1=-1,1
              do j1=-1,1
                 do i1=-1,1
                    xtest(1,1)=pos(1,1)+i1*dx_loc
                    xtest(1,2)=pos(1,2)+j1*dx_loc
                    xtest(1,3)=pos(1,3)+k1*dx_loc
                    call get_cell_index(ind_cell,lev_cell,xtest,nlevelmax,1)

                    if(phi(ind_cell(1)) < phim)then
                       !print*,"offset min by 1c",myid,jj,k1,j1,i1,phim,phi(ind_cell(1))
                       phim=phi(ind_cell(1))
                    end if
                 end do
              end do
           end do

           ! two cells offset
           do k1=-2,2
              do j1=-2,2
                 do i1=-2,2
                    xtest(1,1)=pos(1,1)+i1*dx_loc
                    xtest(1,2)=pos(1,2)+j1*dx_loc
                    xtest(1,3)=pos(1,3)+k1*dx_loc
                    call get_cell_index(ind_cell,lev_cell,xtest,nlevelmax,1)
                    if(phi(ind_cell(1)) < phim)then
                       phim=phi(ind_cell(1))
                       !print*,"offset min by 2c",myid,jj,k1,j1,i1,phim,phi(ind_cell(1))
                       minmatch(jj)=0
                    end if
                 end do
              end do
           end do
           !   phi_min(jj)=phim                                
        end if
        
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
                             phi_min(jj)=min(phi_min(jj),phi(cell_index(1)))                          
                          endif
                       end do                        
                       itest=0
                    endif
                 endif
              end do
           end do
        end do
     
        call cmp_cpumap(xtest_cpu,cc,itest)
        do i=1,nvector
           if(cc(i) == myid) then
              xtest_ind(1,1:3)=xtest_cpu(i,1:3)
              call get_cell_index(cell_index,cell_levl,xtest_ind,nlevelmax,1)
              phi_min(jj)=min(phi_min(jj),phi(cell_index(1)))
           endif
        end do
        
     end if ! end if relevance > than..
  end do ! end loop over clumps



#ifndef WITHOUTMPI     
  call MPI_ALLREDUCE(phi_min,phi_min_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI     
  phi_min_tot=phi_min
#endif
#ifndef WITHOUTMPI     
  call MPI_ALLREDUCE(minmatch,minmatch_tot,npeaks_tot,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI     
  minmatch_tot=minmatch
#endif
end subroutine clump_phi
