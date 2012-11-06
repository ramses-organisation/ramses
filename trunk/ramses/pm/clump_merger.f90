subroutine compute_clump_properties(ntest,ntest_all)
  use amr_commons
  use hydro_commons
  use pm_commons
  use clfind_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  integer::ntest,ntest_all
  !----------------------------------------------------------------------------
  ! this subroutine performs one loop over all particles and collects the 
  ! relevant information from the cells where the particles sit. after a lot
  ! of mpi-communication, all necessary peak-patch properties can be computed
  !----------------------------------------------------------------------------

  ! variables used for the loop over all particles
  integer::igrid,jgrid,ipart,jpart,next_part,ig,ip,npart1,ilevel
  
  integer::info,i,j,jj,peak_nr,nv

  !variables needed temporarily store cell properties
  real(kind=8)::d,vol
  real(kind=8),dimension(1:3)::vd

  ! variables to be used for vector-sweeps
  integer,dimension(1:nvector)::ind_grid,ind_cell,init_ind_cell,init_cell_lev,cell_lev
  integer,dimension(1:nvector)::ind_part,ind_grid_part,indv
  real(dp),dimension(1:nvector,1:ndim)::pos,init_pos

  ! variables related to the size of a cell on a given level
  real(kind=8)::dx,dx_loc,scale,vol_loc
  real(kind=8),dimension(1:nlevelmax)::volume
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  integer::nx_loc,ind
  integer ::ix,iy,iz

  !peak-patch related arrays before sharing information with other cpus
  real(kind=8),dimension(0:npeaks_tot)::max_dens
  real(kind=8),dimension(1:npeaks_tot)::min_dens,av_dens,clump_mass,clump_vol
  real(kind=8),dimension(1:npeaks_tot,1:3)::center_of_mass,clump_momentum
  integer,dimension(1:npeaks_tot)::n_cells

  !units
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2

  min_dens=huge(0.d0);  max_dens=0.d0; n_cells=0
  clump_mass=0.d0; clump_vol=0.d0; clump_momentum=0.d0; center_of_mass=0.d0

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

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

  !---------------------------------------------------------------------------
  ! big loop over all parts to collect information from the cells
  !---------------------------------------------------------------------------
  ! loop over all particles
  do ipart=1,ntest
     nv=1
     ig=1
     ip=1
     ilevel=levp(ipart) !level
     indv(nv)=1+icellp(ipart)/ngridmax
     ind_grid(nv)=icellp(ipart)-ncoarse-(indv(nv)-1)*ngridmax !grid index
     ind_part(nv)=ipart
     ind_grid_part(nv)=1

     dx=0.5D0**ilevel
     nx_loc=(icoarse_max-icoarse_min+1)
     skip_loc(1)=dble(icoarse_min)
     skip_loc(2)=dble(jcoarse_min)
     skip_loc(3)=dble(kcoarse_min)
     scale=boxlen/dble(nx_loc)

     do ind=1,twotondim
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(ind,1)=(dble(ix)-0.5D0)*dx
        xc(ind,2)=(dble(iy)-0.5D0)*dx
        xc(ind,3)=(dble(iz)-0.5D0)*dx
     end do
     init_pos(nv,1)=(xg(ind_grid(nv),1)+xc(indv(nv),1)-skip_loc(1))*scale
     init_pos(nv,2)=(xg(ind_grid(nv),2)+xc(indv(nv),2)-skip_loc(2))*scale
     init_pos(nv,3)=(xg(ind_grid(nv),3)+xc(indv(nv),3)-skip_loc(3))*scale

     ! find peak_nr
     peak_nr=flag2(icellp(ipart))
     d=dble(uold(icellp(ipart),1))
     do i=1,ndim
        vd(i)=uold(icellp(ipart),i+1)
     end do
     vol=volume(levp(ipart))

     if (peak_nr /=0 ) then
        ! number of leaf cells per clump
        n_cells(peak_nr)=n_cells(peak_nr)+1
        
        ! find min density
        if(d<=min_dens(peak_nr))min_dens(peak_nr)=d
        ! find max density
        if(d>=max_dens(peak_nr))max_dens(peak_nr)=d

        ! find clump mass
        clump_mass(peak_nr)=clump_mass(peak_nr)+vol*d
        
        ! center of mass velocity
        do i=1,3
           clump_momentum(peak_nr,i)=clump_momentum(peak_nr,i)+vd(i)*vol
        end do
        
        ! clump size (maybe take center of mass instead of peak as reference point)
        do i=1,ndim
           center_of_mass(peak_nr,i)=center_of_mass(peak_nr,i)+dble(init_pos(nv,i))*vol*d
           
           ! compute second order moments
           do j=1,ndim
              second_moments(peak_nr,i,j)=second_moments(peak_nr,i,j)+init_pos(nv,i)*init_pos(nv,j)*vol*d
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
  call MPI_ALLREDUCE(max_dens,max_dens_tot,npeaks_tot+1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
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

  !---------------------------------------------------------------------------
  ! compute further information related to the clumps
  !---------------------------------------------------------------------------
  tot_mass=0.
  do jj=1,npeaks_tot
     !calculate total mass above threshold
     tot_mass=tot_mass+clump_mass_tot(jj) 
     if (clump_mass_tot(jj) > 0.1d-80) then
        av_dens_tot(jj)=clump_mass_tot(jj)/clump_vol_tot(jj)
        center_of_mass_tot(jj,1:ndim)=center_of_mass_tot(jj,1:ndim)/clump_mass_tot(jj)
     end if
  end do

end subroutine compute_clump_properties
!################################################################
!################################################################
!################################################################
!################################################################
subroutine compute_clump_properties_round2(ntest,ntest_all)
  use poisson_commons
  use amr_commons
  use hydro_commons
  use pm_commons
  use clfind_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  integer::ntest,ntest_all
  !----------------------------------------------------------------------------                
  ! this subroutine performs another loop over all particles and collects the 
  ! information more information like the velocity of a cell realtive to the 
  ! center of mass
  !----------------------------------------------------------------------------

  ! variables used for the loop over all particles
  integer::igrid,jgrid,ipart,jpart,next_part,ig,ip,npart1,ilevel
  
  ! 
  integer::info,i,j,jj,peak_nr,nv

  !variables needed temporarily store cell properties
  real(kind=8)::d,vol,M,ekk,phi_rel,de
  real(kind=8),dimension(1:3)::vd,xcell,xpeak

  ! variables to be used with vector-sweeps
  integer,dimension(1:nvector)::ind_grid,ind_cell,init_ind_cell,init_cell_lev,cell_lev
  integer,dimension(1:nvector)::ind_part,ind_grid_part,indv
  real(dp),dimension(1:nvector,1:ndim)::pos,init_pos

  ! variables related to the size of a cell on a given level
  real(kind=8)::dx,dx_loc,scale,vol_loc
  real(kind=8),dimension(1:nlevelmax)::volume
  integer::nx_loc,ix,iy,iz,ind
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc

  !peak-patch related arrays before sharing information with other cpus
  real(kind=8),dimension(1:npeaks_tot)::e_kin_int,e_bind,e_thermal,e_kin_int4,e_bind4,e_thermal4
  real(kind=8),dimension(1:npeaks_tot,1:3)::clump_size

  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  e_kin_int=0.d0; clump_size=0.d0; e_bind=0.d0; e_thermal=0.d0; e_bind4=0.d0; e_thermal4=0.d0; e_kin_int4=0.d0

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

  !---------------------------------------------------------------------------
  ! big loop over all parts to collect information from the cells
  !---------------------------------------------------------------------------
  ! loop over all particles
  do ipart=1,ntest
     nv=1
     ig=1
     ip=1
     ilevel=levp(ipart) !level
     indv(nv)=1+icellp(ipart)/ngridmax
     ind_grid(nv)=icellp(ipart)-ncoarse-(indv(nv)-1)*ngridmax !grid index
     ind_part(nv)=ipart
     ind_grid_part(nv)=1

     dx=0.5D0**ilevel
     nx_loc=(icoarse_max-icoarse_min+1)
     skip_loc(1)=dble(icoarse_min)
     skip_loc(2)=dble(jcoarse_min)
     skip_loc(3)=dble(kcoarse_min)
     scale=boxlen/dble(nx_loc)

     do ind=1,twotondim
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(ind,1)=(dble(ix)-0.5D0)*dx
        xc(ind,2)=(dble(iy)-0.5D0)*dx
        xc(ind,3)=(dble(iz)-0.5D0)*dx
     end do
     init_pos(nv,1)=(xg(ind_grid(nv),1)+xc(indv(nv),1)-skip_loc(1))*scale
     init_pos(nv,2)=(xg(ind_grid(nv),2)+xc(indv(nv),2)-skip_loc(2))*scale
     init_pos(nv,3)=(xg(ind_grid(nv),3)+xc(indv(nv),3)-skip_loc(3))*scale

     ! use peak number after merge
     peak_nr=new_peak(flag2(icellp(ipart))) 
     ! make sure flag2 is properly updated after merge
     flag2(icellp(ipart))=peak_nr

     pos(nv,1)=peak_pos_tot(peak_nr,1)
     pos(nv,2)=peak_pos_tot(peak_nr,2)
     pos(nv,3)=peak_pos_tot(peak_nr,3)
     call get_cell_index(ind_cell,cell_lev,init_pos,ilevel,ip)

     if(cell_lev(nv) /= ilevel)write(*,*)'alert_problem in get cell index'

     d=dble(uold(ind_cell(nv),1))
     de=uold(ind_cell(nv),ndim+2)
     
     if (peak_nr /=0 ) then
        
        do i=1,ndim
           vd(i)=uold(ind_cell(nv),i+1)
           xcell(i)=init_pos(nv,i)
           xpeak(i)=peak_pos_tot(peak_nr,i)
        end do
        M=clump_mass_tot(peak_nr)
        vol=volume(cell_lev(nv))                  
        phi_rel=(phi(ind_cell(nv))-phi_min_tot(peak_nr))*scale
       
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
        
        ! potential energy of the clump (crude approach)
        ! e_bind(peak_nr)=e_bind(peak_nr)+M*d*vol/(rr+0.5*vol**0.3333333333)*0.4              
        
        ! potential energy using the acutal phi W= 0.5*int phi*rho
        e_bind(peak_nr)=e_bind(peak_nr)+phi_rel*d*vol*5.d-1
        ! add gravitational self energy of every cell
        e_bind(peak_nr)=e_bind(peak_nr)+9.d-1*d**2.d0*vol**(5.d0/3.d0)
        
        
        ! size relative to center of mass 
        do i=1,ndim
           clump_size(peak_nr,i)=clump_size(peak_nr,i)+(xcell(i)-center_of_mass_tot(peak_nr,i))**2.d0*vol
        end do
        
        ! thermal energy
        e_thermal(peak_nr)=e_thermal(peak_nr)+(de-ekk)*vol*(gamma-1)/(g_star-1)
        
        ! repeat same for smaller region if cell is close enough
        if (((xpeak(1)-xcell(1))**2.+(xpeak(2)-xcell(2))**2.+(xpeak(3)-xcell(3))**2.) .LE. 16.*volume(nlevelmax)**(2./3.))then
           do i=1,3
              e_kin_int4(peak_nr)=e_kin_int4(peak_nr)+(vd(i)/d-clump_momentum_tot(peak_nr,i)/M)**2*d*vol*0.5
           end do
           
           e_bind4(peak_nr)=e_bind4(peak_nr)+phi_rel*d*vol*0.5
           e_bind4(peak_nr)=e_bind4(peak_nr)+0.9*d**2.*vol**(5./3.)
           e_thermal4(peak_nr)=e_thermal4(peak_nr)+(de-ekk)*vol*(gamma-1)/(g_star-1)
        end if
        
     end if
  end do

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
  use hydro_commons
  use pm_commons
  use clfind_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  logical::to_file
  !---------------------------------------------------------------------------
  ! this routine writes the clump properties to screen and to file
  !---------------------------------------------------------------------------
  integer::j,jj,ilevel,n_rel,ilun
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
     write(ilun,*)'Cl_N #leaf-cells  peak_x [uu] peak_y [uu] peak_z [uu] size_x [AU] size_y [AU]'//&
          ' size_z [AU] |v|_CM [u.u.] rho- [H/cc] rho+ [H/cc] rho_av [H/cc] M_cl [M_sol] V_cl [AU^3] rel. V/(U+Q) V4/(U4+Q4) m_match'
     do j=npeaks_tot,1,-1
        jj=sort_index(j)
        if (relevance_tot(jj) > 0)then          
           write(ilun,'(I6,X,I10,3(X,F11.5),3(X,F11.5),X,F13.5,3(XE11.2E2),X,F13.5,XE11.2E2,X,F7.3,1X,F6.3,3X,F6.3,4X,I1)')jj&
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
  use hydro_commons
  use pm_commons
  use clfind_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  integer::ntest
  !---------------------------------------------------------------------------
  ! subroutine which creates a npeaks_tot**2 sized array of saddlepoint densities
  !---------------------------------------------------------------------------

  ! variables for particle loop
  integer::igrid,jgrid,ipart,jpart,next_part,ig,ip,npart1,ilevel,nv
  ! dummy variables
  integer::jj,i,j,kk,info,peak_nr
  ! variables to be used for vector-sweeps
  integer,dimension(1:nvector),save::ind_grid,ind_cell,init_ind_cell,init_cell_lev,cell_lev,clump_nr
  integer,dimension(1:nvector),save::ind_part,ind_grid_part,indv
  real(dp),dimension(1:nvector,1:ndim)::pos,init_pos

  ! variables related to the size of a cell on a given level
  real(kind=8)::dx,dx_loc,scale,vol_loc
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  integer::nx_loc,ind,ix,iy,iz

  ! saddle point array for 1 cpu
  real(kind=8),dimension(1:npeaks_tot,1:npeaks_tot)::saddle_dens
  saddle_dens=0.

  !---------------------------------------------------------------------------
  ! big loop over all parts, pass the information of nvector parts to 
  ! subroutine find_best_neighbor
  !---------------------------------------------------------------------------
  do ipart=1,ntest
     nv=1
     ig=1
     ip=1
     ilevel=levp(ipart) !level
     indv(nv)=1+icellp(ipart)/ngridmax
     ind_grid(nv)=icellp(ipart)-ncoarse-(indv(nv)-1)*ngridmax !grid index
     ind_part(nv)=ipart
     ind_grid_part(nv)=1

     dx=0.5D0**ilevel
     nx_loc=(icoarse_max-icoarse_min+1)
     skip_loc(1)=dble(icoarse_min)
     skip_loc(2)=dble(jcoarse_min)
     skip_loc(3)=dble(kcoarse_min)
     scale=boxlen/dble(nx_loc)

     do ind=1,twotondim
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(ind,1)=(dble(ix)-0.5D0)*dx
        xc(ind,2)=(dble(iy)-0.5D0)*dx
        xc(ind,3)=(dble(iz)-0.5D0)*dx
     end do

     ! test particle position
     init_pos(nv,1)=(xg(ind_grid(nv),1)+xc(indv(nv),1)-skip_loc(1))*scale
     init_pos(nv,2)=(xg(ind_grid(nv),2)+xc(indv(nv),2)-skip_loc(2))*scale
     init_pos(nv,3)=(xg(ind_grid(nv),3)+xc(indv(nv),3)-skip_loc(3))*scale

     ! peak number
     peak_nr=flag2(icellp(ipart))

     ! local arrays
     clump_nr(nv)=flag2(icellp(ipart))
     ind_grid(nv)=mod((icellp(ipart)-ncoarse),ngridmax)

     call find_best_neighbor(ind_grid,clump_nr,init_pos,ip,ilevel,saddle_dens)

  end do
  ! end loop over all particles
  
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(saddle_dens,saddle_dens_tot,(npeaks_tot**2),MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  saddle_dens_tot=saddle_dens
#endif

  ! symmetrize saddle_point array
  do i=1,npeaks_tot
     do j=1,i
        saddle_dens_tot(i,j)=0.5*(saddle_dens_tot(i,j)+saddle_dens_tot(j,i))
        saddle_dens_tot(j,i)=saddle_dens_tot(i,j)
     end do
  end do

  ! compute saddle_max value and relevance
  do i=1,npeaks_tot
     do j=1,npeaks_tot
        if (saddle_dens_tot(i,j)>saddle_max_tot(i))then
           saddle_max_tot(i)=saddle_dens_tot(i,j)
        end if
     end do
     if (saddle_max_tot(i)>1.d-40)then
        relevance_tot(i)=max_dens_tot(i)/saddle_max_tot(i)
     else
        relevance_tot(i)=max_dens_tot(i)/min_dens_tot(i)
     end if
  end do

end subroutine saddlepoint_search
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine find_best_neighbor(ind_grid_part,clump_nr,pos,np,ilevel,saddle_dens)
  use amr_commons
  use pm_commons
  use poisson_commons
  use hydro_commons, ONLY: uold
  use clfind_commons
  implicit none
  integer::np,ilevel
  integer,dimension(1:nvector)::ind_grid_part,cell_levl,clump_nr,ind_cell
  ! saddle point array for 1 cpu
  real(kind=8),dimension(1:npeaks_tot,1:npeaks_tot)::saddle_dens
  !------------------------------------------------------------
  ! 
  ! 
  ! 
  ! 
  ! 
  !------------------------------------------------------------
  logical::error
  integer::j,ind,idim,nx_loc
  integer::i1,j1,k1,i2,j2,k2
  real(kind=8)::dx,dx_loc,scale,vol_loc
  integer::i1min,i1max,j1min,j1max,k1min,k1max
  integer::i2min,i2max,j2min,j2max,k2min,k2max
  integer::i3min,i3max,j3min,j3max,k3min,k3max
  ! Grid-based arrays

  real(dp),dimension(1:nvector,1:ndim),save::x0
  ! Particle-based arrays

  real(dp),dimension(1:nvector,1:ndim)::x,xtest,pos
  integer ,dimension(1:nvector,1:ndim),save::ig,id
  integer ,dimension(1:nvector)::cc

  real(kind=8),dimension(1:nvector,1:ndim,1:twotondim),save::xpart

  real(kind=8),dimension(1:3)::skip_loc

  ! Mesh spacing in that level
  dx=0.5D0**ilevel 
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**3

  ! Integer constants
  i1min=0; i1max=0; i2min=0; i2max=0; i3min=1; i3max=1
  j1min=0; j1max=0; j2min=0; j2max=0; j3min=1; j3max=1
  k1min=0; k1max=0; k2min=0; k2max=0; k3min=1; k3max=1
  if(ndim>0)then
     i1max=2; i2max=3; i3max=2
  end if
  if(ndim>1)then
     j1max=2; j2max=3; j3max=2
  end if
  if(ndim>2)then
     k1max=2; k2max=3; k3max=2
  end if

  do j=1,np
     xtest(j,1:ndim)=pos(j,1:ndim)
  end do
  !====================================================
  ! Check for potential new positions at level ilevel-1
  !====================================================
  if(ilevel>levelmin)then

     ! Lower left corner of 3x3x3 grid-cube
     do idim=1,ndim
        do j=1,np
           x0(j,idim)=xg(ind_grid_part(j),idim)-3.0D0*dx
        end do
     end do

     ! Compute parent cell position
     do j=1,np
        ! Particle 1
        xpart(j,1,1)=0.5+ig(j,1)
        xpart(j,2,1)=0.5+ig(j,2)
        xpart(j,3,1)=0.5+ig(j,3)
        ! Particle 2
        xpart(j,1,2)=0.5+id(j,1)
        xpart(j,2,2)=0.5+ig(j,2)
        xpart(j,3,2)=0.5+ig(j,3)
        ! Particle 3
        xpart(j,1,3)=0.5+ig(j,1)
        xpart(j,2,3)=0.5+id(j,2)
        xpart(j,3,3)=0.5+ig(j,3)
        ! Particle 4
        xpart(j,1,4)=0.5+id(j,1)
        xpart(j,2,4)=0.5+id(j,2)
        xpart(j,3,4)=0.5+ig(j,3)
        ! Particle 5
        xpart(j,1,5)=0.5+ig(j,1)
        xpart(j,2,5)=0.5+ig(j,2)
        xpart(j,3,5)=0.5+id(j,3)
        ! Particle 6
        xpart(j,1,6)=0.5+id(j,1)
        xpart(j,2,6)=0.5+ig(j,2)
        xpart(j,3,6)=0.5+id(j,3)
        ! Particle 7
        xpart(j,1,7)=0.5+ig(j,1)
        xpart(j,2,7)=0.5+id(j,2)
        xpart(j,3,7)=0.5+id(j,3)
        ! Particle 8
        xpart(j,1,8)=0.5+id(j,1)
        xpart(j,2,8)=0.5+id(j,2)
        xpart(j,3,8)=0.5+id(j,3)
     end do

     ! Test those particles
     do ind=1,twotondim
        do idim=1,ndim
           do j=1,np
              xtest(j,idim)=xpart(j,idim,ind)*2.*dx+x0(j,idim)
           end do
           do j=1,np
              xtest(j,idim)=(xtest(j,idim)-skip_loc(idim))*scale
           end do
        end do

        call get_cell_index(ind_cell,cell_levl,xtest,ilevel-1,np)

        do j=1,np
           if(son(ind_cell(j))==0 .and. cell_levl(j)==(ilevel-1) .and. flag2(ind_cell(j))/=0 )then
              if(clump_nr(j) /= flag2(ind_cell(j)).and.uold(ind_cell(j),1)>saddle_dens(clump_nr(j),flag2(ind_cell(j))))then   
                 saddle_dens(clump_nr(j),flag2(ind_cell(j)))=uold(ind_cell(j),1)
              end if
           endif
        end do
     end do

  endif

  !====================================================
  ! Check for potential new positions at level ilevel
  !====================================================
  ! Generate 3x3x3 neighboring cells at level ilevel
  do k1=k1min,k1max
     do j1=j1min,j1max
        do i1=i1min,i1max
           do j=1,np
              xtest(j,1)=pos(j,1)+(i1-1)*dx_loc
              xtest(j,2)=pos(j,2)+(j1-1)*dx_loc
              xtest(j,3)=pos(j,3)+(k1-1)*dx_loc
           end do
           
           call get_cell_index(ind_cell,cell_levl,xtest,ilevel,np)

           do j=1,np
              if(son(ind_cell(j))==0 .and. cell_levl(j)==ilevel .and. flag2(ind_cell(j))/=0)then
                 if((clump_nr(j) /= flag2(ind_cell(j))).and.(uold(ind_cell(j),1)>saddle_dens(clump_nr(j),flag2(ind_cell(j)))))then   
                    saddle_dens(clump_nr(j),flag2(ind_cell(j)))=uold(ind_cell(j),1)
                 end if
              end if
           end do
        end do
     end do
  end do

  !====================================================
  ! Check for potential new positions at level ilevel+1
  !====================================================
  if(ilevel<nlevelmax)then

     ! Generate 4x4x4 neighboring cells at level ilevel+1
     do k2=k2min,k2max
        do j2=j2min,j2max
           do i2=i2min,i2max

              do j=1,np
                 xtest(j,1)=pos(j,1)+(i2-1.5)*dx_loc/2.0
                 xtest(j,2)=pos(j,2)+(j2-1.5)*dx_loc/2.0
                 xtest(j,3)=pos(j,3)+(k2-1.5)*dx_loc/2.0
              end do

              call get_cell_index(ind_cell,cell_levl,xtest,ilevel+1,np)

              do j=1,np
                 if(son(ind_cell(j))==0 .and. cell_levl(j)==(ilevel+1) .and. flag2(ind_cell(j))/=0)then
                    if(clump_nr(j) /= flag2(ind_cell(j)).and.uold(ind_cell(j),1)>saddle_dens(clump_nr(j),flag2(ind_cell(j))))then
                       saddle_dens(clump_nr(j),flag2(ind_cell(j)))=uold(ind_cell(j),1)   
                    end if
                 end if
              end do
           end do
        end do
     end do

  endif

end subroutine find_best_neighbor
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine merge_clumps(ntest)
  use amr_commons
  use pm_commons
  use poisson_commons
  use hydro_commons, ONLY: uold
  use clfind_commons
  implicit none

  integer::ntest
  integer::j,jj,i,ii
  integer::igrid,jgrid,ipart,jpart,next_part,npart1,info,ilevel,merge_to
  real(kind=8)::max_val
  logical::merging


  if (verbose)write(*,*)'Now merging clumps'

  ! Sort clumps by peak density in ascending order
  call heapsort_index(max_dens_tot,sort_index,npeaks_tot)

  do i=1,npeaks_tot
     ii=sort_index(i)
     new_peak(ii)=ii

     ! If the relevance is below the threshold -> merge
     if (relevance_tot(ii)<relevance_threshold .and. relevance_tot(ii) > 1.d-80) then
        
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
        if(verbose)write(*,*)'clump ',ii,'merged to ',merge_to
        
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
  
  if (verbose)write(*,*)'Done merging clumps'  
end subroutine merge_clumps
!################################################################
!################################################################
!################################################################
!################################################################
subroutine write_peak_map(ntest)
  use amr_commons
  use hydro_commons
  use pm_commons
  use clfind_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  integer::ntest
  !----------------------------------------------------------------------------                
  ! One loop over all the particles creating a textfile for each cpu
  ! -number of parts on cpu
  ! -xpos for all parts
  ! -ypos for all parts
  ! -zpos for all parts
  ! -clump number for all parts
  !----------------------------------------------------------------------------

  integer::igrid,jgrid,ipart,jpart,next_part,ig,ip,npart1,info,ilevel
  integer::jj,i,kk,nv,peak_nr
  integer,dimension(1:nvector)::ind_grid,ind_cell,init_ind_cell,init_cell_lev,cell_lev
  integer,dimension(1:nvector)::ind_part,ind_grid_part,indv
  real(dp),dimension(1:nvector,1:ndim)::pos,init_pos
  character(LEN=5)::myidstring,nchar

  ! variables related to the size of a cell on a given level                                                                                                                
  real(kind=8)::dx,dx_loc,scale,vol_loc
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  integer::nx_loc,ind,ix,iy,iz


#if NDIM>2
  call title(ifout-1,nchar)
  call title(myid,myidstring)
  open(unit=20,file=TRIM('output_'//TRIM(nchar)//'/clump_map.csv'//myidstring),form='formatted')
  !write(20,*)nparts

  !loop over all particles
  do ipart=1,ntest
     nv=1
     ig=1
     ip=1
     ilevel=levp(ipart) !level                                                                                                                                              
     indv(nv)=1+icellp(ipart)/ngridmax
     ind_grid(nv)=icellp(ipart)-ncoarse-(indv(nv)-1)*ngridmax !grid index                                                                                                   
     ind_part(nv)=ipart
     ind_grid_part(nv)=1

     dx=0.5D0**ilevel
     nx_loc=(icoarse_max-icoarse_min+1)
     skip_loc=(/0.0d0,0.0d0,0.0d0/)
     if(ndim>0)skip_loc(1)=dble(icoarse_min)
     if(ndim>1)skip_loc(2)=dble(jcoarse_min)
     if(ndim>2)skip_loc(3)=dble(kcoarse_min)
     scale=boxlen/dble(nx_loc)
     do ind=1,twotondim
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(ind,1)=(dble(ix)-0.5D0)*dx
        xc(ind,2)=(dble(iy)-0.5D0)*dx
        xc(ind,3)=(dble(iz)-0.5D0)*dx
     end do

     init_pos(nv,1)=(xg(ind_grid(nv),1)+xc(indv(nv),1)-skip_loc(1))*scale
     init_pos(nv,2)=(xg(ind_grid(nv),2)+xc(indv(nv),2)-skip_loc(2))*scale
     init_pos(nv,3)=(xg(ind_grid(nv),3)+xc(indv(nv),3)-skip_loc(3))*scale

     ! peak number                                                                                                                                                         
     peak_nr=flag2(icellp(ipart))

     if (flag2(peak_nr)>0)then
        write(20,'(F11.8,A,F11.8,A,F11.8,A,I8)')init_pos(nv,1),',',init_pos(nv,2),','&
             ,init_pos(nv,3),',',flag2(icellp(ipart))
     end if
  end do
  !end loop over all particles

  close(20)
#endif

end subroutine write_peak_map
!################################################################                 
!################################################################ 
!################################################################                 
!################################################################     
subroutine allocate_peak_patch_arrays
  use amr_commons
  use hydro_commons
  use pm_commons
  use clfind_commons
  implicit none

  ! Allocate peak-patch_properties
  allocate(n_cells_tot(1:npeaks_tot))
  allocate(clump_size_tot(1:npeaks_tot,1:ndim))
  allocate(center_of_mass_tot(1:npeaks_tot,1:ndim))
  allocate(second_moments(1:npeaks_tot,1:ndim,1:ndim)); allocate(second_moments_tot(1:npeaks_tot,1:ndim,1:ndim))
  allocate(min_dens_tot(1:npeaks_tot))
  allocate(av_dens_tot(1:npeaks_tot))
  allocate(max_dens_tot(0:npeaks_tot))
  allocate(clump_mass_tot(1:npeaks_tot))
  allocate(clump_vol_tot(1:npeaks_tot))
  allocate(saddle_max_tot(npeaks_tot))
  allocate(relevance_tot(0:npeaks_tot))
  allocate(sort_index(npeaks_tot))
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
  use amr_commons
  use hydro_commons
  use pm_commons
  use clfind_commons
  implicit none

  deallocate(n_cells_tot)
  deallocate(clump_size_tot)
  deallocate(center_of_mass_tot)
  deallocate(second_moments); deallocate(second_moments_tot)
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

  integer::k1,j1,i1,jj,info
  integer,dimension(1:nvector)::cell_index,cell_levl,ind_cell,lev_cell,cc
  real(dp),dimension(1:nvector,1:3)::pos,xtest
  real(kind=8)::dx,dx_loc,scale,vol_loc,phim
  integer::nx_loc

  real(kind=8),dimension(1:npeaks_tot)::phi_min
  integer,dimension(1:npeaks_tot)::minmatch

  phi_min=0.d0; minmatch=1

  dx=0.5D0**nlevelmax 
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale  
  
  do jj=1,npeaks_tot
     if (relevance_tot(jj)>0)then
        pos(1,1:3)=peak_pos_tot(jj,1:3)
        call cmp_cpumap(pos,cc,1)
        if(cc(1) == myid) then
           call get_cell_index(cell_index,cell_levl,pos,nlevelmax,1)
           phim=phi(cell_index(1))

           !Check for neighbors

           !one cell offset
           do k1=-1,1
              do j1=-1,1
                 do i1=-1,1
                    xtest(1,1)=pos(1,1)+i1*dx_loc
                    xtest(1,2)=pos(1,2)+j1*dx_loc
                    xtest(1,3)=pos(1,3)+k1*dx_loc
                    call get_cell_index(ind_cell,lev_cell,xtest,nlevelmax,1)
                    call cmp_cpumap(xtest,cc,1)

                    if(phi(ind_cell(1)) < 1.0000000001*phim)then
                       !print*,"offset min by 1c",myid,jj,k1,j1,i1,phim,phi(ind_cell(1))
                       phim=phi(ind_cell(1))
                    end if
                 end do
              end do
           end do


           ! !two cells offset
           do k1=-2,2
              do j1=-2,2
                 do i1=-2,2
                    xtest(1,1)=pos(1,1)+i1*dx_loc
                    xtest(1,2)=pos(1,2)+j1*dx_loc
                    xtest(1,3)=pos(1,3)+k1*dx_loc
                    call get_cell_index(ind_cell,lev_cell,xtest,nlevelmax,1)
                    if(phi(ind_cell(1)) < 1.0000000001*phim)then
                       phim=phi(ind_cell(1))
                       !print*,"offset min by 2c",myid,jj,k1,j1,i1,phim,phi(ind_cell(1))
                       minmatch(jj)=0
                    end if
                 end do
              end do
           end do
           phi_min(jj)=phim
        end if
     end if
  end do

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
