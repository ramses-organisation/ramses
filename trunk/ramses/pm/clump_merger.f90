subroutine compute_clump_properties()
  use amr_commons
  use hydro_commons
  use pm_commons
  use clfind_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  !----------------------------------------------------------------------------                
  ! this subroutine performs one loop over all particles and collects the 
  ! relevant information from the cells where the particles sit. after a lot
  ! of mpi-communication, all necessary peak-batch properties can be computed
  !----------------------------------------------------------------------------

  ! variables used for the loop over all particles
  integer::igrid,jgrid,ipart,jpart,next_part,ig,ip,npart1,ilevel
  
  ! 
  integer::info,i,j,jj,peak_nr

  ! variables to be used with get_cell_indices (vector-sweeps)
  integer,dimension(1:nvector)::ind_grid,ind_cell,init_ind_cell,init_cell_lev,cell_lev
  integer,dimension(1:nvector)::ind_part,ind_grid_part
  real(dp),dimension(1:nvector,1:ndim)::pos,init_pos

  ! variables related to the size of a cell on a given level
  real(dp)::dx,dx_loc,scale,vol_loc
  real(dp),dimension(1:nlevelmax)::volume
  integer::nx_loc

  !peak-batch related arrays before sharing information with other cpus
  real(dp),dimension(0:npeaks_tot)::max_dens
  real(dp),dimension(1:npeaks_tot)::min_dens,av_dens,clump_mass,clump_vol
  real(dp),dimension(1:npeaks_tot,1:3)::center_of_mass,clump_momentum

  max_dens=0.; min_dens=1.d99
  clump_mass=0.; clump_vol=0.; clump_momentum=0.; center_of_mass=0.

#if NDIM==3

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

  !---------------------------------
  ! find number of particles on cpus
  !---------------------------------
  nparts=0
  !loop over all particles
  do ilevel=nlevelmax,levelmin,-1
     ! Loop over grids
     igrid=headl(myid,ilevel) 
     do jgrid=1,numbl(myid,ilevel) 
        nparts=nparts+numbp(igrid)
        igrid=next(igrid)   ! Go to next grid
     end do
  end do

#ifndef WITHOUTMPI     
  call MPI_ALLREDUCE(nparts,nparts_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI     
  nparts_tot=nparts
#endif
  if(verbose)write(*,*)'number of particles on ',myid,' after moving back = ',nparts,' of ',nparts_tot


  
  !set flag2 to zero for the reuse as clump number a cell belongs to
  flag2=0. 

     



  !---------------------------------------------------------------------------
  ! big loop over all parts to collect information from the cells
  !---------------------------------------------------------------------------
  !loop over all particles
  do ilevel=nlevelmax,levelmin,-1
     ig=0
     ip=0
     ! Loop over grids
     igrid=headl(myid,ilevel) 
     do jgrid=1,numbl(myid,ilevel) ! Number of grids in the level ilevel on process myid
        npart1=numbp(igrid)  ! Number of particles in the grid
        if(npart1>0)then
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle  <---- Very important !!!
              if (xp(ipart,1)>0.75.or.xp(ipart,1)<0.25)write(*,*)'alert_x',ipart,xp(ipart,1)
              if (xp(ipart,2)>0.75.or.xp(ipart,2)<0.25)write(*,*)'alert_y',ipart,xp(ipart,2)
              if (xp(ipart,3)>0.75.or.xp(ipart,3)<0.25)write(*,*)'alert_z',ipart,xp(ipart,3)
              if (vp(ipart,1)>0.75.or.vp(ipart,1)<0.25)write(*,*)'alert_vx',ipart,vp(ipart,1)
              if (vp(ipart,2)>0.75.or.vp(ipart,2)<0.25)write(*,*)'alert_vy',ipart,vp(ipart,2)
              if (vp(ipart,3)>0.75.or.vp(ipart,3)<0.25)write(*,*)'alert_vz',ipart,vp(ipart,3)
              next_part=nextp(ipart)
              if(ig==0)then
                 ig=1
                 ind_grid(ig)=igrid
              end if
              ip=ip+1
              ind_part(ip)=ipart
              ind_grid_part(ip)=ig   
              if(ip==nvector)then 
                 call get_cell_indices(init_ind_cell,init_cell_lev,ind_cell,cell_lev,ind_part,init_pos,pos,ip,nlevelmax)
                 do jj=1,nvector
                    if(cell_lev(jj) /= ilevel)write(*,*)'alert_problem in get cell index'
                    !find peak_nr
                    peak_nr=int(mp(ind_part(jj)))
                    if (peak_nr /=0 ) then
                       !find min density
                       if(uold(ind_cell(jj),1)<=min_dens(peak_nr))then 
                          min_dens(peak_nr)=uold(ind_cell(jj),1)
                       end if

                       !find max density
                       if(uold(ind_cell(jj),1)>=max_dens(peak_nr))then 
                          max_dens(peak_nr)=uold(ind_cell(jj),1)
                       end if

                       !find clump mass
                       clump_mass(peak_nr)=clump_mass(peak_nr)+ &
                            volume(cell_lev(jj))*uold(ind_cell(jj),1)


                       !center of mass velocity
                       do i=1,ndim
                          clump_momentum(peak_nr,i)=clump_momentum(peak_nr,i)+ &
                               uold(ind_cell(jj),(i+1))*volume(cell_lev(jj))
                       end do

                       !clump size (maybe take center of mass instead of peak as reference point)
                       do i=1,ndim
                          center_of_mass(peak_nr,i)=center_of_mass(peak_nr,i)+&
                               init_pos(jj,i)*volume(cell_lev(jj))*uold(ind_cell(jj),1)

                          !compute second order moments
                          do j=1,ndim
                             second_moments(peak_nr,i,j)=second_moments(peak_nr,i,j)+init_pos(jj,i)*init_pos(jj,j)*&
                                  volume(cell_lev(jj))*uold(ind_cell(jj),1)
                          end do
                       end do

                       !clump volume
                       clump_vol(peak_nr)=clump_vol(peak_nr)+volume(cell_lev(jj))

                       !now change the purpose of flag2 and make it a lookup table for the
                       !clump a certain cell belongs to
                       flag2(ind_cell(jj))=peak_nr
                    end if
                 end do
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if
        igrid=next(igrid)   ! Go to next grid
     end do
     if(ip>0)then 
        call get_cell_indices(init_ind_cell,init_cell_lev,ind_cell,cell_lev,ind_part,init_pos,pos,ip,nlevelmax)
        do jj=1,ip
           !find peak_nr
           peak_nr=int(mp(ind_part(jj)))
           if (peak_nr /=0 ) then
              !find min density
              if(uold(ind_cell(jj),1)<=min_dens(peak_nr))then 
                 min_dens(peak_nr)=uold(ind_cell(jj),1)
              end if

              !find max density
              if(uold(ind_cell(jj),1)>=max_dens(peak_nr))then 
                 max_dens(peak_nr)=uold(ind_cell(jj),1)
              end if

              !find clump mass
              clump_mass(peak_nr)=clump_mass(peak_nr)+ &
                   volume(cell_lev(jj))*uold(ind_cell(jj),1)

              !center of mass velocity
              do i=1,ndim
                 clump_momentum(peak_nr,i)=clump_momentum(peak_nr,i)+ &
                      uold(ind_cell(jj),(1+i))*volume(cell_lev(jj))
              end do

              !clump size (maybe take center of mass instead of peak as reference point)
              do i=1,ndim
                 center_of_mass(peak_nr,i)=center_of_mass(peak_nr,i)+&
                      init_pos(jj,i)*volume(cell_lev(jj))*uold(ind_cell(jj),1)

                 !compute second order moments
                 do j=1,ndim
                    second_moments(peak_nr,i,j)=second_moments(peak_nr,i,j)+init_pos(jj,i)*init_pos(jj,j)*&
                         volume(cell_lev(jj))*uold(ind_cell(jj),1)
                 end do
              end do

              !clump volume
              clump_vol(peak_nr)=clump_vol(peak_nr)+volume(cell_lev(jj))

              !now change the purpose of flag2 and make it a lookup table for the
              !clump a certain cell belongs to
              flag2(ind_cell(jj))=peak_nr
           end if
        end do
     end if
  end do
  !end loop over all particles





!---------------------------------------------------------------------------
! a lot of MPI communication to collect the results from the different cpu's
!---------------------------------------------------------------------------
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

#endif

end subroutine compute_clump_properties

!################################################################
!################################################################
!################################################################
!################################################################

subroutine compute_clump_properties_round2()
  use amr_commons
  use hydro_commons
  use pm_commons
  use clfind_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  !----------------------------------------------------------------------------                
  ! this subroutine performs another loop over all particles and collects the 
  ! information more information like the velocity of a cell realtive to the 
  ! center of mass
  !----------------------------------------------------------------------------

  ! variables used for the loop over all particles
  integer::igrid,jgrid,ipart,jpart,next_part,ig,ip,npart1,ilevel
  
  ! 
  integer::info,i,j,jj,peak_nr

  ! variables to be used with get_cell_indices (vector-sweeps)
  integer,dimension(1:nvector)::ind_grid,ind_cell,init_ind_cell,init_cell_lev,cell_lev
  integer,dimension(1:nvector)::ind_part,ind_grid_part
  real(dp),dimension(1:nvector,1:ndim)::pos
  real(dp),dimension(1:nvector,1:3)::init_pos

  ! variables related to the size of a cell on a given level
  real(dp)::dx,dx_loc,scale,vol_loc
  real(dp),dimension(1:nlevelmax)::volume
  integer::nx_loc

  !peak-batch related arrays before sharing information with other cpus
  real(dp),dimension(1:npeaks_tot)::e_kin_int,e_bind
  real(dp),dimension(1:npeaks_tot,1:3)::clump_size

  e_kin_int=0.; clump_size=0.; e_bind=0.

#if NDIM==3

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

  !---------------------------------
  ! find number of particles on cpus
  !---------------------------------
  nparts=0
  !loop over all particles
  do ilevel=nlevelmax,levelmin,-1
     ! Loop over grids
     igrid=headl(myid,ilevel) 
     do jgrid=1,numbl(myid,ilevel) 
        nparts=nparts+numbp(igrid)
        igrid=next(igrid)   ! Go to next grid
     end do
  end do

#ifndef WITHOUTMPI     
  call MPI_ALLREDUCE(nparts,nparts_tot,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI     
  nparts_tot=nparts
#endif
  if(verbose)write(*,*)'number of particles on ',myid,' after moving back = ',nparts,' of ',nparts_tot


  
  !set flag2 to zero for the reuse as clump number a cell belongs to
  flag2=0. 

     



  !---------------------------------------------------------------------------
  ! big loop over all parts to collect information from the cells
  !---------------------------------------------------------------------------
  !loop over all particles
  do ilevel=nlevelmax,levelmin,-1
     ig=0
     ip=0
     ! Loop over grids
     igrid=headl(myid,ilevel) 
     do jgrid=1,numbl(myid,ilevel) ! Number of grids in the level ilevel on process myid
        npart1=numbp(igrid)  ! Number of particles in the grid
        if(npart1>0)then
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle  <---- Very important !!!
              next_part=nextp(ipart)
              if(ig==0)then
                 ig=1
                 ind_grid(ig)=igrid
              end if
              ip=ip+1
              ind_part(ip)=ipart
              ind_grid_part(ip)=ig   
              if(ip==nvector)then 
                 call get_cell_indices(init_ind_cell,init_cell_lev,ind_cell,cell_lev,ind_part,init_pos,pos,ip,nlevelmax)
                 do jj=1,nvector
                    if(cell_lev(jj) /= ilevel)write(*,*)'alert_problem in get cell index'
                    !find peak_nr
                    peak_nr=int(mp(ind_part(jj)))
                    if (peak_nr /=0 ) then
                       
                       !internal kinetic energy
                       do i=1,ndim
                          e_kin_int(peak_nr)=e_kin_int(peak_nr)+ &
                               ((uold(ind_cell(jj),2)/uold(ind_cell(jj),1)&
                               -clump_momentum_tot(peak_nr,1)/clump_mass_tot(peak_nr))**2&
                               +(uold(ind_cell(jj),3)/uold(ind_cell(jj),1)-&
                               clump_momentum_tot(peak_nr,2)/clump_mass_tot(peak_nr))**2&
                               +(uold(ind_cell(jj),4)/uold(ind_cell(jj),1)-&
                               clump_momentum_tot(peak_nr,3)/clump_mass_tot(peak_nr))**2)*&
                               uold(ind_cell(jj),1)*volume(cell_lev(jj))*0.5
                       end do

                       !potential energy of the clump (crude approach)
                       e_bind(peak_nr)=e_bind(peak_nr)+&
                            clump_mass_tot(peak_nr)*uold(ind_cell(jj),1)*volume(cell_lev(jj))&
                            /(((init_pos(jj,1)-center_of_mass_tot(peak_nr,1))**2&
                            +(init_pos(jj,2)-center_of_mass_tot(peak_nr,2))**2&
                            +(init_pos(jj,3)-center_of_mass_tot(peak_nr,3))**2)**0.5)

                       !size relative to center of mass 
                       do i=1,ndim
                          clump_size(peak_nr,i)=clump_size(peak_nr,i)+ &
                               (init_pos(jj,i)-center_of_mass_tot(peak_nr,i))**2 * volume(cell_lev(jj))
                       end do


                       !now change the purpose of flag2 and make it a lookup table for the
                       !clump a certain cell belongs to
                       flag2(ind_cell(jj))=peak_nr
                    end if
                 end do
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if
        igrid=next(igrid)   ! Go to next grid
     end do
     if(ip>0)then 
        call get_cell_indices(init_ind_cell,init_cell_lev,ind_cell,cell_lev,ind_part,init_pos,pos,ip,nlevelmax)
        do jj=1,ip
           !find peak_nr
           peak_nr=int(mp(ind_part(jj)))
           if (peak_nr /=0 ) then
              
              !internal kinetic energy
              do i=1,ndim
                 e_kin_int(peak_nr)=e_kin_int(peak_nr)+ &
                      ((uold(ind_cell(jj),2)/uold(ind_cell(jj),1)&
                      -clump_momentum_tot(peak_nr,1)/clump_mass_tot(peak_nr))**2&
                      +(uold(ind_cell(jj),3)/uold(ind_cell(jj),1)-&
                      clump_momentum_tot(peak_nr,2)/clump_mass_tot(peak_nr))**2&
                      +(uold(ind_cell(jj),4)/uold(ind_cell(jj),1)-&
                      clump_momentum_tot(peak_nr,3)/clump_mass_tot(peak_nr))**2)*&
                      uold(ind_cell(jj),1)*volume(cell_lev(jj))*0.5
              end do

              !potential energy of the clump (crude approach)
              e_bind(peak_nr)=e_bind(peak_nr)+&
                   clump_mass_tot(peak_nr)*uold(ind_cell(jj),1)*volume(cell_lev(jj))/&
                   (((init_pos(jj,1)-center_of_mass_tot(peak_nr,1))**2+&
                   (init_pos(jj,2)-center_of_mass_tot(peak_nr,2))**2+&
                   (init_pos(jj,3)-center_of_mass_tot(peak_nr,3))**2)**0.5)              

              !size relative to center of mass 
              do i=1,ndim
                 clump_size(peak_nr,i)=clump_size(peak_nr,i)+ &
                      (init_pos(jj,i)-center_of_mass_tot(peak_nr,i))**2 * volume(cell_lev(jj))
              end do

              !now change the purpose of flag2 and make it a lookup table for the
              !clump a certain cell belongs to
              flag2(ind_cell(jj))=peak_nr
           end if
        end do
     end if
  end do
  !end loop over all particles





!---------------------------------------------------------------------------
! a lot of MPI communication to collect the results from the different cpu's
!---------------------------------------------------------------------------
#ifndef WITHOUTMPI     
  call MPI_ALLREDUCE(e_kin_int,e_kin_int_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI     
  e_kin_int_tot=e_kin_int
#endif
#ifndef WITHOUTMPI     
  call MPI_ALLREDUCE(e_bind,e_bind_tot,npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI     
  e_bind_tot=e_bind
#endif
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(clump_size,clump_size_tot,3*npeaks_tot,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  clump_size_tot=clump_size
#endif


#endif

end subroutine compute_clump_properties_round2

!################################################################
!################################################################
!################################################################
!################################################################


subroutine write_clump_properties(to_file,output_threshold)
  use amr_commons
  use hydro_commons
  use pm_commons
  use clfind_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  logical::to_file
  real(dp)::output_threshold
  !---------------------------------------------------------------------------
  ! this routine writes the clump properties to screen and to file
  !---------------------------------------------------------------------------
  integer::j,jj,ilevel,n_rel
  real(dp)::rel_mass
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v

  character(LEN=5)::nchar

#if NDIM==3

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  !sort clumps by peak density in ascending order
  call heapsort_index(max_dens_tot,sort_index,npeaks_tot)

  !print results in descending order to screen
  if(myid==1) then 
     rel_mass=0.
     n_rel=0
     write(*,*)'Cl_N  peak_x [uu] peak_y [uu] peak_z [uu] size_x [AU] size_y [AU]'//&
          & ' size_z [AU] |v|_CM [u.u.] rho- [g/cc] rho+ [g/cc] rho_av [g/cc] M_cl [M_sol] V_cl [AU^3] relevance E_bind/E_kin'
     do j=npeaks_tot,1,-1
        jj=sort_index(j)
        if (relevance_tot(jj) > output_threshold)then
           write(*,'(I6,3(X,F11.5),3(X,F11.2),X,F13.5,3(XE11.2E2),X,F13.5,XE11.2E2,X,F9.3,X,F7.3)'),jj&
                ,peak_pos_tot(jj,1),peak_pos_tot(jj,2),peak_pos_tot(jj,3)&
                ,(5.*clump_size_tot(jj,1)/clump_vol_tot(jj))**0.5*scale_l/1.496d13&
                ,(5.*clump_size_tot(jj,2)/clump_vol_tot(jj))**0.5*scale_l/1.496d13&
                ,(5.*clump_size_tot(jj,3)/clump_vol_tot(jj))**0.5*scale_l/1.496d13&
                ,(clump_momentum_tot(jj,1)**2+clump_momentum_tot(jj,2)**2+ &
                clump_momentum_tot(jj,3)**2)**0.5/clump_mass_tot(jj)&
                ,min_dens_tot(jj)*scale_d,max_dens_tot(jj)*scale_d&
                ,clump_mass_tot(jj)/clump_vol_tot(jj)*scale_d&
                ,clump_mass_tot(jj)*scale_d*scale_l**3/1.98892d33&
                ,clump_vol_tot(jj)*(scale_l/1.496d13)**3&
                ,relevance_tot(jj)&
                ,e_bind_tot(jj)/e_kin_int_tot(jj)
           
           rel_mass=rel_mass+clump_mass_tot(jj)*scale_d*scale_l**3/1.98892d33
           n_rel=n_rel+1
        end if
     end do
     write(*,*)'total mass above threshold =',tot_mass*scale_d*scale_l**3/1.98892d33
     write(*,*)'total mass in ',n_rel,' listed clumps =',rel_mass


     !write clump_info to file
     if (to_file .eqv. .true.) then
        call title(ifout-1,nchar)
        open(unit=20,file=TRIM('output_'//TRIM(nchar)//'/clump_info.txt'),form='formatted')
        open(unit=21,file=TRIM('output_'//TRIM(nchar)//'/clump_masses.txt'),form='formatted')
        write(20,*),' Cl_N peak_x [uu] peak_y [uu] peak_z [uu] size_x [AU] size_y [AU] size_z [AU]'//& 
             & '|v|_CM [u.u.] rho- [g/cc] rho+ [g/cc] rho_av [g/cc] M_cl [M_sol] V_cl [AU^3] relevance E_bind/E_kin'
        do j=npeaks_tot,1,-1
           jj=sort_index(j)
           if (relevance_tot(jj) > output_threshold)then
              write(20,'(I6,3(X,F11.5),3(X,F11.2),3(XE11.2E2),X,F13.5,XE11.2E2,X,F9.3,X,F7.3)'),jj&
                   ,peak_pos_tot(jj,1),peak_pos_tot(jj,2),peak_pos_tot(jj,3)&
                   ,(5.*clump_size_tot(jj,1)/clump_vol_tot(jj))**0.5*scale_l/1.496d13&
                   ,(5.*clump_size_tot(jj,2)/clump_vol_tot(jj))**0.5*scale_l/1.496d13&
                   ,(5.*clump_size_tot(jj,3)/clump_vol_tot(jj))**0.5*scale_l/1.496d13&
                   ,(clump_momentum_tot(jj,1)**2+clump_momentum_tot(jj,2)**2+ &
                   clump_momentum_tot(jj,3)**2)**0.5/clump_mass_tot(jj)&
                   ,min_dens_tot(jj)*scale_d,max_dens_tot(jj)*scale_d&
                   ,clump_mass_tot(jj)/clump_vol_tot(jj)*scale_d&
                   ,clump_mass_tot(jj)*scale_d*scale_l**3/1.98892d33&
                   ,clump_vol_tot(jj)*(scale_l/1.496d13)**3&
                   ,relevance_tot(jj)&
                   ,e_bind_tot(jj)/e_kin_int_tot(jj)
           end if
        end do
        write(21,*)n_rel
        do j=npeaks_tot,1,-1
           jj=sort_index(j)
           if (relevance_tot(jj)>output_threshold)write(21,*)clump_mass_tot(jj)*scale_d*scale_l**3/1.98892d33
        end do
        write(20,'(A,F9.3)')'total mass above threshold =',tot_mass*scale_d*scale_l**3/1.98892d33
        write(20,'(A,I6,A,F9.3)')'total mass in',n_rel,' listed clumps =',rel_mass
        close(20)
        close(21)
     end if
  end if

#endif

end subroutine write_clump_properties

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine saddlepoint_search()
  use amr_commons
  use hydro_commons
  use pm_commons
  use clfind_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  !---------------------------------------------------------------------------
  ! subroutine which creates a npeaks_tot**2 sized array of saddlepoint densities
  !---------------------------------------------------------------------------

  ! variables for particle loop
  integer::igrid,jgrid,ipart,jpart,next_part,ig,ip,npart1,ilevel
  ! dummy variables
  integer::jj,i,j,kk,info
  ! variables to be used with get_cell_indices (vector-sweeps)
  integer,dimension(1:nvector),save::ind_grid,ind_cell,init_ind_cell,init_cell_lev,cell_lev,clump_nr
  integer,dimension(1:nvector),save::ind_part,ind_grid_part
  real(dp),dimension(1:nvector,1:ndim)::pos
  real(dp),dimension(1:nvector,1:3)::init_pos
  ! saddle point array for 1 cpu
  real(kind=8),dimension(1:npeaks_tot,1:npeaks_tot)::saddle_dens
  saddle_dens=0.



  !---------------------------------------------------------------------------
  ! big loop over all parts, pass the information of nvector parts to 
  ! subroutine find_best_neighbor
  !---------------------------------------------------------------------------
  do ilevel=nlevelmax,levelmin,-1
     ig=0
     ip=0
     ! Loop over grids
     igrid=headl(myid,ilevel) 
     do jgrid=1,numbl(myid,ilevel) ! Number of grids in the level ilevel on process myid
        npart1=numbp(igrid)  ! Number of particles in the grid
        if(npart1>0)then
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle  <---- Very important !!!
              next_part=nextp(ipart)
              if(ig==0)then
                 ig=1
                 ind_grid(ig)=igrid
              end if
              ip=ip+1
              ind_part(ip)=ipart
              ind_grid_part(ip)=ig   
              if(ip==nvector)then 
                 call get_cell_indices(init_ind_cell,init_cell_lev,ind_cell,cell_lev,ind_part,init_pos,pos,ip,nlevelmax)
                 do jj=1,nvector
                    clump_nr(jj)=int(mp(ind_part(jj)))
                    ind_grid(jj)=mod((ind_cell(jj)-ncoarse),ngridmax)
                 end do
                 call find_best_neighbor(ind_grid,clump_nr,init_pos,ip,ilevel,saddle_dens)
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if
        igrid=next(igrid)   ! Go to next grid
     end do
     if(ip>0)then 
        call get_cell_indices(init_ind_cell,init_cell_lev,ind_cell,cell_lev,ind_part,init_pos,pos,ip,nlevelmax)
        do jj=1,ip
           clump_nr(jj)=int(mp(ind_part(jj)))
           ind_grid(jj)=mod((ind_cell(jj)-ncoarse),ngridmax)
        end do
        call find_best_neighbor(ind_grid,clump_nr,init_pos,ip,ilevel,saddle_dens)
     end if
  end do
  !end loop over all particles

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(saddle_dens,saddle_dens_tot,(npeaks_tot**2),MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  saddle_dens_tot=saddle_dens
#endif

  !symmetrize saddle_point array
  do i=1,npeaks_tot
     do j=1,i
        saddle_dens_tot(i,j)=0.5*(saddle_dens_tot(i,j)+saddle_dens_tot(j,i))
        saddle_dens_tot(j,i)=saddle_dens_tot(i,j)
     end do
  end do

  !compute saddle_max value and relevance
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
  ! This routine computes the force on each particle by
  ! inverse CIC and computes new positions for all particles.
  ! If particle sits entirely in fine level, then CIC is performed
  ! at level ilevel. Otherwise, it is performed at level ilevel-1.
  ! This routine is called by move_fine.
  !------------------------------------------------------------
  logical::error
  integer::j,ind,idim,nx_loc
  integer::i1,j1,k1,i2,j2,k2
  real(dp)::dx,dx_loc,scale,vol_loc
  integer::i1min,i1max,j1min,j1max,k1min,k1max
  integer::i2min,i2max,j2min,j2max,k2min,k2max
  integer::i3min,i3max,j3min,j3max,k3min,k3max
  ! Grid-based arrays

  real(dp),dimension(1:nvector,1:ndim),save::x0
  ! Particle-based arrays

  real(dp),dimension(1:nvector,1:ndim)::x,xtest,pos
  integer ,dimension(1:nvector,1:ndim),save::ig,id


  real(dp),dimension(1:nvector,1:ndim,1:twotondim),save::xpart

  real(dp),dimension(1:3)::skip_loc


  ! Meshspacing in that level
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

     ! Rescale particle position at level ilevel-1
     do idim=1,ndim
        do j=1,np
           x(j,idim)=pos(j,idim)/scale+skip_loc(idim)
        end do
     end do
     do idim=1,ndim
        do j=1,np
           x(j,idim)=x(j,idim)-x0(j,idim)
        end do
     end do
     do idim=1,ndim
        do j=1,np
           x(j,idim)=x(j,idim)/dx
        end do
     end do
     do idim=1,ndim
        do j=1,np
           x(j,idim)=x(j,idim)/2.0D0
        end do
     end do

     ! Check for illegal moves
     error=.false.
     do idim=1,ndim
        do j=1,np
           if(x(j,idim)<0.5D0.or.x(j,idim)>2.5D0)error=.true.
        end do
     end do
     if(error)then
        write(*,*)'problem in neighbor-check'
        do idim=1,ndim
           do j=1,np
              if(x(j,idim)<0.5D0.or.x(j,idim)>2.5D0)then
                 write(*,*)x(j,1:ndim)
              endif
           end do
        end do
        stop
     end if

     ! Compute parent cell position
#if NDIM==1
     do j=1,np
        xpart(j,1,1)=0.5+ig(j,1)
        xpart(j,1,2)=0.5+id(j,1)
     end do
#endif
#if NDIM==2
     do j=1,np
        ! Particle 1
        xpart(j,1,1)=0.5+ig(j,1)
        xpart(j,2,1)=0.5+ig(j,2)
        ! Particle 2
        xpart(j,1,2)=0.5+id(j,1)
        xpart(j,2,2)=0.5+ig(j,2)
        ! Particle 3
        xpart(j,1,3)=0.5+ig(j,1)
        xpart(j,2,3)=0.5+id(j,2)
        ! Particle 4
        xpart(j,1,4)=0.5+id(j,1)
        xpart(j,2,4)=0.5+id(j,2)
     end do
#endif
#if NDIM==3
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
#endif

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
#if NDIM>1
              xtest(j,2)=pos(j,2)+(j1-1)*dx_loc
#endif     
#if NDIM>2
              xtest(j,3)=pos(j,3)+(k1-1)*dx_loc
#endif     
           end do

           call get_cell_index(ind_cell,cell_levl,xtest,ilevel,np)
           do j=1,np
              if(son(ind_cell(j))==0 .and. flag2(ind_cell(j))/=0)then
                 if(clump_nr(j) /= flag2(ind_cell(j)).and.uold(ind_cell(j),1)>saddle_dens(clump_nr(j),flag2(ind_cell(j))))then   
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
#if NDIM>1
                 xtest(j,2)=pos(j,2)+(j2-1.5)*dx_loc/2.0
#endif     
#if NDIM>2
                 xtest(j,3)=pos(j,3)+(k2-1.5)*dx_loc/2.0
#endif     
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
subroutine merge_clumps()
  use amr_commons
  use pm_commons
  use poisson_commons
  use hydro_commons, ONLY: uold
  use clfind_commons
  implicit none

  integer::j,jj,i,ii
  integer::igrid,jgrid,ipart,jpart,next_part,npart1,info,ilevel,merge_to
  real(dp)::max_val
  logical::merging



  if (verbose)write(*,*)'Now merging clumps'


  merging=.true.
  do while(merging .eqv. .true.)
     merging=.false.

     !sort clumps by peak density ind ascending order
     call heapsort_index(relevance_tot,sort_index,npeaks_tot)

     i=1
     ii=sort_index(i)
     
     !find the first clump which has not yet been merged to another one
     do while (relevance_tot(ii) < 1.d-80 .and. i<npeaks_tot)
        i=i+1
        ii=sort_index(i)
     end do
     
     !if the relevance is below the threshold -> merge
     !if the mass is below the threshold -> merge
     !if the relevance is above the threshold and the mass is above the threshold -> done merging
     !if there is no relevant clump -> done merging
     if ((relevance_tot(ii)<relevance_threshold .or. clump_mass_tot(ii) < mass_threshold) .and. relevance_tot(ii) > 1.d-80) then
        merging=.true.

        !go through the ii-th line in the saddle point array to find the neighbor to merge to
        merge_to=0; max_val=0.
        do j=1,npeaks_tot
           if (saddle_dens_tot(ii,j)>max_val)then
              merge_to=j
              max_val=saddle_dens_tot(ii,j)
           end if
        end do

        !loop over all particles
        do ilevel=nlevelmax,levelmin,-1
           ! Loop over grids
           igrid=headl(myid,ilevel) 
           do jgrid=1,numbl(myid,ilevel) ! Number of grids in the level ilevel on process myid
              npart1=numbp(igrid)  ! Number of particles in the grid
              if(npart1>0)then
                 ipart=headp(igrid)
                 ! Loop over particles
                 do jpart=1,npart1
                    ! Save next particle  <---- Very important !!!
                    next_part=nextp(ipart)
                    if (int(mp(ipart))==ii) then
                       mp(ipart)=merge_to*1.0
                    end if
                    ipart=next_part  ! Go to next particle
                 end do
                 ! End loop over particles
              end if
              igrid=next(igrid)   ! Go to next grid
           end do
        end do
        !end loop over all particles
        if(myid==1)write(*,*)'clump ',ii,'merged to ',merge_to

        !update clump properties
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
           clump_vol_tot(merge_to)=clump_vol_tot(ii)+clump_vol_tot(merge_to)
           max_dens_tot(merge_to)=max(max_dens_tot(merge_to),max_dens_tot(ii))
           min_dens_tot(merge_to)=min(min_dens_tot(merge_to),min_dens_tot(ii))
           clump_mass_tot(merge_to)=clump_mass_tot(merge_to)+clump_mass_tot(ii)
        end if
        clump_vol_tot(ii)=0.
        max_dens_tot(ii)=0.
        min_dens_tot(ii)=0.
        clump_mass_tot(ii)=0.


        !update saddle point array
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

        !update saddle_max value
        if (merge_to>0)then
           saddle_max_tot(merge_to)=0
           do j=1,npeaks_tot
              if (saddle_dens_tot(merge_to,j)>saddle_max_tot(merge_to))then
                 saddle_max_tot(merge_to)=saddle_dens_tot(merge_to,j)
              end if
           end do
        end if

        !update relevance of clumps
        if (merge_to>0)then
           if (saddle_max_tot(merge_to)>1.d-40)then
              relevance_tot(merge_to)=max_dens_tot(merge_to)/saddle_max_tot(merge_to)
           else 
              relevance_tot(merge_to)=max_dens_tot(merge_to)/min_dens_tot(merge_to)
           end if
        end if
        relevance_tot(ii)=0.
     end if
     
     !if done merging, look for clumps which are too light and divide their relevance by factor of 10
     if (merging .eqv. .false.) then
        do j=1,npeaks_tot
           if (clump_mass_tot(j) < mass_threshold .and. clump_mass_tot(j)>1.d-80) then
              relevance_tot(j)=relevance_tot(j)/10.
              merging=.true.
           end if
        end do
     end if
     
  end do
  if (verbose)write(*,*)'Done merging clumps'  
end subroutine merge_clumps

!################################################################
!################################################################
!################################################################
!################################################################
subroutine write_peak_map
  use amr_commons
  use hydro_commons
  use pm_commons
  use clfind_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  !----------------------------------------------------------------------------                
  ! One loop over all the particles creating a textfile for each cpu
  ! -number of parts on cpu
  ! -xpos for all parts
  ! -ypos for all parts
  ! -zpos for all parts
  ! -clump number for all parts
  !----------------------------------------------------------------------------

  integer::igrid,jgrid,ipart,jpart,next_part,ig,ip,npart1,info,ilevel
  integer::jj,i,kk 
  integer,dimension(1:nvector)::ind_grid,ind_cell,init_ind_cell,init_cell_lev,cell_lev
  integer,dimension(1:nvector)::ind_part,ind_grid_part
  real(dp),dimension(1:nvector,1:ndim)::pos
  real(dp),dimension(1:nvector,1:3)::init_pos
  character(LEN=5)::myidstring,nchar

#if NDIM==3

  call title(ifout-1,nchar)
  call title(myid,myidstring)
  open(unit=20,file=TRIM('output_'//TRIM(nchar)//'/clump_map.csv'//myidstring),form='formatted')
  !write(20,*)nparts

  !loop over all particles
  do ilevel=nlevelmax,levelmin,-1
     ig=0
     ip=0
     ! Loop over grids
     igrid=headl(myid,ilevel) 
     do jgrid=1,numbl(myid,ilevel) 
        npart1=numbp(igrid)  ! Number of particles in the grid
        if(npart1>0)then
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle  <---- Very important !!!
              next_part=nextp(ipart)
              if(ig==0)then
                 ig=1
                 ind_grid(ig)=igrid
              end if
              ip=ip+1
              ind_part(ip)=ipart
              ind_grid_part(ip)=ig   
              if(ip==nvector)then 
                 call get_cell_indices(init_ind_cell,init_cell_lev,ind_cell,cell_lev,ind_part,init_pos,pos,ip,nlevelmax)
                 do jj=1,nvector
                    if (int(mp(ind_part(jj)))>0)then
                       write(20,'(F11.8,A,F11.8,A,F11.8,A,I8)')init_pos(jj,1),',',init_pos(jj,2),','&
                            ,init_pos(jj,3),',',int(mp(ind_part(jj)))
                       end if
                 end do
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if
        igrid=next(igrid)   ! Go to next grid
     end do
     if(ip>0)then 
        call get_cell_indices(init_ind_cell,init_cell_lev,ind_cell,cell_lev,ind_part,init_pos,pos,ip,nlevelmax)
        do jj=1,ip
           if (int(mp(ind_part(jj)))>0)then
              write(20,'(F11.8,A,F11.8,A,F11.8,A,I8)')init_pos(jj,1),',',init_pos(jj,2),','&
                   ,init_pos(jj,3),',',int(mp(ind_part(jj)))
           end if
        end do
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
subroutine allocate_peak_batch_arrays
  use amr_commons
  use hydro_commons
  use pm_commons
  use clfind_commons
  implicit none


  !allocate peak-batch_properties
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
  allocate(e_bind_tot(npeaks_tot))
  
  !initialize all peak based arrays                 
  saddle_max_tot=0.
  relevance_tot=1.
  clump_size_tot=0.
  min_dens_tot=1.d99
  max_dens_tot=0.
  av_dens_tot=0.
  clump_mass_tot=0.
  clump_vol_tot=0.
  center_of_mass_tot=0.
  second_moments=0.; second_moments_tot=0.
  saddle_dens_tot=0.
  clump_momentum_tot=0.
  e_kin_int_tot=0.
  e_bind_tot=0.

end subroutine allocate_peak_batch_arrays

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
  deallocate(e_kin_int_tot)
  deallocate(e_bind_tot)
  
  deallocate(saddle_dens_tot)

  deallocate(peak_pos_tot)

end subroutine deallocate_all
