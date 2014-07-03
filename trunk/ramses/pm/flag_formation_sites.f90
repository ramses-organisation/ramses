subroutine flag_formation_sites
  use amr_commons
  use pm_commons
  use clfind_commons
  use hydro_commons, only:uold
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  !=============================================================================
  ! This routine flags (flag 2 = 1 )  the cells where a sink is going to be formed
  !=============================================================================

  real(dp),dimension(1:nvector,1:3)::pos
  integer,dimension(1:nvector)::cell_index,cell_levl,cc
  integer::j,jj,i,nx_loc
  integer::flag_form,flag_form_tot,info
  logical::ok
  real(dp)::dx,dx_min,dist,scale,tff,acc_r
  real(dp)::fourpi,threepi2
  real(dp),dimension(1:npeaks)::peakd
  integer,dimension(1:npeaks)::ind_sort
  !gridspacing and physical scales
  dx=0.5D0**nlevelmax
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_min=dx*scale



  ! loop over sinks and mark all clumps which are already occupied by a sink
  allocate(occupied(1:npeaks_max))
  occupied=0
  pos=0.0
  if(myid==1 .and. clinfo)write(*,*)'looping over ',nsink,' sinks and marking their clumps'

  if (smbh)then 
     !block clumps that contain a sink for formation
     do j=1,nsink
        pos(1,1:3)=xsink(j,1:3)
        call cmp_cpumap(pos,cc,1)
        if (cc(1) .eq. myid)then
           call get_cell_index(cell_index,cell_levl,pos,nlevelmax,1)
           if (flag2(cell_index(1))>0)then
              occupied(flag2(cell_index(1)))=1
              if(clinfo)write(*,*)'CPU # ',myid,'blocked clump # ',flag2(cell_index(1)),' for sink production because of sink # ',idsink(j)
           end if
        end if
     end do
  else
     !block peaks that are closer than R_accretion from existing sinks
     do j=1,nsink
        do i=1,npeaks
           dist=(xsink(j,1)-peak_pos(i,1))**2+&
                (xsink(j,2)-peak_pos(i,2))**2+&
                (xsink(j,3)-peak_pos(i,3))**2
           if (dist<(ir_cloud*dx_min)**2)then
              occupied(i)=1
              if(myid==1 .and. clinfo)write(*,*)'blocked clump # ',i,' for sink production because of sink # ',idsink(j)
           end if
        end do
     end do
  end if


#ifndef WITHOUTMPI
  call virtual_peak_int(occupied,'max')
#endif


  !------------------------------------------------------------------------------
  ! determine whether a peak patch is allowed to form a new sink.
  ! if a new sink has to be created, flag2 is set to the clump number at the peak position
  ! -> criteria to be chosen depend on the physics
  ! -> this routine can be patched
  !------------------------------------------------------------------------------     
  pos=0.0
  flag2=0

  !sort clumps by peak density in ascending order
  do i=1,npeaks
     peakd(i)=max_dens(i)
     ind_sort(i)=i
  end do
  call quick_sort_dp(peakd,ind_sort,npeaks)
  do j=npeaks,1,-1
     jj=ind_sort(j)
     ok=.true.
     ok=ok.and.relevance(jj)>relevance_threshold
     ok=ok.and.occupied(jj)==0
     ok=ok.and.max_dens(jj)>d_sink
     ok=ok.and.contracting(jj)
     ok=ok.and.Icl_dd(jj)<0.
     if (ok)then
        pos(1,1:3)=peak_pos(jj,1:3)
        call cmp_cpumap(pos,cc,1)
        if (cc(1) .eq. myid)then
           call get_cell_index(cell_index,cell_levl,pos,nlevelmax,1)
           flag2(cell_index(1))=jj
           write(*,*)'cpu ',myid,' produces a new sink for clump number ',jj
        end if
     end if
  end do

  deallocate(occupied)


end subroutine flag_formation_sites

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

  integer::ipart,ilevel,info,i,peak_nr,global_peak_id,j,ii,jj
  integer::grid,nx_loc,ix,iy,iz,ind
  real(dp)::d,vol,M,ekk,phi_rel,de,c_sound,d0,v_bulk2,p,cty
  real(dp)::t_larson1,cont_speed=0.
  real(dp)::dx,dx_loc,scale,vol_loc,abs_err,A1=0.,A2=0.,A3=0.
  real(dp),dimension(1:nlevelmax)::volume
  real(dp),dimension(1:3)::vd,xcell,xpeak,v_cl,rrel,vrel,frel,skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:3,1:3)::eigenv,a
  real(dp),dimension(1:npeaks,1:3)::contractions

  call surface_pressure(ntest)
  
  !initialize arrays
  e_kin_int=0.d0; clump_size=0.d0; e_thermal=0.d0
  grav_term=0.d0; Icl_d=0.d0; Icl=0.; Icl_dd=0.
  Icl_3by3=0.;  Icl_d_3by3=0.
  contracting=.false.
  cty=1.

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
  
#ifndef WITHOUTMPI
  call boundary_peak_dp(max_dens)
  call boundary_peak_dp(clump_mass)
  call boundary_peak_dp(peak_pos)
#endif
  !---------------------------------------------------------------------------
  ! loop over all test particles to collect information from the cells
  !---------------------------------------------------------------------------
  do ipart=1,ntest     
     global_peak_id=flag2(icellp(ipart)) 
     if (global_peak_id /=0 ) then
        call get_local_peak_id(global_peak_id,peak_nr)
        
        ! Cell coordinates
        ind=(icellp(ipart)-ncoarse-1)/ngridmax+1 ! cell position
        grid=icellp(ipart)-ncoarse-(ind-1)*ngridmax ! grid index
        dx=0.5D0**levp(ipart)
        xcell(1:ndim)=(xg(grid,1:ndim)+xc(ind,1:ndim)*dx-skip_loc(1:ndim))*scale

        ! gas density and energydensity
        d=uold(icellp(ipart),1)
        de=uold(icellp(ipart),ndim+2) 
        do i=1,ndim
           vd(i)=uold(icellp(ipart),i+1)
           xpeak(i)=peak_pos(peak_nr,i)
        end do
     
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

        ! thermal energy
        ekk=0.
        do i=1,3 
           ekk=ekk+0.5*vd(i)**2/d                          
        end do
        p=(de-ekk)*(gamma-1)
        e_thermal(peak_nr)=e_thermal(peak_nr)+1.5*vol*p

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
  call virtual_peak_dp(e_kin_int,'sum')
  call virtual_peak_dp(clump_size,'sum')
  call virtual_peak_dp(grav_term,'sum')
  call virtual_peak_dp(e_thermal,'sum')
  call virtual_peak_dp(Icl,'sum')
  call virtual_peak_dp(Icl_3by3,'sum')
  call virtual_peak_dp(Icl_d_3by3,'sum')
#endif

  !second time derivative of I
  Icl_dd(1:npeaks)=2.*(grav_term(1:npeaks)-Psurf(1:npeaks)+2*e_kin_int(1:npeaks)+2*e_thermal(1:npeaks))

!  all_bound=.true.

  do j=npeaks,1,-1
     if (relevance(j)>0.)then

        !compute eigenvalues and eigenvectors of Icl_d_3by3
        a=Icl_3by3(j,1:3,1:3)
        abs_err=1.d-8*Icl(j)**2+1.d-40
!        call jacobi(a,eigenv,abs_err)
        A1=a(1,1); A2=a(2,2); A3=a(3,3)

        !compute the contractions along the eigenvectors of Icl
        contractions(j,1:3)=0._dp
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
        
        clump_check(j)=(-1.*grav_term(j)+Psurf(j))/(tiny(0.d0)+2*e_kin_int(j)+2*e_thermal(j))
        
        !update the all_bound property
!        all_bound=all_bound.and.(isodens_check(j)>1.)

     endif
  end do

  !write to the log file some information that could be of interest for debugging etc.
  if(myid==1 .and. clinfo .and. .not. smbh .and. sink)then 
     write(*,'(135A)')'==========================================================================================='
     write(*,'(135A)')'Cl_N     t1[y]      t2[y]      t3[y] |I_d|/I_dd[y] tidal_Fg   Psurf      e_kin      e_therm'
     write(*,'(135A)')'==========================================================================================='
     do j=npeaks,1,-1
        if (relevance(j)>0.)then
           write(*,'(I4,2X,8(E9.2E2,3X))'),j&
                ,A1/(contractions(j,1)+tiny(0.d0))*cty,A2/(contractions(j,2)+tiny(0.d0))*cty,A3/(contractions(j,3)+tiny(0.d0))*cty&
                ,abs(Icl_d(j))/Icl_dd(j)*cty&
                ,grav_term(j),-1.*Psurf(j)&
                ,e_kin_int(j),e_thermal(j)
        end if
     end do
     write(*,'(135A)')'==========================================================================================='
  end if
     
end subroutine compute_clump_properties_round2
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
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
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  integer ::ix,iy,iz,grid,peak_nr,glob_peak_nr

  real(dp),dimension(1:3)::skip_loc,xcell
  real(dp),dimension(1:twotondim,1:3)::xc

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  
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

#ifndef WITHOUTMPI
  call boundary_peak_dp(peak_pos)
#endif

  !update flag 2
  do ipart=1,ntest
     glob_peak_nr=flag2(icellp(ipart))
     if (glob_peak_nr /=0 ) then
        call get_local_peak_id(glob_peak_nr,peak_nr)
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

  !return if already diagonal "enough"
  if (b2 <= err2) then
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
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine surface_int(ind_cell,np,ilevel)
  use amr_commons
  use clfind_commons, ONLY: center_of_mass,Psurf
  use hydro_commons, ONLY: uold,gamma
  implicit none
  integer::np,ilevel
  integer,dimension(1:nvector)::ind_grid,ind_cell

  !------------------------------------------------------------
  ! This routine constructs all neighboring leaf cells that 
  ! have a common cell surface at levels 
  ! ilevel-1, ilevel, ilevel+1. Then, it computes the pressure
  ! pressure onto these surfaces and integrates over the surface
  ! of the clumps.
  !------------------------------------------------------------

  integer::j,ind,nx_loc,i2,j2,k2,ix,iy,iz,idim,jdim,i3,j3,k3
  real(dp)::dx,dx_loc,scale,vol_loc
  integer ,dimension(1:nvector)::cell_index,cell_levl,clump_nr,indv,neigh_cl
  integer ,dimension(1:nvector)::loc_clump_nr
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim)::xtest,r
  real(dp),dimension(1:nvector)::ekk_cell,ekk_neigh,P_cell,P_neigh,r_dot_n
  real(dp),dimension(1:3)::skip_loc,n
  logical ,dimension(1:nvector)::ok

#if NDIM==3

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

  ! Cells center position relative to grid center position
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     xc(ind,1)=(dble(ix)-0.5D0)*dx
     xc(ind,2)=(dble(iy)-0.5D0)*dx
     xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do
  
  ekk_cell=0.; P_neigh=0; P_cell=0
  ! some preliminary action...
  do j=1,np
     indv(j)=(ind_cell(j)-ncoarse-1)/ngridmax+1 ! cell position in grid
     ind_grid(j)=ind_cell(j)-ncoarse-(indv(j)-1)*ngridmax ! grid index
     clump_nr(j)=flag2(ind_cell(j)) ! save clump number
     do jdim=1,ndim
        ekk_cell(j)=ekk_cell(j)+0.5*uold(ind_cell(j),jdim+1)**2
     end do
     ekk_cell(j)=ekk_cell(j)/uold(ind_cell(j),1)
     P_cell(j)=(gamma-1.0)*(uold(ind_cell(j),ndim+2)-ekk_cell(j))
  end do

  do j=1,np
     if (clump_nr(j) .ne. 0)call get_local_peak_id(clump_nr(j),loc_clump_nr(j))
  end do


  
  !================================
  ! generate neighbors at level ilevel (and ilevel -1)
  !================================
  ! Generate 3x3 neighboring cells at level ilevel
  do k2=0,2
     do j2=0,2
        do i2=0,2
           if((k2-1.)**2+(j2-1.)**2+(i2-1.)**2==1)then !check whether common face exists 
              
              !construct outward facing normal vector
              n=0.
              if (k2==0)n(3)=-1.
              if (k2==2)n(3)=1.
              if (j2==0)n(2)=-1.
              if (j2==2)n(2)=1.
              if (i2==0)n(1)=-1.
              if (i2==2)n(1)=1.
              if (n(1)**2+n(2)**2+n(3)**2/=1)print*,'n has wrong lenght'
              
              
              r=0.
              do j=1,np                 
                 xtest(j,1)=(xg(ind_grid(j),1)+xc(indv(j),1)-skip_loc(1))*scale+(i2-1)*dx_loc
                 xtest(j,2)=(xg(ind_grid(j),2)+xc(indv(j),2)-skip_loc(2))*scale+(j2-1)*dx_loc
                 xtest(j,3)=(xg(ind_grid(j),3)+xc(indv(j),3)-skip_loc(3))*scale+(k2-1)*dx_loc

                 if (clump_nr(j)>0)then                    
                    r(j,1)=(xg(ind_grid(j),1)+xc(indv(j),1)-skip_loc(1))*scale+(i2-1)*dx_loc*0.5&
                         -center_of_mass(clump_nr(j),1)
                    r(j,2)=(xg(ind_grid(j),2)+xc(indv(j),2)-skip_loc(2))*scale+(j2-1)*dx_loc*0.5&
                         -center_of_mass(clump_nr(j),2)
                    r(j,3)=(xg(ind_grid(j),3)+xc(indv(j),3)-skip_loc(3))*scale+(k2-1)*dx_loc*0.5&
                         -center_of_mass(clump_nr(j),3)
                 endif                 
              end do
              
              call get_cell_index(cell_index,cell_levl,xtest,ilevel,np)
              do j=1,np           
                 ok(j)=(son(cell_index(j))==0)
              end do
              do j=1,np
                 neigh_cl(j)=flag2(cell_index(j))!nuber of clump the neighboring cell is in 
                 ok(j)=ok(j).and. neigh_cl(j)/=clump_nr(j) !neighboring cell is in another clump
                 ok(j)=ok(j).and. 0/=clump_nr(j) !clump number is not zero
              end do
              
              r_dot_n=0.
              do j=1,np
                 do idim=1,3
                    r_dot_n(j)=r_dot_n(j)+n(idim)*r(j,idim)
                 end do
              end do
              
              ekk_neigh=0.
              do j=1,np
                 if (ok(j))then 
                    do jdim=1,ndim
                       ekk_neigh(j)=ekk_neigh(j)+0.5*uold(cell_index(j),jdim+1)**2
                    end do
                    ekk_neigh(j)=ekk_neigh(j)/uold(cell_index(j),1)
                    P_neigh(j)=(gamma-1.0)*(uold(cell_index(j),ndim+2)-ekk_neigh(j))
                    Psurf(loc_clump_nr(j))=Psurf(loc_clump_nr(j))+r_dot_n(j)*dx_loc**2*0.5*(P_neigh(j)+P_cell(j))
                 endif
              end do
           endif
        end do
     end do
  end do
  

  !===================================
  ! generate neighbors at level ilevel+1
  !====================================  
  if(ilevel<nlevelmax)then  
     ! Generate 4x4x4 neighboring cells at level ilevel+1 
     do k3=0,3
        do j3=0,3
           do i3=0,3
              if((k3-1.5)**2+(j3-1.5)**2+(i3-1.5)**2==2.75)then !check whether common face exists

                 n=0.
                 if (k3==0)n(3)=-1. 
                 if (k3==3)n(3)=1.
                 if (j3==0)n(2)=-1. 
                 if (j3==3)n(2)=1.
                 if (i3==0)n(1)=-1. 
                 if (i3==3)n(1)=1.
                 if (n(1)**2+n(2)**2+n(3)**2/=1)print*,'n has wrong lenght'

                 r=0.
                 do j=1,np 

                    xtest(j,1)=(xg(ind_grid(j),1)+xc(indv(j),1)-skip_loc(1))*scale+(i3-1.5)*dx_loc/2.0
                    xtest(j,2)=(xg(ind_grid(j),2)+xc(indv(j),2)-skip_loc(2))*scale+(j3-1.5)*dx_loc/2.0
                    xtest(j,3)=(xg(ind_grid(j),3)+xc(indv(j),3)-skip_loc(3))*scale+(k3-1.5)*dx_loc/2.0
                    
                    if (clump_nr(j)>0)then                       
                       r(j,1)=(xg(ind_grid(j),1)+xc(indv(j),1)-skip_loc(1))*scale+(i3-1.5)*dx_loc/2.0*0.5&
                            -center_of_mass(clump_nr(j),1)
                       r(j,2)=(xg(ind_grid(j),2)+xc(indv(j),2)-skip_loc(2))*scale+(j3-1.5)*dx_loc/2.0*0.5&
                            -center_of_mass(clump_nr(j),2)
                       r(j,3)=(xg(ind_grid(j),3)+xc(indv(j),3)-skip_loc(3))*scale+(k3-1.5)*dx_loc/2.0*0.5&
                            -center_of_mass(clump_nr(j),3)
                    endif
                 end do
                 call get_cell_index(cell_index,cell_levl,xtest,ilevel+1,np)

                 ok=.false.
                 do j=1,np
                    !check wether neighbor is in a leaf cell at the right level
                    if(son(cell_index(j))==0.and.cell_levl(j)==(ilevel+1))ok(j)=.true.
                 end do                 

                 do j=1,np
                    neigh_cl(j)=flag2(cell_index(j))!nuber of clump the neighboring cell is in 
                    ok(j)=ok(j).and. neigh_cl(j)/=clump_nr(j) !neighboring cell is in another clump
                    ok(j)=ok(j).and. 0/=clump_nr(j) !clump number is not zero 
                 end do
                 
                 r_dot_n=0.
                 do j=1,np
                    do idim=1,3
                       r_dot_n(j)=r_dot_n(j)+n(idim)*r(j,idim)
                    end do
                 end do

                 do j=1,np
                    if (ok(j))then
                       do jdim=1,ndim
                          ekk_neigh(j)=ekk_neigh(j)+0.5*uold(cell_index(j),jdim+1)**2
                       end do
                       ekk_neigh(j)=ekk_neigh(j)/uold(cell_index(j),1)
                       P_neigh(j)=(gamma-1.0)*(uold(cell_index(j),ndim+2)-ekk_neigh(j))
                       Psurf(loc_clump_nr(j))=Psurf(loc_clump_nr(j))+r_dot_n(j)*0.25*dx_loc**2*0.5*(P_neigh(j)+P_cell(j))
                       if(debug.and.((P_neigh(j)-P_cell(j))/P_cell(j))**2>4.)print*,'caution, very high p contrast',(((P_neigh(j)-P_cell(j))/P_cell(j))**2)**0.5
                    endif
                 end do                 
              endif
           end do
        end do
     end do
  endif
#endif
     
end subroutine surface_int
!################################################################                 
!################################################################ 
!################################################################                 
!################################################################     
subroutine surface_pressure(ntest)
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
  ! This loop is also used to compute the surface pressure term.
  !---------------------------------------------------------------

  integer::info   
  integer::ipart,ip,ilevel,next_level
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
        call surface_int(ind_cell,ip,ilevel)
        ip=0
     endif
  end do
  if (ip>0)then 
     call surface_int(ind_cell,ip,ilevel)
  endif
   
#ifndef WITHOUTMPI     
  call virtual_peak_dp(Psurf,'sum')
#endif

end subroutine surface_pressure
