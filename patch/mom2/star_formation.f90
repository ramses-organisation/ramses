#if NDIM==3
subroutine star_formation(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  use cooling_module, ONLY: XH=>X
  use constants, only: Myr2sec, Gyr2sec, mH, pi, rhoc, twopi
  use random
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info,info2,dummy_io
  integer,parameter::tag=1120
#endif
  integer::ilevel
  !----------------------------------------------------------------------
  ! Description: This subroutine spawns star-particle of constant mass
  ! using a Poisson probability law if some gas condition are fulfilled.
  ! It modifies hydrodynamic variables according to mass conservation
  ! and assumes an isothermal transformation...
  ! On exit, the gas velocity and sound speed are unchanged.
  ! New star particles are synchronized with other collisionless particles.
  ! Array flag2 is used as temporary work space.
  ! Yann Rasera  10/2002-01/2003
  !----------------------------------------------------------------------
  ! local constants
  real(dp)::d0,mgas,mcell,t0
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:twotondim,1:3)::xc
  ! other variables
  integer ::ncache,nnew,ivar,ngrid,icpu,index_star,ndebris_tot,ilun=10
  integer ::igrid,ix,iy,iz,ind,i,n,iskip,nx_loc,idim
  integer ::ntot,ntot_all,nstar_corrected
  logical ::ok_free
  real(dp)::d,x,y,z,u,v,w,e,tg,zg
  real(dp)::mstar,dstar,tstar,nISM,nCOM,phi_t,phi_x,theta,sigs,scrit,b_turb
  real(dp)::T2,nH,T_poly,cs2,cs2_poly,trel,t_dyn,t_ff,uvar
  real(dp)::alpha0
  real(dp)::sigma2
  real(dp)::birth_epoch,factG,M2
  real(kind=8)::mlost_all,mtot_all
#ifndef WITHOUTMPI
  real(kind=8)::mlost,mtot
#endif
  real(kind=8)::PoissMean
  real(dp),dimension(1:3)::skip_loc
  real(dp)::dx,dx_loc,scale,vol_loc,dx_min,vol_min
  real(dp)::mdebris
  real(dp),dimension(1:nvector)::sfr_ff
  integer ,dimension(1:ncpu,1:IRandNumSize)::allseed
  integer ,dimension(1:nvector),save::ind_grid,ind_cell,nstar
  integer ,dimension(1:nvector),save::ind_grid_new,ind_cell_new,ind_part
  integer ,dimension(1:nvector),save::ind_debris

  logical ,dimension(1:nvector),save::ok,ok_new=.true.
  integer ,dimension(1:ncpu)::ntot_star_cpu,ntot_star_all
  character(LEN=80)::filename,filedir,fileloc,filedirini
  character(LEN=5)::nchar,ncharcpu
  logical::file_exist
#ifdef SOLVERmhd
  real(dp)::bx1,bx2,by1,by2,bz1,bz2,A,B,C,emag,beta,fbeta
#endif
#if NENER>0
  integer::irad
#endif

  if(numbtot(1,ilevel)==0) return
  if(.not. hydro)return
  if(ndim.ne.3)return
  if(static)return

  if(verbose)write(*,*)' Entering star_formation'

  if(sf_log_properties.and.ifout.gt.1) then
     call title(ifout-1,nchar)
     if(IOGROUPSIZEREP>0) then
        call title(((myid-1)/IOGROUPSIZEREP)+1,ncharcpu)
        filedirini='output_'//TRIM(nchar)//'/'
        filedir='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/'
     else
        filedir='output_'//TRIM(nchar)//'/'
     endif
     filename=TRIM(filedir)//'stars_'//TRIM(nchar)//'.out'
     ilun=myid+103
     call title(myid,nchar)
     fileloc=TRIM(filename)//TRIM(nchar)
     ! Wait for the token
#ifndef WITHOUTMPI
     if(IOGROUPSIZE>0) then
        if (mod(myid-1,IOGROUPSIZE)/=0) then
           call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
                & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
        end if
     endif
#endif

     inquire(file=fileloc,exist=file_exist)
     if((.not.file_exist).or.(abs(t-trestart).lt.dtnew(ilevel))) then
        open(ilun, file=fileloc, form='formatted')
        write(ilun,'(A24)',advance='no') '# event id  ilevel  mp  '
        do idim=1,ndim
           write(ilun,'(A2,I1,A2)',advance='no') 'xp',idim,'  '
        enddo
        do idim=1,ndim
           write(ilun,'(A2,I1,A2)',advance='no') 'vp',idim,'  '
        enddo
        do ivar=1,nvar
           if(ivar.ge.10) then
              write(ilun,'(A1,I2,A2)',advance='no') 'u',ivar,'  '
           else
              write(ilun,'(A1,I1,A2)',advance='no') 'u',ivar,'  '
           endif
        enddo
        write(ilun,'(A5)',advance='no') 'tag  '
        write(ilun,'(A1)') ' '
     else
        open(ilun, file=fileloc, status='old', position='append', action='write', form='formatted')
     endif
  endif

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim
  dx_min=(0.5D0**nlevelmax)*scale
  vol_min=dx_min**ndim

  trel=sf_trelax*Myr2sec/scale_t ! relaxation timescale

  ! ISM density threshold from H/cc to code units
  nISM = n_star
  if(cosmo)then
     nCOM = del_star*omega_b*rhoc*(h0/100)**2/aexp**3*XH/mH
     nISM = MAX(nCOM,nISM)
  endif
  d0   = nISM/scale_nH

  ! Initial star particle mass
  if(m_star < 0d0)then
     mstar=n_star/(scale_nH*aexp**3)*vol_min
  else
     mstar=m_star*mass_sph
  endif
  dstar=mstar/vol_loc

  factG = 1d0
  if(cosmo) factG = 3d0/4d0/twopi*omega_m*aexp

  ! Birth epoch as proper time
  if(use_proper_time)then
     birth_epoch=texp
  else
     birth_epoch=t
  endif

  ! Cells center position relative to grid center position
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     xc(ind,1)=(dble(ix)-0.5D0)*dx
     xc(ind,2)=(dble(iy)-0.5D0)*dx
     xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  ! If necessary, initialize random number generator
  if(localseed(1)==-1)then
     call rans(ncpu,iseed,allseed)
     localseed=allseed(myid,1:IRandNumSize)
  end if

  !------------------------------------------------
  ! Convert hydro variables to primitive variables
  !------------------------------------------------
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        do i=1,ngrid
           d=uold(ind_cell(i),1)
           u=uold(ind_cell(i),2)/d
           v=uold(ind_cell(i),3)/d
           w=uold(ind_cell(i),4)/d
           e=uold(ind_cell(i),5)
#ifdef SOLVERmhd
           bx1=uold(ind_cell(i),6)
           by1=uold(ind_cell(i),7)
           bz1=uold(ind_cell(i),8)
           bx2=uold(ind_cell(i),nvar+1)
           by2=uold(ind_cell(i),nvar+2)
           bz2=uold(ind_cell(i),nvar+3)
           e=e-0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)
#endif
           e=e-0.5d0*d*(u**2+v**2+w**2)
#if NENER>0
           do irad=0,nener-1
              e=e-uold(ind_cell(i),inener+irad)
           end do
#endif
           uold(ind_cell(i),1)=d
           uold(ind_cell(i),2)=u
           uold(ind_cell(i),3)=v
           uold(ind_cell(i),4)=w
           uold(ind_cell(i),5)=e/d
        end do
        do ivar=imetal,nvar
           do i=1,ngrid
              d=uold(ind_cell(i),1)
              w=uold(ind_cell(i),ivar)/d
              uold(ind_cell(i),ivar)=w
           end do
        end do
     end do
  end do

! get values of uold for density and velocities in virtual boundaries
#ifndef WITHOUTMPI
  do ivar=1,4
     call make_virtual_fine_dp(uold(1,ivar),ilevel)
  end do
#endif

  !------------------------------------------------
  ! Compute number of new stars in each cell
  !------------------------------------------------
  ntot=0
  ndebris_tot=0
  ! Loop over grids
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     ! Star formation criterion ---> logical array ok(i)
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        ! Flag leaf cells
        do i=1,ngrid
           ok(i)=son(ind_cell(i))==0
        end do
        if(sf_model > 0)then
           do i=1,ngrid
              ! if cell is a leaf cell
              if (ok(i)) then
                 ! Gas density
                 d         = uold(ind_cell(i),1)
                 ! Compute temperature in K/mu
                 T2        = (gamma-1.0d0)*uold(ind_cell(i),5)*scale_T2
                 ! Correct from polytrope
                 T_poly    = T2_star*(uold(ind_cell(i),1)*scale_nH/nISM)**(g_star-1.0d0)
                 T2        = T2-T_poly
                 ! Compute sound speed squared
                 cs2       = (gamma-1.0d0)*uold(ind_cell(i),5)
                 ! prevent numerical crash due to negative temperature
                 cs2       = max(cs2,smallc**2)
                 ! Correct from polytrope
                 cs2_poly  = (T2_star/scale_T2)*(uold(ind_cell(i),1)*scale_nH/nISM)**(g_star-1.0d0)
                 cs2       = cs2-cs2_poly
                 ! Turbulence 1D velocity dispersion
                 sigma2 = uold(ind_cell(i),ivirial1)*2.0/3.0

                 ! Density criterion
                 if(d<=d0) ok(i)=.false.
                 if(ok(i)) then
                    SELECT CASE (sf_model)

                       ! Multi-ff KM model
                       CASE (1)
                          ! Virial parameter
                          alpha0    = (5.0*(sigma2+cs2))/(pi*factG*d*dx_loc**2)
                          M2        = max(sigma2/cs2,1.0)
                          ! Turbulent forcing parameter (Federrath 2008 & 2010)
                          b_turb    = 0.4
                          ! Fudge for alpha dependence (KM 2005).
                          phi_x     = 1.12
                          ! The prefered value for eps_star = 1.0,
                          ! which represents the theoretical maximum efficiency.
                          sigs      = log(1.0+(b_turb**2)*(M2))
                          scrit     = log(((pi**2)/5)*(phi_x**2)*alpha0*(M2))
                          sfr_ff(i) = (eps_star/2.0)*exp(3.0/8.0*sigs)*(2.0-erfc((sigs-scrit)/sqrt(2.0*sigs)))

                       ! Multi-ff PN model
                       CASE (2)
                          ! Virial parameter
                          alpha0    = (5.0*(sigma2+cs2))/(pi*factG*d*dx_loc**2)
                          b_turb    = 0.4
                          ! Best fit values to the Multi-ff PN model (Hydro)
                          phi_t     = 1.0/0.49
                          theta     = 0.97
                          sigs      = log(1.0+(b_turb**2)*(sigma2/cs2))
                          scrit     = log(0.067/(theta**2)*alpha0*(sigma2/cs2))
                          sfr_ff(i) = (eps_star/(2.0*phi_t))*exp(3.0/8.0*sigs)*(2.0-erfc((sigs-scrit)/sqrt(2.0*sigs)))

                       ! Virial criterion threshold a la Hopkins
                       CASE (3)
                          alpha0       = (5.0*(sigma2+cs2))/(pi*factG*d*dx_loc**2)
                          if(alpha0<1.0) then
                             sfr_ff(i) = eps_star
                          else
                             sfr_ff(i) = 0.0
                             ok(i)     = .false.
                          endif

                       ! Padoan 2012 a la Semenov
                       CASE (4)
                          ! Feedback efficiency
                          t_dyn     = dx_loc/(2.0*sqrt(sigma2+cs2))
                          t_ff      = 0.5427*sqrt(1.0/(factG*max(d,smallr)))
                          sfr_ff(i) = eps_star*exp(-1.6*t_ff/t_dyn)

                       CASE (5)
                          ! Virial parameter
                          alpha0    = (5.0*(sigma2+cs2))/(pi*factG*d*dx_loc**2)
                          M2        = max(sigma2/cs2,smallr)
                          ! Turbulent forcing parameter (Federrath 2008 & 2010)
                          b_turb    = 0.4
                          ! Fudge for alpha dependence (KM 2005).
                          ! phi_x     = 1.12
                          ! The prefered value for eps_star = 1.0,
                          ! which represents the theoretical maximum efficiency.
                          sigs      = log(1.0+(b_turb**2)*(M2))
                          ! scrit     = log(alpha0*(1.0+(M2*(pi**2)/5.0)*(phi_x**2)))
                          scrit     = log(alpha0*(1.0+(2*(M2**2)/(1.0+M2))))
                          sfr_ff(i) = (eps_star/2.0)*exp(3.0/8.0*sigs)*(2.0-erfc((sigs-scrit)/sqrt(2.0*sigs)))

                       END SELECT
                 endif
              endif
           end do
        else
           ! Density criterion
           do i=1,ngrid
              sfr_ff(i) = eps_star
              d=uold(ind_cell(i),1)
              if(d<=d0)ok(i)=.false.
           end do
           ! Temperature criterion
           do i=1,ngrid
              T2=uold(ind_cell(i),5)*scale_T2*(gamma-1.0d0)
              nH=max(uold(ind_cell(i),1),smallr)*scale_nH
              T_poly=T2_star*(nH/nISM)**(g_star-1.0d0)
              T2=T2-T_poly
              if(T2>2d4)ok(i)=.false.
           end do
        endif
        ! Geometrical criterion
        if(ivar_refine>0)then
           do i=1,ngrid
              d=uold(ind_cell(i),ivar_refine)
              if(d<=var_cut_refine)ok(i)=.false.
           end do
        endif
        ! Calculate number of new stars in each cell using Poisson statistics
        do i=1,ngrid
           nstar(i)=0
           if(ok(i))then
              ! Compute mean number of events
              d=uold(ind_cell(i),1)
              mcell=d*vol_loc
              ! Free fall time of an homogeneous sphere
              tstar= .5427d0*sqrt(1.0d0/(factG*max(d,smallr)))
              ! Gas mass to be converted into stars
              mgas=dtnew(ilevel)*(sfr_ff(i)/tstar)*mcell
              ! Poisson mean
              PoissMean=mgas/mstar
              if((trel>0.).and.(.not.cosmo)) PoissMean = PoissMean*min((t/trel), 1.0d0)
              ! Compute Poisson realisation
              call poissdev(localseed,PoissMean,nstar(i))
              ! Compute depleted gas mass
              mgas=nstar(i)*mstar
              ! Security to prevent more than 90% of gas depletion
              if (mgas > 0.9d0*mcell) then
                 nstar_corrected=int(0.9d0*mcell/mstar)
                 mstar_lost=mstar_lost+(nstar(i)-nstar_corrected)*mstar
                 nstar(i)=nstar_corrected
              endif
              ! Compute new stars local statistics
              mstar_tot=mstar_tot+nstar(i)*mstar
              if(nstar(i)>0)then
                 ntot=ntot+1
                 if(f_w>0)ndebris_tot=ndebris_tot+1
              endif
           endif
        enddo
        ! Store nstar in array flag2
        do i=1,ngrid
           flag2(ind_cell(i))=nstar(i)
        end do
     end do
  end do

  !---------------------------------
  ! Check for free particle memory
  !---------------------------------
  ok_free=(numbp_free-ntot-ndebris_tot)>=0
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(numbp_free,numbp_free_tot,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  numbp_free_tot=numbp_free
#endif
  if(.not. ok_free)then
     write(*,*)'No more free memory for particles'
     write(*,*)'Increase npartmax'
#ifndef WITHOUTMPI
    call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
    stop
#endif
  end if

  !---------------------------------
  ! Compute global stars statistics
  !---------------------------------
#ifndef WITHOUTMPI
  mlost=mstar_lost; mtot=mstar_tot
  call MPI_ALLREDUCE(ntot,ntot_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(mtot,mtot_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(mlost,mlost_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  ntot_all=ntot
  mtot_all=mstar_tot
  mlost_all=mstar_lost
#endif
  ntot_star_cpu=0; ntot_star_all=0
  ntot_star_cpu(myid)=ntot
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(ntot_star_cpu,ntot_star_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  ntot_star_cpu(1)=ntot_star_all(1)
#endif
  do icpu=2,ncpu
     ntot_star_cpu(icpu)=ntot_star_cpu(icpu-1)+ntot_star_all(icpu)
  end do
  nstar_tot=nstar_tot+ntot_all
  if(myid==1)then
     if(ntot_all.gt.0)then
        write(*,'(" Level=",I6," New star=",I6," Tot=",I10," Mass=",1PE10.3," Lost=",0PF5.1,"%")')&
             & ilevel,ntot_all,nstar_tot,mtot_all,mlost_all/(mlost_all+mtot_all)*100.
     endif
  end if

  !------------------------------
  ! Create new star particles
  !------------------------------
  ! Starting identity number
  if(myid==1)then
     index_star=nstar_tot-ntot_all
  else
     index_star=nstar_tot-ntot_all+ntot_star_cpu(myid-1)
  end if

  ! Loop over grids
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do

     ! Loop over cells
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do

        ! Flag cells with at least one new star
        do i=1,ngrid
           ok(i)=flag2(ind_cell(i))>0
        end do

        ! Gather new star arrays
        nnew=0
        do i=1,ngrid
           if (ok(i))then
              nnew=nnew+1
              ind_grid_new(nnew)=ind_grid(i)
              ind_cell_new(nnew)=ind_cell(i)
           end if
        end do

        ! Update linked list for stars
        call remove_free(ind_part,nnew)
        call add_list(ind_part,ind_grid_new,ok_new,nnew)

        ! Update linked list for debris
        if(f_w>0)then
           call remove_free(ind_debris,nnew)
           call add_list(ind_debris,ind_grid_new,ok_new,nnew)
        endif

        ! Calculate new star particle and modify gas density
        do i=1,nnew
           index_star=index_star+1

           ! Get gas variables
           n=flag2(ind_cell_new(i))
           d=uold(ind_cell_new(i),1)
           u=uold(ind_cell_new(i),2)
           v=uold(ind_cell_new(i),3)
           w=uold(ind_cell_new(i),4)
           x=(xg(ind_grid_new(i),1)+xc(ind,1)-skip_loc(1))*scale
           y=(xg(ind_grid_new(i),2)+xc(ind,2)-skip_loc(2))*scale
           z=(xg(ind_grid_new(i),3)+xc(ind,3)-skip_loc(3))*scale
           tg=uold(ind_cell_new(i),5)*(gamma-1)*scale_T2
           if(metal)zg=uold(ind_cell_new(i),imetal)

           ! Set star particle variables
           tp(ind_part(i)) = birth_epoch  ! Birth epoch
           mp(ind_part(i)) = n*mstar      ! Mass
           levelp(ind_part(i)) = ilevel   ! Level
           idp(ind_part(i)) = index_star  ! Star identity
           typep(ind_part(i))%family = FAM_STAR
           typep(ind_part(i))%tag = 0
           xp(ind_part(i),1) = x
           xp(ind_part(i),2) = y
           xp(ind_part(i),3) = z
           vp(ind_part(i),1) = u
           vp(ind_part(i),2) = v
           vp(ind_part(i),3) = w
           if(metal)zp(ind_part(i)) = zg  ! Initial star metallicity

           ! Set GMC particle variables
           if(f_w>0)then
              ! Compute GMC mass without more than 50% of gas depletion
              mdebris=min(f_w*n*mstar,0.5d0*(d*vol_loc-n*mstar))
              ! Add supernova ejecta
              mdebris=mdebris+eta_sn*n*mstar
              ! Remove ejecta from the long lived star mass
              mp(ind_part(i))=n*mstar-eta_sn*n*mstar
              ! Set GMC particle variables
              tp(ind_debris(i))=birth_epoch  ! Birth epoch
              mp(ind_debris(i))=mdebris      ! Mass
              levelp(ind_debris(i))=ilevel   ! Level
              idp(ind_debris(i))=-n          ! Number of individual stars
              typep(ind_debris(i))%family = FAM_DEBRIS
              typep(ind_debris(i))%tag = 0
              xp(ind_debris(i),1)=x
              xp(ind_debris(i),2)=y
              xp(ind_debris(i),3)=z
              vp(ind_debris(i),1)=u
              vp(ind_debris(i),2)=v
              vp(ind_debris(i),3)=w
              ! GMC metallicity + yield from ejecta
              if(metal)zp(ind_debris(i))=zg+eta_sn*yield*(1-zg)*n*mstar/mdebris
           endif

           if(sf_log_properties) then
              write(ilun,'(I10)',advance='no') 0
              write(ilun,'(2I10,E24.12)',advance='no') idp(ind_part(i)),ilevel,mp(ind_part(i))
              do idim=1,ndim
                 write(ilun,'(E24.12)',advance='no') xp(ind_part(i),idim)
              enddo
              do idim=1,ndim
                 write(ilun,'(E24.12)',advance='no') vp(ind_part(i),idim)
              enddo
              write(ilun,'(E24.12)',advance='no') uold(ind_cell_new(i),1)
              do ivar=2,nvar
                 if(ivar.eq.ndim+2)then
                    ! Temperature
                    uvar=(gamma-1.0d0)*(uold(ind_cell_new(i),ndim+2))*scale_T2
                 else
                    uvar=uold(ind_cell_new(i),ivar)
                 endif
                 write(ilun,'(E24.12)',advance='no') uvar
              enddo
              write(ilun,'(I10)',advance='no') typep(ind_part(i))%tag
              write(ilun,'(A1)') ' '
           endif


        end do
        ! End loop over new star particles

        ! Modify gas density according to mass depletion
        do i=1,ngrid
           if(flag2(ind_cell(i))>0)then
              n=flag2(ind_cell(i))
              d=uold(ind_cell(i),1)
              uold(ind_cell(i),1)=max(d-n*dstar,0.5d0*d)
           endif
        end do

     end do
     ! End loop over cells
  end do
  ! End loop over grids

  !---------------------------------------------------------
  ! Convert hydro variables back to conservative variables
  !---------------------------------------------------------
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        do i=1,ngrid
           d=uold(ind_cell(i),1)
           u=uold(ind_cell(i),2)
           v=uold(ind_cell(i),3)
           w=uold(ind_cell(i),4)
           e=uold(ind_cell(i),5)*d
#ifdef SOLVERmhd
           bx1=uold(ind_cell(i),6)
           by1=uold(ind_cell(i),7)
           bz1=uold(ind_cell(i),8)
           bx2=uold(ind_cell(i),nvar+1)
           by2=uold(ind_cell(i),nvar+2)
           bz2=uold(ind_cell(i),nvar+3)
           e=e+0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)
#endif
           e=e+0.5d0*d*(u**2+v**2+w**2)
#if NENER>0
           do irad=0,nener-1
              e=e+uold(ind_cell(i),inener+irad)
           end do
#endif
           uold(ind_cell(i),1)=d
           uold(ind_cell(i),2)=d*u
           uold(ind_cell(i),3)=d*v
           uold(ind_cell(i),4)=d*w
           uold(ind_cell(i),5)=e
        end do
        do ivar=imetal,nvar
           do i=1,ngrid
              d=uold(ind_cell(i),1)
              w=uold(ind_cell(i),ivar)
              uold(ind_cell(i),ivar)=d*w
           end do
        end do
     end do
  end do

  if(sf_log_properties) close(ilun)

end subroutine star_formation
#endif
!################################################################
!################################################################
!################################################################
!################################################################
subroutine getnbor(ind_cell,ind_father,ncell,ilevel)
  use amr_commons
  implicit none
  integer::ncell,ilevel
  integer,dimension(1:nvector)::ind_cell
  integer,dimension(1:nvector,0:twondim)::ind_father
  !-----------------------------------------------------------------
  ! This subroutine determines the 2*ndim neighboring cells
  ! cells of the input cell (ind_cell).
  ! If for some reasons they don't exist, the routine returns
  ! the input cell.
  !-----------------------------------------------------------------
  integer::i,j,iok,ind
  integer,dimension(1:nvector),save::ind_grid_father,pos
  integer,dimension(1:nvector,0:twondim),save::igridn,igridn_ok
  integer,dimension(1:nvector,1:twondim),save::icelln_ok


  if(ilevel==1)then
     write(*,*) 'Warning: attempting to form stars on level 1 --> this is not allowed ...'
     return
  endif

  ! Get father cell
  do i=1,ncell
     ind_father(i,0)=ind_cell(i)
  end do

  ! Get father cell position in the grid
  do i=1,ncell
     pos(i)=(ind_father(i,0)-ncoarse-1)/ngridmax+1
  end do

  ! Get father grid
  do i=1,ncell
     ind_grid_father(i)=ind_father(i,0)-ncoarse-(pos(i)-1)*ngridmax
  end do

  ! Get neighboring father grids
  call getnborgrids(ind_grid_father,igridn,ncell)

  ! Loop over position
  do ind=1,twotondim

     ! Select father cells that sit at position ind
     do j=0,twondim
        iok=0
        do i=1,ncell
           if(pos(i)==ind)then
              iok=iok+1
              igridn_ok(iok,j)=igridn(i,j)
           end if
        end do
     end do

     ! Get neighboring cells for selected cells
     if(iok>0)call getnborcells(igridn_ok,ind,icelln_ok,iok)

     ! Update neighboring father cells for selected cells
     do j=1,twondim
        iok=0
        do i=1,ncell
           if(pos(i)==ind)then
              iok=iok+1
              if(icelln_ok(iok,j)>0)then
                 ind_father(i,j)=icelln_ok(iok,j)
              else
                 ind_father(i,j)=ind_cell(i)
              end if
           end if
        end do
     end do

  end do


end subroutine getnbor
!##############################################################
!##############################################################
!##############################################################
!##############################################################
function erfc(x)

! complementary error function
  use amr_commons, ONLY: dp
  implicit none
  real(dp) erfc
  real(dp) x, y
  real(kind=8) pv, ph
  real(kind=8) q0, q1, q2, q3, q4, q5, q6, q7
  real(kind=8) p0, p1, p2, p3, p4, p5, p6, p7
  parameter(pv= 1.26974899965115684d+01, ph= 6.10399733098688199d+00)
  parameter(p0= 2.96316885199227378d-01, p1= 1.81581125134637070d-01)
  parameter(p2= 6.81866451424939493d-02, p3= 1.56907543161966709d-02)
  parameter(p4= 2.21290116681517573d-03, p5= 1.91395813098742864d-04)
  parameter(p6= 9.71013284010551623d-06, p7= 1.66642447174307753d-07)
  parameter(q0= 6.12158644495538758d-02, q1= 5.50942780056002085d-01)
  parameter(q2= 1.53039662058770397d+00, q3= 2.99957952311300634d+00)
  parameter(q4= 4.95867777128246701d+00, q5= 7.41471251099335407d+00)
  parameter(q6= 1.04765104356545238d+01, q7= 1.48455557345597957d+01)

  y = x*x
  y = exp(-y)*x*(p7/(y+q7)+p6/(y+q6) + p5/(y+q5)+p4/(y+q4)+p3/(y+q3) &
       &       + p2/(y+q2)+p1/(y+q1)+p0/(y+q0))
  if (x < ph) y = y+2d0/(exp(pv*x)+1.0)
  erfc = y

  return

end function erfc
