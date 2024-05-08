subroutine init_time
  use amr_commons
  use hydro_commons
  use pm_commons
  use cooling_module
#ifdef grackle
  use grackle_parameters
#endif
#ifdef RT
  use rt_cooling_module
#endif
  use mpi_mod
  implicit none
  integer::i,Nmodel
  real(kind=8)::T2_sim
#ifdef grackle
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  logical::file_exists
#ifndef WITHOUTMPI
  integer::info
#endif
#endif

  if(nrestart==0)then
     if(cosmo)then
        ! Get cosmological parameters from input files
        call init_cosmo
     else
        ! Get parameters from input files
        if(initfile(levelmin).ne.' '.and.filetype.eq.'grafic')then
           call init_file
        endif
        t=0
        aexp=1
     end if
  end if

  if(cosmo)then

     ! Allocate look-up tables
     n_frw=1000
     allocate(aexp_frw(0:n_frw),hexp_frw(0:n_frw))
     allocate(tau_frw(0:n_frw),t_frw(0:n_frw))

     ! Compute Friedman model look up table
     if(myid==1)write(*,*)'Computing Friedman model'
     call friedman(dble(omega_m),dble(omega_l),dble(omega_k), &
          & 1d-6,dble(aexp_ini), &
          & aexp_frw,hexp_frw,tau_frw,t_frw,n_frw)

     ! Compute initial conformal time
     ! Find neighboring expansion factors
     i=1
     do while(aexp_frw(i)>aexp.and.i<n_frw)
        i=i+1
     end do
     ! Interploate time
     if(nrestart==0)then
        t=tau_frw(i)*(aexp-aexp_frw(i-1))/(aexp_frw(i)-aexp_frw(i-1))+ &
             & tau_frw(i-1)*(aexp-aexp_frw(i))/(aexp_frw(i-1)-aexp_frw(i))
        aexp=aexp_frw(i)*(t-tau_frw(i-1))/(tau_frw(i)-tau_frw(i-1))+ &
             & aexp_frw(i-1)*(t-tau_frw(i))/(tau_frw(i-1)-tau_frw(i))
        hexp=hexp_frw(i)*(t-tau_frw(i-1))/(tau_frw(i)-tau_frw(i-1))+ &
             & hexp_frw(i-1)*(t-tau_frw(i))/(tau_frw(i-1)-tau_frw(i))
     end if
     texp=t_frw(i)*(t-tau_frw(i-1))/(tau_frw(i)-tau_frw(i-1))+ &
          & t_frw(i-1)*(t-tau_frw(i))/(tau_frw(i-1)-tau_frw(i))
  else
     texp=t
  end if

  ! Initialize cooling model
#ifdef grackle
  if(use_grackle==1)then
     if(myid==1)then
        write(*,'(A50)')"__________________________________________________"
        write(*,*)'Grackle - Computing cooling model'
        write(*,*)'Grackle - Loading ',TRIM(grackle_data_file)
        write(*,'(A50)')"__________________________________________________"
     endif
     INQUIRE(FILE=grackle_data_file,EXIST=file_exists)
     if(.not.file_exists) then
        if(myid==1) write(*,*) grackle_data_file," not found"
        call clean_stop
     endif

     iresult = set_default_chemistry_parameters(my_grackle_data)
     if(iresult.eq.0)then
         write(*,*) 'Grackle - error in initialize_chemistry_data'
#ifndef WITHOUTMPI
         call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
         stop
#endif
     endif
     my_grackle_data%use_grackle = use_grackle
     my_grackle_data%with_radiative_cooling = grackle_with_radiative_cooling
     my_grackle_data%primordial_chemistry = grackle_primordial_chemistry
     my_grackle_data%metal_cooling = grackle_metal_cooling
     my_grackle_data%UVbackground = grackle_UVbackground
     my_grackle_data%cmb_temperature_floor = grackle_cmb_temperature_floor
     my_grackle_data%h2_on_dust = grackle_h2_on_dust
     my_grackle_data%photoelectric_heating = grackle_photoelectric_heating
     my_grackle_data%use_volumetric_heating_rate = grackle_use_volumetric_heating_rate
     my_grackle_data%use_specific_heating_rate = grackle_use_specific_heating_rate
     my_grackle_data%three_body_rate = grackle_three_body_rate
     my_grackle_data%cie_cooling = grackle_cie_cooling
     my_grackle_data%h2_optical_depth_approximation = grackle_h2_optical_depth_approximation
     my_grackle_data%ih2co = grackle_ih2co
     my_grackle_data%ipiht = grackle_ipiht
     my_grackle_data%NumberOfTemperatureBins = grackle_NumberOfTemperatureBins
     my_grackle_data%CaseBRecombination = grackle_CaseBRecombination
     my_grackle_data%Compton_xray_heating = grackle_Compton_xray_heating
     my_grackle_data%LWbackground_sawtooth_suppression = grackle_LWbackground_sawtooth_suppression
     my_grackle_data%NumberOfDustTemperatureBins = grackle_NumberOfDustTemperatureBins
     my_grackle_data%use_radiative_transfer = grackle_use_radiative_transfer
     my_grackle_data%radiative_transfer_coupled_rate_solver = grackle_radiative_transfer_coupled_rate_solver
     my_grackle_data%radiative_transfer_intermediate_step = grackle_radiative_transfer_intermediate_step
     my_grackle_data%radiative_transfer_hydrogen_only = grackle_radiative_transfer_hydrogen_only
     my_grackle_data%self_shielding_method = grackle_self_shielding_method
     my_grackle_data%Gamma = grackle_Gamma
     my_grackle_data%photoelectric_heating_rate = grackle_photoelectric_heating_rate
     my_grackle_data%HydrogenFractionByMass = grackle_HydrogenFractionByMass
     my_grackle_data%DeuteriumToHydrogenRatio = grackle_DeuteriumToHydrogenRatio
     my_grackle_data%SolarMetalFractionByMass = grackle_SolarMetalFractionByMass
     my_grackle_data%TemperatureStart = grackle_TemperatureStart
     my_grackle_data%TemperatureEnd = grackle_TemperatureEnd
     my_grackle_data%DustTemperatureStart = grackle_DustTemperatureStart
     my_grackle_data%DustTemperatureEnd = grackle_DustTemperatureEnd
     my_grackle_data%LWbackground_intensity = grackle_LWbackground_intensity
     my_grackle_data%UVbackground_redshift_on = grackle_UVbackground_redshift_on
     my_grackle_data%UVbackground_redshift_off = grackle_UVbackground_redshift_off
     my_grackle_data%UVbackground_redshift_fullon = grackle_UVbackground_redshift_fullon
     my_grackle_data%UVbackground_redshift_drop = grackle_UVbackground_redshift_drop
     my_grackle_data%cloudy_electron_fraction_factor = grackle_cloudy_electron_fraction_factor
     grackle_data_file = TRIM(grackle_data_file)//C_NULL_CHAR
     my_grackle_data%grackle_data_file = C_LOC(grackle_data_file(1:1))

     ! Grackle units
     call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
     my_grackle_units%comoving_coordinates = 0
     my_grackle_units%density_units = scale_d
     my_grackle_units%length_units = scale_l
     my_grackle_units%time_units = scale_t
     my_grackle_units%velocity_units = scale_v
     my_grackle_units%a_units = 1.0d0
     !Set initial expansion factor (for internal units).
     !Set expansion factor to 1 for non-cosmological simulation
     ! Safety for GRACKLE initialisation
     my_grackle_units%a_value = aexp_ini

     if(cosmo) then
        my_grackle_units%comoving_coordinates = 1
        ! Reonization redshift has to be later than starting redshift
        z_reion=min(1d0/(1.1d0*aexp_ini)-1d0,z_reion)
        ! Approximate initial temperature
        T2_start=1.356d-2/aexp_ini**2
        if(nrestart==0)then
           if(myid==1)write(*,*)'Starting with T/mu (K) = ',T2_start
        end if
     endif

     ! Initialize the Grackle data
     iresult = initialize_chemistry_data(my_grackle_units)
     ! Enforce UVbackground starting redshift in grackle
     my_grackle_data%UVbackground_redshift_on = grackle_UVbackground_redshift_on
     my_grackle_data%UVbackground_redshift_fullon = grackle_UVbackground_redshift_fullon
     if(iresult.eq.0)then
         write(*,*) 'Grackle - error in initialize_chemistry_data'
#ifndef WITHOUTMPI
         call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
         stop
#endif
     endif
     my_grackle_fields%grid_dimension = C_LOC(gr_dimension)
     my_grackle_fields%grid_start = C_LOC(gr_start)
     my_grackle_fields%grid_end = C_LOC(gr_end)
     ! Point to grackle fields
     my_grackle_fields%density = C_LOC(gr_density)
     my_grackle_fields%HI_density = C_LOC(gr_HI_density)
     my_grackle_fields%HII_density = C_LOC(gr_HII_density)
     my_grackle_fields%HM_density = C_LOC(gr_HM_density)
     my_grackle_fields%HeI_density = C_LOC(gr_HeI_density)
     my_grackle_fields%HeII_density = C_LOC(gr_HeII_density)
     my_grackle_fields%HeIII_density = C_LOC(gr_HeIII_density)
     my_grackle_fields%H2I_density = C_LOC(gr_H2I_density)
     my_grackle_fields%H2II_density = C_LOC(gr_H2II_density)
     my_grackle_fields%DI_density = C_LOC(gr_DI_density)
     my_grackle_fields%DII_density = C_LOC(gr_DII_density)
     my_grackle_fields%HDI_density = C_LOC(gr_HDI_density)
     my_grackle_fields%e_density = C_LOC(gr_e_density)
     my_grackle_fields%metal_density = C_LOC(gr_metal_density)
     my_grackle_fields%internal_energy = C_LOC(gr_energy)
     my_grackle_fields%x_velocity = C_LOC(gr_x_velocity)
     my_grackle_fields%y_velocity = C_LOC(gr_y_velocity)
     my_grackle_fields%z_velocity = C_LOC(gr_z_velocity)
     my_grackle_fields%volumetric_heating_rate = C_LOC(gr_volumetric_heating_rate)
     my_grackle_fields%specific_heating_rate = C_LOC(gr_specific_heating_rate)
     my_grackle_fields%RT_HI_ionization_rate = C_LOC(gr_RT_HI_ionization_rate)
     my_grackle_fields%RT_HeI_ionization_rate = C_LOC(gr_RT_HeI_ionization_rate)
     my_grackle_fields%RT_HeII_ionization_rate = C_LOC(gr_RT_HeII_ionization_rate)
     my_grackle_fields%RT_H2_dissociation_rate = C_LOC(gr_RT_H2_dissociation_rate)
     my_grackle_fields%RT_heating_rate = C_LOC(gr_RT_heating_rate)
     do i=1,nvector
        gr_x_velocity(i) = 0d0
        gr_y_velocity(i) = 0d0
        gr_z_velocity(i) = 0d0
        gr_HII_density(i) = 0d0
        gr_HM_density(i) = 0d0
        gr_HeII_density(i) = 0d0
        gr_HeIII_density(i) = 0d0
        gr_H2I_density(i) = 0d0
        gr_H2II_density(i) = 0d0
        gr_DII_density(i) = 0d0
        gr_HDI_density(i) = 0d0
        gr_e_density(i) = 0d0
        gr_volumetric_heating_rate(i) = 0d0
        gr_specific_heating_rate(i) = 0d0
        gr_RT_HI_ionization_rate(i) = 0d0
        gr_RT_HeI_ionization_rate(i) = 0d0
        gr_RT_HeII_ionization_rate(i) = 0d0
        gr_RT_H2_dissociation_rate(i) = 0d0
        gr_RT_heating_rate(i) = 0d0
     enddo
  else
     if(cooling.and..not.(neq_chem.or.rt).and..not.cooling_ism)then
        if(myid==1)write(*,*)'Computing cooling model'
        Nmodel=-1
        if(.not. haardt_madau)then
           Nmodel=2
        endif
        if(cosmo)then
           ! Reonization redshift has to be later than starting redshift
           z_reion=min(1d0/(1.1d0*aexp_ini)-1d0,z_reion)
           call set_model(Nmodel,dble(J21*1d-21),-1.0d0,dble(a_spec),-1.0d0,dble(z_reion), &
                & -1,2, &
                & dble(h0/100),dble(omega_b),dble(omega_m),dble(omega_l), &
                & dble(aexp_ini),T2_sim)
           T2_start=T2_sim
           if(nrestart==0)then
              if(myid==1)write(*,*)'Starting with T/mu (K) = ',T2_start
           end if
        else
           call set_model(Nmodel,dble(J21*1d-21),-1.0d0,dble(a_spec),-1.0d0,dble(z_reion), &
                & -1,2, &
                & dble(70d0/100d0),0.04d0,0.3d0,0.7d0, &
                & dble(aexp_ini),T2_sim)
        endif
     end if
  endif
#else
  if(cooling.and..not.(neq_chem.or.rt).and..not.cooling_ism) then
     if(myid==1)write(*,*)'Computing cooling model'
     Nmodel=-1
     if(.not. haardt_madau)then
        Nmodel=2
     endif
     if(cosmo)then
        ! Reonization redshift has to be later than starting redshift
        z_reion=min(1d0/(1.1d0*aexp_ini)-1d0,z_reion)
        call set_model(Nmodel,dble(J21*1d-21),-1.0d0,dble(a_spec),-1.0d0,dble(z_reion), &
             & -1,2, &
             & dble(h0/100.),dble(omega_b),dble(omega_m),dble(omega_l), &
             & dble(aexp_ini),T2_sim)
        T2_start=T2_sim
        if(nrestart==0)then
           if(myid==1)write(*,*)'Starting with T/mu (K) = ',T2_start
        end if
     else
        call set_model(Nmodel,dble(J21*1d-21),-1.0d0,dble(a_spec),-1.0d0,dble(z_reion), &
             & -1,2, &
             & dble(70./100.),dble(0.04),dble(0.3),dble(0.7), &
             & dble(aexp_ini),T2_sim)
     endif
  end if
#endif
#ifdef RT
  if(neq_chem.or.rt) then
     if(myid==1)write(*,*)'Computing thermochemistry model'
     Nmodel=-1
     if(.not. haardt_madau)then
        Nmodel=2
     endif
     if(cosmo)then
        ! Reonization redshift has to be later than starting redshift
        z_reion=min(1d0/(1.1d0*aexp_ini)-1d0,z_reion)
        call rt_set_model(dble(h0/100.),dble(omega_b),dble(omega_m),dble(omega_l), &
             & dble(aexp_ini),T2_sim)
        T2_start=T2_sim
        if(nrestart==0)then
           if(myid==1)write(*,*)'Starting with T/mu (K) = ',T2_start
        end if
     else
        call rt_set_model(dble(70./100.),dble(0.04),dble(0.3),dble(0.7), &
             & dble(aexp_ini),T2_sim)
     endif
  end if
#endif

end subroutine init_time

subroutine init_file
  use amr_commons
  use hydro_commons
  use pm_commons
  use mpi_mod
  implicit none
  !------------------------------------------------------
  ! Read geometrical parameters in the initial condition files.
  ! Initial conditions are supposed to be made by
  ! Bertschinger's grafic version 2.0 code.
  !------------------------------------------------------
  integer:: ilevel,nx_loc,ny_loc,nz_loc
  real(sp)::dxini0,xoff10,xoff20,xoff30,astart0,omega_m0,omega_l0,h00
  character(LEN=80)::filename
  logical::ok
  integer,parameter::tag=1116
#ifndef WITHOUTMPI
  integer::dummy_io,info2
#endif

  if(verbose)write(*,*)'Entering init_file'

  ! Reading initial conditions parameters only


  nlevelmax_part=levelmin-1
  do ilevel=levelmin,nlevelmax
     if(initfile(ilevel).ne.' ')then
        filename=TRIM(initfile(ilevel))//'/ic_d'

        ! Wait for the token
#ifndef WITHOUTMPI
        if(IOGROUPSIZE>0) then
           if (mod(myid-1,IOGROUPSIZE)/=0) then
              call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
                   & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
           end if
        endif
#endif
        INQUIRE(file=filename,exist=ok)
        if(.not.ok)then
           if(myid==1)then
              write(*,*)'File '//TRIM(filename)//' does not exist'
           end if
           call clean_stop
        end if
        open(10,file=filename,form='unformatted')
        if(myid==1)write(*,*)'Reading file '//TRIM(filename)
        rewind 10
        read(10)n1(ilevel),n2(ilevel),n3(ilevel),dxini0 &
             & ,xoff10,xoff20,xoff30 &
             & ,astart0,omega_m0,omega_l0,h00
        close(10)

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

        dxini(ilevel)=dxini0
        xoff1(ilevel)=xoff10
        xoff2(ilevel)=xoff20
        xoff3(ilevel)=xoff30
        nlevelmax_part=nlevelmax_part+1
     endif
  end do


  ! Check compatibility with run parameters
  nx_loc=icoarse_max-icoarse_min+1
  ny_loc=jcoarse_max-jcoarse_min+1
  nz_loc=kcoarse_max-kcoarse_min+1
  if(         nx_loc.ne.n1(levelmin)/2**levelmin &
       & .or. ny_loc.ne.n2(levelmin)/2**levelmin &
       & .or. nz_loc.ne.n3(levelmin)/2**levelmin) then
     write(*,*)'coarser grid is not compatible with initial conditions file'
     write(*,*)'Found    n1=',n1(levelmin),&
          &            ' n2=',n2(levelmin),&
          &            ' n3=',n3(levelmin)
     write(*,*)'Expected n1=',nx_loc*2**levelmin &
          &           ,' n2=',ny_loc*2**levelmin &
          &           ,' n3=',nz_loc*2**levelmin
     call clean_stop
  end if

  ! Write initial conditions parameters
  if(myid==1)then
     do ilevel=levelmin,nlevelmax_part
        write(*,'(" Initial conditions for level =",I4)')ilevel
        write(*,'(" n1=",I4," n2=",I4," n3=",I4)') &
             & n1(ilevel),&
             & n2(ilevel),&
             & n3(ilevel)
        write(*,'(" dx=",1pe10.3)')dxini(ilevel)
        write(*,'(" xoff=",1pe10.3," yoff=",1pe10.3," zoff=",&
             & 1pe10.3)') &
             & xoff1(ilevel),&
             & xoff2(ilevel),&
             & xoff3(ilevel)
     end do
  end if

end subroutine init_file


subroutine init_cosmo
  use amr_commons
  use hydro_commons
  use pm_commons
  use gadgetreadfilemod
  use mpi_mod
  implicit none
  !------------------------------------------------------
  ! Read cosmological and geometrical parameters
  ! in the initial condition files.
  ! Initial conditions are supposed to be made by
  ! Bertschinger's grafic version 2.0 code.
  !------------------------------------------------------
  integer:: ilevel
  real(sp)::dxini0,xoff10,xoff20,xoff30,astart0,omega_m0,omega_l0,h00
  character(LEN=80)::filename
  character(LEN=5)::nchar
  logical::ok
  TYPE(gadgetheadertype) :: gadgetheader
  integer::i
  integer,parameter::tag=1117
#ifndef WITHOUTMPI
  integer::dummy_io,info2
#endif

  if(verbose)write(*,*)'Entering init_cosmo'

  if(initfile(levelmin)==' ')then
     write(*,*)'You need to specifiy at least one level of initial condition'
     call clean_stop
  end if

  SELECT CASE (filetype)
  case ('grafic', 'ascii')
     ! Reading initial conditions parameters only
     aexp=2
     nlevelmax_part=levelmin-1
     do ilevel=levelmin,nlevelmax
        if(initfile(ilevel).ne.' ')then
           if(multiple)then
              call title(myid,nchar)
              filename=TRIM(initfile(ilevel))//'/dir_deltab/ic_deltab.'//TRIM(nchar)
           else
              filename=TRIM(initfile(ilevel))//'/ic_deltab'
           endif

           ! Wait for the token
#ifndef WITHOUTMPI
           if(IOGROUPSIZE>0) then
              if (mod(myid-1,IOGROUPSIZE)/=0) then
                 call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
                      & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
              end if
           endif
#endif

           INQUIRE(file=filename,exist=ok)
           if(.not.ok)then
              if(myid==1)then
                 write(*,*)'File '//TRIM(filename)//' does not exist'
              end if
              call clean_stop
           end if
           open(10,file=filename,form='unformatted')
           if(myid==1)write(*,*)'Reading file '//TRIM(filename)
           rewind 10
           read(10)n1(ilevel),n2(ilevel),n3(ilevel),dxini0 &
                & ,xoff10,xoff20,xoff30 &
                & ,astart0,omega_m0,omega_l0,h00
           close(10)

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

           dxini(ilevel)=dxini0
           xoff1(ilevel)=xoff10
           xoff2(ilevel)=xoff20
           xoff3(ilevel)=xoff30
           astart(ilevel)=astart0
           omega_m=omega_m0
           omega_l=omega_l0
           h0=h00
           aexp=MIN(aexp,astart(ilevel))
           nlevelmax_part=nlevelmax_part+1
           ! Compute SPH equivalent mass (initial gas mass resolution)
           mass_sph=omega_b/omega_m*0.5d0**(ndim*ilevel)

        endif
     end do

     ! Compute initial expansion factor
     if(aexp_ini.lt.1.0)then
        aexp=aexp_ini
     else
        aexp_ini=aexp
     endif

     ! Check compatibility with run parameters
     if(.not. multiple) then
        if(         nx.ne.n1(levelmin)/2**levelmin &
             & .or. ny.ne.n2(levelmin)/2**levelmin &
             & .or. nz.ne.n3(levelmin)/2**levelmin) then
           write(*,*)'coarser grid is not compatible with initial conditions file'
           write(*,*)'Found    n1=',n1(levelmin),&
                &            ' n2=',n2(levelmin),&
                &            ' n3=',n3(levelmin)
           write(*,*)'Expected n1=',nx*2**levelmin &
                &           ,' n2=',ny*2**levelmin &
                &           ,' n3=',nz*2**levelmin
           call clean_stop
        endif
     end if

     ! Compute box length in the initial conditions in units of h-1 Mpc
     boxlen_ini=dble(nx)*2**levelmin*dxini(levelmin)*(h0/100)

  CASE ('gadget')
     if (verbose) write(*,*)'Reading in gadget format from '//TRIM(initfile(levelmin))
     call gadgetreadheader(TRIM(initfile(levelmin)), 0, gadgetheader, ok)
     if(.not.ok) call clean_stop
     do i=1,6
        if (i .ne. 2) then
           if (gadgetheader%nparttotal(i) .ne. 0) then
              write(*,*) 'Non DM particles present in bin ', i
              call clean_stop
           endif
        endif
     enddo
     if (gadgetheader%mass(2) == 0) then
        write(*,*) 'Particles have different masses, not supported'
        call clean_stop
     endif
     omega_m = gadgetheader%omega0
     omega_l = gadgetheader%omegalambda
     h0 = gadgetheader%hubbleparam * 100d0
     boxlen_ini = gadgetheader%boxsize
     aexp = gadgetheader%time
     aexp_ini = aexp
     ! Compute SPH equivalent mass (initial gas mass resolution)
     mass_sph=omega_b/omega_m*0.5d0**(ndim*levelmin)
     nlevelmax_part = levelmin
     astart(levelmin) = aexp
     xoff1(levelmin)=0
     xoff2(levelmin)=0
     xoff3(levelmin)=0
     dxini(levelmin) = boxlen_ini/(nx*2**levelmin*(h0/100))

  CASE DEFAULT
     write(*,*) 'Unsupported input format '//filetype
     call clean_stop
  END SELECT

  ! Write cosmological parameters
  if(myid==1)then
     write(*,'(" Cosmological parameters:")')
     write(*,'(" aexp=",1pe10.3," H0=",1pe10.3," km s-1 Mpc-1")')aexp,h0
     write(*,'(" omega_m=",F7.3," omega_l=",F7.3," omega_b=",F7.3)')omega_m,omega_l,omega_b
     write(*,'(" box size=",1pe10.3," h-1 Mpc")')boxlen_ini
  end if
  omega_k=1d0-omega_l-omega_m

  ! Compute linear scaling factor between aexp and astart(ilevel)
  do ilevel=levelmin,nlevelmax_part
     dfact(ilevel)=d1a(aexp)/d1a(astart(ilevel))
     vfact(ilevel)=astart(ilevel)*fpeebl(astart(ilevel)) & ! Same scale factor as in grafic1
          & *sqrt(omega_m/astart(ilevel)+omega_l*astart(ilevel)*astart(ilevel)+omega_k) &
          & /astart(ilevel)*h0
  end do

  ! Write initial conditions parameters
  do ilevel=levelmin,nlevelmax_part
     if(myid==1)then
        write(*,'(" Initial conditions for level =",I4)')ilevel
        write(*,'(" dx=",1pe10.3," h-1 Mpc")')dxini(ilevel)*h0/100.
     endif
     if(.not.multiple)then
        if(myid==1)then
           write(*,'(" n1=",I4," n2=",I4," n3=",I4)') &
                & n1(ilevel),&
                & n2(ilevel),&
                & n3(ilevel)
           write(*,'(" xoff=",1pe10.3," yoff=",1pe10.3," zoff=",&
                & 1pe10.3," h-1 Mpc")') &
                & xoff1(ilevel)*h0/100.,&
                & xoff2(ilevel)*h0/100.,&
                & xoff3(ilevel)*h0/100.
        endif
     else
        write(*,'(" myid=",I4," n1=",I4," n2=",I4," n3=",I4)') &
             & myid,n1(ilevel),n2(ilevel),n3(ilevel)
        write(*,'(" myid=",I4," xoff=",1pe10.3," yoff=",1pe10.3," zoff=",&
             & 1pe10.3," h-1 Mpc")') &
             & myid,xoff1(ilevel)*h0/100.,&
             & xoff2(ilevel)*h0/100.,&
             & xoff3(ilevel)*h0/100.
     endif
  end do

  ! Scale displacement in Mpc to code velocity (v=dx/dtau)
  ! in coarse cell units per conformal time
  vfact(1)=aexp*fpeebl(aexp)*sqrt(omega_m/aexp+omega_l*aexp*aexp+omega_k)
  ! This scale factor is different from vfact in grafic by h0/aexp

contains

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function fy(a)
    implicit none
    !      Computes the integrand
    real(dp)::fy
    real(dp)::y,a

    y=omega_m*(1d0/a-1d0) + omega_l*(a*a-1d0) + 1d0
    fy=1d0/y**1.5d0

    return
  end function fy
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function d1a(a)
    implicit none
    real(dp)::d1a
    !     Computes the linear growing mode D1 in a Friedmann-Robertson-Walker
    !     universe. See Peebles LSSU sections 11 and 14.
    real(dp)::a,y12,y,eps

    eps=1.0d-6
    if(a .le. 0.0d0)then
       write(*,*)'a=',a
       call clean_stop
    end if
    y=omega_m*(1d0/a-1d0) + omega_l*(a*a-1d0) + 1d0
    if(y .lt. 0.0D0)then
       write(*,*)'y=',y
       call clean_stop
    end if
    y12=y**0.5d0
    d1a=y12/a*rombint(eps,a,eps)

    return
  end function d1a
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!$  function ad1(d1)
!!$    implicit none
!!$    real(dp)::ad1
!!$    real(dp)::a,d1,da
!!$    integer::niter
!!$    ! Inverts the relation d1(a) given by function d1a(a) using
!!$    ! Newton-Raphson.
!!$    if (d1.eq.0.0) stop 'ad1 undefined for d1=0!'
!!$    ! Initial guess for Newton-Raphson iteration, good for Omega near 1.
!!$    a=1.e-7
!!$    niter=0
!!$10  niter=niter+1
!!$    da=(d1/d1a(a)-1d0)/fpeebl(a)*a
!!$    a=a+da
!!$    if (abs(da).gt.1.0e-8.and.niter.lt.10) go to 10
!!$    ad1=a
!!$    return
!!$  end function ad1
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function fpeebl(a)
    implicit none
    real(dp) :: fpeebl,a
    !     Computes the growth factor f=d\log D1/d\log a.
    real(dp) :: fact,y,eps

    eps=1.0d-6
    y=omega_m*(1d0/a-1d0) + omega_l*(a*a-1d0) + 1d0
    fact=rombint(eps,a,eps)
    fpeebl=(omega_l*a*a-0.5d0*omega_m/a)/y - 1d0 + a*fy(a)/fact
    return
  end function fpeebl
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function rombint(a,b,tol)
    implicit none
    real(dp)::rombint
    !
    !     Rombint returns the integral from a to b of f(x)dx using Romberg
    !     integration. The method converges provided that f(x) is continuous
    !     in (a,b). The function f must be double precision and must be
    !     declared external in the calling routine.
    !     tol indicates the desired relative accuracy in the integral.
    !
    integer::maxiter=16,maxj=5
    real(dp),dimension(100):: g
    real(dp)::a,b,tol,fourj
    real(dp)::h,error,gmax,g0,g1
    integer::nint,i,j,k,jmax

    h=0.5d0*(b-a)
    gmax=h*(fy(a)+fy(b))
    g(1)=gmax
    nint=1
    error=1.0d20
    i=0
10  i=i+1
    if(.not.  (i>maxiter.or.(i>5.and.abs(error)<tol)))then
       !  Calculate next trapezoidal rule approximation to integral.

       g0=0.0d0
       do k=1,nint
          g0=g0+fy(a+(k+k-1)*h)
       end do
       g0=0.5d0*g(1)+h*g0
       h=0.5d0*h
       nint=nint+nint
       jmax=min(i,maxj)
       fourj=1.0d0

       do j=1,jmax
          ! Use Richardson extrapolation.
          fourj=4.0d0*fourj
          g1=g0+(g0-g(j))/(fourj-1.0d0)
          g(j)=g0
          g0=g1
       enddo
       if (abs(g0).gt.tol) then
          error=1.0d0-gmax/g0
       else
          error=gmax
       end if
       gmax=g0
       g(jmax+1)=g0
       go to 10
    end if
    rombint=g0
    if (i>maxiter.and.abs(error)>tol) &
         &    write(*,*) 'Rombint failed to converge; integral, error=', &
         &    rombint,error
    return
  end function rombint

end subroutine init_cosmo

subroutine friedman(O_mat_0,O_vac_0,O_k_0,alpha,axp_min, &
     & axp_out,hexp_out,tau_out,t_out,ntable)
  use amr_parameters
  implicit none
  integer::ntable
  real(kind=8)::O_mat_0, O_vac_0, O_k_0
  real(kind=8)::alpha,axp_min
  real(dp),dimension(0:ntable)::axp_out,hexp_out,tau_out,t_out
  ! ######################################################!
  ! This subroutine assumes that axp = 1 at z = 0 (today) !
  ! and that t and tau = 0 at z = 0 (today).              !
  ! axp is the expansion factor, hexp the Hubble constant !
  ! defined as hexp=1/axp*daxp/dtau, tau the conformal    !
  ! time, and t the look-back time, both in unit of 1/H0. !
  ! alpha is the required accuracy and axp_min is the     !
  ! starting expansion factor of the look-up table.       !
  ! ntable is the required size of the look-up table.     !
  ! ######################################################!
  real(kind=8)::axp_tau, axp_t
  real(kind=8)::axp_tau_pre, axp_t_pre
  real(kind=8)::dadtau, dadt
  real(kind=8)::dtau,dt
  real(kind=8)::tau,t
  integer::nstep,nout,nskip

  if( (O_mat_0+O_vac_0+O_k_0) .ne. 1.0D0 )then
     write(*,*)'Error: non-physical cosmological constants'
     write(*,*)'O_mat_0,O_vac_0,O_k_0=',O_mat_0,O_vac_0,O_k_0
     write(*,*)'The sum must be equal to 1.0, but '
     write(*,*)'O_mat_0+O_vac_0+O_k_0=',O_mat_0+O_vac_0+O_k_0
     call clean_stop
  end if

  axp_tau = 1.0D0
  axp_t = 1.0D0
  tau = 0.0D0
  t = 0.0D0
  nstep = 0

  do while ( (axp_tau .ge. axp_min) .or. (axp_t .ge. axp_min) )

     nstep = nstep + 1
     dtau = alpha * axp_tau / dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)
     axp_tau_pre = axp_tau - dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)*dtau/2d0
     axp_tau = axp_tau - dadtau(axp_tau_pre,O_mat_0,O_vac_0,O_k_0)*dtau
     tau = tau - dtau

     dt = alpha * axp_t / dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
     axp_t_pre = axp_t - dadt(axp_t,O_mat_0,O_vac_0,O_k_0)*dt/2d0
     axp_t = axp_t - dadt(axp_t_pre,O_mat_0,O_vac_0,O_k_0)*dt
     t = t - dt

  end do

  if(debug)then
     write(*,666)-t
  end if
  666 format(' Age of the Universe (in unit of 1/H0)=',1pe10.3)

  nskip=nstep/ntable

  axp_t = 1d0
  t = 0d0
  axp_tau = 1d0
  tau = 0d0
  nstep = 0
  nout=0
  t_out(nout)=t
  tau_out(nout)=tau
  axp_out(nout)=axp_tau
  hexp_out(nout)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau

  do while ( (axp_tau .ge. axp_min) .or. (axp_t .ge. axp_min) )

     nstep = nstep + 1
     dtau = alpha * axp_tau / dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)
     axp_tau_pre = axp_tau - dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)*dtau/2d0
     axp_tau = axp_tau - dadtau(axp_tau_pre,O_mat_0,O_vac_0,O_k_0)*dtau
     tau = tau - dtau

     dt = alpha * axp_t / dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
     axp_t_pre = axp_t - dadt(axp_t,O_mat_0,O_vac_0,O_k_0)*dt/2d0
     axp_t = axp_t - dadt(axp_t_pre,O_mat_0,O_vac_0,O_k_0)*dt
     t = t - dt

     if(mod(nstep,nskip)==0)then
        nout=nout+1
        t_out(nout)=t
        tau_out(nout)=tau
        axp_out(nout)=axp_tau
        hexp_out(nout)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau
     end if

  end do
  t_out(ntable)=t
  tau_out(ntable)=tau
  axp_out(ntable)=axp_tau
  hexp_out(ntable)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau

end subroutine friedman

function dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)
  use amr_parameters
  real(kind=8)::dadtau,axp_tau,O_mat_0,O_vac_0,O_k_0
  dadtau = axp_tau*axp_tau*axp_tau *  &
       &   ( O_mat_0 + &
       &     O_vac_0 * axp_tau*axp_tau*axp_tau + &
       &     O_k_0   * axp_tau )
  dadtau = sqrt(dadtau)
  return
end function dadtau

function dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
  use amr_parameters
  real(kind=8)::dadt,axp_t,O_mat_0,O_vac_0,O_k_0
  dadt   = (1.0D0/axp_t)* &
       &   ( O_mat_0 + &
       &     O_vac_0 * axp_t*axp_t*axp_t + &
       &     O_k_0   * axp_t )
  dadt = sqrt(dadt)
  return
end function dadt
