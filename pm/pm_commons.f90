module pm_commons
  use amr_parameters
  use pm_parameters
  use random
  ! Sink particle related arrays
  real(dp),allocatable,dimension(:)::msink,c2sink,oksink_new,oksink_all
  real(dp),allocatable,dimension(:)::tsink,tsink_new,tsink_all
  real(dp),allocatable,dimension(:)::msink_new,msink_all
  real(dp),allocatable,dimension(:)::mseed,mseed_new,mseed_all
  real(dp),allocatable,dimension(:)::xmsink
  real(dp),allocatable,dimension(:)::dMsink_overdt,dMBHoverdt
  real(dp),allocatable,dimension(:)::rho_gas,volume_gas,eps_sink
  real(dp),allocatable,dimension(:,:)::vel_gas
  real(dp),allocatable,dimension(:)::delta_mass,delta_mass_new,delta_mass_all
  real(dp),allocatable,dimension(:)::wden,weth,wvol,wdiv,wden_new,weth_new,wvol_new,wdiv_new
  real(dp),allocatable,dimension(:,:)::wmom,wmom_new
  real(dp),allocatable,dimension(:,:)::vsink,vsink_new,vsink_all
  real(dp),allocatable,dimension(:,:)::fsink,fsink_new,fsink_all
  real(dp),allocatable,dimension(:,:,:)::vsnew,vsold
  real(dp),allocatable,dimension(:,:,:)::fsink_partial,sink_jump
  real(dp),allocatable,dimension(:,:)::lsink,lsink_new,lsink_all!sink angular momentum
  real(dp),allocatable,dimension(:,:)::xsink,xsink_new,xsink_all
  real(dp),allocatable,dimension(:,:)::weighted_density,weighted_volume,weighted_ethermal,weighted_divergence
  real(dp),allocatable,dimension(:,:,:)::weighted_momentum
  real(dp),allocatable,dimension(:)::dt_acc                ! maximum timestep allowed by the sink
  real(dp),allocatable,dimension(:)::rho_sink_tff
  integer,allocatable,dimension(:)::idsink,idsink_new,idsink_old,idsink_all
  logical,allocatable,dimension(:,:)::level_sink,level_sink_new
  logical,allocatable,dimension(:)::ok_blast_agn,ok_blast_agn_all,direct_force_sink
  logical,allocatable,dimension(:)::new_born,new_born_all,new_born_new
  integer,allocatable,dimension(:)::idsink_sort
  integer::ncloud_sink,ncloud_sink_massive
  integer::nindsink=0
  integer::sinkint_level=0         ! maximum level currently active is where the global sink variables are updated
  real(dp)::ssoft                  ! sink softening lenght in code units


  ! Particles related arrays
  real(dp),allocatable,dimension(:,:)::xp       ! Positions
  real(dp),allocatable,dimension(:,:)::vp       ! Velocities
  real(dp),allocatable,dimension(:)  ::mp       ! Masses
#ifdef OUTPUT_PARTICLE_POTENTIAL
  real(dp),allocatable,dimension(:)  ::ptcl_phi ! Potential of particle added by AP for output purposes 
#endif
  real(dp),allocatable,dimension(:)  ::tp       ! Birth epoch
  real(dp),allocatable,dimension(:,:)::weightp  ! weight of cloud parts for sink accretion only
  real(dp),allocatable,dimension(:)  ::zp       ! Birth metallicity
  integer ,allocatable,dimension(:)  ::nextp    ! Next particle in list
  integer ,allocatable,dimension(:)  ::prevp    ! Previous particle in list
  integer ,allocatable,dimension(:)  ::levelp   ! Current level of particle
  integer(i8b),allocatable,dimension(:)::idp    ! Identity of particle
  ! Tree related arrays
  integer ,allocatable,dimension(:)  ::headp    ! Head particle in grid
  integer ,allocatable,dimension(:)  ::tailp    ! Tail particle in grid
  integer ,allocatable,dimension(:)  ::numbp    ! Number of particles in grid
  ! Global particle linked lists
  integer::headp_free,tailp_free,numbp_free=0,numbp_free_tot=0
  ! Local and current seed for random number generator
  integer,dimension(IRandNumSize) :: localseed=-1


  ! Add particle type
  integer(kind=2), allocatable, dimension(:) :: typep  ! Particle type
  integer, parameter :: TYPE_DM=0
  integer, parameter :: TYPE_STAR=100, TYPE_DEAD_STAR=101
  integer, parameter :: TYPE_SINK=200
  integer, parameter :: TYPE_TRACER_DM=1000, TYPE_TRACER_STAR=1100, TYPE_TRACER_DEAD_STAR=1101, TYPE_TRACER_SINK=1200

contains
  function cross(a,b)
    use amr_parameters, only:dp
    real(dp),dimension(1:3)::a,b
    real(dp),dimension(1:3)::cross
    !computes the cross product c= a x b
    cross(1)=a(2)*b(3)-a(3)*b(2)
    cross(2)=a(3)*b(1)-a(1)*b(3)
    cross(3)=a(1)*b(2)-a(2)*b(1)
  end function cross

  logical pure function is_DM(typep)
    integer(kind=2), intent(in) :: typep
    is_DM = (typep >= 0) .and. (typep < 100)
  end function is_DM

  logical pure function is_star(typep)
    integer(kind=2), intent(in) :: typep
    is_star = (typep >= 100) .and. (typep < 200)
  end function is_star

  logical pure function is_sink(typep)
    integer(kind=2), intent(in) :: typep

    is_sink = (typep >= 200) .and. (typep < 300)
  end function is_sink

  logical pure function is_tracer(typep)
    integer(kind=2), intent(in) :: typep
    is_tracer = is_particle_tracer(typep) .or. is_gas_tracer(typep)
  end function is_tracer

  logical pure function is_particle_tracer(typep)
    integer(kind=2), intent(in) :: typep
    is_particle_tracer = (typep >= 1000) .and. (typep < 2000)
  end function is_particle_tracer

  logical pure function is_gas_tracer(typep)
    integer(kind=2), intent(in) :: typep
    is_gas_tracer = (typep < 0)
  end function is_gas_tracer

end module pm_commons
