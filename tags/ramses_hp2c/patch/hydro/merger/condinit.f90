!=============================================================================
! hydro ic for a Toomre/exponential radial density profile merger (03/03/09)
!=============================================================================
subroutine condinit(x,u,dx,nn)
  use amr_commons
  use hydro_commons
  use poisson_parameters, ONLY: gravity_params
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:ndim+1): d.u,d.v,d.w and U(i,ndim+2): E.
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:ndim+1):u,v,w and Q(i,ndim+2): P.
  ! If nvar >= ndim+3, remaining variables are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
  ! Galactic merger IC
  real(dp), dimension(3)::gal_center1 = 0.0D0
  real(dp), dimension(3)::gal_center2 = 0.0D0
  real(dp)::Mgaz_disk1 = 1.0D2
  real(dp)::Mgaz_disk2 = 1.0D2
  real(dp)::typ_radius1 = 3.0D0
  real(dp)::typ_radius2 = 3.0D0
  real(dp)::cut_radius1 = 10.0D0
  real(dp)::cut_radius2 = 10.0D0
  real(dp)::typ_height1 = 1.5D2
  real(dp)::typ_height2 = 1.5D2
  real(dp)::cut_height1 = 4.0D2
  real(dp)::cut_height2 = 4.0D2
  character(len=16)::rad_profile='exponential'
  real(dp)::IG_density_factor = 1.0D-5
  character(len=256)::Vcirc_dat_file1=''
  character(len=256)::Vcirc_dat_file2=''
  real(dp), dimension(3)::gal_axis1= (/ 0.0D0, 0.0D0, 1.0D0 /)
  real(dp), dimension(3)::gal_axis2= (/ 0.0D0, 0.0D0, 1.0D0 /)
  real(dp), dimension(3)::Vgal1 = 0.0D0
  real(dp), dimension(3)::Vgal2 = 0.0D0
  
  real(dp), dimension(:,:), allocatable,save::Vcirc_dat1, Vcirc_dat2


  integer::ivar,i, ind_gal
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables
  real(dp)::v,M,rho,dzz,zint, HH, rdisk, dpdr,dmax
  real(dp)::r, rr, rr1, rr2, abs_z, hcar1, hcar2, hcut1, hcut2
  real(dp), dimension(3)::vgal, axe_rot, xx1, xx2, xx, xx_rad, xc1, xc2, vg1, vg2
  real(dp)::rgal, sum,sum2,dmin,zmin,zmax,pi,tol
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2,M_b,az,eps
  real(dp)::rmin,rmax,a2,aa,Vcirc, HH_max, rcar1, rcar2, rcut1, rcut2
  real(dp)::rho_0_1, rho_0_2, rho_0, weight, da1, Vrot
  logical, save:: init_nml=.false.

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Read user-defined merger parameters from the namelist
  if (.not. init_nml) then
      call read_merger_params()
      init_nml = .true.
  end if

  ! Gaseous disks masses (given in  1.0E9 Msun unit)
  ! !!!!!! The galaxy #1 must be the more heavy !!!!!
  Mgaz_disk1 = Mgaz_disk1 * 1.0D9 * 1.9891D33  / (scale_d * scale_l**3)
  Mgaz_disk2 = Mgaz_disk2 * 1.0D9 * 1.9891D33  / (scale_d * scale_l**3)
  ! Galactic centers (given in kpc unit)
  xc1 = gal_center1 * 3.085677581282D21 / scale_l
  xc2 = gal_center2 * 3.085677581282D21 / scale_l
  ! Galaxy global velocities (given in km/s unit)
  vg1 = Vgal1 * 1.0D5 / scale_v
  vg2 = Vgal2 * 1.0D5 / scale_v
  ! Gaseous disk typical radii (a) (given in kpc unit)
  rcar1 = typ_radius1 * 3.085677581282D21 / scale_l
  rcar2 = typ_radius2 * 3.085677581282D21 / scale_l
  ! Gaseous disk max radii (given in kpc unit)
  rcut1 = cut_radius1 * 3.085677581282D21 / scale_l
  rcut2 = cut_radius2 * 3.085677581282D21 / scale_l
  ! Gaseous disk typical thicknesses (h) (given in pc unit)
  hcar1 = typ_height1 * 3.085677581282D18 / scale_l
  hcar2 = typ_height2 * 3.085677581282D18 / scale_l
  ! Gaseous disk max thicknesses(zmax) (given in pc unit)
  hcut1 = cut_height1 * 3.085677581282D18 / scale_l
  hcut2 = cut_height2 * 3.085677581282D18 / scale_l


  if(maxval(abs(xc1)) .GT. (boxlen/2.0D0))then
     write(*,*)'Error: galactic center (1) coordinates must be in the box [', (-boxlen/2.0D0), ';', (boxlen/2.0D0), ']^3'
     call clean_stop
  endif
  if((maxval(abs(xc2)) .GT. (boxlen/2.0D0)) .AND. (Mgaz_disk2 .NE. 0.0D0))then
     write(*,*)'Error: galactic center (2) coordinates must be in the box [', (-boxlen/2.0D0), ';', (boxlen/2.0D0), ']^3'
     call clean_stop
  endif

  a2=T2_star / scale_T2 ! sound speed squared
  aa=sqrt(a2)

  ! Galactic central gas densities
  pi=dacos(-1.0D0)
  select case (rad_profile)
      case ('exponential')
          rho_0_1 = (1.0D0 - exp(-rcut1 / rcar1) * (1.0D0 + rcut1 / rcar1) ) * (1.0D0 - exp(-hcut1 / hcar1))
          rho_0_1 = Mgaz_disk1 / (4.0D0 * pi * rcar1**2 * hcar1 * rho_0_1)
          rho_0_2 = (1.0D0 - exp(-rcut2 / rcar2) * (1.0D0 + rcut2 / rcar2) ) * (1.0D0 - exp(-hcut2 / hcar2)) 
          rho_0_2 = Mgaz_disk2 / (4.0D0 * pi * rcar2**2 * hcar2 * rho_0_2)
      case ('Toomre')
          rho_0_1 = (sqrt(1.0D0 + rcut1**2/rcar1**2) - 1.0D0) * (1.0D0 - exp(-hcut1 / hcar1))
          rho_0_1 = Mgaz_disk1 / (4.0D0 * pi * rcar1**2 * hcar1 * rho_0_1)
          rho_0_2 = (sqrt(1.0D0 + rcut2**2/rcar2**2) - 1.0D0) * (1.0D0 - exp(-hcut2 / hcar2)) 
          rho_0_2 = Mgaz_disk2 / (4.0D0 * pi * rcar2**2 * hcar2 * rho_0_2)
  end select

  ! Intergalactic gas density
  select case (rad_profile)
      case ('exponential')
          dmin = rho_0_1 * exp(-hcut1 / hcar1) * exp(-rcut1 / rcar1)
      case ('Toomre')
          dmin = rho_0_1 * exp(-hcut1 / hcar1) / sqrt(1.0D0 + rcut1**2/rcar1**2)
  end select
  if(Mgaz_disk2 .NE. 0.0D0) then
      select case (rad_profile)
          case ('exponential')
              dmin = min(dmin, rho_0_2 * exp(-hcut2 / hcar2) * exp(-rcut2 / rcar2))
          case ('Toomre')
              dmin = min(dmin, rho_0_2 * exp(-hcut2 / hcar2) / sqrt(1.0D0 + rcut2**2/rcar2**2))
      end select
  end if
  dmin = IG_density_factor * dmin

  
  ! Loop over cells
  do i=1,nn
     xx1=x(i,:)-(xc1(:)+boxlen/2.0D0)
     xx2=x(i,:)-(xc2(:)+boxlen/2.0D0)

    ! Compute angular velocity

    ! Distance between cell and both galactic centers
     rr1 = norm2(xx1)
     rr2 = norm2(xx2)

    ! Projected cell position over galactic centers axis
     da1 = dot_product(xc2 - xc1, xx1) / norm2(xc2 - xc1)**2

     if(da1 .LT. (Mgaz_disk1 / (Mgaz_disk1 + Mgaz_disk2))) then ! cell belongs to galaxy #1
         ind_gal = 1
         rr = rr1
         xx = xx1
         axe_rot = gal_axis1
         vgal = vg1
         rgal = rcar1
         rdisk = rcut1
         HH = hcar1
         HH_max = hcut1
         rho_0 = rho_0_1
     else ! cell belongs to galaxy #2
         ind_gal = 2
         rr = rr2
         xx = xx2
         axe_rot = gal_axis2
         vgal = vg2
         rgal = rcar2
         rdisk = rcut2
         HH = hcar2
         HH_max = hcut2
         rho_0 = rho_0_2
     end if

     ! Cylindric radius : distance between the cell and the galactic rotation axis
     xx_rad = xx - dot_product(xx,axe_rot) * axe_rot
     r = norm2(xx_rad)

     ! vertical position absolute value
     abs_z = sqrt(rr**2 - r**2)

     if(((r-dx/2.0D0).lt.rdisk) .and. ((abs_z-dx/2.0D0) .lt. HH_max))then ! Cell in the disk : analytical density profile + rotation velocity
        weight = (min(r+dx/2.0D0,rdisk)-(r-dx/2.0D0))/dx
        if (weight .NE. 1.0D0) then
            r = r + (weight-1.0D0)*dx/2.0D0
        end if
        ! Circular velocity
        Vcirc= find_Vcirc(r, ind_gal)
        weight = weight*(min(abs_z+dx/2.0D0,HH_max)-(abs_z-dx/2.0D0))/dx
        ! Density
        select case (rad_profile)
            case ('exponential')
                q(i,1)= rho_0 * exp(-abs_z / HH) * exp(-r / rgal)
            case ('Toomre')
                q(i,1)= rho_0 * exp(-abs_z / HH) / sqrt(1.0D0 + r**2/rgal**2)
        end select
        q(i,1) = max(weight * q(i,1), dmin)
        ! P = rho * a² = rho * Cs²
        q(i,ndim+2)=a2*q(i,1)
        ! V = Vrot * (u_rot^xx_rad)/r + Vx_gal        
        !  -> Vrot = sqrt(Vcirc**2 - 3*Cs² + r/rho * grad(rho) * Cs²)
        select case (rad_profile)
            case ('exponential')
                Vrot = sqrt(max(Vcirc**2 - 3.0D0*a2 - r/rgal * a2,0.0D0))
            case ('Toomre')
                Vrot = sqrt(max(Vcirc**2 - 3.0D0*a2 - r**2/(r**2+rgal**2)*a2,0.0D0))
        end select
        Vrot = weight * Vrot
        q(i,ndim-1:ndim+1) = Vrot * vect_prod(axe_rot,xx_rad)/r + vgal
    else ! Cell out of the gaseous disk : density = peanut and velocity = V_gal
        q(i,1)=dmin
        q(i,ndim+2)=a2*q(i,1)
        ! V = Vgal
        q(i,ndim-1:ndim+1)= vgal
     endif
  enddo


!  deallocate(Vcirc_dat1, Vcirc_dat2)


  ! Convert primitive to conservative variables
  ! density -> density
  u(1:nn,1)=q(1:nn,1)
  ! velocity -> momentum : Omega = rho * V
  u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
#if NDIM>1
  u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
#endif
#if NDIM>2
  u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
#endif

  ! kinetic energy
  ! Total system global velocity : 0
  u(1:nn,ndim+2)=0.0D0
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5D0*q(1:nn,1)*q(1:nn,2)**2
#if NDIM>1
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5D0*q(1:nn,1)*q(1:nn,3)**2
#endif
#if NDIM>2
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5D0*q(1:nn,1)*q(1:nn,4)**2
#endif
  ! pressure -> total fluid energy
  ! E = Ec + P / (gamma - 1)
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+q(1:nn,ndim+2)/(gamma-1.0d0)
  ! passive scalars
  do ivar=ndim+3,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do


contains

  subroutine read_merger_params
    implicit none
#ifndef WITHOUTMPI
    include 'mpif.h'
#endif
    logical::nml_ok=.true.
    character(LEN=80)::infile
    
    !--------------------------------------------------
    ! Local variables  
    !--------------------------------------------------
    real(dp)::norm_u
    logical::vcirc_file1_exists, vcirc_file2_exists
    
    !--------------------------------------------------
    ! Namelist definitions
    !--------------------------------------------------
    namelist/merger_params/ IG_density_factor &
         & ,gal_center1, gal_center2, Mgaz_disk1, Mgaz_disk2 &
         & ,typ_radius1, typ_radius2, cut_radius1, cut_radius2 &
         & ,typ_height1, typ_height2, cut_height1, cut_height2 &
         & ,rad_profile, Vcirc_dat_file1, Vcirc_dat_file2 &
         & ,gal_axis1, gal_axis2, Vgal1, Vgal2
    
    
    CALL getarg(1,infile)
    open(1,file=infile)
    rewind(1)
    read(1,NML=merger_params,END=106)
    goto 107
106 write(*,*)' You need to set up namelist &MERGER_PARAMS in parameter file'
    call clean_stop
107 continue
    close(1)
    
    !-------------------------------------------------
    ! This section deals with the galactic merger initial conditions
    !-------------------------------------------------
    select case (rad_profile)
    case ('Toomre')
       if(myid==1) write(*,*) "Chosen hydro radial density profile :'Toomre'"
    case ('exponential')
       if(myid==1) write(*,*) "Chosen hydro radial density profile :'exponential'"
    case default
       if(myid==1) write(*,*) "Chosen hydro radial density profile :'exponential'"
       rad_profile='exponential'
    end select
    if(Mgaz_disk2 .GT. Mgaz_disk1)then
       if(myid==1)write(*,*)'Error: The galaxy #1 must be bigger than #2'
       nml_ok=.false.
    endif
    inquire(file=trim(Vcirc_dat_file1), exist=vcirc_file1_exists)
    if(.NOT. vcirc_file1_exists) then
       if(myid==1)write(*,*)'Error: Vcirc_dat_file1 ''', trim(Vcirc_dat_file1), ''' doesn''t exist '
       nml_ok=.false.
    end if
    inquire(file=trim(Vcirc_dat_file2), exist=vcirc_file2_exists)
    if(.NOT. vcirc_file2_exists) then
       if(myid==1)write(*,*)'Error: Vcirc_dat_file2 ''', trim(Vcirc_dat_file2), ''' doesn''t exist '
       nml_ok=.false.
    end if
    if ((vcirc_file1_exists == .true.) .AND. (vcirc_file2_exists == .true.)) then
       write(*,*)'Reading vcirc files'
       call read_vcirc_files()
    end if
    
    norm_u = sqrt(gal_axis1(1)**2 + gal_axis1(2)**2 + gal_axis1(3)**2)
    if(norm_u .EQ. 0.0D0) then
       if(myid==1)write(*,*)'Error: Galactic axis(1) is zero '
       nml_ok=.false.
    else
       if(norm_u .NE. 1.0D0) then
          gal_axis1 = gal_axis1 / norm_u
       end if
    end if
    norm_u = sqrt(gal_axis2(1)**2 + gal_axis2(2)**2 + gal_axis2(3)**2)
    if(norm_u .EQ. 0.0D0) then
       if(myid==1)write(*,*)'Error: Galactic axis(2) is zero '
       nml_ok=.false.
    else
       if(norm_u .NE. 1.0D0) then
          gal_axis2 = gal_axis2 / norm_u
       end if
    end if
    
    if(.not. nml_ok)then
       if(myid==1)write(*,*)'Too many errors in the namelist'
       if(myid==1)write(*,*)'Aborting...'
       call clean_stop
    end if
    
  end subroutine read_merger_params
  
  function find_Vcirc(rayon, indice)
    implicit none
    real(dp), intent(in)		:: rayon
    integer, intent(in)		:: indice
    real(dp)					:: find_Vcirc
    real(dp)					:: vitesse, rayon_bin, vitesse_old, rayon_bin_old
    integer					:: k, indmax
    
    k=2
    if (indice .EQ. 1) then
       indmax = size(Vcirc_dat1,1)
       rayon_bin = Vcirc_dat1(k,1)
       rayon_bin_old = Vcirc_dat1(k-1,1)
       vitesse = Vcirc_dat1(k,2)
       vitesse_old = Vcirc_dat1(k-1,2)
    else
       indmax = size(Vcirc_dat2,1)
       rayon_bin = Vcirc_dat2(k,1)
       rayon_bin_old = Vcirc_dat2(k-1,1)
       vitesse = Vcirc_dat2(k,2)
       vitesse_old = Vcirc_dat2(k-1,2)
    end if
    do while (rayon .GT. rayon_bin)
       if(k .GE. indmax) then
          write(*,*) "Hydro IC error : Radius out of rotation curve !!!"
          call clean_stop
       end if
       k = k + 1
       if (indice .EQ. 1) then
          vitesse_old = vitesse
          vitesse = Vcirc_dat1(k,2)
          rayon_bin_old = rayon_bin
          rayon_bin = Vcirc_dat1(k,1)
       else
          vitesse_old = vitesse
          vitesse = Vcirc_dat2(k,2)
          rayon_bin_old = rayon_bin
          rayon_bin = Vcirc_dat2(k,1)
       end if
    end do
    
    find_Vcirc = vitesse_old + (rayon - rayon_bin_old) * (vitesse - vitesse_old) / (rayon_bin - rayon_bin_old)
    
    return
    
  end function find_Vcirc
  
  
  function vect_prod(a,b)
    implicit none
    real(dp), dimension(3), intent(in)::a,b
    real(dp), dimension(3)::vect_prod
    
    vect_prod(1) = a(2) * b(3) - a(3) * b(2)
    vect_prod(2) = a(3) * b(1) - a(1) * b(3)
    vect_prod(3) = a(1) * b(2) - a(2) * b(1)
    
  end function vect_prod
  
  
  function norm2(x)
    implicit none
    real(dp), dimension(3), intent(in)::x
    real(dp) :: norm2
    
    norm2 = sqrt(dot_product(x,x))
    
  end function norm2
  
  !------------------------------------------------------------------------------------- 
  ! Circular velocity files reading
  ! {{{
  subroutine read_vcirc_files
    implicit none
    integer:: nvitesses, ierr, i
    real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
    
    call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ WARNING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    ! Those circular velocity files must contain :                              !
    ! - Column #1 : radius (in pc)                                              !
    ! - Column #2 : circular velocity (in km/s)                                 !
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    
    ! Galaxy #1
    nvitesses = 0
    open(unit=123, file=trim(Vcirc_dat_file1), iostat=ierr)
    do while (ierr==0)
       read(123,*,iostat=ierr)
       if(ierr==0) then
          nvitesses = nvitesses + 1  ! Number of samples
       end if
    end do
    allocate(Vcirc_dat1(nvitesses,2))
    Vcirc_dat1 = 0.0D0
    rewind(123)
    do i=1,nvitesses
       read(123,*) Vcirc_dat1(i,:)
    end do
    close(123)
    ! Unit conversion pc -> code unit and km/s -> code unit
    Vcirc_dat1(:,1) = Vcirc_dat1(:,1) * 3.085677581282D18 / scale_l
    Vcirc_dat1(:,2) = Vcirc_dat1(:,2) * 1.0D5 / scale_v
    
    ! Galaxy #2
    nvitesses = 0
    open(unit=123, file=trim(Vcirc_dat_file2), iostat=ierr)
    do while (ierr==0)
       read(123,*,iostat=ierr)
       if(ierr==0) then
          nvitesses = nvitesses + 1 ! Number of samples
       end if
    end do
    allocate(Vcirc_dat2(nvitesses,2))
    Vcirc_dat2 = 0.0D0
    rewind(123)
    do i=1,nvitesses
       read(123,*) Vcirc_dat2(i,:)
    end do
    close(123)
    ! Unit conversion pc -> code unit and km/s -> code unit
    Vcirc_dat2(:,1) = Vcirc_dat2(:,1) * 3.085677581282D18 / scale_l
    Vcirc_dat2(:,2) = Vcirc_dat2(:,2) * 1.0D5 / scale_v
    
    
  end subroutine read_vcirc_files
  ! }}}
  !--------------------------------------------------------------------------------------
end subroutine condinit


