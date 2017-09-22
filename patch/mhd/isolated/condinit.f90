!==================================================================================
!=== Hydro ic for Toomre/exponential radial density profile galaxies.           ===
!=== Sets up a galactic binary merger.                                          ===
!=== Author : D.Chapon                                                          ===
!=== Date : 2010/09/01                                                          ===
!==================================================================================
module merger_parameters!{{{
  use amr_commons
  !--------------------!
  ! Galactic merger IC !
  !--------------------!

  ! Gas disk masses, given in GMsun in namelist,
  ! then converted in user unit.
  ! !!!!!! The galaxy #1 must be the heaviest !!!!!
  real(dp)::Mgas_disk1 = 1.0D2
  real(dp)::Mgas_disk2 = 1.0D2 ! Set to 0.0 for isolated galaxy runs
  ! Galactic centers, 0-centered, given in user unit
  real(dp), dimension(3)::gal_center1 = 0.0D0
  ! Set gal_center2 to values larger than 2*Lbox for isolated galaxy runs
  real(dp), dimension(3)::gal_center2 = 0.0D0
  ! Galactic disks rotation axis
  real(dp), dimension(3)::gal_axis1= (/ 0.0D0, 0.0D0, 1.0D0 /)
  real(dp), dimension(3)::gal_axis2= (/ 0.0D0, 0.0D0, 1.0D0 /)
  ! Particle ic ascii file names for the galaxies.
  !~~~~~~~~~~~~~~~~~~~~~~~WARNING~~~~~~~~~~~~~~~~~~~~~~~!
  ! Assumed to be J=z-axis, 0-centered galaxy ic files. !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  character(len=256)::ic_part_file_gal1='ic_part1'
  character(len=256)::ic_part_file_gal2='ic_part2'
  ! Rotation curve files for the galaxies
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ WARNING ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  ! Those circular velocity files must contain :                              !
  ! - Column #1 : radius (in pc)                                              !
  ! - Column #2 : circular velocity (in km/s)                                 !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  character(len=256)::Vcirc_dat_file1='Vcirc1.dat'
  character(len=256)::Vcirc_dat_file2='Vcirc2.dat'
  ! Galactic global velocities, given in km/s in namelist,
  ! then converted in user unit.
  real(dp), dimension(3)::Vgal1 = 0.0D0
  real(dp), dimension(3)::Vgal2 = 0.0D0
  ! Gas disk typical/truncation radius/height, given in kpc in namelist,
  ! then converted in user unit.
  real(dp)::typ_radius1 = 3.0D0
  real(dp)::typ_radius2 = 3.0D0
  real(dp)::cut_radius1 = 10.0D0
  real(dp)::cut_radius2 = 10.0D0
  real(dp)::typ_height1 = 1.5D-1
  real(dp)::typ_height2 = 1.5D-1
  real(dp)::cut_height1 = 4.0D-1
  real(dp)::cut_height2 = 4.0D-1
  ! Gas density profile : radial =>'exponential' (default) or 'Toomre'
  !                       vertical => 'exponential' or 'gaussian'
  character(len=16)::rad_profile='exponential'
  character(len=16)::z_profile='gaussian'
  ! Inter galactic gas density contrast factor
  real(dp)::IG_density_factor = 1.0D-5

  ! Magnetic Field Setup
  real(dp)::typ_mag_radius = 3.0D0           ! radial scale length
  real(dp)::typ_mag_height = 1.5D-1          ! vertical scale length
  character(len=16)::mag_topology='constant' ! magnetic topology:
                                             ! 'constant' (along x-axis)
                                             ! 'toroidal'
                                             ! 'dipole'
                                             ! 'quadrupole'

  !~~~~~~~~~~~~~~~~~~~~~~~~ NOTE ~~~~~~~~~~~~~~~~~~~~~~~~!
  ! For isolated galaxy runs :                           !
  ! --------------------------                           !
  ! - set 'Mgas_disk2' to 0.0                            !
  ! - set 'gal_center2' to values larger than Lbox*2     !
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
end module merger_parameters!}}}

module merger_commons

  use merger_parameters
  real(dp), dimension(:,:), allocatable::Vcirc_dat1, Vcirc_dat2

end module merger_commons

subroutine read_merger_params! {{{
  use merger_parameters
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
  logical::ic_part_file1_exists, ic_part_file2_exists
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2

  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------
  namelist/merger_params/ IG_density_factor &
       & ,gal_center1, gal_center2, Mgas_disk1, Mgas_disk2 &
       & ,typ_radius1, typ_radius2, cut_radius1, cut_radius2 &
       & ,typ_height1, typ_height2, cut_height1, cut_height2 &
       & ,rad_profile, z_profile, Vcirc_dat_file1, Vcirc_dat_file2 &
       & ,ic_part_file_gal1, ic_part_file_gal2 &
       & ,gal_axis1, gal_axis2, Vgal1, Vgal2 &
       & ,typ_mag_radius, typ_mag_height, mag_topology


  CALL getarg(1,infile)
  open(1,file=infile)
  rewind(1)
  read(1,NML=merger_params,END=106)
  goto 107
106 write(*,*)' You need to set up namelist &MERGER_PARAMS in parameter file'
  call clean_stop
107 continue
  close(1)

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  !-------------------------------------------------
  ! This section deals with the galactic merger initial conditions
  !-------------------------------------------------
  ! Gaseous disks masses : GMsun => user unit
  Mgas_disk1 = Mgas_disk1 * 1.0D9 * 1.9891D33  / (scale_d * scale_l**3)
  Mgas_disk2 = Mgas_disk2 * 1.0D9 * 1.9891D33  / (scale_d * scale_l**3)
  ! Galaxy global velocities : km/s => user unit
  Vgal1 = Vgal1 * 1.0D5 / scale_v
  Vgal2 = Vgal2 * 1.0D5 / scale_v
  ! Gas disk typical radii (a) : kpc => user unit
  typ_radius1 = typ_radius1 * 3.085677581282D21 / scale_l
  typ_radius2 = typ_radius2 * 3.085677581282D21 / scale_l
  typ_mag_radius = typ_mag_radius * 3.085677581282D21 / scale_l
  ! Gas disk max radii : kpc => user unit
  cut_radius1 = cut_radius1 * 3.085677581282D21 / scale_l
  cut_radius2 = cut_radius2 * 3.085677581282D21 / scale_l
  ! Gas disk typical thicknesses (h) : kpc => user unit
  typ_height1 = typ_height1 * 3.085677581282D21 / scale_l
  typ_height2 = typ_height2 * 3.085677581282D21 / scale_l
  typ_mag_height = typ_mag_height * 3.085677581282D21 / scale_l
  ! Gas disk max thicknesses(zmax) : kpc => user unit
  cut_height1 = cut_height1 * 3.085677581282D21 / scale_l
  cut_height2 = cut_height2 * 3.085677581282D21 / scale_l

  select case (rad_profile)
      case ('Toomre')
          if(myid==1) write(*,*) "Chosen hydro radial density profile :'Toomre'"
      case ('exponential')
          if(myid==1) write(*,*) "Chosen hydro radial density profile :'exponential'"
      case default
          if(myid==1) write(*,*) "Chosen hydro radial density profile :'exponential'"
          rad_profile='exponential'
  end select
  select case (z_profile)
      case ('gaussian')
          if(myid==1) write(*,*) "Chosen hydro vertical density profile :'gaussian'"
      case ('exponential')
          if(myid==1) write(*,*) "Chosen hydro vertical density profile :'exponential'"
      case default
          if(myid==1) write(*,*) "Chosen hydro vertical density profile :'exponential'"
          z_profile='exponential'
  end select
  if(Mgas_disk2 .GT. Mgas_disk1)then
     if(myid==1)write(*,*)'Error: The galaxy #1 must be bigger than #2'
     nml_ok=.false.
  endif
  ! Check whether velocity radial profile files exist or not.
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
  if ((vcirc_file1_exists) .AND. (vcirc_file2_exists)) then
      call read_vcirc_files()
  end if

  ! Check whether ic_part files exist or not.
  inquire(file=trim(initfile(levelmin))//'/'//trim(ic_part_file_gal1), exist=ic_part_file1_exists)
  if(.NOT. ic_part_file1_exists) then
     if(myid==1)write(*,*)'Error: ic_part_file1 ''', trim(ic_part_file_gal1), ''' doesn''t exist in '''&
     & , trim(initfile(levelmin))
     nml_ok=.false.
  end if
  inquire(file=trim(initfile(levelmin))//'/'//trim(ic_part_file_gal2), exist=ic_part_file2_exists)
  if(.NOT. ic_part_file2_exists) then
     if(myid==1)write(*,*)'Error: ic_part_file2 ''', trim(ic_part_file_gal2), ''' doesn''t exist in '''&
     & , trim(initfile(levelmin))
     nml_ok=.false.
  end if

  ! galactic rotation axis
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

  ! Magnetic Topology
  select case (mag_topology)
      case ('constant')
          if(myid==1) write(*,*) "Chosen magnetic topology: 'constant'"
      case ('toroidal')
          if(myid==1) write(*,*) "Chosen magnetic topology: 'toroidal'"
      case ('dipole')
          if(myid==1) write(*,*) "Chosen magnetic topology: 'dipole'"
      case ('quadrupole')
          if(myid==1) write(*,*) "Chosen magnetic topology: 'quadrupole'"
      case default
          if(myid==1) write(*,*) "Magnetic field topology not recognised."
          nml_ok=.false.
  end select

  if(.not. nml_ok)then
     if(myid==1)write(*,*)'Too many errors in the namelist'
     if(myid==1)write(*,*)'Aborting...'
     call clean_stop
  end if

end subroutine read_merger_params
! }}}


subroutine condinit(x,u,dx,nn)
  use amr_parameters
  use hydro_parameters
  use merger_commons
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar+3)::u ! Conservative variables
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
  integer::ivar,i,ind_gal
  real(dp),dimension(1:nvector,1:nvar+3),save::q   ! Primitive variables
  real(dp)::v,M,rho,dzz,zint, HH, rdisk, dpdr,dmax
  real(dp)::r, rr, rr1, rr2, abs_z
  real(dp), dimension(3)::vgal, axe_rot, xx1, xx2, xx, xx_rad
  real(dp)::rgal, sum,sum2,dmin,dmin1,dmin2,zmin,zmax,pi,tol
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp)::M_b,az,eps
  real(dp)::rmin,rmax,a2,aa,Vcirc, HH_max
  real(dp)::rho_0_1, rho_0_2, rho_0, weight, da1, Vrot
  logical, save:: init_nml=.false.

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Read user-defined merger parameters from the namelist
  if (.not. init_nml) then
      call read_merger_params
      init_nml = .true.
  end if

  if(maxval(abs(gal_center1)) .GT. (boxlen/2.0D0))then
     write(*,*)'Error: galactic center (1) coordinates must be in the box [', &
     & (-boxlen/2.0D0), ';', (boxlen/2.0D0), ']^3'
     call clean_stop
  endif
  if((maxval(abs(gal_center2)) .GT. (boxlen/2.0D0)) .AND. (Mgas_disk2 .NE. 0.0D0))then
     write(*,*)'Error: galactic center (2) coordinates must be in the box [', &
     & (-boxlen/2.0D0), ';', (boxlen/2.0D0), ']^3'
     call clean_stop
  endif


!  a2=T2_star / scale_T2 ! sound speed squared
  a2=1d4 / scale_T2 ! sound speed squared
  aa=sqrt(a2)

  ! Galactic central gas densities
  pi=dacos(-1.0D0)
  select case (rad_profile)
      case ('exponential')
          rho_0_1 = 1.0D0 - exp(-cut_radius1 / typ_radius1) * (1.0D0 + cut_radius1 / typ_radius1)
          rho_0_2 = 1.0D0 - exp(-cut_radius2 / typ_radius2) * (1.0D0 + cut_radius2 / typ_radius2)
      case ('Toomre')
          rho_0_1 = sqrt(1.0D0 + cut_radius1**2/typ_radius1**2) - 1.0D0
          rho_0_2 = sqrt(1.0D0 + cut_radius2**2/typ_radius2**2) - 1.0D0
  end select
  select case (z_profile)
      case ('exponential')
          rho_0_1 = rho_0_1 * (1.0D0 - exp(-cut_height1 / typ_height1))
          rho_0_2 = rho_0_2 * (1.0D0 - exp(-cut_height2 / typ_height2))
      case ('gaussian')
          rho_0_1 = rho_0_1 * (dsqrt(pi/2.0D0) * erf(cut_height1 / (dsqrt(2.0D0)*typ_height1)))
          rho_0_2 = rho_0_2 * (dsqrt(pi/2.0D0) * erf(cut_height2 / (dsqrt(2.0D0)*typ_height2)))
  end select
  rho_0_1 = Mgas_disk1 / (4.0D0 * pi * typ_radius1**2 * typ_height1 * rho_0_1)
  rho_0_2 = Mgas_disk2 / (4.0D0 * pi * typ_radius2**2 * typ_height2 * rho_0_2)

  ! Intergalactic gas density
  select case (rad_profile)
      case ('exponential')
          dmin1 = rho_0_1 * exp(-cut_radius1 / typ_radius1)
          dmin2 = rho_0_2 * exp(-cut_radius2 / typ_radius2)
      case ('Toomre')
          dmin1 = rho_0_1 / sqrt(1.0D0 + cut_radius1**2/typ_radius1**2)
          dmin2 = rho_0_2 / sqrt(1.0D0 + cut_radius2**2/typ_radius2**2)
  end select
  select case (z_profile)
      case ('exponential')
          dmin1 = dmin1 * exp(-cut_height1 / typ_height1)
          dmin2 = dmin2 * exp(-cut_height2 / typ_height2)
      case ('gaussian')
          dmin1 = dmin1 * exp(-0.5D0 * (cut_height1 / typ_height1)**2)
          dmin2 = dmin2 * exp(-0.5D0 * (cut_height2 / typ_height2)**2)
  end select

  if(Mgas_disk2 .NE. 0.0D0) then
      dmin = min(dmin1, dmin2)
  else
      dmin = dmin1
  end if
  dmin = IG_density_factor * dmin


  ! Loop over cells
  do i=1,nn
     xx1=x(i,:)-(gal_center1(:)+boxlen/2.0D0)
     xx2=x(i,:)-(gal_center2(:)+boxlen/2.0D0)

    ! Compute angular velocity

    ! Distance between cell and both galactic centers
     rr1 = norm2(xx1)
     rr2 = norm2(xx2)

    ! Projected cell position over galactic centers axis
     da1 = dot_product(gal_center2 - gal_center1, xx1) / norm2(gal_center2 - gal_center1)**2

     if(da1 .LT. (Mgas_disk1 / (Mgas_disk1 + Mgas_disk2))) then ! cell belongs to galaxy #1
         ind_gal = 1
         rr = rr1
         xx = xx1
         axe_rot = gal_axis1
         vgal = Vgal1
         rgal = typ_radius1
         rdisk = cut_radius1
         HH = typ_height1
         HH_max = cut_height1
         rho_0 = rho_0_1
     else ! cell belongs to galaxy #2
         ind_gal = 2
         rr = rr2
         xx = xx2
         axe_rot = gal_axis2
         vgal = Vgal2
         rgal = typ_radius2
         rdisk = cut_radius2
         HH = typ_height2
         HH_max = cut_height2
         rho_0 = rho_0_2
     end if

     ! Cylindric radius : distance between the cell and the galactic rotation axis
     xx_rad = xx - dot_product(xx,axe_rot) * axe_rot
     r = norm2(xx_rad)

     ! vertical position absolute value
     abs_z = sqrt(rr**2 - r**2)

     if(((r-dx/2.0D0).lt.rdisk) .and. ((abs_z-dx/2.0D0) .lt. HH_max)) then
        ! Cell in the disk : analytical density profile + rotation velocity
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
                q(i,1)= rho_0 * exp(-r / rgal)
            case ('Toomre')
                q(i,1)= rho_0 / sqrt(1.0D0 + r**2/rgal**2)
        end select
        select case (z_profile)
            case ('exponential')
                q(i,1)= q(i,1) * exp(-abs_z / HH)
            case ('gaussian')
                q(i,1)= q(i,1) * exp(-0.5D0 * (abs_z / HH)**2)
        end select
        q(i,1) = max(weight * q(i,1), 1d-7/scale_nH)
        ! P=rho*cs^2
        q(i,5) = a2*q(i,1)
        ! V = Vrot * (u_rot^xx_rad)/r + Vx_gal
        !  -> Vrot = sqrt(Vcirc**2 - 3*Cs^2 + r/rho * grad(rho) * Cs)
        select case (rad_profile)
            case ('exponential')
                Vrot = sqrt(max(Vcirc**2 - 3.0D0*a2 - r/rgal * a2,0.0D0))
            case ('Toomre')
                Vrot = sqrt(max(Vcirc**2 - 3.0D0*a2 - r**2/(r**2+rgal**2)*a2,0.0D0))
        end select
        Vrot = weight * Vrot
        q(i,2:4) = Vrot * vect_prod(axe_rot,xx_rad)/r + vgal
        if(metal)q(i,imetal)=z_ave*0.02 ! z_ave is in solar units

    else ! Cell out of the gaseous disk : density = peanut and velocity = V_gal
        !q(i,1) = dmin
        !q(i,5) = 100.*a2*q(i,1)
        q(i,1)=1d-7/scale_nH
        q(i,ndim+2)=1d7/scale_T2*q(i,1)
        ! V = Vgal
        q(i,ndim-1:ndim+1)= vgal
        if(metal)q(i,imetal)=0.0

     endif
  enddo

  ! Magnetic Topology
  select case (mag_topology)
    case ('constant')
      call mag_constant(q,nn,B_ave)
    case ('toroidal')
      call mag_toroidal(x,q,dx,nn,B_ave,typ_mag_radius,typ_mag_height)
    case ('dipole')
      call mag_dipole(x,q,dx,nn,B_ave,typ_mag_radius,typ_mag_height)
    case ('quadrupole')
      call mag_quadrupole(x,q,dx,nn,B_ave,typ_mag_radius,typ_mag_height)
  end select


  ! Convert primitive to conservative variables
  ! density -> density
  u(1:nn,1)=q(1:nn,1)
  ! velocity -> momentum : Omega = rho * V
  u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
  u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
  u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
  ! kinetic energy
  ! Total system global velocity : 0
  u(1:nn,5)=0.0D0
  u(1:nn,5)=u(1:nn,5)+0.5D0*q(1:nn,1)*q(1:nn,2)**2
  u(1:nn,5)=u(1:nn,5)+0.5D0*q(1:nn,1)*q(1:nn,3)**2
  u(1:nn,5)=u(1:nn,5)+0.5D0*q(1:nn,1)*q(1:nn,4)**2
  ! pressure -> total fluid energy
  ! E = Ec + P / (gamma - 1)
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+q(1:nn,ndim+2)/(gamma-1.0d0)
  ! magnetic energy -> total fluid energy
  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,6)+q(1:nn,nvar+1))**2
  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,7)+q(1:nn,nvar+2))**2
  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,8)+q(1:nn,nvar+3))**2
  ! Magnetic field components
  u(1:nn,6:8)=q(1:nn,6:8)
  u(1:nn,nvar+1:nvar+3)=q(1:nn,nvar+1:nvar+3)
  ! passive scalars
  do ivar=9,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do

contains
!{{{
function find_Vcirc(rayon, indice)
implicit none
real(dp), intent(in)    :: rayon
integer, intent(in)    :: indice
real(dp)          :: find_Vcirc
real(dp)          :: vitesse, rayon_bin, vitesse_old, rayon_bin_old
integer          :: k, indmax

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
!}}}

end subroutine condinit

subroutine mag_constant(q,nn,B_0)
  use amr_parameters
  use hydro_parameters, only: nvar
  implicit none

  real(dp),dimension(1:nvector,1:ndim)::x    ! Cell center position
  real(dp),dimension(1:nvector,1:nvar+3)::q  ! Primitive variables
  integer ::nn                               ! Number of cells
  real(dp)::dx                               ! Cell size
  real(dp)::B_0,mag_radius,mag_height        ! Default B-strength

  ! B-left
  q(1:nn,6) = B_0
  q(1:nn,7) = 0.0
  q(1:nn,8) = 0.0

  ! B-right
  q(1:nn,nvar+1) = B_0
  q(1:nn,nvar+2) = 0.0
  q(1:nn,nvar+3) = 0.0

end subroutine mag_constant

subroutine mag_toroidal(x,q,dx,nn,B_0,mag_radius,mag_height)
  use amr_parameters
  use hydro_parameters, only: nvar
  use const
  implicit none

  real(dp),dimension(1:nvector,1:ndim)::x    ! Cell center position
  real(dp),dimension(1:nvector,1:nvar+3)::q  ! Primitive variables
  integer::nn                                ! Number of cells
  real(dp)::dx                               ! Cell size
  real(dp)::B_0,mag_radius,mag_height        ! Default B-strength
  real(dp)::xx,yy,zz
  integer::i,it,nticks
  real(dp)::dxmin,Al,Ar,ztick,rltick,rrtick,zavg,Bscale,r_scale,z_scale

  ! Toroidal field, prop. to rho**(2/3)

  dxmin=boxlen*0.5d0**nlevelmax
  nticks=nint(dx/dxmin)

  r_scale = two3rd / mag_radius
  z_scale = two3rd / mag_height

  do i=1,nn
    ! box-centered coordinates
    xx=x(i,1) - boxlen * 0.5d0
    yy=x(i,2) - boxlen * 0.5d0
    zz=x(i,3) - boxlen * 0.5d0

    ! average exp(-h) over the z axis
    ztick = zz + 0.5*(dxmin-dx)
    zavg = 0.0
    do it=1,nticks
      zavg=zavg + exp(-ABS(ztick)*z_scale)
      ztick=ztick + dxmin
    end do
    Bscale = B_0 * zavg/DBLE(nticks) / r_scale / dx

    ! B left
    ! X direction
    rltick = SQRT( (xx - 0.5*dx)**2 + (yy - 0.5*dx)**2 )
    rrtick = SQRT( (xx - 0.5*dx)**2 + (yy + 0.5*dx)**2 )

    Al= exp(-rltick*r_scale)
    Ar= exp(-rrtick*r_scale)

    q(i,6)=Bscale * (Ar-Al)

    ! Y direction
    rltick = SQRT( (xx - 0.5*dx)**2 + (yy - 0.5*dx)**2 )
    rrtick = SQRT( (xx + 0.5*dx)**2 + (yy - 0.5*dx)**2 )

    Al= exp(-rltick*r_scale)
    Ar= exp(-rrtick*r_scale)

    q(i,7)=Bscale * (Al-Ar)

    ! Z direction
    q(i,8)=0.0

    ! B right
    ! X direction
    rltick = SQRT( (xx + 0.5*dx)**2 + (yy - 0.5*dx)**2 )
    rrtick = SQRT( (xx + 0.5*dx)**2 + (yy + 0.5*dx)**2 )

    Al= exp(-rltick*r_scale)
    Ar= exp(-rrtick*r_scale)

    q(i,nvar+1)=Bscale * (Ar-Al)

    ! Y direction
    rltick = SQRT( (xx - 0.5*dx)**2 + (yy + 0.5*dx)**2 )
    rrtick = SQRT( (xx + 0.5*dx)**2 + (yy + 0.5*dx)**2 )

    Al= exp(-rltick*r_scale)
    Ar= exp(-rrtick*r_scale)

    q(i,nvar+2)=Bscale * (Al-Ar)

    ! Z direction
    q(i,nvar+3)=0.0
  end do

end subroutine mag_toroidal

subroutine mag_dipole(x,q,dx,nn,B_0,mag_radius,mag_height)
  use amr_parameters
  use hydro_parameters, only: nvar
  use const
  implicit none

  real(dp),dimension(1:nvector,1:ndim)::x    ! Cell center position
  real(dp),dimension(1:nvector,1:nvar+3)::q  ! Primitive variables
  integer::nn                                ! Number of cells
  real(dp)::dx                               ! Cell size
  real(dp)::B_0,mag_radius,mag_height        ! Default B-strength
  real(dp)::xx,yy,zz
  integer::i,it,nticks
  real(dp)::dxmin,zfl,zceil,rtick,x_edge,y_edge,r_scale,z_scale
  real(dp),dimension(1:2,1:2,1:2)::A_mag     ! x/y,left/right,up/down

  ! Dipole field, prop. to rho**(2/3)

  dxmin=boxlen*0.5d0**nlevelmax
  nticks=nint(dx/dxmin)

  r_scale = two3rd / mag_radius
  z_scale = two3rd / mag_height

  do i=1,nn
    ! box-centered coordinates
    xx=x(i,1) - boxlen * 0.5d0
    yy=x(i,2) - boxlen * 0.5d0
    zz=x(i,3) - boxlen * 0.5d0

    zfl   = ABS( zz - 0.5*dx )
    zceil = ABS( zz + 0.5*dx )

    A_mag = 0.0

    ! A_(x,l)
    x_edge = xx + 0.5*(dxmin-dx)
    y_edge = yy - 0.5*dx
    do it=1,nticks
      rtick = SQRT( x_edge**2 + y_edge**2 )

      A_mag(1,1,1) = A_mag(1,1,1) - exp(-r_scale*rtick)*exp(-z_scale*zfl)  *y_edge
      A_mag(1,1,2) = A_mag(1,1,2) - exp(-r_scale*rtick)*exp(-z_scale*zceil)*y_edge
      x_edge=x_edge + dxmin
    end do

    ! A_(x,r)
    x_edge = xx + 0.5*(dxmin-dx)
    y_edge = yy + 0.5*dx
    do it=1,nticks
      rtick = SQRT( x_edge**2 + y_edge**2 )

      A_mag(1,2,1) = A_mag(1,2,1) - exp(-r_scale*rtick)*exp(-z_scale*zfl)  *y_edge
      A_mag(1,2,2) = A_mag(1,2,2) - exp(-r_scale*rtick)*exp(-z_scale*zceil)*y_edge
      x_edge=x_edge + dxmin
    end do

    ! A_(y,l)
    x_edge = xx - 0.5*dx
    y_edge = yy + 0.5*(dxmin-dx)
    do it=1,nticks
      rtick = SQRT( x_edge**2 + y_edge**2 )

      A_mag(2,1,1) = A_mag(2,1,1) + exp(-r_scale*rtick)*exp(-z_scale*zfl)  *x_edge
      A_mag(2,1,2) = A_mag(2,1,2) + exp(-r_scale*rtick)*exp(-z_scale*zceil)*x_edge
      y_edge=y_edge + dxmin
    end do

    ! A_(y,r)
    x_edge = xx + 0.5*dx
    y_edge = yy + 0.5*(dxmin-dx)
    do it=1,nticks
      rtick = SQRT( x_edge**2 + y_edge**2 )

      A_mag(2,2,1) = A_mag(2,2,1) + exp(-r_scale*rtick)*exp(-z_scale*zfl)  *x_edge
      A_mag(2,2,2) = A_mag(2,2,2) + exp(-r_scale*rtick)*exp(-z_scale*zceil)*x_edge
      y_edge=y_edge + dxmin
    end do

    ! average value
    A_mag = A_mag / DBLE(nticks)

    ! B left
    ! X direction
    q(i,6)= - B_0 * (A_mag(2,1,2)-A_mag(2,1,1)) / dx ! [-d/dz A_y]
    ! Y direction
    q(i,7)=   B_0 * (A_mag(1,1,2)-A_mag(1,1,1)) / dx ! [ d/dz A_x]
    ! Z direction  [d/dx A_y - d/dy A_x]
    q(i,8)=   B_0 * ( A_mag(2,2,1)-A_mag(2,1,1) - A_mag(1,2,1)+A_mag(1,1,1) ) / dx

    ! B right
    ! X direction
    q(i,nvar+1)= - B_0 * (A_mag(2,2,2)-A_mag(2,2,1)) / dx ! [-d/dz A_y]
    ! Y direction
    q(i,nvar+2)=   B_0 * (A_mag(1,2,2)-A_mag(1,2,1)) / dx ! [ d/dz A_x]
    ! Z direction  [d/dx A_y - d/dy A_x]
    q(i,nvar+3)=   B_0 * ( A_mag(2,2,2)-A_mag(2,1,2) - A_mag(1,2,2)+A_mag(1,1,2) ) / dx
  end do

end subroutine mag_dipole

subroutine mag_quadrupole(x,q,dx,nn,B_0,mag_radius,mag_height)
  use amr_parameters
  use hydro_parameters, only: nvar
  use const
  implicit none

  real(dp),dimension(1:nvector,1:ndim)::x    ! Cell center position
  real(dp),dimension(1:nvector,1:nvar+3)::q  ! Primitive variables
  integer::nn                                ! Number of cells
  real(dp)::dx                               ! Cell size
  real(dp)::B_0,mag_radius,mag_height        ! Default B-strength
  real(dp)::xx,yy,zz
  integer::i,it,nticks
  real(dp)::dxmin,zfl,zceil,rtick,x_edge,y_edge,r_scale,z_scale
  real(dp),dimension(1:2,1:2,1:2)::A_mag     ! x/y,left/right,up/down

  ! Quadrupole field, prop. to rho**(2/3)

  dxmin=boxlen*0.5d0**nlevelmax
  nticks=nint(dx/dxmin)

  r_scale = two3rd / mag_radius
  z_scale = two3rd / mag_height

  do i=1,nn
    ! box-centered coordinates
    xx=x(i,1) - boxlen * 0.5d0
    yy=x(i,2) - boxlen * 0.5d0
    zz=x(i,3) - boxlen * 0.5d0

    zfl   = zz - 0.5*dx
    zceil = zz + 0.5*dx

    A_mag = 0.0

    ! A_(x,l)
    x_edge = xx + 0.5*(dxmin-dx)
    y_edge = yy - 0.5*dx
    do it=1,nticks
      rtick = SQRT( x_edge**2 + y_edge**2 )

      A_mag(1,1,1) = A_mag(1,1,1) - zfl*exp(-r_scale*rtick)*exp(-z_scale*ABS(zfl))  *y_edge
      A_mag(1,1,2) = A_mag(1,1,2) - zceil*exp(-r_scale*rtick)*exp(-z_scale*ABS(zceil))*y_edge
      x_edge=x_edge + dxmin
    end do

    ! A_(x,r)
    x_edge = xx + 0.5*(dxmin-dx)
    y_edge = yy + 0.5*dx
    do it=1,nticks
      rtick = SQRT( x_edge**2 + y_edge**2 )

      A_mag(1,2,1) = A_mag(1,2,1) - zfl*exp(-r_scale*rtick)*exp(-z_scale*ABS(zfl))  *y_edge
      A_mag(1,2,2) = A_mag(1,2,2) - zceil*exp(-r_scale*rtick)*exp(-z_scale*ABS(zceil))*y_edge
      x_edge=x_edge + dxmin
    end do

    ! A_(y,l)
    x_edge = xx - 0.5*dx
    y_edge = yy + 0.5*(dxmin-dx)
    do it=1,nticks
      rtick = SQRT( x_edge**2 + y_edge**2 )

      A_mag(2,1,1) = A_mag(2,1,1) + zfl*exp(-r_scale*rtick)*exp(-z_scale*ABS(zfl))  *x_edge
      A_mag(2,1,2) = A_mag(2,1,2) + zceil*exp(-r_scale*rtick)*exp(-z_scale*ABS(zceil))*x_edge
      y_edge=y_edge + dxmin
    end do

    ! A_(y,r)
    x_edge = xx + 0.5*dx
    y_edge = yy + 0.5*(dxmin-dx)
    do it=1,nticks
      rtick = SQRT( x_edge**2 + y_edge**2 )

      A_mag(2,2,1) = A_mag(2,2,1) + zfl*exp(-r_scale*rtick)*exp(-z_scale*ABS(zfl))  *x_edge
      A_mag(2,2,2) = A_mag(2,2,2) + zceil*exp(-r_scale*rtick)*exp(-z_scale*ABS(zceil))*x_edge
      y_edge=y_edge + dxmin
    end do

    ! average value
    A_mag = A_mag / DBLE(nticks)

    ! B left
    ! X direction
    q(i,6)= - B_0 * (A_mag(2,1,2)-A_mag(2,1,1)) / dx ! [-d/dz A_y]
    ! Y direction
    q(i,7)=   B_0 * (A_mag(1,1,2)-A_mag(1,1,1)) / dx ! [ d/dz A_x]
    ! Z direction  [d/dx A_y - d/dy A_x]
    q(i,8)=   B_0 * ( A_mag(2,2,1)-A_mag(2,1,1) - A_mag(1,2,1)+A_mag(1,1,1) ) / dx

    ! B right
    ! X direction
    q(i,nvar+1)= - B_0 * (A_mag(2,2,2)-A_mag(2,2,1)) / dx ! [-d/dz A_y]
    ! Y direction
    q(i,nvar+2)=   B_0 * (A_mag(1,2,2)-A_mag(1,2,1)) / dx ! [ d/dz A_x]
    ! Z direction  [d/dx A_y - d/dy A_x]
    q(i,nvar+3)=   B_0 * ( A_mag(2,2,2)-A_mag(2,1,2) - A_mag(1,2,2)+A_mag(1,1,2) ) / dx
  end do

end subroutine mag_quadrupole

subroutine velana(x,v,dx,t,ncell)
  use amr_parameters
  use hydro_parameters
  implicit none
  integer ::ncell                         ! Size of input arrays
  real(dp)::dx                            ! Cell size
  real(dp)::t                             ! Current time
  real(dp),dimension(1:nvector,1:3)::v    ! Velocity field
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine computes the user defined velocity fields.
  ! x(i,1:ndim) are cell center position in [0,boxlen] (user units).
  ! v(i,1:3) is the imposed 3-velocity in user units.
  !================================================================
  integer::i
  real(dp)::xx,yy,zz,vx,vy,vz,rr,tt,omega,aa,twopi

  ! Add here, if you wish, some user-defined initial conditions
  aa=1.0
  twopi=2d0*ACOS(-1d0)
  do i=1,ncell

     xx=x(i,1)
     yy=x(i,2)
     zz=x(i,3)

     ! ABC
     vx=aa*(cos(twopi*yy)+sin(twopi*zz))
     vy=aa*(sin(twopi*xx)+cos(twopi*zz))
     vz=aa*(cos(twopi*xx)+sin(twopi*yy))

     v(i,1)=vx
     v(i,2)=vy
     v(i,3)=vz

  end do


end subroutine velana


!-------------------------------------------------------------------------------------
! Circular velocity files reading
! {{{
subroutine read_vcirc_files
  use merger_commons
  implicit none
  integer:: nvitesses, ierr, i
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

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
  ! Unit conversion kpc -> code unit and km/s -> code unit
  Vcirc_dat1(:,1) = Vcirc_dat1(:,1) * 3.085677581282D21 / scale_l
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
  ! Unit conversion kpc -> code unit and km/s -> code unit
  Vcirc_dat2(:,1) = Vcirc_dat2(:,1) * 3.085677581282D21 / scale_l
  Vcirc_dat2(:,2) = Vcirc_dat2(:,2) * 1.0D5 / scale_v


end subroutine read_vcirc_files
! }}}
!--------------------------------------------------------------------------------------
