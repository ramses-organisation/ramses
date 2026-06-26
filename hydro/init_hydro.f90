subroutine init_hydro
  use amr_commons
  use hydro_commons
#ifdef RT
  use rt_parameters,only: convert_birth_times
#endif
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info,info2,dummy_io
#endif
  integer::ncell,ncache,iskip,igrid,i,ilevel,ind,ivar
  integer::nvar2,ilevel2,numbl2,ilun,ibound,istart
  integer::ncpu2,ndim2,nlevelmax2,nboundary2
  integer ,dimension(:),allocatable::ind_grid
  real(dp),dimension(:),allocatable::xx
  real(dp)::gamma2
  character(LEN=80)::fileloc
  character(LEN=5)::nchar,ncharcpu
  integer,parameter::tag=1108
#if NENER>0
  integer::irad
#endif

  if(verbose)write(*,*)'Entering init_hydro'

  !------------------------------------------------------
  ! Allocate conservative, cell-centered variables arrays
  !------------------------------------------------------
  ncell=ncoarse+twotondim*ngridmax
  allocate(uold(1:ncell,1:nvar_all))
  allocate(unew(1:ncell,1:nvar_all))
  uold=0.0d0; unew=0.0d0
  if(MC_tracer) then
     allocate(fluxes(1:ncell,1:twondim))
     fluxes(1:ncell,1:twondim)=0.0d0
  end if
  if(momentum_feedback>0)then
     allocate(pstarold(1:ncell))
     allocate(pstarnew(1:ncell))
     pstarold=0.0d0; pstarnew=0.0d0
  endif
  if(pressure_fix)then
     allocate(divu(1:ncell))
     allocate(enew(1:ncell))
     divu=0.0d0; enew=0.0d0
  end if
  if(strict_equilibrium>0)then
     allocate(rho_eq(1:ncell))
     allocate(p_eq(1:ncell))
     rho_eq=0.0d0; p_eq=0.0d0
  endif

  !--------------------------------
  ! For a restart, read hydro file
  !--------------------------------
  if(nrestart>0)then
     ilun=ncpu+myid+103
     call title(nrestart,nchar)

     if(IOGROUPSIZEREP>0)then
        call title(((myid-1)/IOGROUPSIZEREP)+1,ncharcpu)
        fileloc='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/hydro_'//TRIM(nchar)//'.out'
     else
        fileloc='output_'//TRIM(nchar)//'/hydro_'//TRIM(nchar)//'.out'
     endif

     call title(myid,nchar)
     fileloc=TRIM(fileloc)//TRIM(nchar)

     ! Wait for the token
#ifndef WITHOUTMPI
     if(IOGROUPSIZE>0) then
        if (mod(myid-1,IOGROUPSIZE)/=0) then
           call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
                & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
        end if
     endif
#endif

     open(unit=ilun,file=fileloc,form='unformatted')
     read(ilun)ncpu2
     read(ilun)nvar2
     if(strict_equilibrium>0)nvar2=nvar2-2
#ifdef SOLVERmhd
     nvar2=nvar2-3
#endif
     read(ilun)ndim2
     read(ilun)nlevelmax2
     read(ilun)nboundary2
     read(ilun)gamma2
     if(myid==1)then
        write(*,*)'Restart - Non-thermal pressure / Passive scalar mapping'
        write(*,'(A50)')"__________________________________________________"
        do i=1,nvar2-nhydro
            if(remap_pscalar(i).gt.0) then
               write(*,'(A,I3,A,I3)') ' Restart var',i+nhydro,' loaded in var',remap_pscalar(i)
            else if(remap_pscalar(i).gt.-1)then
               write(*,'(A,I3,A)') ' Restart var',i+nhydro,' read but not loaded'
            else
               write(*,'(A,I3,A)') ' Restart var',i+nhydro,' not read'
            endif
        enddo
        write(*,'(A50)')"__________________________________________________"
     endif
#ifdef RT
     if((neq_chem.or.rt).and.nvar2.lt.nvar)then ! OK to add ionization fraction vars
        ! Convert birth times for RT postprocessing:
        if(rt.and.static) convert_birth_times=.true.
        if(myid==1) write(*,*)'File hydro.tmp is not compatible'
        if(myid==1) write(*,*)'Found nvar2  =',nvar2
        if(myid==1) write(*,*)'Expected=',nvar
        if(myid==1) write(*,*)'..so only reading available variables and setting the rest to zero'
     end if
     if((neq_chem.or.rt).and.nvar2.gt.nvar)then ! Not OK to drop variables
#else
     if(nvar2.ne.(nvar))then
#endif
        if(myid==1) write(*,*)'File hydro.tmp is not compatible'
        if(myid==1) write(*,*)'Found   =',nvar2
        if(myid==1) write(*,*)'Expected=',nvar
        call clean_stop
     end if
     do ilevel=1,nlevelmax2
        do ibound=1,nboundary+ncpu
           if(ibound<=ncpu)then
              ncache=numbl(ibound,ilevel)
              istart=headl(ibound,ilevel)
           else
              ncache=numbb(ibound-ncpu,ilevel)
              istart=headb(ibound-ncpu,ilevel)
           end if
           read(ilun)ilevel2
           read(ilun)numbl2
           if(numbl2.ne.ncache)then
              write(*,*)'File hydro.tmp is not compatible'
              write(*,*)'Found   =',numbl2,' for level ',ilevel2
              write(*,*)'Expected=',ncache,' for level ',ilevel
           end if
           if(ncache>0)then
              allocate(ind_grid(1:ncache))
              allocate(xx(1:ncache))
              ! Loop over level grids
              igrid=istart
              do i=1,ncache
                 ind_grid(i)=igrid
                 igrid=next(igrid)
              end do
              ! Loop over cells
              do ind=1,twotondim
                 iskip=ncoarse+(ind-1)*ngridmax

                 ! Loop over conservative variables
                 ! Read density (no conversion needed)
                 read(ilun)xx
                 call scatter_conservative_to_uold(ind_grid, iskip, 1, xx, ncache)
                 ! Read velocities --> momenta
                 do ivar=2,neul-1
                    read(ilun)xx
                    if (read_conservative) then
                       call scatter_conservative_to_uold(ind_grid, iskip, ivar, xx, ncache)
                    else
                       call scatter_primitive_to_uold(ind_grid, iskip, ivar, xx, ncache)
                    endif
                 end do
#ifdef SOLVERmhd
                 ! Read left magnetic field (no conversion needed)
                 do ivar=6,8
                    read(ilun)xx
                    call scatter_conservative_to_uold(ind_grid, iskip, ivar, xx, ncache)
                 end do
                 ! Read right magnetic field (no conversion needed)
                 do ivar=nvar+1,nvar+3
                    read(ilun)xx
                    call scatter_conservative_to_uold(ind_grid, iskip, ivar, xx, ncache)
                 end do
#endif
#if NENER>0
                 ! Read non-thermal pressures --> non-thermal energies
                 do ivar=nhydro+1,nhydro+nener
                    if(remap_pscalar(ivar-nhydro).gt.-1) read(ilun)xx
                    if(remap_pscalar(ivar-nhydro).gt.0) then
                       if (read_conservative) then
                          call scatter_conservative_to_uold(ind_grid, iskip, remap_pscalar(ivar-nhydro), xx, ncache)
                       else
                          call scatter_primitive_to_uold(ind_grid, iskip, remap_pscalar(ivar-nhydro), xx, ncache)
                       endif
                    else if(remap_pscalar(ivar-nhydro).lt.0) then
                       do i=1,ncache
                          uold(ind_grid(i)+iskip,abs(remap_pscalar(ivar-nhydro)))=0d0
                       end do
                    endif
                 end do
#endif
                 ! Read thermal pressure --> total fluid energy
                 read(ilun)xx
                 if (read_conservative) then
                    call scatter_conservative_to_uold(ind_grid, iskip, neul, xx, ncache)
                 else
                    call calc_total_energy_from_thermal_pressure(ind_grid, iskip, xx, ncache)
                 endif
#if NVAR>NHYDRO+NENER
                 ! Read passive scalars if any
                 do ivar=nhydro+1+nener,max(nvar2,nvar)
                    if(remap_pscalar(ivar-nhydro).gt.-1) read(ilun)xx
                    if(ivar.gt.nvar)then
                       continue
                    endif
                    if(remap_pscalar(ivar-nhydro).gt.0)then
                       if (read_conservative) then
                          call scatter_conservative_to_uold(ind_grid, iskip, remap_pscalar(ivar-nhydro), xx, ncache)
                       else
                          call scatter_primitive_to_uold(ind_grid, iskip, remap_pscalar(ivar-nhydro), xx, ncache)
                       endif
                    else if(remap_pscalar(ivar-nhydro).lt.0) then
                       do i=1,ncache
                          uold(ind_grid(i)+iskip,abs(remap_pscalar(ivar-nhydro)))=0d0
                       end do
                    endif
                 end do
#endif
                 ! Read equilibrium density and pressure profiles
                 if(strict_equilibrium>0)then
                    read(ilun)xx
                    do i=1,ncache
                       rho_eq(ind_grid(i)+iskip)=xx(i)
                    end do
                    read(ilun)xx
                    do i=1,ncache
                       p_eq(ind_grid(i)+iskip)=xx(i)
                    end do
                 endif

              end do
              deallocate(ind_grid,xx)
           end if
        end do
     end do
     close(ilun)

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

#ifndef WITHOUTMPI
     if(debug)write(*,*)'hydro.tmp read for processor ',myid
     call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
     if(verbose)write(*,*)'HYDRO backup files read completed'
  end if

end subroutine init_hydro
!#####################################################################
!#####################################################################
!#####################################################################
subroutine scatter_conservative_to_uold(ind_grid, iskip, ivar, xx, ncache)
   use amr_parameters, only:dp
   use hydro_commons
   implicit none
   integer,intent(in)::ivar,iskip,ncache
   integer, dimension(1:ncache),intent(in)::ind_grid
   real(dp), dimension(1:ncache),intent(in)::xx
   !----------------------------------------------------------------------
   ! Scatter a variable from the array xx directly into uold at index ivar
   ! (e.g. density, B_left, B_right)
   !----------------------------------------------------------------------
   integer::i

   do i = 1, ncache
      uold(ind_grid(i)+iskip, ivar) = xx(i)
   end do

end subroutine scatter_conservative_to_uold
!#####################################################################
!#####################################################################
!#####################################################################
subroutine scatter_primitive_to_uold(ind_grid, iskip, ivar, xx, ncache)
   use amr_parameters, only:dp
   use hydro_commons
   implicit none
   integer,intent(in)::ivar,iskip,ncache
   integer, dimension(1:ncache),intent(in)::ind_grid
   real(dp), dimension(1:ncache),intent(in)::xx
   !-----------------------------------------------------------------------
   ! Scatter a primitive variable from the array xx into index ivar of uold
   ! (which contains conservative quantities)
   ! (e.g. velocity to momentum)
   !-----------------------------------------------------------------------
   integer::i

   select case(ivar)

#if NENER > 0
   case(nhydro+1:nhydro+nener)
      ! non-thermal pressures
      do i = 1, ncache
         uold(ind_grid(i)+iskip,ivar)=xx(i)/(gamma_rad(ivar-nhydro)-1d0)
      end do
#endif

   case default
      do i = 1, ncache
         uold(ind_grid(i)+iskip,ivar)=xx(i)*max(uold(ind_grid(i)+iskip,1),smallr)
      end do

   end select

end subroutine scatter_primitive_to_uold
!#####################################################################
!#####################################################################
!#####################################################################
subroutine calc_total_energy_from_thermal_pressure(ind_grid, iskip, pressure, ncache)
   use amr_parameters, only:dp
   use hydro_commons
   implicit none
   integer,intent(in)::iskip,ncache
   integer, dimension(1:ncache),intent(in)::ind_grid
   real(dp), dimension(1:ncache),intent(in)::pressure
   !--------------------------------------------------------------------------------------
   ! Calculate the total energy from the thermal and store it in uold(:,neul)
   !--------------------------------------------------------------------------------------
   integer::i
#if NENER > 0
   integer :: irad
#endif
   real(dp) :: d,energy
#ifdef SOLVERmhd
   real(dp) :: A, B, C
#endif

   do i=1,ncache
      d = max(uold(ind_grid(i)+iskip,1),smallr)
      ! convert thermal pressure to thermal energy
      energy=pressure(i)/(gamma-1d0)
      ! add kinetic energy
      if (uold(ind_grid(i)+iskip,1)>0.)then
         energy = energy + 0.5d0*uold(ind_grid(i)+iskip,2)**2/d
#if NDIM>1 || SOLVERmhd
         energy = energy + 0.5d0*uold(ind_grid(i)+iskip,3)**2/d
#endif
#if NDIM>2 || SOLVERmhd
         energy = energy + 0.5d0*uold(ind_grid(i)+iskip,4)**2/d
#endif
#ifdef SOLVERmhd
         ! add magnetic energy
         A = 0.5d0*(uold(ind_grid(i)+iskip, 6)+uold(ind_grid(i)+iskip, nvar+1))
         B = 0.5d0*(uold(ind_grid(i)+iskip, 7)+uold(ind_grid(i)+iskip, nvar+2))
         C = 0.5d0*(uold(ind_grid(i)+iskip, 8)+uold(ind_grid(i)+iskip, nvar+3))
         energy = energy + 0.5*(A**2+B**2+C**2)
#endif
#if NENER>0
         ! add non-thermal energies
         do irad=1,nener
            energy = energy + uold(ind_grid(i)+iskip,nhydro+irad)
         end do
#endif
      else
          energy=0
      end if
      ! Store the total energy
      uold(ind_grid(i)+iskip, neul) = energy
   end do

end subroutine calc_total_energy_from_thermal_pressure
