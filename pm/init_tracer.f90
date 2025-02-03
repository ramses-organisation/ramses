subroutine init_tracer
    use amr_commons
    use pm_commons
    use amr_parameters, only: i8b
    use random
    use mpi_mod
    implicit none
    integer, dimension(1:ncpu,1:IRandNumSize)::allseed
    integer(i8b),dimension(1:ncpu)::npart_cpu,npart_all
    integer(i8b)::indglob, npart_tot
#ifndef WITHOUTMPI
    integer :: info
#endif

    if(tracer_seed(1)==-1)then
        call rans(ncpu, tseed, allseed)
        tracer_seed = allseed(myid, 1:IRandNumSize)
    end if

#ifndef WITHOUTMPI
    ! Broadcast the number of particles for the id of the tracers
#ifdef LONGINT
    call MPI_ALLREDUCE(npart, npart_tot, 1, MPI_INTEGER8, MPI_SUM, MPI_COMM_WORLD, info)
#else
    call MPI_ALLREDUCE(npart, npart_tot, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, info)
#endif
#else
    npart_tot = npart
#endif
    indglob = npart_tot

    if (trim(tracer_feed_fmt) == 'binary') then
        call load_tracers_bin(1)
    else if (trim(tracer_feed_fmt) == 'binary2') then
        call load_tracers_bin(2)
    else if (trim(tracer_feed_fmt) == 'inplace') then
        call load_tracers_inplace
    else if (trim(tracer_feed_fmt) == 'ascii') then
        call load_tracers
    else
        write(*, '(a,a,a)')'Data input format not understood: "', (tracer_feed_fmt), '"'
        stop
    end if

    ! Reset first balance flags
    tracer_first_balance_levelmin = nlevelmax + 1

contains

!################################################################
!################################################################
!################################################################
!################################################################
subroutine load_tracers
    use amr_commons
    use pm_commons
    use mpi_mod
    implicit none

    real(dp):: xx1, xx2, xx3
    integer(1), dimension(1:nvector) :: ixx
    real(kind=8),dimension(1:nvector,1:3)::xx,vv
#ifndef WITHOUTMPI
    integer,dimension(1:nvector)::cc
#endif
    integer::i,ipart,jpart,icpu,buf_count
    logical::eof
    integer ,dimension(1:nvector)::ii

    ! The tracers are loaded after all the other particles, so the first tracer
    ! is the particle number 'npart'
    ipart=npart

    if(myid==1)then
       open(10, file=trim(tracer_feed), form='formatted', status='old')
       write(*, *) 'Reading initial tracers from ', trim(tracer_feed)
    end if
    eof=.false.

    ! Reading line by line until end of file
    do while (.not.eof)
       xx=0.0
       if(myid==1)then
          jpart=0
          do i=1,nvector
             read(10,*,end=100)xx1,xx2,xx3
             jpart=jpart+1
             indglob=indglob+1
             xx(i,1)=xx1*boxlen !+boxlen/2.0
             xx(i,2)=xx2*boxlen !+boxlen/2.0
             xx(i,3)=xx3*boxlen !+boxlen/2.0
             ii(i  )=indglob
             ixx(i )=FAM_TRACER_GAS
          end do
100       continue
          if(jpart<nvector)eof=.true.
       endif
       buf_count=nvector*3
#ifndef WITHOUTMPI
       call MPI_BCAST(xx,buf_count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
       call MPI_BCAST(ii,nvector  ,MPI_INTEGER         ,0,MPI_COMM_WORLD,info)
       call MPI_BCAST(ixx,nvector ,MPI_INTEGER1        ,0,MPI_COMM_WORLD,info)
       call MPI_BCAST(eof,1       ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,info)
       call MPI_BCAST(jpart,1     ,MPI_INTEGER         ,0,MPI_COMM_WORLD,info)
       call cmp_cpumap(xx,cc,jpart)
#endif

       do i=1,jpart
#ifndef WITHOUTMPI
          if(cc(i)==myid)then
#endif
             ipart=ipart+1
             if(ipart>npartmax)then
                write(*,*)'Maximum number of particles incorrect'
                write(*,*)'npartmax should be greater than',ipart, 'got', npartmax
                stop
             endif
             xp(ipart,:)  = xx(i,:)
             vp(ipart,:)  = vv(i,:)
             mp(ipart)    = tracer_mass
             levelp(ipart)= levelmin
             idp(ipart)   = ii(i)
             typep(ipart)%family  = ixx(i)
#ifndef WITHOUTMPI
          endif
#endif
       enddo

    end do
    if(myid==1)close(10)
    ! end if
    npart=ipart

    ! Compute total number of particle
    npart_cpu=0; npart_all=0
    npart_cpu(myid)=count(is_tracer(typep(:)) .and. (levelp(:) > 0))
#ifndef WITHOUTMPI
#ifndef LONGINT
    call MPI_ALLREDUCE(npart_cpu,npart_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#else
    call MPI_ALLREDUCE(npart_cpu,npart_all,ncpu,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
#endif
    npart_cpu(1)=npart_all(1)
#endif
    write(*,*)'npart=',npart_cpu(myid),'/',sum(npart_cpu), '(tracers)'

    do icpu=2,ncpu
       npart_cpu(icpu)=npart_cpu(icpu-1)+npart_all(icpu)
    end do

end subroutine load_tracers

!=================================================================
! Loads the tracer mass and decide whether we're using v1 or v2 of
! initial tracers
!=================================================================
subroutine load_tracers_bin(iversion)
    use amr_commons
    use pm_commons
    use mpi_mod
    implicit none
    integer, intent(in) :: iversion
    integer :: unit_record, ntot, ipos
    real(dp) :: tmp_tracer_mass

    if (myid == 1) then
       open(newunit=unit_record, file=trim(tracer_feed), &
            form='unformatted', status='old')
       read(unit_record) ntot
       read(unit_record) tmp_tracer_mass

       if (tracer_mass > 0) then
          if (tmp_tracer_mass /= tracer_mass) then
             write(*, *) 'WARNING: the tracer mass from file differs from the one from namelist. Keeping latter.'
          end if
          write(*, *) 'Using a tracer mass of ', tracer_mass
       else
          if (tmp_tracer_mass > 0) then
             tracer_mass = tmp_tracer_mass
             write(*, *) 'Using a tracer mass of ', tracer_mass
          end if
       end if

       call ftell(unit_record, ipos)
       close(unit_record)
    end if

#ifndef WITHOUTMPI
    call MPI_BCAST(ntot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, info)
    call MPI_BCAST(tracer_mass, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, info)
#endif

    if (iversion == 2) then
       call load_tracers_bin_v2(ntot)
    else
       if (myid == 1)  write(*, *) 'Reading initial tracers (binary v1) from ', trim(tracer_feed)

       call load_tracers_bin_v1(ntot)
    end if

end subroutine load_tracers_bin

!------------------------------------------------------------
! Create the tracer inplace by looping on the amr grid
subroutine load_tracers_inplace
    use amr_commons
    use pm_commons
    use hydro_commons, only : uold
    use mpi_mod
    implicit none
    integer :: nx_loc, icpu, jgrid, igrid, j, icell, iskip
#if NDIM > 1
    integer :: ix, iy
#if NDIM > 2
    integer :: iz
#endif
#endif
    real(dp) :: scale, dx, dx_loc, vol_loc, d

    real(dp) :: xcell(ndim), skip_loc(ndim)
    integer(i8b) :: itracer_start
    integer(i8b) :: ntracer_loc, ntracer_cpu(ncpu), idp_start

    real(dp) :: dx_cell(twotondim, ndim)

    real(dp) :: npart_loc_real, rand
    integer :: npart_loc
    integer::ind,ipart,ilevel,idim

    ! Build the positions of the cells w.r.t. their grid in dx unit
    do ind = 1, twotondim
#if NDIM == 3
       iz = (ind-1)/4
       iy = (ind-1-4*iz)/2
       ix = (ind-1-4*iz-2*iy)
       dx_cell(ind, 1) = (real(ix, dp)-0.5_dp)
       dx_cell(ind, 2) = (real(iy, dp)-0.5_dp)
       dx_cell(ind, 3) = (real(iz, dp)-0.5_dp)
#elif NDIM == 2
       iy = (ind-1) / 2
       ix = (ind-1-2*iy)
       dx_cell(ind, 1) = (real(ix, dp)-0.5_dp)
       dx_cell(ind, 2) = (real(iy, dp)-0.5_dp)
#elif NDIM == 1
       dx_cell(ind, 1) = real(ind, dp)-0.5_dp
#endif
    end do

    nx_loc = (icoarse_max - icoarse_min + 1)
    scale = boxlen / dble(nx_loc)

    skip_loc = 0

    skip_loc(1) = dble(icoarse_min)
#if NDIM > 1
    skip_loc(2) = dble(jcoarse_min)
#if NDIM > 2
    skip_loc(3) = dble(kcoarse_min)
#endif
#endif

    ! Store index of first tracer
    itracer_start = npart
    ipart = npart

    npart_loc_real = 0

    ! Loop over levels
    do ilevel = levelmin, nlevelmax
       dx = 0.5_dp**(ilevel)
       dx_loc = dx * scale
       vol_loc = dx_loc**ndim

       do jgrid = 1, active(ilevel)%ngrid
          igrid = active(ilevel)%igrid(jgrid)
          ! Loop on cells
          do ind = 1, twotondim
             iskip = ncoarse + (ind-1) * ngridmax
             icell = iskip + igrid

             ! Select leaf cells
             if (son(icell) == 0) then
                ! In zoomed region (if any)
                if (ivar_refine > 0) then
                   if (uold(icell, ivar_refine) / uold(icell, 1) < var_cut_refine) then
                      cycle
                   end if
                end if

                ! Compute number of tracers to create
                d = uold(icell, 1) * vol_loc
                npart_loc_real = d / tracer_mass
                npart_loc = int(npart_loc_real)

                ! The number of tracer is real, so we have to decide
                ! whether the number is the floor or ceiling of the
                ! real number.
                call ranf(tracer_seed, rand)

                if (rand < npart_loc_real-npart_loc) then
                   npart_loc = npart_loc + 1
                end if

                ! Get cell position
                xcell(:) = (xg(igrid, :) - skip_loc(:) + dx_cell(ind, :) * dx) * scale

                ! Now create the right number of tracers
                !
                ! Note: we don't create the idp of the tracers here. See below.
                do j = 1, npart_loc
                   ipart = ipart+1
                   if (ipart > npartmax) then
                      write(*,*) 'Maximum number of particles incorrect'
                      write(*,*) 'npartmax should be greater than', ipart, 'got', npartmax
                      stop
                   end if
                   do idim = 1, ndim
                     xp(ipart, idim) = xcell(idim)
                   end do

                   vp(ipart, :) = 0._dp
                   mp(ipart) = tracer_mass
                   levelp(ipart) = ilevel
                   typep(ipart)%family = FAM_TRACER_GAS
                end do
             end if
          end do
          ! Get next grid
       end do
       ! End loop over active grids
    end do ! End loop over levels

    ! Store total number of particules
    npart = ipart

    ! Count tracers and scatter to other CPUs
    ntracer_loc = npart - itracer_start
    ntracer_cpu(myid) = ntracer_loc

#ifndef WITHOUTMPI
#ifndef LONGINT
    call MPI_ALLGATHER(ntracer_loc, 1, MPI_INTEGER, ntracer_cpu, 1, MPI_INTEGER, MPI_COMM_WORLD, info)
#else
    call MPI_ALLGATHER(ntracer_loc, 1, MPI_INTEGER8, ntracer_cpu, 1, MPI_INTEGER8, MPI_COMM_WORLD, info)
#endif
#endif

    ! Compute number of tracer in CPUs of lesser rank
    do icpu = 2, ncpu
       ntracer_cpu(icpu) = ntracer_cpu(icpu-1) + ntracer_cpu(icpu)
    end do

    ! Get first available index: this is the total number of
    ! particules + the number of tracers in CPUs with smaller ranks
    if (myid == 1) then
       idp_start = npart_tot
    else
       idp_start = npart_tot + ntracer_cpu(myid-1)
    end if

    ! Now loop on the created particles and give them an id
    do ipart = itracer_start+1, npart+1
       idp_start = idp_start + 1
       idp(ipart) = idp_start
    end do

    ! Update the global counter (useless if nothing is loaded after the tracers)
    if (myid == 1) then
       indglob = indglob + ntracer_cpu(ncpu)
    end if

    if (myid == 1 .and. ntracer_cpu(ncpu) == 0) then
       write(*,*) '______ NO TRACER CREATED! ______'
    end if
    if (ntracer_loc > 0) &
         write(*,'(a,i15,a,i15,a,i7)') 'ntracer=', ntracer_loc, '/', ntracer_cpu(ncpu), &
         '(tracers) for PE=', myid
  end subroutine load_tracers_inplace

!------------------------------------------------------------
! Read the tracer in version 1 of the format
!
! The version 1 is a record based format with the following structure
! * ntracer [integer]
! * mtracer [float64]
! * x[ntracer] [float64]
! * y[ntracer] [float64]
! * z[ntracer] [float64]
!
subroutine load_tracers_bin_v1(ntot)
    use amr_commons
    use pm_commons
    implicit none
    integer, intent(in) :: ntot
#if NDIM < 3
    write(*,*) "Can only initialize tracer particles in 3D with method 'binary'"
    call clean_end()
#else

    integer :: unit_in

    real(dp), dimension(nvector) :: xx1, xx2, xx3
    integer(1), dimension(1:nvector) :: ixx

    integer, dimension(nvector) :: ii, cmap
    real(dp), dimension(1:nvector, 1:ndim) :: tmpxx
    integer :: icpu

    integer :: j, jj, nbuffer
    integer :: ix, iy, iz, ipart

    ! The tracers are loaded after all the other particles, so the first tracer
    ! is the particle number 'npart'
    ipart=npart

    ! Because we read the file in stream mode, we need to know where
    ! each element is using this map:
    !
    !  Length | Starting position | Comment
    ! --------+-------------------+-----------------------
    !       4 | 1                 | record length
    !       4 | 5                 | number of particles N
    !       4 | 9                 | record end
    ! --------+-------------------+-----------------------
    !       4 | 14                | record length
    !       4 | 17                | tracer mass
    !       4 | 21                | record end
    ! --------+-------------------+-----------------------
    !       4 | 25                | record length
    !      8N | ix=29             | x
    !       4 | ix+8N             | record end
    ! --------+-------------------+-----------------------
    !       4 | ix+8N+4           | record length
    !      8N | iy=ix+8N+8        | y
    !       4 | iy+8N             | record end
    ! --------+-------------------+-----------------------
    !       4 | iy+8N+4           | record length
    !      8N | iz=iy+8N+8        | z
    !       4 | iz+8N             | record end

    if (myid == 1) then
       ix = 1 + 28
       iy = ix + ntot * 8 + 8
       iz = iy + ntot * 8 + 8

       open(newunit=unit_in, file=trim(tracer_feed), access='stream', status='old')
    end if


    do jj = 1, ntot, nvector
       if (jj + nvector > ntot) then
          nbuffer = ntot - jj + 1
       else
          nbuffer = nvector
       end if

       ! Read nbuffer elements from file
       if(myid == 1) then
          do j = 1, nbuffer
             read(unit_in, pos=ix + 8*(j+jj-1)) xx1(j)
          end do

          do j = 1, nbuffer
             read(unit_in, pos=iy + 8*(j+jj-1)) xx2(j)
          end do

          do j = 1, nbuffer
             read(unit_in, pos=iz + 8*(j+jj-1)) xx3(j)
          end do

          do j = 1, nbuffer
             indglob = indglob + 1
             ii(j) = indglob
             ixx(j) = FAM_TRACER_GAS
          end do
       end if

#ifndef WITHOUTMPI
       call MPI_BCAST(xx1, nbuffer, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, info)
       call MPI_BCAST(xx2, nbuffer, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, info)
       call MPI_BCAST(xx3, nbuffer, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, info)
       call MPI_BCAST(ii,  nbuffer, MPI_INTEGER         , 0, MPI_COMM_WORLD, info)
       call MPI_BCAST(ixx, nbuffer, MPI_INTEGER         , 0, MPI_COMM_WORLD, info)
#endif
       tmpxx(1:nbuffer, 1) = xx1(1:nbuffer) * boxlen
       tmpxx(1:nbuffer, 2) = xx2(1:nbuffer) * boxlen
       tmpxx(1:nbuffer, 3) = xx3(1:nbuffer) * boxlen

       call cmp_cpumap(tmpxx, cmap, nbuffer)

       do j = 1, nbuffer
#ifndef WITHOUTMPI
          if (cmap(j)==myid) then
#endif
             ipart=ipart+1
             if(ipart>npartmax)then
                write(*,*)'Maximum number of particles incorrect'
                write(*,*)'npartmax should be greater than',ipart, 'got', npartmax
                stop
             end if
             xp(ipart, 1)  = xx1(j)
             xp(ipart, 2)  = xx2(j)
             xp(ipart, 3)  = xx3(j)

             vp(ipart,:)  = 0._dp
             mp(ipart)    = tracer_mass
             levelp(ipart)= levelmin
             idp(ipart)   = ii(j)
             typep(ipart)%family = int(ixx(j), 1)
#ifndef WITHOUTMPI
          endif
#endif
       end do ! End loop on buffer
    end do ! End loop on particle number

    if (myid == 1) close(unit_in)

    ! end if
    npart=ipart

    ! Compute total number of particle
    npart_cpu=0; npart_all=0
    npart_cpu(myid)=npart
#ifndef WITHOUTMPI
#ifndef LONGINT
    call MPI_ALLREDUCE(npart_cpu,npart_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#else
    call MPI_ALLREDUCE(npart_cpu,npart_all,ncpu,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
#endif
    npart_cpu(1)=npart_all(1)
#endif
    do icpu=2,ncpu
       npart_cpu(icpu)=npart_cpu(icpu-1)+npart_all(icpu)
    end do
    write(*,*)'npart=',npart,'/',npart_cpu(ncpu), '(tracers)'
#endif
  end subroutine load_tracers_bin_v1

!------------------------------------------------------------
! Read the tracer in version 2 of the format
!
! This version starts by 2 records containing:
! * ntracer [integer]
! * mtracer [float64]
!
! The data is then written directly in binary, encoded as float64 in
! C order. If you read 3 contiguous parts, you'll get x, y and z of
! a given particle.
subroutine load_tracers_bin_v2(ntot)
    use amr_commons
    use pm_commons
    use mpi_mod
    implicit none
    integer, intent(in) :: ntot

#if NDIM < 3
    write(*,*) "Can only initialize tracer particles in 3D with method 'binary2'"
    call clean_end()
#else
    integer :: unit_in
    real(dp), dimension(:, :), allocatable :: allpos
    integer, dimension(:), allocatable :: allcmap
    integer, dimension(nvector) :: cmap
    real(dp), dimension(1:nvector, 1:ndim) :: tmpxx
    real(dp):: xx1, xx2, xx3
    integer :: ipart

    type(communicator), dimension(:), allocatable :: sender  ! To send data
    type(communicator) :: receiver                           ! To receive data

#ifndef WITHOUTMPI
    integer, dimension(MPI_STATUS_SIZE) :: status
    integer, dimension(MPI_STATUS_SIZE, 2:ncpu) :: statuses
    integer :: ierror
    integer, dimension(2:ncpu) :: reqsend1, reqsend2
#endif

    integer :: icpu
    integer :: iwrite

    integer(i8b) :: pos
    integer :: j, jj, nbuffer, ibuff
    integer, dimension(ncpu) :: cpu_count
    integer :: cpu_count_loc

    ! The tracers are loaded after all the other particles, so the first tracer
    ! is the particle number 'npart'
    ipart=npart

    if (myid == 1) then
       allocate(allpos(ntot, 3), allcmap(ntot))
    end if

#ifndef WITHOUTMPI
    call MPI_BCAST(ntot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, info)
    call MPI_BCAST(tracer_mass, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, info)
#endif

    cpu_count(:) = 0
    if (myid == 1) then
       write(*, 321, advance='no') trim(tracer_feed)
       unit_in = 11
       open(unit=unit_in, file=trim(tracer_feed), access='stream', status='old')

       iwrite = 0
       do jj = 1, ntot, nvector
          if (iwrite < int((100. * jj) / ntot)) then
             iwrite = int((100. * jj) / ntot)
             write(*, '(".")', advance='no')
          end if

          nbuffer = min(nvector, ntot-jj+1)

          ! Read data from file
          do j = 1, nbuffer
             pos = 33 + int(jj+j-2, i8b)*24

             read(unit_in, pos=pos) xx1, xx2, xx3

             tmpxx(j, 1) = xx1 * boxlen
             tmpxx(j, 2) = xx2 * boxlen
             tmpxx(j, 3) = xx3 * boxlen
          end do

          ! Compute cpu_map
          call cmp_cpumap(tmpxx, cmap, nbuffer)

          ! Fill the sender
          do j = 1, nbuffer
             cpu_count(cmap(j)) = cpu_count(cmap(j)) + 1
             allcmap(jj+j-1) = cmap(j)
             allpos(jj+j-1, :) = tmpxx(j, :)
          end do
       end do

       close(unit_in)

       ! Force termination of line
       write (*,*) ''
    end if

#ifndef WITHOUTMPI
    ! Send count of particles to each CPU
    call MPI_SCATTER(&
         cpu_count, 1, MPI_INTEGER,     &
         cpu_count_loc, 1, MPI_INTEGER, &
         0, MPI_COMM_WORLD, ierror)
#endif
    write(*,'(a,i10,a,i10,a,i6)') 'ntracer=',cpu_count_loc,' /',ntot, ' for PE=', myid

    ! Allocate sender/receiver
    if (myid == 1) then
       allocate(sender(1:ncpu))
       do icpu = 1, ncpu
          allocate(&
               sender(icpu)%f(cpu_count(icpu), 1:1), &
               sender(icpu)%up(cpu_count(icpu), 1:3))
       end do
    end if
    allocate(&
         receiver%f(cpu_count_loc, 1:1), &
         receiver%up(cpu_count_loc, 1:3))

    ! Fill the send buffer
    if (myid == 1) then
       cpu_count(:) = 0
       do j = 1, ntot
          ibuff = cpu_count(allcmap(j)) + 1

          cpu_count(allcmap(j)) = ibuff
          indglob = indglob + 1

          sender(allcmap(j))%up(ibuff, :) = allpos(j, :)
          sender(allcmap(j))%f(ibuff, :) = indglob
       end do
       deallocate(allpos, allcmap)

#ifndef WITHOUTMPI
       do icpu = 2, ncpu
          call MPI_ISEND(sender(icpu)%up, 3*cpu_count(icpu), MPI_DOUBLE_PRECISION, &
               icpu-1, 0, MPI_COMM_WORLD, reqsend1(icpu), ierror)
          call MPI_ISEND(sender(icpu)%f,    cpu_count(icpu), MPI_INTEGER, &
               icpu-1, 0, MPI_COMM_WORLD, reqsend2(icpu), ierror)
       end do
       receiver%up = sender(myid)%up
       receiver%f = sender(myid)%f
#endif
    else
#ifndef WITHOUTMPI
       call MPI_RECV(receiver%up, 3*cpu_count_loc, MPI_DOUBLE_PRECISION, &
            0, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierror)
       call MPI_RECV(receiver%f,    cpu_count_loc, MPI_INTEGER, &
            0, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierror)
#endif
    end if

    ! Save the particles
    do j = 1, cpu_count_loc
       ipart=ipart+1
       if(ipart>npartmax)then
          write(*,*)'Maximum number of particles incorrect'
          write(*,*)'npartmax should be greater than', ipart, 'got', npartmax, 'for PE=', myid
          stop
       end if
       xp(ipart, 1)  = receiver%up(j, 1)
       xp(ipart, 2)  = receiver%up(j, 2)
       xp(ipart, 3)  = receiver%up(j, 3)

       vp(ipart,:)  = 0._dp
       mp(ipart)    = tracer_mass
       levelp(ipart)= levelmin
       idp(ipart)   = receiver%f(j, 1)
       typep(ipart)%family = FAM_TRACER_GAS
    end do

    npart=ipart

    ! Compute total number of particle
    npart_cpu=0; npart_all=0
    npart_cpu(myid)=npart

#ifndef WITHOUTMPI
    ! Wait for transmission end
    if (myid == 1) then
       call MPI_WAITALL(ncpu-1, reqsend1, statuses, ierror)
       call MPI_WAITALL(ncpu-1, reqsend2, statuses, ierror)
       do icpu = 2, ncpu
          deallocate(sender(icpu)%f, sender(icpu)%up)
       end do
       deallocate(sender)
    end if
    deallocate(receiver%f, receiver%up)
#endif
321 format('Reading initial tracers (binary v2) from ', A)
#endif
end subroutine load_tracers_bin_v2

end subroutine

!################################################################
!################################################################
!################################################################
!################################################################
subroutine read_tracer_mass
    use amr_commons
    use pm_commons
    use mpi_mod
    implicit none
#ifndef WITHOUTMPI
    integer :: info
#endif
    ! Attempt to read mass from binary file
    if (myid == 1) then
        if (trim(tracer_feed_fmt) == 'binary' .and. tracer_mass < 0) then
            open(unit=10, file=trim(tracer_feed), form='unformatted', status='old')
            read(10) ! ntot
            read(10) tracer_mass
            close(10)
        end if
    end if
    ! Broadcast to all CPUs the value of the tracer mass
#ifndef WITHOUTMPI
    call MPI_BCAST(tracer_mass, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, info)
#endif
    if (myid == 1) write(*, *) 'Using a tracer mass of ', tracer_mass
end subroutine read_tracer_mass

!################################################################
!################################################################
!################################################################
!################################################################
! On load, we need to convert the ids of the star to their local index
! in the linked list for usage for the tracer particles
subroutine convert_global_index_to_local_index(npart)
    use pm_commons, only: is_star, is_star_tracer, partp, typep, idp
    use amr_parameters, only: i8b
    implicit none
    integer, intent(in) :: npart

    integer :: i, j
    integer(i8b) :: minidp, maxidp
    integer(i8b), allocatable, dimension(:) :: isp8
    logical :: ok

    minidp = npart
    maxidp = 0
    do i = 1, npart
       if (is_star(typep(i))) then
          minidp = min(minidp, idp(i))
          maxidp = max(maxidp, idp(i))
       end if
    end do

    ! We now need to convert partp for star tracers (idp -> local index)
    ! Either the difference between the smallest star id and the largest is
    ! small enough (here, less than 50,000,000 -- that's 50M in memory)...
    if (maxidp - minidp < 50000000) then
       allocate(isp8(minidp:maxidp))
       isp8(:) = -1
       do i = 1, npart
          if (is_star(typep(i))) then
             isp8(idp(i)) = i
          end if
       end do

       ok = .true.
       do i = 1, npart
          if (is_star_tracer(typep(i))) then
             if (isp8(partp(i)) == -1) then
                write(*, *) 'An error occured while loading star tracers. Aborting.'
                stop 1
             end if
             partp(i) = isp8(partp(i))
          end if
       end do

       deallocate(isp8)
    else
    ! ... or for each tracers we loop on *all* the particles.
    ! It however costs 0 memory and is only runned once.
       do i = 1, npart
          ! Get star tracers
          if (is_star_tracer(typep(i))) then
             star_loop: do j = 1, npart
                if (is_star(typep(j))) then
                   ! Check that star's id == tracer partp
                   if (partp(i) == idp(j)) then
                      partp(i) = j
                      exit star_loop
                   end if
                end if
             end do star_loop
             if (.not. is_star(typep(partp(i)))) then
                write(*, *) 'An error occured while loading star tracers. Aborting.'
                stop 1
             end if
          end if
       end do
    end if

end subroutine convert_global_index_to_local_index
