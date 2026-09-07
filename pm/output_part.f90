subroutine backup_part(filename, filename_desc)
  use amr_commons
  use pm_commons
  use dump_utils, only : generic_dump, dump_header_info, dim_keys
  use iso_fortran_env
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer :: dummy_io, info2
  integer, parameter :: tag = 1122
#endif
  character(len=80) :: filename, filename_desc

  integer :: i, idim, unit_out, ipart
  integer :: ilevel, ibound, igrid, jgrid, jpart, npart1, ncache, counter
  character(len=80) :: fileloc
  character(len=5) :: nchar
  real(dp), allocatable, dimension(:) :: xdp
  integer(i8b), allocatable, dimension(:) :: ii8
  integer, allocatable, dimension(:) :: ll
  integer(int8), allocatable, dimension(:) :: ii1
  integer, allocatable, dimension(:) :: ind_dump

  integer :: unit_info, ivar
  logical :: dump_info

  if (verbose) write(*,*) 'Entering backup_part'

  ! Set ivar to 1 for first variable
  ivar = 1

  ! Wait for the token
#ifndef WITHOUTMPI
  if (IOGROUPSIZE > 0) then
     if (mod(myid-1, IOGROUPSIZE) /= 0) then
        call MPI_RECV(dummy_io, 1, MPI_INTEGER, myid-1-1, tag, &
             & MPI_COMM_WORLD, MPI_STATUS_IGNORE, info2)
     end if
  end if
#endif


  call title(myid, nchar)
  fileloc = TRIM(filename) // TRIM(nchar)
  open(newunit=unit_out, file=TRIM(fileloc), form='unformatted')
  if (myid == 1) then
     open(newunit=unit_info, file=trim(filename_desc), form='formatted')
     call dump_header_info(unit_info)
     dump_info = .true.
  else
     dump_info = .false.
  end if

  rewind(unit_out)
  ! Write header
  write(unit_out) ncpu
  write(unit_out) ndim
  write(unit_out) npart
  if (MC_tracer) then
     write(unit_out) localseed, tracer_seed
  else
     write(unit_out) localseed
  end if
  write(unit_out) nstar_tot
  write(unit_out) mstar_tot
  write(unit_out) mstar_lost
  write(unit_out) nsink

  ! Gather the order in which particles should be outputted. Dumping in tree
  ! order preserves the particle linked-list order upon restart. This is
  ! needed for making processes that assign random numbers in tree order
  ! (e.g. MC tracers) reproducable.
  ! Remark that merge_tree_fine does not clean up the linked list, so
  ! particles are duplicated accross levels. We track which particle is
  ! already gathered with the dumpedp array.

  allocate(ind_dump(1:npart))  !particle order, to be gathered
  dumpedp = .false.            !particle index already gathered?
  counter=0
  do ilevel = 1, nlevelmax
     ! Loop over cpus and boundary regions:
     ! - The dump happens between make_tree_fine and virtual_tree_fine, so
     !   particles can sit in the ghostzones waiting for MPI transfer.
     ! - Some particles could have gotten stuck in non-periodic boundaries.
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           igrid=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           igrid=headb(ibound-ncpu,ilevel)
        end if
        do jgrid=1,ncache
           npart1=numbp(igrid)
           if(npart1>0)then
              ipart=headp(igrid)
              do jpart=1,npart1
                 if (levelp(ipart)>0 .and. .not.dumpedp(ipart)) then
                    counter=counter+1
                    ind_dump(counter)=ipart
                    dumpedp(ipart)=.true.
                 end if
                 ipart=nextp(ipart)
              end do
           end if
           igrid=next(igrid)
        end do
     end do
     ! exit when all particles have been accounted for
     if (counter == npart) exit
  end do

  ! Safety net: a particle attached to no grid at all would be missing from
  ! the walk. Append it in array order rather than lose it from the dump.
  ! Normally the code should never get here.
  if (counter < npart) then
     write(*,*) 'Warning in backup_part: ', npart-counter, ' particles are not attached to any grid'
     do i = 1, npartmax
        if (levelp(i) > 0 .and. .not. dumpedp(i)) then
           counter = counter+1
           if (counter > npart) then
              write(*,*) 'Error in backup_part: more particles than npart=', npart
              call clean_stop
           end if
           ind_dump(counter) = i
        end if
     end do
  end if

  ! Write position
  allocate(xdp(1:npart))
  do idim = 1, ndim
     do i = 1, npart
        xdp(i) = xp(ind_dump(i), idim)
     end do
     call generic_dump("position_"//dim_keys(idim), ivar, xdp, unit_out, dump_info, unit_info)
  end do
  ! Write velocity
  do  idim = 1, ndim
     do i = 1, npart
        xdp(i) = vp(ind_dump(i), idim)
     end do
     call generic_dump("velocity_"//dim_keys(idim), ivar, xdp, unit_out, dump_info, unit_info)
  end do
  ! Write mass
  do i = 1, npart
     xdp(i) = mp(ind_dump(i))
  end do
  call generic_dump("mass", ivar, xdp, unit_out, dump_info, unit_info)
  deallocate(xdp)
  ! Write identity
  allocate(ii8(1:npart))
  do i = 1, npart
     ii8(i) = idp(ind_dump(i))
  end do
  call generic_dump("identity", ivar, ii8, unit_out, dump_info, unit_info)
  deallocate(ii8)

  ! Write level
  allocate(ll(1:npart))
  do i = 1, npart
     ll(i) = levelp(ind_dump(i))
  end do
  call generic_dump("levelp", ivar, ll, unit_out, dump_info, unit_info)

  deallocate(ll)

  ! Write family
  allocate(ii1(1:npart))
  do i = 1, npart
     ii1(i) = int(typep(ind_dump(i))%family, 1)
  end do
  call generic_dump("family", ivar, ii1, unit_out, dump_info, unit_info)

  ! Write tag
  do i = 1, npart
     ii1(i) = int(typep(ind_dump(i))%tag, 1)
  end do
  call generic_dump("tag", ivar, ii1, unit_out, dump_info, unit_info)
  deallocate(ii1)

#ifdef OUTPUT_PARTICLE_POTENTIAL
  ! Write potential (added by AP)
  allocate(xdp(1:npart))
  do i = 1, npart
     xdp(i) = ptcl_phi(ind_dump(i))
  end do
  call generic_dump("potential", ivar, xdp, unit_out, dump_info, unit_info)

  deallocate(xdp)
#endif

  ! Write birth epoch
  if (star .or. sink) then
     allocate(xdp(1:npart))
     do i = 1, npart
        xdp(i) = tp(ind_dump(i))
     end do
     call generic_dump("birth_time", ivar, xdp, unit_out, dump_info, unit_info)
     ! Write metallicity
     if (metal) then
        do i = 1, npart
           xdp(i) = zp(ind_dump(i))
        end do
        call generic_dump("metallicity", ivar, xdp, unit_out, dump_info, unit_info)
     end if
     deallocate(xdp)
  end if

  if (MC_tracer) then
     ! Dump particle pointer
     allocate(ll(1:npart))
     ! Get the idp of the stars on which tracers are attached
     do i = 1, npart
        ipart = ind_dump(i)
        ! For star tracers, store the id of the star instead of local index
        if (is_star_tracer(typep(ipart))) then
           ll(i) = idp(partp(ipart))
        else ! store the relative location
           ll(i) = partp(ipart)
        end if
     end do

     call generic_dump("partp", ivar, ll, unit_out, dump_info, unit_info)
     deallocate(ll)
  end if

  deallocate(ind_dump)

  !------------!
  close(unit_out)
  if (myid == 1) close(unit_info)

  ! Send the token
#ifndef WITHOUTMPI
  if (IOGROUPSIZE > 0) then
     if (mod(myid, IOGROUPSIZE) /= 0 .and. (myid .lt. ncpu)) then
        dummy_io = 1
        call MPI_SEND(dummy_io, 1, MPI_INTEGER, myid-1+1, tag, &
             & MPI_COMM_WORLD, info2)
     end if
  end if
#endif

end subroutine backup_part
