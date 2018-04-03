!--------------------------------------------------------
! Merger Tree patch. See README for more information.
!
!
! Contains:
! subroutine make_merger_tree()
! subroutine process_progenitor_data()
! subroutine create_prog_desc_links()
! subroutine make_trees()
!   contains subroutine find_main_desc()
!            subroutine find_main_prog()
!            subroutine add_new_pmprog()
!            subroutine find_prog_in_older_snapshots()
! subroutine read_progenitor_data()
! subroutine write_trees()
! subroutine write_progenitor_data()
! subroutine get_local_prog_id()
! subroutine fill_matrix()
! subroutine deallocate_mergertree()
! subroutine mark_tracer_particles()
! subroutine dissolve_small_clumps()
!  countains subroutine get_exclusive_clump_mass()
!--------------------------------------------------------



#if NDIM == 3

!================================
subroutine make_merger_tree()
!================================

  !-----------------------------------------
  ! This subroutine is the main routine for
  ! creating merger trees that calls all
  ! others.
  !-----------------------------------------
  use clfind_commons 
  use amr_commons

  implicit none


  !----------------------------
  ! Create mergertrees
  !----------------------------

  if (myid==1)write(*,'(A31,x,I5)') " Calling merger tree for output", ifout
  ! create trees only if progenitors might exist
  if (ifout > 1) then 

    ! first allocate some global arrays for current peaks
    allocate(main_prog(1:npeaks_max))
    main_prog = 0

    ! Read in progenitor files
    call read_progenitor_data()

    if (nprogs > 0) then
      ! Sort out progenitor data.
      call process_progenitor_data()
      
      ! link progenitors and descendants
      call create_prog_desc_links()

      ! create trees
      call make_trees()

    endif

    ! write tree to file
    call write_trees()

  endif




  !----------------------------
  ! Prepare for next round
  !----------------------------

  ! Mark tracer particles
  call mark_tracer_particles()

  ! write output
  call write_progenitor_data()

  ! Finish
  call deallocate_mergertree()



  !--------------------
  ! Say good bye.
  !--------------------

  if(verbose) write(*,*) "Finished making merger tree."

end subroutine make_merger_tree






!=====================================
subroutine process_progenitor_data()
!=====================================

  !-----------------------------------------------------
  ! This subroutine processes the progenitor data:
  !   - It finds which tracer particles are currently 
  !     on this CPU
  !   - It finds which CPU is owner of which progenitor, 
  !     and communicates that info across all CPUs.
  !   - Counts how many tracer particles of any
  !     progenitor are on this processor, create cleaned
  !     up list
  !------------------------------------------------------

  use clfind_commons
  use amr_parameters, only: i8b
  use pm_commons, only: idp, npartmax

  implicit none

#ifndef WITHOUTMPI
  include 'mpif.h'
  integer, dimension(:), allocatable :: local_owners_info
#endif

  integer(i8b), allocatable, dimension(:) :: idp_copy, galaxy_tracers_copy
  integer, allocatable, dimension(:)      :: part_local_ind, sort_ind_past, dummy
  real(dp), allocatable, dimension(:)     :: dummy_real
  integer :: itrace, ipart, igalaxy, ipastprog, i, iprog

  integer, dimension(:), allocatable :: tracers_local_pid_long, tracer_local_ids_long


  if (verbose) write(*,*) " processing progenitor data."

  

  !-------------------
  ! Allocate stuff
  !-------------------

  ! Create copies of arrays you mustn't modify
  allocate(idp_copy(1:npartmax))
  idp_copy = idp

  allocate(galaxy_tracers_copy(1:nprogs))
  galaxy_tracers_copy = galaxy_tracers

  ! allocate others
  if (nprogs > npastprogs) then
    allocate(dummy(1:nprogs))
  else
    allocate(dummy(1:npastprogs))
  endif
  dummy = 1
  allocate(dummy_real(1:npastprogs))

  allocate(part_local_ind(1:npartmax))
  part_local_ind = [(i, i = 1, npartmax)]

  allocate(sort_ind_past(1:npastprogs))
  sort_ind_past = [(i, i=1, npastprogs)]

  allocate(tracers_local_pid_long(1:nprogs*nmost_bound)) ! local particle   id for tracers (room enough for all tracers)
  allocate(tracer_local_ids_long(1:nprogs*nmost_bound))  ! local progenitor id for tracers (room enough for all tracers)




  !---------------------------------
  ! Sort arrays for quick matching
  !---------------------------------

  ! Sort arrays for matching
  call quick_sort_int_int(idp_copy, part_local_ind, npartmax)
  call quick_sort_int_int(tracers_all, tracer_loc_progids_all, nprogs*nmost_bound)
  call quick_sort_int_int(galaxy_tracers_copy, dummy, nprogs)

  ! sort past progenitors in-place by galaxy particle ID;
  ! will be needed later for multi-snapshot progenitor search
  call quick_sort_int_int(pmprogs_galaxy, sort_ind_past, npastprogs)

  dummy(1:npastprogs) = pmprogs(1:npastprogs)
  do i = 1, npastprogs
    pmprogs(i) = dummy(sort_ind_past(i))
  enddo

  dummy(1:npastprogs) = pmprogs_t(1:npastprogs)
  do i = 1, npastprogs
    pmprogs_t(i) = dummy(sort_ind_past(i))
  enddo

  dummy_real(1:npastprogs) = pmprogs_mass(1:npastprogs)

  do i = 1, npastprogs
    pmprogs_mass(i) = dummy_real(sort_ind_past(i))
  enddo

  deallocate(dummy, dummy_real, sort_ind_past)




  !------------------------
  ! find starting indices
  !------------------------

  itrace = 0; ipart = 0; igalaxy = 0; iprog = 0;

  do i = 1, npartmax
    if (idp_copy(i) > 0) then
      ipart = i
      exit
    endif
  enddo

  do i = 1, nprogs*nmost_bound
    if (tracers_all(i) > 0) then
      itrace = i
      exit
    endif
  enddo

  igalaxy = 1
  ipastprog = 1

  ntracers = 0 




  !--------------------------------------------------------------
  ! Identify which progenitor tracer particles are on this CPU.
  ! loop over all local particles and all progenitor particles.
  ! at this point, both arrays are sorted.
  ! Raise the array index of array which has lower particle ID 
  ! to find matches.
  !--------------------------------------------------------------

  do while ( ipart <= npartmax)

    !-------------------------------------------------------------------------------------
    ! Check for tracers and past progenitor galaxies while itrace <= nprogs*nmost_bound 
    !-------------------------------------------------------------------------------------

    do while (itrace <= nprogs * nmost_bound .and. ipart <= npartmax) 
      ! if particles aren't a match, raise the index where lower ID is
      if (tracers_all(itrace) < idp_copy(ipart)) then
        itrace = itrace + 1

      else if (tracers_all(itrace) > idp_copy(ipart)) then
        ! before raising local particle index, check whether you own a past progenitor
        ! no need to check whether npastprogs > 0: arrays are allocated
        ! (1:npastprogs+nprogs) to have extra space in case new progs
        ! need to be added to the list
        do while (ipastprog <= npastprogs)
          if (pmprogs_galaxy(ipastprog) < idp_copy(ipart)) then
              ipastprog = ipastprog + 1
          else if (pmprogs_galaxy(ipastprog) == idp_copy(ipart)) then 
            ! you found a match!
            pmprogs_owner(ipastprog) = myid
            ipastprog = ipastprog + 1
          else  
            exit
          endif
        enddo
        ipart = ipart + 1


      else
        !----------------
        ! found a match!
        !----------------
        iprog = tracer_loc_progids_all(itrace)
        ! count the tracer for corresponding prog
        ! add to list of this CPU
        ntracers = ntracers + 1                                   ! count one more tracer on this cpu
        tracers_local_pid_long(ntracers) = part_local_ind(ipart)  ! save local particle ID
        tracer_local_ids_long(ntracers) = iprog                   ! save local prog ID it belongs to

        ! check if found tracer is also a galaxy:
        do while (igalaxy <= nprogs)
          if (galaxy_tracers_copy(igalaxy) < tracers_all(itrace)) then
            igalaxy = igalaxy + 1
          else if (galaxy_tracers_copy(igalaxy) == tracers_all(itrace)) then
            prog_owner(iprog) = myid
            igalaxy = igalaxy + 1
            exit
          else
            exit
          endif
        enddo

        ipart = ipart + 1
        itrace = itrace + 1
      endif

    enddo
      

    ! reset ipart in case it reached the limit in the previous loop
    if (ipart>npartmax) ipart = npartmax

    !---------------------------------------------------------------------------
    ! If there are no more tracers to check for, check only for past galaxies
    !---------------------------------------------------------------------------
    do while (ipastprog <= npastprogs)
      if (pmprogs_galaxy(ipastprog) < idp_copy(ipart)) then
          ipastprog = ipastprog + 1
      else if (pmprogs_galaxy(ipastprog) == idp_copy(ipart)) then 
        ! you found a match!
        pmprogs_owner(ipastprog) = myid
        ipastprog = ipastprog + 1
      else  
        exit
      endif
    enddo

    if (ipastprog > npastprogs) exit
    ipart = ipart + 1

  enddo


  deallocate(idp_copy, galaxy_tracers_copy, part_local_ind)




  !---------------------------------------------------
  ! Clean up local tracers information for this CPU
  !---------------------------------------------------

  allocate(tracers_loc_pid(1:ntracers))
  allocate(tracer_loc_progids(1:ntracers))

  tracers_loc_pid(1:ntracers) = tracers_local_pid_long(1:ntracers)
  tracer_loc_progids(1:ntracers) = tracer_local_ids_long(1:ntracers)

  ! Deallocate auxilliary arrays
  deallocate(tracer_local_ids_long, tracers_local_pid_long)

  ! Deallocate global arrays that aren't used anymore
  deallocate(tracers_all) 
  deallocate(tracer_loc_progids_all)
  




#ifndef WITHOUTMPI
  !--------------------------------
  ! communicate progenitor owners
  !--------------------------------

  ! first the actual progenitors
  allocate(local_owners_info(1:nprogs))     ! progenitor owners array, local to each cpu
  local_owners_info = prog_owner
  call MPI_ALLREDUCE(local_owners_info, prog_owner, nprogs, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, i)
  deallocate(local_owners_info)

  ! then the past progenitors
  allocate(local_owners_info(1:npastprogs)) ! progenitor owners array, local to each cpu
  local_owners_info = pmprogs_owner(1:npastprogs)
  call MPI_ALLREDUCE(local_owners_info, pmprogs_owner, npastprogs, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, i)
  deallocate(local_owners_info)
#endif


  return

end subroutine process_progenitor_data







!=====================================
subroutine create_prog_desc_links()
!=====================================

  !--------------------------------------------------------------------
  ! Establishes connections between progenitors and descendants via
  ! tracer particles, i.e. in which descendant they ended up.
  ! Initialises p2d_links matrix, then all data are combined at the
  ! progenitor owner's CPU.
  !--------------------------------------------------------------------

  use clfind_commons

  implicit none

  integer :: iprog, ipart, ind, idesc, i

#ifndef WITHOUTMPI
  include 'mpif.h'
  integer, dimension(:), allocatable :: sendcount, receivecount
  integer, dimension(:), allocatable :: sendcount2, receivecount2
  integer, dimension(:), allocatable :: sendbuf, recvbuf
  integer, dimension(:), allocatable :: sendbuf2, recvbuf2
  integer, dimension(:), allocatable :: send_displ, rec_displ
  integer :: icpu, ind2, recndesc
#endif

  if (verbose) write(*,*) "linking progenitors and descendants."


  !------------------------------------------------------------------
  ! Initialise sparse matrix for progenitors->descendant (p2d) links
  !------------------------------------------------------------------

  allocate(p2d_links%first(1:nprogs))
  p2d_links%first = 0
  allocate(p2d_links%cnt(1:nprogs))
  p2d_links%cnt = 0
  allocate(p2d_links%clmp_id(1:10*npeaks_max))
  p2d_links%clmp_id = 0
  allocate(p2d_links%ntrace(1:10*npeaks_max))
  p2d_links%ntrace = 0
  allocate(p2d_links%next(1:10*npeaks_max))
  p2d_links%next = 0






  !----------------------------
  ! Fill up matrix with values 
  !----------------------------

  i = 0
  do ipart = 1, ntracers
    ! If tracer is currently in a descendant clump, add it to list
    if (clmpidp(tracers_loc_pid(ipart)) /= 0) then
      call fill_matrix(p2d_links, tracer_loc_progids(ipart), clmpidp(tracers_loc_pid(ipart)), 1, 'add')
    else
      i = i + 1 !count how many zeros 
    endif
  enddo

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(MPI_IN_PLACE, i, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD, ipart)
#endif

  if (myid == 1) write(*, '(A6,x,I9,x,A58)') " Found", i, "progenitor tracer particles that are not in clumps anymore."
  
  deallocate(tracers_loc_pid, tracer_loc_progids)





#ifndef WITHOUTMPI

  !=====================================================
  ! PART 1: 
  ! Send local progenitor data to progenitor owners
  !=====================================================

  !------------------------------------------------------------------------
  ! Find out how much stuff you need to send where.
  ! You'll send an array of the form: ("#" = "number of")
  ! prog1_local_id, # descendants, descendant 1 ID, 
  !   #tracers in desc1, descendant 2 ID, #tracers in desc 2, ...,
  ! prog2_local_id, #descendants, descendant 1 ID, #tracers in desc1, ...
  !
  ! prog_local_id should be the same on every CPU
  !-------------------------------------------------------------------------

  allocate(sendcount(1:ncpu))
  sendcount = 0
  allocate(receivecount(1:ncpu))
  receivecount = 0


  do iprog = 1, nprogs
    ! if progenitor has tracers on this cpu
    if ( p2d_links%cnt(iprog)> 0 .and. prog_owner(iprog) /= myid) then
      sendcount(prog_owner(iprog)) = sendcount(prog_owner(iprog)) + 2 + 2 * p2d_links%cnt(iprog)
    endif
  enddo



  !---------------------------------
  ! Share infos across processors
  !---------------------------------

  call MPI_ALLTOALL(sendcount, 1, MPI_INT, receivecount, 1, MPI_INT, MPI_COMM_WORLD, i)


 
  !-------------------------
  ! Send actual data
  !-------------------------

  allocate(sendbuf(1:sum(sendcount)))
  allocate(recvbuf(1:sum(receivecount)))

  ! Fill up sendbuffer
  ind = 1
  do icpu = 1, ncpu
    if (icpu /= myid) then
      do iprog = 1, nprogs

        if (p2d_links%cnt(iprog) > 0 .and. prog_owner(iprog) == icpu) then

          sendbuf(ind) = iprog
          sendbuf(ind+1) = p2d_links%cnt(iprog)
          ind = ind + 2

          idesc = p2d_links%first(iprog)

          !for each descendant:
          do i = 1, p2d_links%cnt(iprog)
            sendbuf(ind) = p2d_links%clmp_id(idesc)
            sendbuf(ind+1) = p2d_links%ntrace(idesc)
            idesc = p2d_links%next(idesc)
            ind = ind + 2
          enddo
        endif
      enddo
    endif
  enddo


  allocate(send_displ(1:ncpu), rec_displ(1:ncpu))
  send_displ = 0; rec_displ = 0;

  do i = 2, ncpu
    send_displ(i) = send_displ(i-1) + sendcount(i-1)
    rec_displ(i) = rec_displ(i-1) + receivecount(i-1)
  enddo

  call MPI_ALLTOALLV(sendbuf, sendcount, send_displ, MPI_INT, &
    recvbuf, receivecount, rec_displ, MPI_INT, MPI_COMM_WORLD, i)



  !----------------------------------------------
  ! Sum up newly arrived values if you are owner
  !----------------------------------------------

  ind = 1
  do while (ind <= sum(receivecount))
    iprog = recvbuf(ind)
    recndesc = recvbuf(ind + 1)
    ind = ind + 2
    do idesc = 1, recndesc
      call fill_matrix(p2d_links, iprog, recvbuf(ind), recvbuf(ind+1), 'add')
      ind = ind + 2 ! skip 2: first is descendant ID, second is # particles
    enddo
  enddo






  !=========================================================
  ! PART 2: 
  ! Broadcast progenitor data from owners to not-owner CPUs
  !=========================================================

  !-----------------------------------------------
  ! Find out how much to send where
  ! Must be done second: Depends on how much data
  ! you received. Same data structure to send.
  !-----------------------------------------------

  allocate(sendcount2(1:ncpu))
  sendcount2 = 0
  allocate(receivecount2(1:ncpu))
  receivecount2 = 0

  ind = 1
  icpu = 1
  do while (ind <= sum(receivecount))
    if (icpu < ncpu) then
      do while(ind > rec_displ(icpu+1))
        icpu = icpu + 1
        if (icpu == ncpu) exit
      enddo
    endif
    iprog = recvbuf(ind)
    recndesc = recvbuf(ind + 1)
    sendcount2(icpu) = sendcount2(icpu) + 2 + 2*p2d_links%cnt(iprog)
    ind = ind + 2 + 2 * recndesc
  enddo



  !---------------------------------
  ! Share infos across processors
  !---------------------------------

  call MPI_ALLTOALL(sendcount2, 1, MPI_INT, receivecount2, 1, MPI_INT, MPI_COMM_WORLD, i)



  !-------------------------
  ! Send actual data
  !-------------------------

  allocate(sendbuf2(1:sum(sendcount2)))
  allocate(recvbuf2(1:sum(receivecount2)))

  ind = 1  ! index for previously received array; 
           ! needed to identify which progenitor to send where
  ind2 = 1 ! index for new array to send stuff
  do while (ind <= sum(receivecount))
    iprog = recvbuf(ind)
    recndesc = recvbuf(ind+1) !needed to skip correctly to next iprog
    ind = ind + 2 + 2 * recndesc

    sendbuf2(ind2) = iprog
    sendbuf2(ind2 + 1) = p2d_links%cnt(iprog)
    ind2 = ind2 + 2

    idesc = p2d_links%first(iprog)
    ! write each descendant
    do i = 1, p2d_links%cnt(iprog)
      sendbuf2(ind2) = p2d_links%clmp_id(idesc)
      sendbuf2(ind2+1) = p2d_links%ntrace(idesc)
      idesc = p2d_links%next(idesc)
      ind2 = ind2 + 2
    enddo
  enddo


  send_displ = 0; rec_displ = 0;

  do i = 2, ncpu
    send_displ(i) = send_displ(i-1) + sendcount2(i-1)
    rec_displ(i) = rec_displ(i-1) + receivecount2(i-1)
  enddo

  call MPI_ALLTOALLV(sendbuf2, sendcount2, send_displ, MPI_INT, &
    recvbuf2, receivecount2, rec_displ, MPI_INT, MPI_COMM_WORLD, i)





  !-----------------------------
  ! Sum up newly arrived values
  !-----------------------------

  ind = 1
  do while (ind <= sum(receivecount2))
    iprog = recvbuf2(ind)
    recndesc = recvbuf2(ind + 1)
    ind = ind + 2
    do idesc = 1, recndesc
      call fill_matrix(p2d_links, iprog, recvbuf2(ind), recvbuf2(ind+1), 'set')
      ind = ind + 2 ! skip 2: first is descendant ID, second is # particles
    enddo
  enddo




  !-------------------
  ! Cleanup
  !-------------------

  deallocate(sendcount, receivecount)
  deallocate(sendcount2, receivecount2)
  deallocate(sendbuf, recvbuf)
  deallocate(sendbuf2, recvbuf2)
  deallocate(send_displ, rec_displ)


  ! end infdef WITHOUTMPI: No sending at all necessary.
#endif 




  !-----------------------------------------------
  ! Create descendants-to-progs ( = d2p) matrix
  !-----------------------------------------------

  allocate(d2p_links%first(1:npeaks_max))
  d2p_links%first = 0
  allocate(d2p_links%cnt(1:npeaks_max))
  d2p_links%cnt = 0
  allocate(d2p_links%clmp_id(1:10*npeaks_max))
  d2p_links%clmp_id = 0
  allocate(d2p_links%ntrace(1:10*npeaks_max))
  d2p_links%ntrace = 0
  allocate(d2p_links%next(1:10*npeaks_max))
  d2p_links%next = 0


  do iprog = 1, nprogs
    ind = p2d_links%first(iprog)
    do i = 1, p2d_links%cnt(iprog)
      call get_local_peak_id(p2d_links%clmp_id(ind), idesc)
      call fill_matrix(d2p_links, idesc, iprog, p2d_links%ntrace(ind), 'add')
      ind = p2d_links%next(ind)
    enddo
  enddo





  call build_peak_communicator()
  call boundary_peak_dp(clmp_mass_exclusive(:))

  return

end subroutine create_prog_desc_links








!===========================
subroutine make_trees()
!===========================

  !---------------------------------------------------------------------------
  ! This subroutine establishes the tree.
  ! First we obtain initial guesses for main progenitors of descendants and
  ! main descendants of progenitors.
  ! Secondly, a loop is performed until each descendant has found a main
  ! progenitor or has no more suitable candidates. In each step of the loop,
  ! all descendant candidates of all progenitors which still haven't found a
  ! matching direct descendant are checkt for a match.
  ! When the loop is over, all descendants that haven't got a main progenitor
  ! are checked for the possibility of containing a past merged progenitor.
  !---------------------------------------------------------------------------

  use clfind_commons

  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer :: i
#endif

  real(dp), dimension(:), allocatable ::  merit_desc
  real(dp), dimension(:), allocatable ::  merit_desc_copy
  logical, dimension(:), allocatable  :: to_iter_prog

  integer :: iprog, ipeak, idl
  integer :: store_id
  integer :: peakshift

  real(dp):: r_null

  logical :: is_first, found, reiter

  ! ! For debug
  ! character(len=80) :: filename
  ! character(len=5)  :: id_to_string, output_to_string



  if (verbose) write(*,*) "ID", myid, "making trees."


  !============================================
  ! PART 1: Preparation and initial guesses
  !============================================

  allocate(main_desc(1:nprogs))
  main_desc = 0
  allocate(merit_desc(1:npeaks_max))
  merit_desc = 0

  allocate(merit_desc_copy(1:npeaks_max))
  allocate(to_iter_prog(1:nprogs))


  !-----------------------------------
  ! Make sure you have necessary
  ! clump data available on this cpu
  !-----------------------------------
  call build_peak_communicator()
  call boundary_peak_dp(clmp_mass_exclusive)


  
  !-------------------------------------------------------------
  ! Find initial guess for main descendant for each progenitor
  !-------------------------------------------------------------

  do iprog = 1, nprogs
    if (p2d_links%cnt(iprog)>0) call find_main_desc(iprog, found)
  enddo

  ! no need to communicate results; data for every instance of prog is identical,
  ! so it should give identical results



  
  !-------------------------------------------------------------
  ! Find initial guess for main progenitor for each descendant
  !-------------------------------------------------------------

  do ipeak = 1, hfree-1
    if (d2p_links%cnt(ipeak) > 0) then
      call find_main_prog(ipeak, merit_desc, found)
    endif
  enddo


#ifndef WITHOUTMPI
  ! Check whether there are better progenitor candidates 
  ! on other processors
  merit_desc_copy = merit_desc
  call virtual_peak_dp(merit_desc(:), 'max')
  call boundary_peak_dp(merit_desc(:))

  do ipeak = 1, hfree-1
    if ( merit_desc_copy(ipeak) /= merit_desc(ipeak)) then
      main_prog(ipeak) = 0
    endif
  enddo
  
  call virtual_peak_int(main_prog(:), 'max')
  call boundary_peak_int(main_prog(:))
#endif



  !------------------------------------------------------------------
  ! Introduce peak shift:
  ! Shift the peak of special cases like past merged progenitors as
  ! main progenitors and mergers by a high number so you still can
  ! use MAX reduction and know that they're special cases
  !------------------------------------------------------------------

  peakshift = 10*(ipeak_start(ncpu)+npeaks_max)

 

    ! if (debug) then
    ! call title(ifout, output_to_string)
    ! call title(myid, id_to_string)
    ! filename = "output_"//output_to_string//"/MATRIXCHECK_PROG_BEFORE_TREE"//id_to_string//".txt"
    ! open(unit=666, file=filename, form='formatted')
    ! write(666,'(A30,x,I9)') "MATRIXCHECK BEFORE TREE ID", myid
    ! do iprog = 1, nprogs
    !   if (p2d_links%cnt(iprog) > 0) then
    !     write(666, '(4(A10,x,I9x),A7,x,E14.6,x,A10,x,I9,8x)', advance='no') &
    !       "Prog:", prog_id(iprog), "local id: ", iprog, "# desc:", p2d_links%cnt(iprog), &
    !       "owner:", prog_owner(iprog), "mass:", prog_mass(iprog), "main desc:", main_desc(iprog)
    !
    !     ipeak = p2d_links%first(iprog)
    !     do i = 1, p2d_links%cnt(iprog)
    !       call get_local_peak_id(p2d_links%clmp_id(ipeak), idl)
    !       write(666, '(A3,x,2(I9,x),A9,x,I9)', advance='no') "D:", p2d_links%clmp_id(ipeak), idl,"tracers:", p2d_links%ntrace(ipeak)
    !       ipeak = p2d_links%next(ipeak)
    !     enddo
    !     write(666,*)
    !   endif
    !
    ! enddo
    ! close(666)
    !
    ! filename = "output_"//output_to_string//"/MATRIXCHECK_DESC_BEFORE_TREE"//id_to_string//".txt"
    ! open(unit=666, file=filename, form='formatted')
    ! write(666,'(A30,x,I9)') "MATRIXCHECK BEFORE TREE ID", myid
    ! do ipeak = 1, hfree-1
    !   if (d2p_links%cnt(ipeak) > 0) then
    !     write(666, '(2(A10,x,I9x),A7,x,E14.6,x,A10,x,I9,8x)', advance='no') &
    !       "Desc:", ipeak, "# progs:", d2p_links%cnt(ipeak), &
    !        "mass:", clmp_mass_exclusive(ipeak), "main prog:", main_prog(ipeak)
    !
    !     iprog = d2p_links%first(ipeak)
    !     do i = 1, d2p_links%cnt(ipeak)
    !       write(666, '(A3,x,I9,x,A9,x,I9)', advance='no') "P:", d2p_links%clmp_id(iprog), "tracers:", d2p_links%ntrace(iprog)
    !       iprog = d2p_links%next(iprog)
    !     enddo
    !     write(666,*)
    !   endif
    !
    ! enddo
    ! close(666)
    ! endif


  !==================================
  ! PART 2: TREEMAKING LOOP
  !==================================

  !------------------------------------------------------------------------------
  ! To establish a connection over snapshots, the main descendant of each
  ! progenitor must have said progenitor as the main progenitor.
  ! In each iteration, progenitors which haven't found such a match go through
  ! all possible descendant candidates, trying to find a match. If no such 
  ! match is found, it is assumed that the progenitor merged into the initial
  ! best guess.
  ! Then descendants are checked for having the correct main progenitor
  ! identified. If this isn't the case, the guess is moved to the next best
  ! candidate (doesn't loop through all of them!), after which the entire loop
  ! is restarted.
  !------------------------------------------------------------------------------

  ! initialise what needs to be checked. Both are arrays!
  to_iter = (main_prog > 0)
  to_iter_prog = (main_desc > 0)

  reiter = .true.

  do while (reiter)

    reiter = .false.

    !----------------------------------------------------------------------
    ! update to_iter_prog array for every loop: Check whether they still
    ! need to be updated
    !----------------------------------------------------------------------

    do iprog = 1, nprogs
      if (to_iter_prog(iprog)) then
        if (main_desc(iprog)>0) then
          ! to_iter_prog(iprog): only re-check those that aren't finished
          ! main_desc(iprog) > 0: only check those that aren't considered to have merged
          ! (mergers need to be re-checked every time)
          call get_local_peak_id(main_desc(iprog), idl) 
          if (iprog==main_prog(idl)) then
            to_iter_prog(iprog) = .false.
          endif
        endif
      endif
    enddo


    !--------------------------------------------------------------
    ! Check all candidates of still available progenitors. 
    ! If no match is found, mark the progenitor as merged into the
    ! best candidate.
    !--------------------------------------------------------------

    is_first = .true.
    iprog = 1

    do while (iprog <= nprogs)

      if (to_iter_prog(iprog)) then
        if (is_first) then
          ! store id of best initial match; only for first iteration of each prog per loop
          store_id = abs(main_desc(iprog))

          ! reset values of all candidates for loop
          ! candidates that have been checked already will get negative tracer numbers
          if (p2d_links%cnt(iprog) > 0) then
            ipeak = p2d_links%first(iprog)
            do i = 1, p2d_links%cnt(iprog)
              if (p2d_links%ntrace(ipeak) < 0) p2d_links%ntrace(ipeak)=-p2d_links%ntrace(ipeak)
              ipeak = p2d_links%next(ipeak)
            enddo
          endif
        endif

        is_first = .true.
        to_iter_prog(iprog) = .false.

        ! take abs(main_desc) as key: mergers will have negative main_desc, and
        ! fill_matrix() will add a new descendant instead of inverting the correct one
        call fill_matrix(p2d_links, iprog, abs(main_desc(iprog)), 0, 'inv')
        call find_main_desc(iprog, found)

        ! if you found another candidate:
        if (found) then
          call get_local_peak_id(main_desc(iprog), idl)
          if (iprog/=main_prog(idl)) then
            ! if the new candidate still doesn't match:
            ! need to re-iterate
            ! else: to_iter_prog is set to = .false. anyways
            to_iter_prog(iprog) = .true.
            is_first = .false.
          endif
          iprog = iprog - 1 ! check again
        else
          ! if nothing found, assume progenitor merged into descendant
          main_desc(iprog) = store_id + peakshift
          ! check this one again next round
          to_iter_prog(iprog) = .true.
        endif

      endif ! to_iter_prog
      iprog = iprog + 1

    enddo



#ifndef WITHOUTMPI
    !---------------------------------------------------------------------------------
    ! Communicate results. All CPUs that have tracer particles of any progenitor
    ! have full data of that progenitor; But it might be that there is a descendant,
    ! split among multiple CPU's, where some CPU's have missing progenitor data  
    ! because they don't have any progenitor's particles on their domain. Such
    ! descendants will need to know main_desc values which will otherwise be unknown.
    ! First reset values: If you had a main_desc with a higher ID previously, you'll
    ! get junk results. Unlike with descendants, keeping only main_desc of the owner 
    ! CPU is safe, as all CPUs that have tracer particles of a progenitor have full 
    ! data of that progenitor.
    !---------------------------------------------------------------------------------
    do iprog = 1, nprogs
      if (prog_owner(iprog)/=myid) main_desc(iprog) = 0
    enddo
    call MPI_ALLREDUCE(MPI_IN_PLACE, main_desc, nprogs, MPI_INT, MPI_MAX, MPI_COMM_WORLD, i)
#endif

    ! revert peakshift
    do iprog = 1, nprogs
      if (main_desc(iprog)>peakshift) main_desc(iprog) = -(main_desc(iprog)-peakshift) ! make it negative!
    enddo




    !------------------------------------------------------------
    ! Check descendants. If no match found, try next best 
    ! candidate.
    !------------------------------------------------------------

    merit_desc = 0 ! is array!

    do ipeak = 1, hfree-1
      if ( to_iter(ipeak) ) then ! if there is something to check for

        to_iter(ipeak) = .false.
        
        if (main_desc(main_prog(ipeak))/=0) then
          ! abs needed here: mergers are signified by a negative main descendant ID
          call get_local_peak_id(abs(main_desc(main_prog(ipeak))), idl)
        else
          idl=0
        endif
        
        if (ipeak /= idl) then
          ! if this descendant is not main descendant
          ! of its own main progenitor, look for next best candidate
          ! first reset value for current best fit
          call fill_matrix(d2p_links, ipeak, main_prog(ipeak), 0, 'inv')
          ! now look for better candidate
          call find_main_prog(ipeak, merit_desc, found)
          if (.not.found) then 
            ! if you run out of candidates:
            main_prog(ipeak) = 0
          else 
            ! check this peak for matches again
            to_iter(ipeak) = .true.
            reiter = .true.
          endif

        endif !ipeak /= idl
      endif
    enddo



#ifndef WITHOUTMPI
    !---------------------------------------------
    ! Communicate results
    !---------------------------------------------

    ! communicate merits, find max
    merit_desc_copy = merit_desc
    call build_peak_communicator()
    call virtual_peak_dp(merit_desc(:), 'max')
    call boundary_peak_dp(merit_desc(:))

    do ipeak = 1, hfree-1
      ! if you didn't have the max, reset stuff
      if (merit_desc_copy(ipeak) /= merit_desc(ipeak)) then
        ! reset mark in matrix for later use
        call fill_matrix(d2p_links, ipeak, main_prog(ipeak), 0, 'inv')
        ! reset value in array
        main_prog(ipeak) = 0
      endif
    enddo

    ! call build_peak_communicator()
    call virtual_peak_int(main_prog, 'max')
    call boundary_peak_int(main_prog)

    do ipeak = 1, hfree-1
      ! check whether you need to reiterate peak first, while you have data
      ! synchronized globally
      if (to_iter(ipeak)) then
        if (main_desc(main_prog(ipeak))/=0) then
          ! abs needed here: mergers are signified by a negative main descendant ID
          call get_local_peak_id(abs(main_desc(main_prog(ipeak))), idl)
        else
          idl=0
        endif
        if (ipeak == idl) to_iter(ipeak) = .false.
      endif
    enddo

    ! check globally whether you need to reiterate treebuilding loop
    call MPI_ALLREDUCE(MPI_IN_PLACE, reiter, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, i)
#endif
  enddo
  !-------------------------
  ! End of treemaking loop
  !-------------------------
 




  !===================================================
  ! PART 3: Look for progenitors in older snapshots
  !===================================================

  !---------------------------------------------------------------------------
  ! After tree is made, add merged progenitors to past merged progenitors
  !---------------------------------------------------------------------------
  do iprog = 1, nprogs
    if (main_desc(iprog)<0) then
      if (prog_owner(iprog) == myid) call add_new_pmprog(iprog)
    endif
  enddo



  !---------------------------------------------------------------------------------
  ! update all necessary unbinding arrays in case you introduced new virtual peaks
  !---------------------------------------------------------------------------------
  do i = 1, nmassbins
    call boundary_peak_dp(cmp(1,i))
  enddo

  do i = 1, 3
    call boundary_peak_dp(peak_pos(1,i))
    call boundary_peak_dp(clmp_vel_pb(1,i))
  enddo

  call boundary_peak_dp(cmp_distances(1,nmassbins))

  ! recompute cumulative mass profile bin distances for virtual peaks
  do ipeak=npeaks+1, hfree-1

    ! set up distances
    if (cmp_distances(ipeak, nmassbins)>0) then
      if (logbins) then
        do i=1, nmassbins-1
          ! rmin is declared in clfind_commons
          cmp_distances(ipeak,i)=rmin*(cmp_distances(ipeak,nmassbins)/rmin)**(real(i)/real(nmassbins))
        enddo
      else ! linear binnings
        r_null=cmp_distances(ipeak,nmassbins)/real(nmassbins)
        do i=0, nmassbins-1
          cmp_distances(ipeak,i)=r_null*i
        enddo
      endif
    endif
  enddo



  !------------------------------------------------------------------
  ! If you still haven't found a progenitor for a descendant, 
  ! try looking in previous, non-adjacent snapshots.
  !------------------------------------------------------------------

  ! reset merit and to_iter arrays
  merit_desc = 0
  to_iter = (clmp_mass_exclusive > 0) .and. (main_prog == 0)

  do ipeak = 1, hfree-1
    if ( to_iter(ipeak) ) then !if there is something to check for
        call find_prog_in_older_snapshots(ipeak, peakshift, merit_desc)
    endif
  enddo




  !-------------------------------------------------------
  ! Check that you found best candidate across processors
  !-------------------------------------------------------
  merit_desc_copy = merit_desc
  ! Here the best candidate has the lowest merit value! 
  call virtual_peak_dp(merit_desc(:), 'min')
  call boundary_peak_dp(merit_desc(:))

  do ipeak = 1, hfree-1
    if ( to_iter(ipeak) .and. (merit_desc_copy(ipeak) /= merit_desc(ipeak)) ) main_prog(ipeak) = 0
  enddo


  ! Communicate results
  call virtual_peak_int(main_prog, 'max')
  call boundary_peak_int(main_prog)

  ! Mark progenitors from earlier snapshots as such
  do ipeak = 1, hfree -1
    if (main_prog(ipeak) > peakshift ) then
      main_prog(ipeak) = - (main_prog(ipeak)-peakshift)
    endif
  enddo




  !-------------------------
  ! Cleanup before finish
  !-------------------------

  deallocate(merit_desc_copy)
  deallocate(merit_desc)
  deallocate(to_iter_prog)



    ! if (debug) then
    !   call title(ifout, output_to_string)
    !   call title(myid, id_to_string)
    !   filename = "output_"//output_to_string//"/MATRIXCHECK_PROG_AFTER_TREE"//id_to_string//".txt"
    !   open(unit=666, file=filename, form='formatted')
    !   write(666,'(A30,x,I9)') "MATRIXCHECK AFTER TREE ID", myid
    !
    !   do iprog = 1, nprogs
    !     if (p2d_links%cnt(iprog) > 0) then
    !       write(666, '(4(A10,x,I9x),A7,x,E14.6,x,A10,x,I9,8x)', advance='no') &
    !         "Prog:", prog_id(iprog), "local id: ", iprog, "# desc:", p2d_links%cnt(iprog), &
    !         "owner:", prog_owner(iprog), "mass:", prog_mass(iprog), "main desc:", main_desc(iprog)
    !
    !       ipeak = p2d_links%first(iprog)
    !       do i = 1, p2d_links%cnt(iprog)
    !         call get_local_peak_id(p2d_links%clmp_id(ipeak), idl)
    !         write(666, '(A3,x,2(I9,x),A9,x,I9)', advance='no') "D:", p2d_links%clmp_id(ipeak), idl,"tracers:", p2d_links%ntrace(ipeak)
    !         ipeak = p2d_links%next(ipeak)
    !       enddo
    !       write(666,*)
    !     endif
    !   enddo
    !
    !   close(666)
    !
    !   filename = "output_"//output_to_string//"/MATRIXCHECK_DESC_AFTER_TREE"//id_to_string//".txt"
    !   open(unit=666, file=filename, form='formatted')
    !   write(666,'(A30,x,I9)') "MATRIXCHECK AFTER TREE ID", myid
    !
    !   do ipeak = 1, hfree-1
    !     if (d2p_links%cnt(ipeak) > 0) then
    !       write(666, '(2(A10,x,I9x),A7,x,E14.6,x,A10,x,I9,8x)', advance='no') &
    !         "Desc:", ipeak, "# progs:", d2p_links%cnt(ipeak), &
    !          "mass:", clmp_mass_exclusive(ipeak), "main prog", main_prog(ipeak)
    !
    !       iprog = d2p_links%first(ipeak)
    !       do i = 1, d2p_links%cnt(ipeak)
    !         write(666, '(A3,x,I9,x,A9,x,I9)', advance='no') "P:", d2p_links%clmp_id(iprog), "tracers:", d2p_links%ntrace(iprog)
    !         iprog = d2p_links%next(iprog)
    !       enddo
    !       write(666,*)
    !     endif
    !   enddo
    !
    !   close(666)
    ! endif
    !





  contains 

    !==============================================
    subroutine find_main_desc(iprog, found_one)
    !==============================================

      !-------------------------------------------------
      ! Find the main descendant of progenitor iprog
      !-------------------------------------------------

      use clfind_commons
      
      implicit none
      integer, intent(in) :: iprog
      logical, intent(out) :: found_one
      real(dp):: merit_max, merit_calc, a 
      integer :: merit_max_id

      integer :: ind, i, idesc, idl
      

      merit_max = 0
      merit_max_id = 0

      ind = p2d_links%first(iprog)

      ! loop over all descendants
      do i = 1, p2d_links%cnt(iprog)
        
        ! if there are particles left 
        ! (might have been removed during matching)
        if (p2d_links%ntrace(ind) > 0) then
          idesc = p2d_links%clmp_id(ind)
          call get_local_peak_id(idesc, idl)

          ! calculate merits
          if (clmp_mass_exclusive(idl) > prog_mass(iprog)) then
            a = abs(1 - clmp_mass_exclusive(idl)/prog_mass(iprog))
          else
            a = abs(1 - prog_mass(iprog)/clmp_mass_exclusive(idl))
          endif
          if (a<1d-100) a = 1d-100 ! don't devide by zero
          merit_calc =  real(p2d_links%ntrace(ind)) / a**2

          if (merit_calc > merit_max) then
            merit_max = merit_calc
            merit_max_id = idesc
          endif
        endif

        ind = p2d_links%next(ind)

      enddo

      if (merit_max_id > 0) then
        found_one = .true.
        main_desc(iprog) = merit_max_id
      else
        found_one = .false.
        main_desc(iprog) = -7
        ! essentially this means there is no descendant for this progenitor.
        ! so.. disappearing clumps? Might be possible if you had very small
        ! clumps that dissolved over time
      endif

      return

    end subroutine find_main_desc





    !===========================================================
    subroutine find_main_prog(ipeak, merit_desc, found_one)
    !===========================================================

      !-------------------------------------------------
      ! Find the main progenitor of descendant idesc
      !-------------------------------------------------

      use clfind_commons
      
      implicit none
      integer, intent(in) :: ipeak
      real(dp), dimension(1:npeaks_max), intent(inout):: merit_desc
      logical, intent(out) :: found_one

      real(dp):: merit_max, merit_calc, a 
      integer :: merit_max_id

      integer :: ind, i, iprog

      merit_max = 0
      merit_max_id = 0

      ind = d2p_links%first(ipeak)

      do i = 1, d2p_links%cnt(ipeak)
        
        iprog = d2p_links%clmp_id(ind)

        if (d2p_links%ntrace(ind) > 0) then
          ! calculate merit
          if (clmp_mass_exclusive(ipeak) > prog_mass(iprog)) then
            a = abs(1 - clmp_mass_exclusive(ipeak)/prog_mass(iprog))
          else
            a = abs(1 - prog_mass(iprog)/clmp_mass_exclusive(ipeak))
          endif
          if (a<1d-100) a = 1d-100 ! don't devide by zero
          merit_calc = real(d2p_links%ntrace(ind))/a**2

          if (merit_calc > merit_max) then
            merit_max = merit_calc
            merit_max_id = iprog
          endif
        endif

        ind = d2p_links%next(ind)

      enddo

      if (merit_max_id > 0) then
        ! save results
        main_prog(ipeak) = merit_max_id
        merit_desc(ipeak) = merit_max
        found_one = .true.
      else
        main_prog(ipeak) = 0
        found_one = .false.
      endif

      return

    end subroutine find_main_prog






    !=====================================
    subroutine add_new_pmprog(iprog)
    !=====================================

      !-----------------------------------------------------------
      ! Add new progenitor to the list of past merged progenitors
      !-----------------------------------------------------------

      use clfind_commons

      integer, intent(in) :: iprog ! local prog index to add

      pmprogs(pmprog_free) = prog_id(iprog)
      pmprogs_galaxy(pmprog_free) = galaxy_tracers(iprog)
      pmprogs_t(pmprog_free) = ifout-1 ! prog was active clump for the last time at this timestep 
      pmprogs_owner(pmprog_free) = myid
      pmprogs_mass(pmprog_free) = prog_mass(iprog)
      pmprog_free = pmprog_free + 1

      return

    end subroutine add_new_pmprog


    




    !=====================================================================
    subroutine find_prog_in_older_snapshots(ipeak, peakshift, merit_desc)
    !=====================================================================

      !--------------------------------------------------------
      ! If no progenitor from adjacent snapshot has been found,
      ! this subroutine looks for a possible progenitor in an
      ! earlier snapshot.
      ! Here, the merit is to be minimized: It is the binding
      ! energy of possible candidates
      !--------------------------------------------------------

      use clfind_commons
      use pm_commons, only: idp

      implicit none

      integer, intent(in) :: ipeak
      integer, intent(in) :: peakshift
      real(dp), intent(inout), dimension(1:npeaks_max) :: merit_desc

      integer, dimension(:), allocatable :: particlelist   ! list of particle IDs of this clump
      integer, dimension(:), allocatable :: canddts        ! list of progenitor candidates
      real(dp),dimension(:), allocatable :: merit          ! merit of progenitor candidates
      integer, dimension(:), allocatable :: part_local_ind ! local particle index for clumpparticles
      integer :: ncand                                     ! number of candidates

      integer :: ipart, thispart, iclump, ipastprog, lpcid, merit_min_id
      real(dp) :: merit_min


      if (nclmppart(ipeak) > 0) then
        !-------------
        ! Set up
        !-------------

        ncand = 0

        allocate(canddts(1:npastprogs_max))
        allocate(merit(1:npastprogs_max))

        allocate(particlelist(1:nclmppart(ipeak)))
        particlelist = 0
        allocate(part_local_ind(1:nclmppart(ipeak)))


        ! Get particle ID list for this clump
        thispart = clmppart_first(ipeak)
        do ipart = 1, nclmppart(ipeak)
          if (clmpidp(thispart) > 0) then
            call get_local_peak_id(clmpidp(thispart), lpcid)
            if (ipeak == lpcid) then
              particlelist(ipart) = idp(thispart)
              part_local_ind(ipart) = thispart
            endif
          endif
          thispart = clmppart_next(thispart)
        enddo
        
        ! Sort the list for matching with galaxy particles
        call quick_sort_int_int(particlelist, part_local_ind, nclmppart(ipeak))

        ! Find starting values
        iclump = 1
        do ipart = 1, nclmppart(ipeak)
          if (particlelist(ipart) > 0) then
            iclump = ipart
            exit
          endif
        enddo

        ipastprog = 1



        !-------------------------------------------------------------------
        ! Find candidates by checking if any particle in clump is a galaxy
        ! of a past merged progenitor
        !-------------------------------------------------------------------

        call compute_phi(ipeak)

        do while (iclump <= nclmppart(ipeak) .and. ipastprog <= npastprogs)
          
          if (particlelist(iclump) < pmprogs_galaxy(ipastprog)) then
            iclump = iclump + 1
          else if (particlelist(iclump) > pmprogs_galaxy(ipastprog)) then
            ipastprog = ipastprog + 1
          else
            ! Found a match!
            ncand = ncand + 1
            canddts(ncand) = ipastprog ! store local ID
            call eparttot(ipeak, part_local_ind(iclump), merit(ncand))
            ipastprog = ipastprog + 1
            iclump = iclump + 1
          endif
        enddo
        


        ! Find the one with actual minimal merit
        merit_min_id = 0
        merit_min = HUGE(0.d0)
        do ipastprog = 1, ncand
          if (merit(ipastprog) < merit_min) then
            merit_min = merit(ipastprog)
            merit_min_id = canddts(ipastprog)
          endif
        enddo

        ! If there is one, store values
        if (merit_min_id > 0) then
          main_prog(ipeak) = pmprogs(merit_min_id) + peakshift
          merit_desc(ipeak) = merit_min
        endif

        deallocate(particlelist, canddts, merit, part_local_ind)
      endif

      return

    end subroutine find_prog_in_older_snapshots

end subroutine make_trees








!====================================
subroutine read_progenitor_data()
!====================================
  
  !---------------------------------
  ! Reads in all progenitor data
  !---------------------------------

  use clfind_commons
  use amr_commons

  implicit none

#ifndef WITHOUTMPI
  include 'mpif.h'
#endif


  integer :: iprog, i

  logical :: exists
  integer :: prog_read, prog_read_local, startind, tracer_free, nprogs_to_read, progcount_to_read, np
  integer, allocatable, dimension(:) :: read_buffer   ! temporary array for reading in data
  real(dp),allocatable, dimension(:) :: read_buffer_2 ! temporary array for reading in data

  character(LEN=80)     :: fileloc
  character(LEN=5)      :: output_to_string


#ifndef WITHOUTMPI
  integer, dimension (1:MPI_STATUS_SIZE):: state
  integer :: mpi_err, filehandle
  integer, dimension(1:4) :: buf
#endif


  if (verbose) write(*,*) " Calling read progenitor data."

  call title(ifout-1, output_to_string)
  ! ifout -1: read from previous output!
  nprogs = 0
  nprogs_to_read = 0
  progcount_to_read = 0






  !========================
  ! Read progenitor counts
  !========================

  if (myid == 1) then ! read in stuff
    
    !-------------------------------------------------------------
    ! Current progenitor counts
    ! Both of these count files need to be present in any case.
    !-------------------------------------------------------------

    fileloc=TRIM('output_'//TRIM(output_to_string)//'/progenitorcount.dat')

    ! check that file exists
    inquire(file=fileloc, exist=exists)
    if (.not.exists) then
      write(*, *) "ID", myid, "didn't find file ", fileloc
      stop
    endif


    open(unit=666,file=fileloc,form='unformatted')
    read(666) nprogs, nprogs_to_read, progcount_to_read, npastprogs
    ! open(unit=666,file=fileloc,form='formatted')
    ! read(666, '(2(I7,x))') nprogs, nprogs_to_read
    close(666)


  endif


#ifndef WITHOUTMPI
    buf = (/nprogs, progcount_to_read, nprogs_to_read, npastprogs/)
    call MPI_BCAST(buf, 4, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
    nprogs = buf(1)
    progcount_to_read = buf(2)
    nprogs_to_read = buf(3)
    npastprogs = buf(4)
#endif
  




  !==================================
  ! READ CURRENT PROGENITOR DATA
  !==================================

  !---------------------------
  ! Allocate arrays
  !---------------------------

  allocate(prog_id(1:nprogs))
  prog_id = 0       ! list of progenitor global IDs
  prog_free = 1     ! first free local progenitor ID

  allocate(prog_owner(1:nprogs))
  prog_owner = -1

  allocate(prog_mass(1:nprogs))
  prog_mass = 0

  allocate(tracers_all(1:nprogs*nmost_bound))
  tracers_all = 0

  allocate(galaxy_tracers(1:nprogs))

  allocate(tracer_loc_progids_all(1:nprogs*nmost_bound))
  tracer_loc_progids_all = 0
  tracer_free = 1   ! first free local tracer index




  if (nprogs > 0) then

    !--------------------------------
    ! Read progenitor particles
    !--------------------------------

    fileloc=TRIM('output_'//TRIM(output_to_string)//'/progenitor_data.dat')

    inquire(file=fileloc, exist=exists)
    if (.not.exists) then
      write(*,*) "ID", myid, "didn't find file ", fileloc
      stop
    endif


    allocate(read_buffer(1:progcount_to_read))

#ifndef WITHOUTMPI
    call MPI_FILE_OPEN(MPI_COMM_WORLD, fileloc, &
      MPI_MODE_RDONLY, MPI_INFO_NULL,filehandle, mpi_err)
    call MPI_FILE_READ(filehandle, read_buffer, &
      progcount_to_read, MPI_INTEGER, state, mpi_err)
    call MPI_FILE_CLOSE(filehandle, mpi_err)
#else
    open(unit=666,file=fileloc,form='unformatted')
    read(666) read_buffer
    ! open(unit=666,file=fileloc,form='formatted')
    ! i = 1
    ! do iprog = 1, nprogs_to_read
    !   write(*,*) "read ", i, i+250, iprog
    !   read(666, '(251(I7,x))') read_buffer(i:i+250)
    !   i = i + 251
    ! enddo
    close(666)
#endif





    !--------------------------------
    ! Read progenitor masses
    !--------------------------------

    fileloc=TRIM('output_'//TRIM(output_to_string)//'/progenitor_mass.dat')

    inquire(file=fileloc, exist=exists)
    if (.not.exists) then
      write(*, *) "ID", myid, "didn't find file ", fileloc
      stop
    endif


    allocate(read_buffer_2(1:nprogs_to_read))

#ifndef WITHOUTMPI
    call MPI_FILE_OPEN(MPI_COMM_WORLD, fileloc, MPI_MODE_RDONLY, MPI_INFO_NULL, filehandle, mpi_err)
    call MPI_FILE_READ(filehandle, read_buffer_2, nprogs_to_read, MPI_DOUBLE_PRECISION, state, mpi_err)
    call MPI_FILE_CLOSE(filehandle, mpi_err)
#else
    open(unit=666,file=fileloc,form='unformatted')
    read(666) read_buffer_2
    ! open(unit=666,file=fileloc,form='formatted')
    ! do iprog = 1, nprogs_to_read
    !   read(666, '(E14.6,x)') read_buffer_2(iprog)
    ! enddo
    close(666)
#endif






    !----------------------------------
    ! Sort out the data you just read
    !----------------------------------

    tracer_free = 1

    iprog = 1
    startind = 1

    do while (startind <= progcount_to_read)

      prog_read = read_buffer(startind)
      np = read_buffer(startind + 1)

      ! get local instead global ID in prog_read (past tense "read")
      call get_local_prog_id(prog_read, prog_read_local)

      prog_mass(prog_read_local) = read_buffer_2(iprog)

      do i = startind+2, startind+1+np
        if (read_buffer(i) > 0) then
          tracers_all(tracer_free) = read_buffer(i)               ! add new tracer particle
          tracer_loc_progids_all(tracer_free) = prog_read_local   ! write which progenitor tracer belongs to
          tracer_free = tracer_free + 1                           ! raise index for next tracer
        else 
          ! found a galaxy particle
          tracers_all(tracer_free) = -read_buffer(i)              ! add new tracer particle
          galaxy_tracers(prog_read_local) = -read_buffer(i)       ! add new galaxy tracer
          tracer_loc_progids_all(tracer_free) = prog_read_local   ! write which progenitor tracer belongs to
          tracer_free = tracer_free + 1                           ! raise index for next tracer
        endif
      enddo

      iprog = iprog + 1
      startind = startind + 2 + np

    enddo

    deallocate(read_buffer, read_buffer_2)

  endif ! nprogs > 0










  !========================================
  ! READ PAST PROGENITOR DATA
  !========================================
  
  !-------------------------
  ! Allocate arrays
  !-------------------------

  ! overestimate size to fit new ones if necessary
  npastprogs_max = npastprogs + nprogs
  pmprog_free = npastprogs + 1


  ! Past Merged Progenitors for multi-snapshot matching
  allocate(pmprogs(1:npastprogs_max))
  pmprogs = 0 

  ! Current owner Merged Progenitors' galaxy particles
  allocate(pmprogs_owner(1:npastprogs_max))
  pmprogs_owner = 0

  ! Past Merged Progenitors' galaxy particles
  allocate(pmprogs_galaxy(1:npastprogs_max))
  pmprogs_galaxy = 0

  ! Time at which past progenitors have been merged (= ifout at merging time)
  allocate(pmprogs_t(1:npastprogs_max))
  pmprogs_t = 0

  ! past merged progenitor mass
  allocate(pmprogs_mass(1:npastprogs_max))
  pmprogs_mass = 0


  

  if (npastprogs > 0) then

    !-----------------------------
    ! Read in data
    !-----------------------------

    allocate(read_buffer(1:3*npastprogs))

      fileloc=TRIM('output_'//TRIM(output_to_string)//'/past_merged_progenitors.dat')

#ifndef WITHOUTMPI
    call MPI_FILE_OPEN(MPI_COMM_WORLD, fileloc, &
      MPI_MODE_RDONLY, MPI_INFO_NULL,filehandle, mpi_err)
    call MPI_FILE_READ(filehandle, read_buffer, &
      npastprogs*3, MPI_INTEGER, state, mpi_err)
    call MPI_FILE_CLOSE(filehandle, mpi_err)
#else
    open(unit=666,file=fileloc,form='unformatted')
    read(666) read_buffer
    ! open(unit=666,file=fileloc,form='formatted')
    ! i = 1
    ! do iprog = 1, npastprogs
    !   read(666, '(3(I7,x))') read_buffer(i:i+3)
    !   write(*,*) "Read past prog", read_buffer(i)
    !   i = i + 3
    ! enddo
    close(666)
#endif




    !---------------------------------
    ! Read past progenitor's masses
    !---------------------------------

    fileloc=TRIM('output_'//TRIM(output_to_string)//'/past_merged_progenitor_mass.dat')

    inquire(file=fileloc, exist=exists)
    if (.not.exists) then
      write(*, *) "ID", myid, "didn't find file ", fileloc
      stop
    endif



#ifndef WITHOUTMPI
    call MPI_FILE_OPEN(MPI_COMM_WORLD, fileloc, MPI_MODE_RDONLY, MPI_INFO_NULL, filehandle, mpi_err)
    call MPI_FILE_READ(filehandle, pmprogs_mass(1:npastprogs), npastprogs, MPI_DOUBLE_PRECISION, state, mpi_err)
    call MPI_FILE_CLOSE(filehandle, mpi_err)
#else
    open(unit=666,file=fileloc,form='unformatted')
    read(666) pmprogs_mass(1:npastprogs)
    ! open(unit=666,file=fileloc,form='formatted')
    ! do iprog = 1, nprogs_to_read
    !   read(666, '(E14.6,x)') read_buffer_2(iprog)
    ! enddo
    close(666)
#endif





    !----------------------------------
    ! Sort out the data you just read
    !----------------------------------

    pmprog_free = 1
    iprog = 1
    do while (iprog <= 3*npastprogs)
      pmprogs(pmprog_free) = read_buffer(iprog)
      pmprogs_galaxy(pmprog_free) = read_buffer(iprog + 1)
      pmprogs_t(pmprog_free) = read_buffer(iprog + 2)
      iprog = iprog + 3
      pmprog_free = pmprog_free + 1
    enddo


    deallocate(read_buffer)


  endif ! npastprogs > 0

end subroutine read_progenitor_data







!=============================
subroutine write_trees()
!=============================

  !-------------------------------
  ! Write the tree to file
  !-------------------------------

  use clfind_commons

  implicit none

#ifndef WITHOUTMPI
  include 'mpif.h'
  integer :: err
#endif

  character (len=5)  :: dir, idnr
  character (len=80) :: fileloc
  integer:: ipeak, iprog, ipastprog, snapshot_nr

  integer, dimension(1:nprogs) :: printed

  if (verbose) write(*,*) " writing trees to file."

  printed = 0

  call title(ifout, dir)
  call title(myid, idnr) 
  fileloc=TRIM('output_'//TRIM(dir)//'/mergertree.txt'//TRIM(idnr))

  open(unit=666,file=fileloc,form='formatted')
  write(666, '(4(A15))') "clump", "progenitor", "prog outputnr", "case"
  !----------------------------------
  ! Possible cases:
  ! 1: adjacent link found
  ! 2: no link or new clump found
  ! 3: non-adjacent link found
  ! 4: merging detected
  !----------------------------------


  do ipeak = 1, npeaks
    if (clmp_mass_exclusive(ipeak) > 0) then
      
      !----------------------
      ! Adjacent link found
      !----------------------
      if (main_prog(ipeak) > 0 ) then
        write(666,'(4(I15))') &
          ipeak+ipeak_start(myid), prog_id(main_prog(ipeak)), &
          ifout-1, 1
        printed(main_prog(ipeak)) = 1


      !------------------------------
      ! No link or new clump found
      !------------------------------
      else if (main_prog(ipeak) == 0) then
        write(666,'(4(I15))') &
          ipeak+ipeak_start(myid), 0, &
          ifout-1, 2


      !----------------------------------------------
      ! Progenitor from non-adjacent snapshot found
      !----------------------------------------------
      else 
        ! First find its snapshot number
        snapshot_nr = 0
        do ipastprog = 1, npastprogs
          if (pmprogs(ipastprog) == -main_prog(ipeak)) then
            snapshot_nr = pmprogs_t(ipastprog)
            pmprogs_owner(ipastprog) = 0
            exit
          endif
        enddo

        write(666,'(4(I15),x,E14.6)') &
          ipeak+ipeak_start(myid), main_prog(ipeak), &
          snapshot_nr, 3, clmp_mass_exclusive(ipeak)

      endif
    endif
  enddo


#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(MPI_IN_PLACE, printed, nprogs, MPI_INT, MPI_SUM, MPI_COMM_WORLD, err)
#endif



  !-----------------------------
  ! Merged progenitors
  !-----------------------------

  ! Print whatever hasn't been printed yet.
  ! If a progenitor merged with another into a composite descendant, the non-main progenitor
  ! will have a negative main descendant ID.
  do iprog = 1, nprogs
    if ( printed(iprog) == 0 .and. main_desc(iprog) /= 0 .and. prog_owner(iprog) == myid ) then
      write(666, '(4(I15))') main_desc(iprog), prog_id(iprog), ifout-1, 4
    endif
  enddo


  close(666)

end subroutine write_trees






!======================================
subroutine write_progenitor_data()
!======================================

  !--------------------------------------------
  ! This subroutine writes the formatted and
  ! unformatted output needed for mergertrees.
  !--------------------------------------------

  use clfind_commons
  use amr_commons
  use pm_commons, only: idp

  implicit none

  integer :: ipeak, ipart, pind, startind, first_bound, partcount

  integer :: ihalo, haloid
  integer, allocatable, dimension(:)  :: particlelist, pastproglist
  real(dp), allocatable, dimension(:) :: masslist, pastprogmasslist

  character(LEN=80)     :: fileloc
  character(LEN=5)      :: output_to_string, id_to_string 

#ifndef WITHOUTMPI
  include 'mpif.h'
  integer, dimension (1:MPI_STATUS_SIZE):: state
  integer :: mpi_err, filehandle, npastprogs_all
  integer, dimension(1:4) :: buf
#endif



  !=======================================================
  ! Write current clumps as progenitors for next snapshot 
  !=======================================================

  call title(ifout, output_to_string)
  call title(myid, id_to_string)

  !------------------------------------------------------------------------------
  ! Prepare particle list
  ! Format: haloID, number of tracers written, [specified nr of particle IDs]
  !------------------------------------------------------------------------------

  ihalo = 0           ! count how many clumps you're actually writing 
  progenitorcount = 0
  pind = 1            ! index where to write

  ! progenitorcount_written is local to each CPU at this point
  ! and overestimated to surely have enough array length

  allocate(particlelist(1:progenitorcount_written*(nmost_bound+2)))
  allocate(masslist(1:progenitorcount_written))

  if (progenitorcount_written > 0) then
    particlelist = 0
    masslist = 0

    do ipeak = 1, hfree-1
      !write only peaks that have most bound particles
      if (clmp_mass_exclusive(ipeak) > 0 ) then

        haloid = 0
        first_bound = 0
        partcount = 0
        if (ipeak <= npeaks) progenitorcount = progenitorcount + 1 ! count only non-virtuals

        do ipart = 1, nmost_bound
          ! check if there are mostbound particles for this clump on this CPU 
          if (most_bound_pid(ipeak, ipart) > 0) then
            first_bound = ipart
            haloid = -clmpidp(most_bound_pid(ipeak, ipart))
            exit
          endif
        enddo

        if (first_bound > 0) then
          ihalo = ihalo + 1

          startind = pind
          particlelist(startind) = haloid
          pind = pind + 2
          masslist(ihalo) = clmp_mass_exclusive(ipeak)

          ! loop over mostbound particle list

          do ipart = first_bound, nmost_bound
            ! write only halos that have mbp on this proc
            if (most_bound_pid(ipeak, ipart) > 0) then
              partcount = partcount + 1
              if (ipart == 1) then ! if galaxy particle: mark it
                particlelist(pind) = -idp(most_bound_pid(ipeak, ipart))
              else !if regular most bound particle
                particlelist(pind) = idp(most_bound_pid(ipeak, ipart))
              endif
              pind = pind + 1

            endif
          enddo

          particlelist(startind + 1) = partcount
        endif
      endif !clump mass > 0
    enddo ! loop over clumps


    ! reset written progenitorcount: Too small clumps might've 
    ! been dissolved after being counted.
    progenitorcount_written = pind-1

  endif ! if there is potentially stuff to write




  !--------------------------------
  ! write mostbound particle list
  !--------------------------------

  fileloc=TRIM('output_'//TRIM(output_to_string)//'/progenitor_data.dat')

#ifndef WITHOUTMPI
  ! Need to call MPI routines even if this CPU has nothing to write!
  call MPI_FILE_OPEN(MPI_COMM_WORLD, fileloc, MPI_MODE_WRONLY + MPI_MODE_CREATE, &
    MPI_INFO_NULL, filehandle, mpi_err)
  call MPI_FILE_WRITE_ORDERED(filehandle, particlelist, & 
    progenitorcount_written, MPI_INTEGER, state, mpi_err) 
  call MPI_FILE_CLOSE(filehandle, mpi_err)
#else
  open(unit=666,file=fileloc,form='unformatted')
  write(666) particlelist
  close(666)
  ! open(unit=666, file=fileloc, form='formatted')
  ! do ipeak = 1, ihalo
  !   do ipart=(ipeak-1)*(nmost_bound+1) + 1, (ipeak)*(nmost_bound+1)
  !     write(666, '(I7,x)', advance='no') particlelist(ipart)
  !   enddo
  !   write(666,*)
  ! enddo
  ! close(666)
#endif




  !--------------------------------
  ! write progenitor mass list
  !--------------------------------

  fileloc=TRIM('output_'//TRIM(output_to_string)//'/progenitor_mass.dat')


#ifndef WITHOUTMPI
  call MPI_FILE_OPEN(MPI_COMM_WORLD, fileloc, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, filehandle, mpi_err)
  call MPI_FILE_WRITE_ORDERED(filehandle, masslist, ihalo, MPI_DOUBLE_PRECISION, state, mpi_err) 
  call MPI_FILE_CLOSE(filehandle, mpi_err)
#else
  open(unit=666,file=fileloc,form='unformatted')
  write(666) masslist
  close(666)
  ! open(unit=666, file=fileloc, form='formatted')
  ! do ipeak = 1, ihalo
  !   write(666, '(E14.6,x)') masslist(ipeak)
  ! enddo
  ! close(666)
#endif

  deallocate(particlelist)
  deallocate(masslist)






  !========================================
  ! Write past merged progenitors data
  !========================================

  !---------------------------------------------------------
  ! Write past merged progenitor IDs, their galaxy particle
  ! and their merging time
  ! Format:
  ! ID prog1, galaxy prog1, time prog1, ID prog2, ...
  !---------------------------------------------------------

  allocate(pastproglist(1:3*(pmprog_free - 1)))
  pastproglist = 0
  allocate(pastprogmasslist(1:3*(pmprog_free - 1)))
  pastprogmasslist = 0

  npastprogs_all = 0 ! count how many pmprogs you write. will be communicated later.
  pind = 0

  ! If the past merged progenitor was used, the owner was overwritten to 0.
  ! Share that info before continuing.
#ifndef WHITOUTMPI
  if (npastprogs > 0) then
    call MPI_ALLREDUCE(MPI_IN_PLACE, pmprogs_owner, npastprogs, &
      MPI_INT, MPI_MIN, MPI_COMM_WORLD, ipeak)
  endif
#endif


  if (pmprog_free > 1) then ! if you have stuff to write
    do ipeak = 1, pmprog_free - 1

      ! Only write if you own this progenitor
      if (pmprogs_owner(ipeak)==myid) then 

        ! Only write if the progenitor hasn't merged too long ago
        if (max_past_snapshots > 0) then
          if ((ifout - pmprogs_t(ipeak)) > max_past_snapshots) then
            pastproglist(pind+1) = pmprogs(ipeak)
            pastproglist(pind+2) = pmprogs_galaxy(ipeak)
            pastproglist(pind+3) = pmprogs_t(ipeak)
            pind = pind + 3
            npastprogs_all = npastprogs_all + 1
            pastprogmasslist(npastprogs_all) = pmprogs_mass(ipeak)
          endif
        else
          pastproglist(pind+1) = pmprogs(ipeak)
          pastproglist(pind+2) = pmprogs_galaxy(ipeak)
          pastproglist(pind+3) = pmprogs_t(ipeak)
          pind = pind + 3
          npastprogs_all = npastprogs_all + 1
          pastprogmasslist(npastprogs_all) = pmprogs_mass(ipeak)
        endif
      endif
    enddo
  endif

  ! henceforth, pind is the number of elements in the array to be written





  !-------------------------------------
  ! Write past merged progenitors list
  !-------------------------------------

  fileloc=TRIM('output_'//TRIM(output_to_string)//'/past_merged_progenitors.dat')

#ifndef WITHOUTMPI
  call MPI_FILE_OPEN(MPI_COMM_WORLD, fileloc, MPI_MODE_WRONLY + MPI_MODE_CREATE, &
    MPI_INFO_NULL, filehandle, mpi_err)
  call MPI_FILE_WRITE_ORDERED(filehandle, pastproglist, & 
    pind, MPI_INTEGER, state, mpi_err) 
  call MPI_FILE_CLOSE(filehandle, mpi_err)
#else
  open(unit=666,file=fileloc,form='unformatted')
  write(666) pastproglist
  close(666)
  ! open(unit=666, file=fileloc, form='formatted')
  ! ipart = 1
  ! do while (ipart< pind)
  !   write(666, '(3(I7,x))') pastproglist(ipart), pastproglist(ipart+1), pastproglist(ipart+2)
  !   ipart = ipart + 3
  ! enddo
  ! close(666)
#endif  




  !-----------------------------------------
  ! Write past merged progenitors mass list
  !-----------------------------------------
  fileloc=TRIM('output_'//TRIM(output_to_string)//'/past_merged_progenitor_mass.dat')

#ifndef WITHOUTMPI
  call MPI_FILE_OPEN(MPI_COMM_WORLD, fileloc, MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, filehandle, mpi_err)
  call MPI_FILE_WRITE_ORDERED(filehandle, pastprogmasslist, npastprogs_all, MPI_DOUBLE_PRECISION, state, mpi_err) 
  call MPI_FILE_CLOSE(filehandle, mpi_err)
#else
  open(unit=666,file=fileloc,form='unformatted')
  write(666) pastprogmasslist
  close(666)
  ! open(unit=666, file=fileloc, form='formatted')
  ! do ipeak = 1, ihalo
  !   write(666, '(E14.6,x)') pastprogmasslist(ipeak)
  ! enddo
  ! close(666)
#endif

  deallocate(pastproglist, pastprogmasslist)





  !======================================
  ! Write number of progenitors to file
  ! (both current and past)
  !======================================

#ifndef WITHOUTMPI
  buf = (/progenitorcount, ihalo, progenitorcount_written, npastprogs_all/)
  if (myid == 1) then
    call MPI_REDUCE(MPI_IN_PLACE, buf, 4, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
    progenitorcount = buf(1)
    ihalo = buf(2)
    progenitorcount_written = buf(3)
    npastprogs_all = buf(4)
  else
    call MPI_REDUCE(buf, buf, 4, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, mpi_err)
  endif
#endif


  if (myid == 1) then 
    fileloc=TRIM('output_'//TRIM(output_to_string)//'/progenitorcount.dat')
    open(unit=666,file=fileloc,form='unformatted')
    write(666) progenitorcount, ihalo, progenitorcount_written, npastprogs_all
    ! open(unit=666,file=fileloc,form='formatted')
    ! write(666, '(2(I7,x))') progenitorcount, progenitorcount_written
    close(666)
  endif

end subroutine write_progenitor_data











!=================================================
subroutine get_local_prog_id(global_id, local_id)
!=================================================
  
  use clfind_commons
  implicit none
  integer, intent(in) :: global_id 
  integer, intent(out) :: local_id 
  integer :: i

  !-------------------------------
  ! Gets the local progenitor ID
  ! given the global progenitor ID.
  !-------------------------------



  ! Check if prog is already included
  do i = 1, prog_free
    if (prog_id(i) == global_id) then
      local_id = i
      return
    endif
  enddo

  !progenitor is not in there; Give him a new index
  i = prog_free !save for later
  prog_id(i) = global_id
  prog_free = prog_free + 1
  local_id = i
  return

end subroutine get_local_prog_id 





!==================================================
subroutine fill_matrix(mat, key, tar, np, act)
!==================================================

  !----------------------------------------------------
  ! Counts that np tracer particles for the
  ! clump key has been found in the clump target
  ! and adds it to the matrix mat.
  !
  ! clump key:    local id for progenitor,
  !               local id for descendant
  ! clump target: local id for progenitor,
  !               global id for descendant
  ! act:          how to add the particles to matrix
  !   act = 'add': add up
  !   act = 'set': don't add, just set
  !   act = 'inv': make value negative
  !----------------------------------------------------

  use clfind_commons

  implicit none

  integer, intent(in)                :: key, tar, np
  character(len=3), intent(in)       :: act
  type(prog_desc_mat), intent(inout) :: mat

  integer :: i, itar

  if (mat%cnt(key) == 0) then
    ! if there is no first value, 
    ! add new first value
    mat%first(key) = mat%mat_free_ind
    mat%clmp_id(mat%mat_free_ind) = tar
    mat%cnt(key) = mat%cnt(key) + 1

    if (act=='add') then
      mat%ntrace(mat%mat_free_ind) = mat%ntrace(mat%mat_free_ind) + np
    else if (act=='set') then
      mat%ntrace(mat%mat_free_ind) = np
    else if (act=='inv') then
      mat%ntrace(mat%mat_free_ind) = - mat%ntrace(mat%mat_free_ind)
    endif

    mat%mat_free_ind = mat%mat_free_ind + 1
    return

  else
    ! Try to find a match first
    itar = mat%first(key)
    do i = 1, mat%cnt(key)
      ! If you found a match:
      if (mat%clmp_id(itar) == tar) then
        if (act=='add') then
          mat%ntrace(itar) = mat%ntrace(itar) + np
        else if (act=='set') then
          mat%ntrace(itar) = np
        else if (act=='inv') then
          mat%ntrace(itar) = -mat%ntrace(itar)
        endif
        return
      endif

      if (mat%next(itar) > 0) then
        !cycle
        itar = mat%next(itar)
      else
        exit
      endif
    enddo
  endif




  ! if you didn't find anything, add new value
  mat%next(itar) = mat%mat_free_ind
  mat%clmp_id(mat%mat_free_ind) = tar
  mat%cnt(key) = mat%cnt(key) + 1
  if (act=='add') then
    mat%ntrace(mat%mat_free_ind) = mat%ntrace(mat%mat_free_ind) + np
  else if (act=='set') then
    mat%ntrace(mat%mat_free_ind) = np
  else if (act=='inv') then
    mat%ntrace(mat%mat_free_ind) = - mat%ntrace(mat%mat_free_ind)
  endif

  mat%mat_free_ind = mat%mat_free_ind + 1
  return

end subroutine fill_matrix







!======================================
subroutine deallocate_mergertree()
!======================================

  !---------------------------------
  ! Deallocate arrays for mergertree
  !---------------------------------

  use clfind_commons 
  implicit none
  

  if (allocated(prog_id)) then
    deallocate(prog_id)
    deallocate(prog_owner)
    deallocate(galaxy_tracers)
    deallocate(prog_mass)
  endif


  if (allocated(p2d_links%first)) then
    deallocate(p2d_links%first)
    deallocate(p2d_links%cnt)
    deallocate(p2d_links%ntrace)
    deallocate(p2d_links%clmp_id)
    deallocate(p2d_links%next)

    deallocate(d2p_links%first)
    deallocate(d2p_links%cnt)
    deallocate(d2p_links%ntrace)
    deallocate(d2p_links%clmp_id)
    deallocate(d2p_links%next)

    deallocate(main_desc)
  endif

  if(ifout > 1) deallocate(main_prog)

  if (allocated(pmprogs)) then
    deallocate(pmprogs)
    deallocate(pmprogs_owner)
    deallocate(pmprogs_galaxy)
    deallocate(pmprogs_t)
    deallocate(pmprogs_mass)
  endif


  if (allocated(tracers_all)) then
      deallocate(tracers_all)
      deallocate(tracer_loc_progids_all)
  endif
   

  deallocate(most_bound_energy)
  deallocate(most_bound_pid)
  deallocate(clmp_mass_exclusive)
  ! deallocate(clmp_vel_exclusive)

end subroutine deallocate_mergertree






!=====================================
subroutine mark_tracer_particles()
!=====================================

  !------------------------------------------------------------------
  ! goes through all clumps
  ! first determines most bound particles via mpi communications
  ! then marks all most bound particles by making their clumpid 
  ! negative
  !------------------------------------------------------------------
  use clfind_commons 
  use amr_commons 
  implicit none
  
  integer :: ipeak, ipart, i
  real(dp), dimension(1:npeaks_max) :: temp_energy
  temp_energy = HUGE(0.d0)



  ! for each of nmost_bound particles
  do ipart = 1, nmost_bound 

    ! copy this particle of each clump in an array, communicate and find min
    do ipeak = 1, hfree -1
      temp_energy(ipeak) = most_bound_energy(ipeak, ipart)
    enddo

    call build_peak_communicator()
    call virtual_peak_dp(temp_energy(:), 'min')
    call boundary_peak_dp(temp_energy(:))


    do ipeak = 1, hfree -1
      !need to check whether it is relevant clump. otherwise temp_energy == most_bound_energy = HUGE anyways
      if (clmp_mass_exclusive(ipeak) > 0) then 

        !if this is true minimal value for this peak
        if (temp_energy(ipeak) == most_bound_energy(ipeak, ipart) ) then

          ! mark particle if it is bound
          if (most_bound_energy(ipeak, ipart) < 0) then
            clmpidp(most_bound_pid(ipeak, ipart)) = - clmpidp(most_bound_pid(ipeak, ipart))
          endif

        else
          ! move particle energy list one to the right
          ! so you won't miss true minimum at later checks

          do i = nmost_bound, ipart + 1, -1
            most_bound_energy(ipeak, i) = most_bound_energy(ipeak, i-1)
            most_bound_pid(ipeak, i) = most_bound_pid(ipeak, i-1)
          enddo 

          ! also reset value at this position
          most_bound_energy(ipeak, ipart) = HUGE(0.d0)
          most_bound_pid(ipeak, ipart) = 0

        endif !true minimum
      
      endif !relevant clump

    enddo

  enddo

end subroutine mark_tracer_particles






!=======================================================
subroutine dissolve_small_clumps(ilevel, for_halos)
!=======================================================

  !----------------------------------------------------------------------
  ! Dissolve clumps with too small mass into parents/nothing.
  ! Clump is required to have at least mass_threshold number of
  ! its own particles.
  !----------------------------------------------------------------------

  use amr_commons 
  use clfind_commons 

  implicit none

  integer, intent(in)    :: ilevel
  logical, intent(in)    :: for_halos ! whether to do it for halos or for subhalos

  integer       :: killed, appended
  integer       :: ipeak, ipart, thispart, particle_local_id

#ifndef WITHOUTMPI
  include 'mpif.h'
  integer, dimension(1:2) :: buf
  integer :: info
#endif 



  !------------------------------------
  ! Kill or append too small clumps
  !------------------------------------

  call get_exclusive_clump_mass(ilevel) !subroutine further below

  killed = 0; appended = 0

  do ipeak=1, hfree-1
    if (lev_peak(ipeak) == ilevel) then

      ! if there are too few particles in there, but at least 1 (i.e. don't do it for noise)
      if (relevance(ipeak) > relevance_threshold .and. clmp_mass_exclusive(ipeak) < (mass_threshold * partm_common) ) then

        if (is_namegiver(ipeak) .and. for_halos) then

          !--------------------------------
          ! if clump is namegiver, kill it
          !--------------------------------

          if(ipeak <= npeaks) killed = killed + 1 !count only non-virtuals
          
          ! remove particles from clump
          thispart = clmppart_first(ipeak)
          do ipart = 1, nclmppart(ipeak)

            if (clmpidp(thispart) > 0) then
              call get_local_peak_id(clmpidp(thispart), particle_local_id)
              if (particle_local_id == ipeak) then
                clmpidp(thispart) = 0
              endif
            endif
            thispart = clmppart_next(thispart)

          enddo 

          nclmppart(ipeak) = 0
          clmp_mass_pb(ipeak) = 0
          clmp_mass_exclusive(ipeak) = 0

        elseif (.not.is_namegiver(ipeak) .and. .not.for_halos) then
          !---------------------------------------------------
          ! if clump isn't namegiver, add particles to parent
          !---------------------------------------------------

          if(ipeak <= npeaks) appended = appended + 1 ! count only non-virtuals

          ! remove particles from clump
          thispart = clmppart_first(ipeak)
          do ipart = 1, nclmppart(ipeak)

            if (clmpidp(thispart) > 0) then
              call get_local_peak_id(clmpidp(thispart), particle_local_id)
              if (particle_local_id == ipeak) then
                clmpidp(thispart) = new_peak(ipeak)
              endif
            endif
            thispart = clmppart_next(thispart)

          enddo 

          nclmppart(ipeak) = 0
          clmp_mass_pb(ipeak) = 0
          clmp_mass_exclusive(ipeak) = 0

        endif !namegiver or not
      endif !if too small
    endif ! correct peak level
  enddo !all peaks





  !---------------------------------------
  ! speak to me
  !---------------------------------------


  killed_tot = killed_tot + killed
  appended_tot = appended_tot + appended




  if (for_halos) then

#ifndef WITHOUTMPI
    buf = (/killed_tot, appended_tot/)
    if (myid == 1) then
      call MPI_REDUCE(MPI_IN_PLACE, buf, 2, MPI_INTEGER,MPI_SUM, 0, MPI_COMM_WORLD, info)
      killed_tot = buf(1)
      appended_tot = buf(2)
    else
      call MPI_REDUCE(buf, buf, 2, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, info)
    endif
#endif

    if(myid == 1) then
      write(*,'(A43,I6,A14,I6,A18)') " Handling too small clumps TOTAL: Dissolved ", killed_tot, " halos; Merged ", appended_tot, " to their parents." 
    endif

    !reset values for next output step
    killed_tot = 0
    appended_tot = 0

  endif




  contains 
    !==============================================
    subroutine get_exclusive_clump_mass(ilevel)
    !==============================================

      use clfind_commons
      use pm_commons, only: mp!, vp

      implicit none
      integer, intent(in) :: ilevel
      integer :: ipeak, ipart, thispart!, i

      !----------------------------------------------------
      ! recompute clump properties after unbinding
      !----------------------------------------------------

      do ipeak=1, hfree-1 !loop over all peaks
       
        !reset values for virtual peaks to communicate multiple times
        if (ipeak > npeaks) then
          clmp_mass_exclusive(ipeak) = 0
          ! clmp_vel_exclusive(ipeak,:) = 0
        endif


        if (lev_peak(ipeak) == ilevel ) then

          clmp_mass_exclusive(ipeak) = 0
          ! clmp_vel_exclusive(ipeak,:) = 0

          if (nclmppart(ipeak) > 0 ) then
            ! if there is work to do on this processing unit for this peak
            thispart=clmppart_first(ipeak)
            
            do ipart=1, nclmppart(ipeak)        ! while there is a particle linked list
              if (clmpidp(thispart) > 0) then
                call get_local_peak_id(clmpidp(thispart), particle_local_id) 
                if (particle_local_id == ipeak) then

                  clmp_mass_exclusive(ipeak)=clmp_mass_exclusive(ipeak)+mp(thispart)
                  ! do i=1,3
                  !   clmp_vel_exclusive(ipeak,i)=clmp_vel_exclusive(ipeak,i)+vp(thispart,i)*mp(thispart) !get velocity sum
                  ! enddo
                endif
              endif

              thispart=clmppart_next(thispart) ! go to next particle in linked list
            enddo   ! loop over particles
          endif     ! clump has particles on this processor 
        endif       ! there is work for this peak on this processor
      enddo         ! loop over peaks


      !----------------------------------------------------------
      ! communicate clump mass and velocity across processors
      !----------------------------------------------------------
      call build_peak_communicator
      call virtual_peak_dp(clmp_mass_exclusive,'sum')       !collect
      call boundary_peak_dp(clmp_mass_exclusive)            !scatter
      ! do i=1,3
      !   call virtual_peak_dp(clmp_vel_exclusive(1,i),'sum')  !collect
      !   call boundary_peak_dp(clmp_vel_exclusive(1,i))       !scatter
      ! enddo

    end subroutine get_exclusive_clump_mass 

end subroutine dissolve_small_clumps 

! #endif NDIM == 3
#endif


