!-----------------------------------------------------------------
! This file contains the routines for particle unbinding.
! See wiki for more information.
!
! There are three optional preprocessing definitions for particle
! unbinding and merger trees only:
! -DUNBINDINGCOM
!   use (and iteratively determine) the center of mass as the
!   center of clumps
! -DMTREEDEBUG
!   create a lot of formatted output to help debugging the merger
!   tree routines. Don't use this unless you're fighting bugs, it
!   will create a looooot of otherwise unnecessary output.
! -DMTREE_INDIVIDUAL_FILES
!   instead of collective writes into a single file for progenitor
!   data (not mergertree/galaxy result files), every task writes an
!   individual file. This was added because some MPI implementations
!   had issues with collective writing.
!
! Contains:
!   subroutine unbinding()
!   subroutine get_clumpparticles()
!   subroutine get_clump_properties_pb()
!   subroutine get_cmp_iter()
!   subroutine get_cmp_noiter()
!   subroutine get_closest_border()
!   subroutine unbinding_neighborsearch()
!   subroutine bordercheck()
!   subroutine particle_unbinding()
!   subroutine eparttot()
!   subroutine potential()
!   subroutine compute_phi()
!   subroutine dissolve_small_clumps()
!     contains subroutine get_exclusive_clump_mass()
!   subroutine read_unbinding_params()
!   subroutine allocate_unbinding_arrays()
!   subroutine deallocate_unbinding_arrays()
!-----------------------------------------------------------------

#if NDIM==3
subroutine unbinding()

  use amr_commons    ! MPI stuff
  use pm_commons, only: mp, npart, npartmax, levelp
  use hydro_commons, ONLY:mass_sph
  use clfind_commons ! unbinding stuff
  use mpi_mod

!
  implicit none
#ifndef WITHOUTMPI
  integer :: info
  real(dp):: partm_common_all
  logical :: loop_again_global
#endif

  !------------------------------------------------------------
  ! This subroutine assigns all particles that are in cells
  ! identified to be in a clump by the clump finder the peak
  ! ID of that clump and checks whether they are energetically
  ! bound to that structure. If not, they are passed on to the
  ! clump's parents.
  !------------------------------------------------------------

  integer                     :: ipeak, ilevel, ipart, i, parent_local_id
  integer                     :: loop_counter=0
  integer, dimension(1:npart) :: clump_ids
  character(LEN=80)           :: fileloc, filedir
  character(LEN=5)            :: nchar,nchar2
  logical                     :: is_final_round, check




  !========================
  ! Initial set-up
  !========================

  ! Logging/Announcing stuff
  if(myid==1) write(*,*) "Started unbinding."

  ! update boundary relevance
  call build_peak_communicator
  call boundary_peak_dp(relevance)
  ! call boundary_peak_int(new_peak) ! already done in clump_finder, doesn't
                                     ! need an update
  call boundary_peak_int(lev_peak)
  do i = 1, 3
    call boundary_peak_dp(peak_pos(1,i))
  enddo

#ifdef MTREEDEBUG
  ! if there are no clumps yet, the output directories haven't been made yet.
  ! if MTREEDEBUG, you'll need the directory for file dumps.
  call title(ifout, nchar)
  filedir = 'output_'//TRIM(nchar)
  call create_output_dirs(filedir)
#endif


  ! set up constants and counters
  GravConst=1d0                                       ! Gravitational constant
  if(cosmo) GravConst=3d0/8d0/3.1415926*omega_m*aexp

  periodical=(nx==1)                                  ! true if periodic
  nunbound=0                                          ! count unbound particles
  candidates=0                                        ! count unbinding candidates: Particles of child clumps
                                                      ! that will be tested. If a particle is passed on to the
                                                      ! parent clump, it will be counted twice.

  progenitorcount = 0                                 ! count clumps that will be progenitors in the next snapshot
  progenitorcount_written = 0                         ! count halos that you'll write to file for mergertree

  killed_tot = 0; appended_tot = 0;                   ! count how many too small clumps have been dissolved

  mergelevel_max = max(0, mergelevel_max)             ! make sure you go at least once; if no subhalos exist,
                                                      ! the initial mergelevel_max will be -1




  ! get particle mass (copied from subroutine write_clump_properties)
  if(ivar_clump==0)then
    partm_common=MINVAL(mp, MASK=(mp.GT.0.))
#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(partm_common,partm_common_all,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,info)
    partm_common=partm_common_all
#endif
  else
    if(hydro)then
      partm_common=mass_sph
    endif
  endif



  ! allocate necessary arrays
  call allocate_unbinding_arrays()


  if (npeaks_tot > 0) then

    ! initialise constant part of is_namegiver array
    do ipeak=1, npeaks
      if (new_peak(ipeak)>0) then
        call get_local_peak_id(new_peak(ipeak), parent_local_id)
        if (ipeak == parent_local_id) is_namegiver(ipeak) = .true.
      endif
    enddo


    !===================
    ! Gather particles
    !===================

    ! Get particles in substructure, create linked lists
    call get_clumpparticles()


#ifndef UNBINDINGCOM
    !-------------------------------
    ! get cumulative mass profiles
    !-------------------------------
    call get_cmp_noiter()
#endif
#ifdef MTREEDEBUG
    call mtreedebug_dump_unbinding_data('init')
#endif



    !==================
    ! Unbinding loop
    !==================

    ! go level by level
    do ilevel=0, mergelevel_max

      !---------------------------------
      ! Preparation for iteration loop
      !---------------------------------

      is_final_round=.true.
      if (iter_properties) is_final_round=.false.

      loop_again=.true.
      loop_counter=0

      ! reset values
      ! WARNING: Here to_iter doesn't include whether it is a namegiver!!
      to_iter = (lev_peak==ilevel)
      hasatleastoneptcl=1 ! set array value to 1

      do ipeak=npeaks+1, hfree-1
        if (new_peak(ipeak) > 0) then
          call get_local_peak_id(new_peak(ipeak), parent_local_id)
          if (ipeak == parent_local_id) is_namegiver(ipeak) = .true.
        endif
      enddo


      !-------------------------
      ! iteration per level
      !-------------------------

      do while(loop_again)
        loop_again = .false.  ! set boolean whether to continue to false as default;
                              ! will get reset if conditions are met

        loop_counter=loop_counter+1
        niterunbound=0

        !--------------------------------------------------------
        ! get particle based clump properties :
        ! bulk velocity, particle furthest away from clump center
        ! This subroutine determines whether to loop again
        !--------------------------------------------------------
        call get_clump_properties_pb(loop_counter==1, ilevel)

        ! forcibly stop loop if necessary
        if (loop_counter==repeat_max) loop_again=.false.

#ifndef WITHOUTMPI
        ! sync with other processors whether you need to repeat loop
        call MPI_ALLREDUCE(loop_again, loop_again_global, 1, MPI_LOGICAL, MPI_LOR,MPI_COMM_WORLD, info)
        loop_again=loop_again_global
#endif

        if (.not.loop_again) is_final_round=.true.

#ifdef UNBINDINGCOM
        !-------------------------------
        ! get cumulative mass profiles
        !-------------------------------
        call get_cmp_iter()
#endif


        !-----------------------------------------
        ! get closest border to the peak position
        !-----------------------------------------
        if (saddle_pot) call get_closest_border()



        !---------------
        ! Unbinding
        !---------------

        do ipeak=1, hfree-1
          ! don't apply to_iter here!
          ! needs to be done in final round even if to_iter = .false.
          check = clmp_mass_pb(ipeak)>0.0 .and. lev_peak(ipeak) == ilevel
          if (check) then
            call particle_unbinding(ipeak, is_final_round)
          endif
        enddo



        !------------------------
        ! prepare for next round
        !------------------------

        if (loop_again) then
          ! communicate whether peaks have remaining contributing particles
          call build_peak_communicator()
          call virtual_peak_int(hasatleastoneptcl,'max')
          call boundary_peak_int(hasatleastoneptcl)

          do ipeak=1,hfree-1
            ! if peak has no contributing particles anymore
            if (hasatleastoneptcl(ipeak)==0) then
              to_iter(ipeak)=.false. ! don't loop anymore over this peak
            endif

            if (loop_counter == 1) then ! only do this for first round
              if(is_namegiver(ipeak)) then
                ! Don't iterate over halo-namegivers.
                ! Only set here to false and not earlier, so they'll be considered
                ! every first time the loop over levels starts.
                to_iter(ipeak) = .false. ! don't iterate over halo namegivers
              endif
            endif

          enddo

        else ! if it is final round

          if(make_mergertree) call dissolve_small_clumps(ilevel, .false., .false.)

        endif


        !---------------------
        ! Talk to me.
        !---------------------
        if (clinfo) then
#ifndef WITHOUTMPI
            call MPI_ALLREDUCE(niterunbound, niterunbound_tot, 1, MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD, info)
#else
            niterunbound_tot=niterunbound
#endif
          if (iter_properties .and. myid==1) then
            if (loop_again) then
              write(*,'(A10,I10,A30,I5,A7,I5)') " Unbound", niterunbound_tot, &
                "particles at level", ilevel, "loop", loop_counter

            else ! do not loop again

              if (loop_counter < repeat_max) then
                write(*,'(A10,I10,A30,I5,A7,I5)') " Unbound", niterunbound_tot, &
                  "particles at level", ilevel, "loop", loop_counter
                write(*, '(A7,I5,A35,I5,A12)') "Level ", ilevel, &
                  "clump properties converged after ", loop_counter, "iterations."
              else
                 write(*,'(A15,I5,A20,I5,A35)') "WARNING: Level ", ilevel, &
                   "not converged after ", repeat_max, "iterations. Moving on to next step."
              endif

            endif   ! loop again
          endif     ! iter properties

          if (.not. iter_properties .and. myid==1) then
            write(*,'(A10,I10,A30,I5)') " Unbound", niterunbound_tot, &
              "particles at level", ilevel
          endif

        endif       ! clinfo



      enddo ! loop again for ilevel

    enddo ! loop over levels





    !================================
    ! Talk to me when unbinding done
    !================================


    if (clinfo) then
#ifndef WITHOUTMPI
      call MPI_ALLREDUCE(nunbound, nunbound_tot, 1, MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD, info)
      call MPI_ALLREDUCE(candidates, candidates_tot, 1, MPI_INTEGER, MPI_SUM,MPI_COMM_WORLD, info)
#else
      nunbound_tot=nunbound
      candidates_tot=candidates
#endif
      if (myid==1) then
        write(*,'(A6,I10,A30,I10,A12)') " Found", nunbound_tot, "unbound particles out of ", candidates_tot, " candidates"
      endif
    endif


  ! After the loop: Dissolve too small halos (namegivers)
  call dissolve_small_clumps(0, .true., .false.)


  else ! if npeaks_tot == 0
    if(myid==1) write(*,*) "I have no clumps to work with. Skipping unbinding."
  endif






  !=================
  ! Write output
  !=================

#ifndef MTREEDEBUG
  ! if there are no clumps yet, the output directories haven't been made yet.
  ! if not MTREEDEBUG, the directories haven't been necessary yet.
  call title(ifout, nchar)
  filedir = 'output_'//TRIM(nchar)
  call create_output_dirs(filedir)
#endif

  call title(myid, nchar2)
  fileloc=TRIM(filedir)//'/unbinding_'//TRIM(nchar)//'.out'//TRIM(nchar2)

  open(unit=666,file=fileloc,form='unformatted')

  ipart=0
  do i=1,npartmax
    if(levelp(i)>0)then
      ipart=ipart+1
      clump_ids(ipart)=clmpidp(i)
    endif
  enddo
  write(666) clump_ids
  close(666)

  ! create_output = .true., otherwise unbdinging() wouldn't have been called
  if(npeaks_tot>0)then
    if(myid==1)write(*,*)"Outputing clump properties to disc."
    call write_clump_properties(.true.)
  endif






  !=========================================
  ! After unbinding: Do mergertree stuff
  !=========================================

  if (make_mergertree) then
    ! sum up masses if necessary
    if (.not. use_exclusive_mass .or. make_mock_galaxies) then
      ! recompute clmp_mass_pb if necessary
      ! clmp_mass_pb will only be needed in any of the two cases in the if-condition

      ! reset mass: clmp_mass_pb is not updated after unbinding;
      ! which particles are bound in the end is not included in there.
      ! There might've been some changes. Fix this now.
      clmp_mass_pb(1:npeaks) = clmp_mass_exclusive(1:npeaks)

      do ilevel = 0, mergelevel_max

        ! First reset virtual's mass for comms
        clmp_mass_pb(npeaks+1:hfree-1) = 0d0

        do ipeak = 1, npeaks
          ! only do this part for non-virtuals as "source", otherwise you'll get wrong additions!
          if (lev_peak(ipeak)== ilevel) then
            if (clmp_mass_exclusive(ipeak) > 0 .and. .not. is_namegiver(ipeak)) then
              call get_local_peak_id(new_peak(ipeak), parent_local_id)
              if (.not. is_namegiver(parent_local_id)) then ! namegivers already take all particles in get_clumpproperties()
                clmp_mass_pb(parent_local_id) = clmp_mass_pb(parent_local_id) + clmp_mass_pb(ipeak)
              endif
            endif
          endif
        enddo

        ! communicate: gather only, no need to scatter to virtuals, they will be reset immediately
        call virtual_peak_dp(clmp_mass_pb(:), 'sum')

      enddo
    endif

    ! now scatter to virtuals
    call boundary_peak_dp(clmp_mass_pb(:))

#ifdef MTREEDEBUG
    call mtreedebug_dump_unbinding_data('after')
#endif

    ! Now call mergertree
    call make_merger_tree()

  endif



  !====================
  ! Deallocate arrays
  !====================
  call deallocate_unbinding_arrays()



  !====================
  ! Say good bye.
  !====================
  if(verbose.or.myid==1) write(*,*) "Finished unbinding."


  return

end subroutine unbinding
!######################################
!######################################
!######################################
subroutine get_clumpparticles()

  !---------------------------------------------------------------------------
  ! This subroutine loops over all test cells and assigns all particles in a
  ! testcell the peak ID the testcell has. If the peak is not a namegiver
  ! (= not its own parent), the particles are added to a linked list for
  ! unbinding later.
  !---------------------------------------------------------------------------

  use amr_commons
  use clfind_commons        ! unbinding stuff is all in here
  use pm_commons, only: numbp, headp, nextp, xp
  use amr_parameters
  use hydro_commons         ! using mass_sph
  implicit none

  ! for looping over test cells and getting particle list
  integer   :: itestcell, ipart,this_part, global_peak_id, local_peak_id, prtcls_in_grid

  ! getting particles per peak
  integer   :: ind, grid


  ! getting in which cell of a grid a particle is
  integer   :: part_cell_ind,i,j,k

  ! appending linked lists
  integer   :: ipeak, new_peak_local_id, ilevel


  if(verbose) write(*,*) "Entered get_clumpparticles"


  !-----------------------------------------------------------
  ! Get particles from testcells into linked lists for clumps
  !-----------------------------------------------------------

  do itestcell=1, ntest ! loop over all test cells
    global_peak_id=flag2(icellp(itestcell))

    if (global_peak_id /= 0) then

      ! get local peak id
      call get_local_peak_id(global_peak_id, local_peak_id)

      if (relevance(local_peak_id) > relevance_threshold) then
        ! create linked particle list

        ind=(icellp(itestcell)-ncoarse-1)/ngridmax+1    ! get cell position
        grid=icellp(itestcell)-ncoarse-(ind-1)*ngridmax ! get grid index
        prtcls_in_grid = numbp(grid)                    ! get number of particles in grid
        this_part=headp(grid)                           ! get index of first particle


        ! loop over particles in grid
        do ipart=1, prtcls_in_grid
          !check cell index of particle so you loop only once over each
          i=0
          j=0
          k=0
          if(xg(grid,1)-xp(this_part,1)/boxlen+(nx-1)/2.0 .le. 0) i=1
          if(xg(grid,2)-xp(this_part,2)/boxlen+(ny-1)/2.0 .le. 0) j=1
          if(xg(grid,3)-xp(this_part,3)/boxlen+(nz-1)/2.0 .le. 0) k=1

          part_cell_ind=i+2*j+4*k+1

          ! If index is correct, assign clump id to particle
          if (part_cell_ind==ind) then
            ! assign peak ID
            clmpidp(this_part)=global_peak_id
            ! add particle to linked list of clumpparticles
            ! check if already particles are assigned
            if (nclmppart(local_peak_id)>0) then
              ! append to the last particle of the list
              clmppart_next(clmppart_last(local_peak_id))=this_part
            else
              ! assign particle as first particle
              ! for this peak of linked list
              clmppart_first(local_peak_id)=this_part
            endif
            ! update last particle for this clump
            nclmppart(local_peak_id)=nclmppart(local_peak_id)+1
            clmppart_last(local_peak_id)=this_part
          endif
          ! go to next particle in this grid
          this_part=nextp(this_part)
        enddo
      endif     ! if clump is relevant
    endif       ! global peak /=0
  enddo         ! loop over test cells




  !------------------------------------------------------
  ! Append substructure particles to parents' linked list
  !------------------------------------------------------

  ! first do it for subhalos only to get full lists
  ! must be done level by level!
  do ilevel=0,mergelevel_max
    do  ipeak=1, hfree-1
      ! append substructure linked lists to parent linked lists
      if(lev_peak(ipeak)==ilevel) then
        if (nclmppart(ipeak)>0) then
          ! get local id of parent
          call get_local_peak_id(new_peak(ipeak),new_peak_local_id)
          ! if peak is namegiver, don't append to yourself
          if(ipeak/=new_peak_local_id) then
            ! It might happen that the parent peak doesn't have a
            ! particle linked list yet (on this processor).
            if (nclmppart(new_peak_local_id)>0) then ! particle ll exists
              clmppart_next(clmppart_last(new_peak_local_id))=clmppart_first(ipeak)
            else
              clmppart_first(new_peak_local_id)=clmppart_first(ipeak)
            endif

            clmppart_last(new_peak_local_id)=clmppart_last(ipeak)
            nclmppart(new_peak_local_id)=nclmppart(new_peak_local_id)+nclmppart(ipeak)
          endif
        endif
      endif
    enddo
  enddo

end subroutine get_clumpparticles
!########################################
!########################################
!########################################
subroutine get_clump_properties_pb(first, ilevel)
  use amr_commons
#ifdef UNBINDINGCOM
  use pm_commons, only: mp, vp, xp
#else
  use pm_commons, only: mp, vp
#endif
  use clfind_commons
  implicit none

  logical, intent(in) :: first  ! if it is the first time calculating
  integer, intent(in) :: ilevel ! current clump level

  !--------------------------------------------------------------------------
  ! This subroutine computes the particle-based properties of the clumps:
  ! namely the center of mass and the clump's velocity.
  ! If it's called for the first time, it will compute the properties for
  ! all peak IDs. If not, it will go level by level.
  !--------------------------------------------------------------------------

  ! iterators
  integer :: ipeak, i, ipart
  integer :: thispart

  real(dp) :: vsq
  real(dp),dimension(1:3) :: period
  real(dp),dimension(1:npeaks_max) :: clmp_vel_sq_pb_old
  logical :: check

#ifdef UNBINDINGCOM
  real(dp) :: biggest, distance
#endif

  if (verbose) write(*,*) "Entered get_clump_properties (particle based)"


  !------------------------------------------------------------
  ! If iterative: Store old values, reset virtual peak values
  !------------------------------------------------------------


  period = 0d0
  if (.not. first) clmp_vel_sq_pb_old = 0

  do ipeak=1, hfree-1
    check = (.not. first) .and. to_iter(ipeak)
    if (check) then
      clmp_vel_sq_pb_old(ipeak)=clmp_vel_pb(ipeak,1)**2+clmp_vel_pb(ipeak,2)**2+clmp_vel_pb(ipeak,3)**2
    endif

    if (iter_properties .and. ipeak>npeaks) then
      ! for communication: set virtual peak values=0
      ! so they won't contribute in the communication sum
      ! reset values
      do i=1,3
        clmp_vel_pb(ipeak,i)=0d0
#ifdef UNBINDINGCOM
        clmp_com_pb(ipeak,i)=0d0
#endif
      enddo
      clmp_mass_pb(ipeak)=0d0
    endif
  enddo


  !------------------------------------------------------
  ! GET CENTER OF MASS, CENTER OF MOMENTUM FRAME VELOCITY
  !------------------------------------------------------

  do ipeak=1, hfree-1 ! loop over all peaks

    if (to_iter(ipeak)) then ! if peak has particles and needs to be iterated over

      ! reset values
      do i=1,3
        clmp_vel_pb(ipeak,i)=0d0
#ifdef UNBINDINGCOM
        clmp_com_pb(ipeak,i)=0d0
#endif
      enddo
#ifdef UNBINDINGCOM
        cmp_distances(ipeak,nmassbins)=0d0
#endif
      clmp_mass_pb(ipeak)=0d0


      if (hasatleastoneptcl(ipeak)>0 .and. nclmppart(ipeak)>0) then
        ! if there is work to do on this processing unit for this peak

        thispart=clmppart_first(ipeak)

        do ipart=1, nclmppart(ipeak)       ! while there is a particle linked list
          if (contributes(thispart).or.is_namegiver(ipeak)) then  ! if the particle should be considered
#ifdef UNBINDINGCOM
            if (periodical) then           ! determine periodic correction
              period=0d0
              do i=1, 3
                if (xp(thispart,i)-peak_pos(ipeak,i) > 0.5*boxlen)    period(i)=(-1.0)*boxlen
                if (xp(thispart,i)-peak_pos(ipeak,i) < (-0.5*boxlen)) period(i)=boxlen
              enddo
            endif

            do i=1,3
              clmp_com_pb(ipeak,i)=clmp_com_pb(ipeak,i)+(xp(thispart,i)+period(i))*mp(thispart) ! get center of mass sum
            enddo
#endif

            clmp_mass_pb(ipeak)=clmp_mass_pb(ipeak)+mp(thispart)
            do i=1,3
              clmp_vel_pb(ipeak,i)=clmp_vel_pb(ipeak,i)+vp(thispart,i)*mp(thispart) ! get velocity sum
            enddo
          else
            contributes(thispart)=.true.   ! reset value
          endif
          thispart=clmppart_next(thispart) ! go to next particle in linked list
        enddo   ! loop over particles
      endif     ! there is work for this peak on this processor
    endif       ! peak needs to be looked at
  enddo         ! loop over peaks


  !----------------------------------------------------------------------
  ! communicate clump mass and velocity across processors
  !----------------------------------------------------------------------
  call build_peak_communicator
  call virtual_peak_dp(clmp_mass_pb,'sum')        ! collect
  call boundary_peak_dp(clmp_mass_pb)             ! scatter
  do i=1,3
#ifdef UNBINDINGCOM
    call virtual_peak_dp(clmp_com_pb(1,i),'sum')  ! collect
    call boundary_peak_dp(clmp_com_pb(1,i))       ! scatter
#endif
    call virtual_peak_dp(clmp_vel_pb(1,i),'sum')  ! collect
    call boundary_peak_dp(clmp_vel_pb(1,i))       ! scatter
  enddo


  !----------------------------------------------------------------------
  ! Cleanup: Dissolve clumps that are too small to start with
  !----------------------------------------------------------------------
  if (first) then
    call dissolve_small_clumps(ilevel, .false., .true.)
  endif


  !----------------------------------------------------------------------
  ! Now compute the actual clump properties of the remaining clumps
  !----------------------------------------------------------------------

  do ipeak=1, hfree-1
    check = to_iter(ipeak)
    check = check .and. clmp_mass_pb(ipeak)>0

    if (check) then
      ! calculate actual center of momentum frame velocity and CoM if necessary
      do i=1,3
#ifdef UNBINDINGCOM
        clmp_com_pb(ipeak,i)=clmp_com_pb(ipeak,i)/clmp_mass_pb(ipeak)
#endif
        clmp_vel_pb(ipeak,i)=clmp_vel_pb(ipeak,i)/clmp_mass_pb(ipeak)
      enddo

#ifdef UNBINDINGCOM
      !----------------------------------------
      ! FIND PARTICLE FURTHEST AWAY FROM CoM
      !----------------------------------------
      ! The maximal distance of a particle to the CoM is saved in the last
      ! cmp_distances array for every peak.
      ! NOTE: ifndef UNBINDINGCOM, this part is done only once in the get_cmp_noiter routine
      if(nclmppart(ipeak)>0) then
        biggest=0.0
        thispart=clmppart_first(ipeak)
        do ipart=1, nclmppart(ipeak) ! while there is a particle linked list
            period=0d0
            if (periodical) then
              do i=1, 3
                if (xp(thispart,i)-peak_pos(ipeak,i)>0.5*boxlen) period(i)=(-1.0)*boxlen
                if (xp(thispart,i)-peak_pos(ipeak,i)<(-0.5*boxlen)) period(i)=boxlen
              enddo
            endif

            distance=(xp(thispart,1)+period(1)-clmp_com_pb(ipeak,1))**2 + &
              (xp(thispart,2)+period(2)-clmp_com_pb(ipeak,2))**2 + &
              (xp(thispart,3)+period(3)-clmp_com_pb(ipeak,3))**2

            if(distance>biggest) biggest=distance ! save if it is biggest so far

          thispart=clmppart_next(thispart)

        enddo
        if (biggest>0.0) cmp_distances(ipeak,nmassbins)=sqrt(biggest) ! write if you have a result
      endif ! to iterate
#endif

    endif
  enddo     ! over all peaks

#ifdef UNBINDINGCOM

  !-------------------------------------------------
  ! communicate distance of particle furthest away
  !-------------------------------------------------
  call build_peak_communicator
  call virtual_peak_dp(cmp_distances(1,nmassbins), 'max')
  call boundary_peak_dp(cmp_distances(1,nmassbins))
#endif



  !-------------------------------------------------
  ! If iterative clump properties determination:
  ! Check whether bulk velocity converged
  !-------------------------------------------------

  if (iter_properties) then ! if clump properties will be determined iteratively
    if (first) then
      loop_again = .true.
    else
      do ipeak=1, hfree-1

        check = to_iter(ipeak)
        check = check .and.(.not.is_namegiver(ipeak))
        check = check .and. clmp_mass_pb(ipeak)>0d0

        if (check) then
          vsq=clmp_vel_pb(ipeak,1)**2+clmp_vel_pb(ipeak,2)**2+clmp_vel_pb(ipeak,3)**2

          if ( abs( sqrt(clmp_vel_sq_pb_old(ipeak)/vsq) - 1.0) < conv_limit ) then
            to_iter(ipeak) = .false. ! consider bulk velocity as converged
            ! write(*,'(A8,I3,A15,I8,A6,E15.6E2,A5,E15.6E2,A7,E15.6E2,A9,E15.6E2)') &
            ! & "#####ID", myid, "clump CONVERGED", ipeak+ipeak_start(myid), "old:", &
            ! & clmp_vel_sq_pb_old(ipeak), "new:", vsq, "ratio",  abs( sqrt(clmp_vel_sq_pb_old(ipeak)/vsq) - 1.0),&
            ! & "v_bulk=", sqrt(vsq)
          else
            loop_again=.true. ! repeat
          endif
        endif
      enddo
    endif
  endif


end subroutine get_clump_properties_pb
!###################################
!###################################
!###################################
#ifdef UNBINDINGCOM
subroutine get_cmp_iter()

  use amr_commons
  use pm_commons
  use clfind_commons
  use mpi_mod
  implicit none
  !-------------------------------------------------------------
  ! Get cumulative mass profiles with the center of the clump
  ! set as the center of mass of the clump. The center of mass
  ! is determined iteratively.
  !-------------------------------------------------------------

  integer  :: ipeak, i, ipart, levelmax
  real(dp) :: r_null, distance
  integer  :: thispart
  real(dp),dimension(1:3) :: period
  logical  :: check

#ifndef WITHOUTMPI
  integer  :: levelmax_glob, info
#endif

  if(verbose) write(*,*) "Entered get cumulative mass profiles"

  if (logbins) then
    ! get minimal distance:
    levelmax=0
    do i=1,nlevelmax
       if(numbtot(1,i)>0) levelmax=levelmax+1
    enddo

#ifndef WITHOUTMPI
    ! get system-wide levelmax
    call MPI_ALLREDUCE(levelmax,levelmax_glob, 1, MPI_INTEGER, MPI_MAX,MPI_COMM_WORLD, info)

    levelmax=levelmax_glob
#endif

    rmin=boxlen/2**levelmax
  endif




  do ipeak=1, hfree-1

    ! peak must have need to be reiterated
    check=to_iter(ipeak)                         ! =.true. for namegivers for the first computation
    check=check.and.nclmppart(ipeak)>0           ! peak must have particles on this processor

    ! reset values
    if (check .or. ipeak > npeaks) then
      do i = 1, nmassbins
        cmp(ipeak,i) = 0.0
      enddo
    endif

    if (check) then
      !------------------------------------------
      ! Compute cumulative mass binning distances
      !------------------------------------------
      ! The distances are not communicated later, but computed on each
      ! processor independently, because each processor has all information it needs
      ! with cmp_distances(ipeak,nmassbins) and CoM

      if (logbins) then
        do i=1, nmassbins-1
          cmp_distances(ipeak,i)=rmin*(cmp_distances(ipeak,nmassbins)/rmin)**(real(i)/real(nmassbins))
        enddo
      else ! linear binnings
        r_null=cmp_distances(ipeak,nmassbins)/real(nmassbins)
        do i=0, nmassbins-1
          cmp_distances(ipeak,i)=r_null*i
        enddo
      endif
      ! The last bin must end with precicely with the maximal
      ! Distance of the particle. That is
      ! needed because precision errors. The maximal distance is
      ! computed via particle data and the outermost bin is equal
      ! to the distance of the outermost particle to the CoM.
      ! Precision errors cause the code to crash here.

      !---------------------------------------------
      ! bin particles in cumulative mass profiles:
      ! get mass of each bin
      ! calculate particle distance to CoM
      !---------------------------------------------
      thispart=clmppart_first(ipeak)
      do ipart=1, nclmppart(ipeak)! while there is a particle linked list
        if (contributes(thispart)) then
          period=0d0
          if (periodical) then
            do i=1, 3
              if (xp(thispart,i)-peak_pos(ipeak,i)>0.5*boxlen)  period(i)=(-1.0)*boxlen
              if (xp(thispart,i)-peak_pos(ipeak,i)<(-0.5*boxlen)) period(i)=boxlen
            enddo
          endif
          distance=(xp(thispart,1)+period(1)-clmp_com_pb(ipeak,1))**2 + &
            (xp(thispart,2)+period(2)-clmp_com_pb(ipeak,2))**2 + &
            (xp(thispart,3)+period(3)-clmp_com_pb(ipeak,3))**2
          ! distance=(xp(thispart,1)+period(1)-peak_pos(ipeak,1))**2 + &
          !   (xp(thispart,2)+period(2)-peak_pos(ipeak,2))**2 + &
          !   (xp(thispart,3)+period(3)-peak_pos(ipeak,3))**2
          distance=sqrt(distance)


          i=1
          do
            if (distance<=cmp_distances(ipeak,i)) then
              cmp(ipeak,i) = cmp(ipeak,i) + mp(thispart)
              exit
            else
              i=i+1
            endif
          enddo
        endif

        thispart=clmppart_next(thispart)
      enddo

      ! sum up masses to get profile instead of mass in shell
      do i=0,nmassbins-1
        cmp(ipeak,i+1)=cmp(ipeak,i+1)+cmp(ipeak,i)
      enddo

    endif  ! check
  enddo    ! loop over peaks

  !--------------------------------------
  !communicate cummulative mass profiles
  !--------------------------------------
  call build_peak_communicator()
  do i=1,nmassbins
    call virtual_peak_dp(cmp(1,i), 'sum')
    call boundary_peak_dp(cmp(1,i))
  enddo

end subroutine get_cmp_iter
!########################################
!########################################
!########################################
#else
! =ifndef UNBINDINGCOM
subroutine get_cmp_noiter()

  use amr_commons
  use pm_commons
  use clfind_commons
  use mpi_mod
  implicit none

  !----------------------------------------------------------
  ! Get cumulative mass profiles
  ! Center of the clump is set at the density peak, so
  ! no iteration is necessary.
  !----------------------------------------------------------

  integer  :: ipeak, i, ipart, levelmax
  real(dp) :: r_null, distance, biggest
  integer  :: thispart
  real(dp),dimension(1:3) :: period

#ifndef WITHOUTMPI
  integer  :: levelmax_glob, info
#endif

  !--------------------------
  ! Prepare stuff
  !--------------------------


  if(verbose) write(*,*) "Entered get cumulative mass profiles"

  if (logbins) then
    ! get minimal distance:
    levelmax=0
    do i=1,nlevelmax
       if(numbtot(1,i)>0) levelmax=levelmax+1
    enddo

#ifndef WITHOUTMPI
    ! get system-wide levelmax
    call MPI_ALLREDUCE(levelmax,levelmax_glob, 1, MPI_INTEGER, MPI_MAX,MPI_COMM_WORLD, info)

    levelmax=levelmax_glob
#endif

    rmin=boxlen/2**levelmax
  endif


  !--------------------------------------
  ! FIND PARTICLE FURTHEST AWAY FROM CoM
  !--------------------------------------
  ! The maximal distance of a particle to the CoM is saved in the last
  ! cmp_distances array for every peak.
  do ipeak=1, hfree-1

    if (nclmppart(ipeak) > 0) then
      biggest=0.0
      thispart=clmppart_first(ipeak)
      do ipart=1, nclmppart(ipeak) ! while there is a particle linked list
        period=0d0
        if (periodical) then
          do i=1, 3
            if (xp(thispart,i)-peak_pos(ipeak,i)>0.5*boxlen) period(i)=(-1.0)*boxlen
            if (xp(thispart,i)-peak_pos(ipeak,i)<(-0.5*boxlen)) period(i)=boxlen
          enddo
        endif

        distance=(xp(thispart,1)+period(1)-peak_pos(ipeak,1))**2 + &
          (xp(thispart,2)+period(2)-peak_pos(ipeak,2))**2 + &
          (xp(thispart,3)+period(3)-peak_pos(ipeak,3))**2

        if(distance>biggest) biggest=distance ! save if it is biggest so far

        thispart=clmppart_next(thispart)
      enddo

      if (biggest>0.0) cmp_distances(ipeak,nmassbins)=sqrt(biggest) ! write if you have a result
    endif
  enddo     ! over all peaks

  !-------------------------------------------------
  ! communicate distance of particle furthest away
  !-------------------------------------------------
  call build_peak_communicator
  call virtual_peak_dp(cmp_distances(1,nmassbins), 'max')
  call boundary_peak_dp(cmp_distances(1,nmassbins))




  !-------------------------------------------
  ! Compute cumulative mass binning distances
  !-------------------------------------------
  ! The distances are not communicated later, but computed on each
  ! processor independently, because each processor has all information it needs
  ! with cmp_distances(ipeak,nmassbins) and CoM

  cmp = 0d0
  do ipeak=1, hfree-1

    ! set up distances
    if (cmp_distances(ipeak, nmassbins)>0) then
      if (logbins) then
        do i=1, nmassbins-1
          cmp_distances(ipeak,i)=rmin*(cmp_distances(ipeak,nmassbins)/rmin)**(real(i,dp)/real(nmassbins,dp))
        enddo
      else !linear binnings
        r_null=cmp_distances(ipeak,nmassbins)/real(nmassbins,dp)
        do i=0, nmassbins-1
          cmp_distances(ipeak,i)=r_null*i
        enddo
      endif
      ! The last bin must end with precicely with the maximal
      ! Distance of the particle. That is
      ! needed because precision errors. The maximal distance is
      ! computed via particle data and the outermost bin is equal
      ! to the distance of the outermost particle to the CoM.
      ! Precision errors cause the code to crash here.

      !---------------------------------------------
      ! bin particles in cumulative mass profiles:
      ! get mass of each bin
      ! calculate particle distance to CoM
      !---------------------------------------------
      thispart=clmppart_first(ipeak)
      do ipart=1, nclmppart(ipeak) ! while there is a particle linked list
        period=0d0
        if (periodical) then
          do i=1, 3
            if (xp(thispart,i)-peak_pos(ipeak,i)>0.5*boxlen)  period(i)=(-1.0)*boxlen
            if (xp(thispart,i)-peak_pos(ipeak,i)<(-0.5*boxlen)) period(i)=boxlen
          enddo
        endif
        distance=(xp(thispart,1)+period(1)-peak_pos(ipeak,1))**2 + &
          (xp(thispart,2)+period(2)-peak_pos(ipeak,2))**2 + &
          (xp(thispart,3)+period(3)-peak_pos(ipeak,3))**2
        distance=sqrt(distance)

        i=1
        do
          if (distance<=cmp_distances(ipeak,i)) then
            cmp(ipeak,i) = cmp(ipeak,i) + mp(thispart)
            exit
          else
            i=i+1
          endif
        enddo

        thispart=clmppart_next(thispart)
      enddo

      ! sum up masses to get profile instead of mass in shell
      do i=0,nmassbins-1
        cmp(ipeak,i+1)=cmp(ipeak,i+1)+cmp(ipeak,i)
      enddo

    endif  ! check
  enddo    ! loop over peaks

  !--------------------------------------
  !communicate cummulative mass profiles
  !--------------------------------------
  call build_peak_communicator()
  do i=1,nmassbins
    call virtual_peak_dp(cmp(1,i), 'sum')
    call boundary_peak_dp(cmp(1,i))
  enddo

end subroutine get_cmp_noiter
#endif
!########################################
!########################################
!########################################
subroutine get_closest_border()
  use amr_commons
  use clfind_commons
  implicit none
  !---------------------------------------------------------------------------
  ! Find closest border to centre of mass. Modified subroutine saddlepoint_search
  !---------------------------------------------------------------------------
  integer                         ::  ipart,ipeak,ip,jlevel,next_level
  integer                         ::  local_peak_id,global_peak_id
  integer,dimension(1:nvector)    ::  ind_cell
  logical,dimension(1:npeaks_max) ::  check

  ! character(len=80) :: fileloc
  ! character(len=5)  :: nchar, nchar2

  if(verbose)write(*,*) "Entered get_closest_border"

  check=.false.
  do ipeak=1, hfree-1
    ! hic sunt dracones
    ! for some reason, changing the conditions for check leads to different results
    ! for other ipeaks, so just leave it here I guess
    check(ipeak)=cmp_distances(ipeak,nmassbins)>0d0 ! peak must have particles somewhere
    check(ipeak)=check(ipeak).and.to_iter(ipeak)
    ! save some work: only check for non-namegivers
    check(ipeak)=check(ipeak).and..not.is_namegiver(ipeak)
    if(check(ipeak)) closest_border(ipeak) = HUGE(1d0) !reset value
  enddo



  !-------------------------
  ! Loop over all testcells
  !-------------------------
  ip=0
  do ipart=1,ntest
    jlevel=levp(ipart)  ! level
    next_level=0        ! level of next particle
    if(ipart<ntest)next_level=levp(ipart+1)


    global_peak_id=flag2(icellp(ipart))
    if (global_peak_id/=0) then

      call get_local_peak_id(global_peak_id,local_peak_id)

      if(check(local_peak_id)) then ! if testcell is of interest:
        ip=ip+1
        ind_cell(ip)=icellp(ipart)
        if(ip==nvector .or. next_level /= jlevel)then
          call unbinding_neighborsearch(ind_cell,ip,jlevel)
          ip=0
        endif
      endif
    endif
  enddo
  if (ip>0)call unbinding_neighborsearch(ind_cell,ip,jlevel)

  !------------------------
  ! Communicate results
  !------------------------

  call build_peak_communicator()
  call virtual_peak_dp(closest_border,'min')
  call boundary_peak_dp(closest_border)



end subroutine get_closest_border
!#####################################################
!#####################################################
!#####################################################
!#####################################################
subroutine unbinding_neighborsearch(ind_cell,np,jlevel)
  use amr_commons
  implicit none
  integer,dimension(1:nvector),intent(in) :: ind_cell  ! array of indices of cells that I want to check
  integer,intent(in)                      :: np        ! number of actual cells in ind_cell
  integer,intent(in)                      :: jlevel    ! cell level

  !------------------------------------------------------------
  ! Modified subroutine neighborsearch
  ! This routine constructs all neighboring leaf cells at levels
  ! jlevel-1, jlevel, jlevel+1.
  ! Then performs the check if the neighbors are a border
  ! in order to find the closest border to the center of mass
  !------------------------------------------------------------

  integer::j,ind,nx_loc,i1,j1,k1,i2,j2,k2,i3,j3,k3,ix,iy,iz
  integer::i1min,i1max,j1min,j1max,k1min,k1max
  integer::i2min,i2max,j2min,j2max,k2min,k2max
  integer::i3min,i3max,j3min,j3max,k3min,k3max
  real(dp)::dx,dx_loc,scale
  integer ,dimension(1:nvector)::clump_nr,indv,ind_grid,grid,ind_cell_coarse

  real(dp),dimension(1:twotondim,1:3)::xc
  integer ,dimension(1:99)::neigh_cell_index,cell_levl,test_levl
  real(dp),dimension(1:99,1:ndim)::xtest,xrel
  logical ,dimension(1:99)::ok
  real(dp),dimension(1:3)::skip_loc
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:threetondim)::nbors_father_cells_pass
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  integer::ntestpos,ntp,idim,ipos

  real(dp),dimension(1:3)::this_cellpos



  ! Mesh spacing in that level
  dx=0.5D0**jlevel
  nx_loc=(icoarse_max-icoarse_min+1)
  !skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  ! Integer constants
  i1min=0; i1max=1; i2min=0; i2max=2; i3min=0; i3max=3
  j1min=0; j1max=1; j2min=0; j2max=2; j3min=0; j3max=3
  k1min=0; k1max=1; k2min=0; k2max=2; k3min=0; k3max=3

  ! Cells center position relative to grid center position
  do ind=1,twotondim
    iz=(ind-1)/4
    iy=(ind-1-4*iz)/2
    ix=(ind-1-2*iy-4*iz)
    xc(ind,1)=(dble(ix)-0.5D0)*dx
    xc(ind,2)=(dble(iy)-0.5D0)*dx
    xc(ind,3)=(dble(iz)-0.5D0)*dx
  enddo

  ! some preliminary action...
  do j=1,np
    indv(j)   = (ind_cell(j)-ncoarse-1)/ngridmax+1         ! cell position in grid
    ind_grid(j) = ind_cell(j)-ncoarse-(indv(j)-1)*ngridmax ! grid index
    clump_nr(j) = flag2(ind_cell(j))                       ! save clump number
  enddo




  ntestpos=3**ndim
  if(jlevel>levelmin)  ntestpos=ntestpos+2**ndim
  if(jlevel<nlevelmax) ntestpos=ntestpos+4**ndim

  !===================================
  ! generate neighbors level jlevel-1
  !===================================
  ntp=0
  if(jlevel>levelmin)then
    ! Generate 2x2x2  neighboring cells at level jlevel-1
    do k1=k1min,k1max
      do j1=j1min,j1max
        do i1=i1min,i1max
          ntp=ntp+1
          xrel(ntp,1)=(2*i1-1)*dx_loc
          xrel(ntp,2)=(2*j1-1)*dx_loc
          xrel(ntp,3)=(2*k1-1)*dx_loc
          test_levl(ntp)=jlevel-1
        enddo
      enddo
    enddo
  endif

  !=====================================
  ! generate neighbors at level jlevel
  !=====================================
  ! Generate 3x3x3 neighboring cells at level jlevel
  do k2=k2min,k2max
    do j2=j2min,j2max
      do i2=i2min,i2max
        ntp=ntp+1
        xrel(ntp,1)=(i2-1)*dx_loc
        xrel(ntp,2)=(j2-1)*dx_loc
        xrel(ntp,3)=(k2-1)*dx_loc
        test_levl(ntp)=jlevel
      enddo
    enddo
  enddo

  !=====================================
  ! generate neighbors at level jlevel+1
  !=====================================
  if(jlevel<nlevelmax)then
    ! Generate 4x4x4 neighboring cells at level jlevel+1
    do k3=k3min,k3max
      do j3=j3min,j3max
        do i3=i3min,i3max
          ntp=ntp+1
          xrel(ntp,1)=(i3-1.5)*dx_loc/2.0
          xrel(ntp,2)=(j3-1.5)*dx_loc/2.0
          xrel(ntp,3)=(k3-1.5)*dx_loc/2.0
          test_levl(ntp)=jlevel+1
        enddo
      enddo
    enddo
  endif



  ! Gather 27 neighboring father cells (should be present anytime !)
  do j=1,np
    ind_cell_coarse(j)=father(ind_grid(j))
  enddo
  call get3cubefather(ind_cell_coarse,nbors_father_cells,nbors_father_grids,np,jlevel)


  do j=1,np
    ok=.false.
    do idim=1,ndim
      ! get real coordinates of neighbours
      xtest(1:ntestpos,idim)=(xg(ind_grid(j),idim)+xc(indv(j),idim)-skip_loc(idim))*scale+xrel(1:ntestpos,idim)
      if(jlevel>levelmin)xtest(1:twotondim,idim)=xtest(1:twotondim,idim)+xc(indv(j),idim)*scale
    enddo
    grid(1)=ind_grid(j)
    nbors_father_cells_pass=nbors_father_cells(j,1:threetondim)
    call get_cell_index_fast(neigh_cell_index,cell_levl,xtest,ind_grid(j),nbors_father_cells_pass,ntestpos,jlevel)

    do ipos=1,ntestpos
      ! make sure neighbour is a leaf cell
      if(son(neigh_cell_index(ipos))==0.and.cell_levl(ipos)==test_levl(ipos)) then
        ok(ipos)=.true.
      endif
    enddo

    ! get position of the cell whose neighbours will be tested
    do idim=1,ndim
      this_cellpos(idim)=(xg(ind_grid(j),idim)+xc(indv(j),idim)-skip_loc(idim))*scale
    enddo

    ! check neighbors
    call bordercheck(this_cellpos,clump_nr(j),xtest,neigh_cell_index,ok,ntestpos)
    ! bordercheck (this_cellpos=position of cell to test;
    ! clump_nr(j)=peak ID of cell to test;
    ! xtest=positions of neighbour cells;
    ! neigh_cell_index=index of neighbour cells;
    ! ok = if neighbour cell is leaf cell;
    ! ntestpos = how many neighbour cells there are
  enddo

end subroutine unbinding_neighborsearch
!########################################
!########################################
!########################################
subroutine bordercheck(this_cellpos,clump_nr,xx,neigh_cell_index,ok,np)
  !----------------------------------------------------------------------
  ! routine to check wether neighbor belongs to another clump and is closer to the center
  ! of mass than all others before
  ! modified subroutine saddlecheck
  !----------------------------------------------------------------------
  use amr_commons
  use clfind_commons
  implicit none
  real(dp), dimension(1:np,1:ndim), intent(in)  :: xx               ! positions of neighbour cells
  real(dp), dimension(1:ndim),      intent(in)  :: this_cellpos     ! position of test cell whose neighbours
                                                                    ! are to be tested
  integer,  dimension(1:99),        intent(in)  :: neigh_cell_index ! cell index of neighbours
  integer,                          intent(in)  :: clump_nr         ! global peak ID of cell whose neighbours
                                                                    ! will be tested
  logical,  dimension(1:99),     intent(inout)  :: ok               ! wether cell should be checked
  integer,                          intent(in)  :: np               ! number of neighbours to be looped over


  real(dp), dimension(1:99,1:ndim)  :: pos        ! position of border for each neighbour
  integer,  dimension(1:99)         :: neigh_cl   ! clump number of neighbour,local peak id of neighbour
  real(dp)                          :: newsum
  integer                           :: i,j,ipeak
  real(dp), dimension(1:3)          :: period


  do j=1,np
    neigh_cl(j)=flag2(neigh_cell_index(j)) ! index of the clump the neighboring cell is in

    ok(j)=ok(j).and. clump_nr/=0           ! temporary fix...
    ok(j)=ok(j).and. neigh_cl(j)/=0        ! neighboring cell is in a clump. If neighbour not in clump, clump is still considered isolated.
    ok(j)=ok(j).and. neigh_cl(j)/=clump_nr ! neighboring cell is in another clump
  enddo


  call get_local_peak_id(clump_nr,ipeak)

  do j=1,np
    if(ok(j))then ! if all criteria met, you've found a neighbour cell that belongs to a different clump

      period=0d0
      if (periodical) then
        do i=1, ndim
          if (xx(j,i)-peak_pos(ipeak,i) > 0.5*boxlen)    period(i)=(-1.0)*boxlen
          if (xx(j,i)-peak_pos(ipeak,i) < (-0.5*boxlen)) period(i)=boxlen
        enddo
      endif


      do i=1, ndim
        ! the cells will be nighbours, so no need to compute two different periodic corrections
        pos(j,i)=(xx(j,i)+period(i)+this_cellpos(i)+period(i))*0.5
      enddo

      newsum=0
      do i=1, ndim
#ifdef UNBINDINGCOM
        newsum=newsum+(pos(j,i)-clmp_com_pb(ipeak,i))**2
#else
        newsum=newsum+(pos(j,i)-peak_pos(ipeak,i))**2
#endif
      enddo

      if (newsum<closest_border(ipeak)) then
        closest_border(ipeak)=newsum
      endif
    endif
  enddo


end subroutine bordercheck
!########################################
!########################################
!########################################
subroutine particle_unbinding(ipeak, final_round)
  use amr_commons, only: dp
  use clfind_commons

  implicit none
  integer, intent(in) :: ipeak       ! peak to loop over
  logical, intent(in) :: final_round ! if it is the final round => whether to write
  !--------------------------------------------------------------
  ! This subroutine loops over all particles in the linked list of
  ! peak ipeak and checks if they are bound.
  ! Also identifies the nmost_bound most strongly bound particles
  ! Of each clump.
  !--------------------------------------------------------------

  integer :: thispart, ipeak_test, ipart, n
  real(dp):: phi_border   ! the potential at the border of the peak patch closest
                          ! to the center of mass
  real(dp):: dist_border  ! distance to the border

  real(dp),dimension(:), allocatable:: particle_energy
  integer, dimension(:), allocatable:: particle_energy_id
  real(dp) :: epart

  n=0;


  if (nclmppart(ipeak) > 0) then

    ! compute the potential for this peak on the points of the mass bin distances
    call compute_phi(ipeak)

    !-----------------------------------------------
    ! If not namegiver
    !-----------------------------------------------
    phi_border=0d0
    if (.not. is_namegiver(ipeak)) then

      ! compute potential at the closest border from the center of mass
      if(saddle_pot) then
        dist_border=sqrt(closest_border(ipeak))
        if(dist_border<=cmp_distances(ipeak,nmassbins)) then
          call potential(ipeak, dist_border, phi_border)
        endif
      endif



      !--------------------------------
      ! Not namegiver, final round
      !--------------------------------

      if (final_round) then
        ! store energy and particle clump ID in a list

        if (make_mergertree) then
          allocate(particle_energy(1:nclmppart(ipeak)))
          particle_energy = HUGE(0d0)
          allocate(particle_energy_id(1:nclmppart(ipeak)))
          particle_energy_id = 0
        endif


        !loop through particle list
        thispart=clmppart_first(ipeak)

        do ipart=1, nclmppart(ipeak)    ! loop over particle LL
          if (clmpidp(thispart)>0) then
            call get_local_peak_id(clmpidp(thispart), ipeak_test)
          else
            ipeak_test=0
          endif
          if (ipeak_test==ipeak) then   ! if this particle needs to be checked for unbinding
                                        ! particle may be assigned to child/parent clump
            candidates=candidates+1

            call eparttot(ipeak, thispart, epart)
            epart = epart + phi_border

            !check if unbound
            if(epart >= 0) then
              nunbound=nunbound+1                  ! counter
              niterunbound=niterunbound+1          ! counter 2
              clmpidp(thispart)=new_peak(ipeak)    ! update clump id
            else
              if (make_mergertree) then
                ! store the values for mergertrees
                particle_energy(ipart) = epart
                particle_energy_id(ipart) = thispart
              endif
            endif
          endif
          thispart=clmppart_next(thispart)
        enddo

        if (make_mergertree) then
          ! sort particles by lowest energy
          ! sort the particle ID's accordingly
          call quick_sort_real_int(particle_energy, particle_energy_id, nclmppart(ipeak))
          n = min(nclmppart(ipeak), nmost_bound)

          ! store bound and sorted particles
          do ipart = 1, n
            most_bound_energy(ipeak, ipart) = particle_energy(ipart)
            most_bound_pid(ipeak, ipart) = particle_energy_id(ipart)
          enddo

          ! count output  += 1
          progenitorcount_written = progenitorcount_written + 1

          deallocate(particle_energy_id, particle_energy)
        endif ! if make_mergertree





      !---------------------------------------
      ! still not namegiver, not final round
      ! (= usual iterative unbinding)
      !---------------------------------------

      else ! not final round
        if (to_iter(ipeak)) then
          hasatleastoneptcl(ipeak)=0       ! set to false

          ! loop through particle list
          thispart=clmppart_first(ipeak)

          do ipart=1, nclmppart(ipeak)     ! loop over particle LL
            if (clmpidp(thispart)>0) then
              call get_local_peak_id(clmpidp(thispart), ipeak_test)
            else
              ipeak_test=0
            endif
            if (ipeak_test==ipeak) then
              call eparttot(ipeak, thispart, epart)
              epart = epart + phi_border
              if(epart >= 0.0 ) then

                niterunbound=niterunbound+1
                contributes(thispart) = .false.   ! particle doesn't contribute to
                                                  ! clump properties
              else
                hasatleastoneptcl(ipeak)=1 ! there are contributing particles for this peak
              endif
            else
              contributes(thispart) = .false.
            endif
            thispart=clmppart_next(thispart)
          enddo

        endif
      endif

    !-----------------------------------------------
    ! If namegiver
    !-----------------------------------------------
    else ! is namegiver; only find most bound for mergertrees
      if (final_round .and. make_mergertree) then
        ! store energy and particle clump ID in a list

        allocate(particle_energy(1:nclmppart(ipeak)))
        particle_energy = HUGE(0d0)
        allocate(particle_energy_id(1:nclmppart(ipeak)))
        particle_energy_id = 0

        ! loop through particle list
        thispart=clmppart_first(ipeak)

        do ipart=1, nclmppart(ipeak)    ! loop over particle LL
          if (clmpidp(thispart)>0) then
            call get_local_peak_id(clmpidp(thispart),ipeak_test)
          else
            ipeak_test = 0
          endif
          if (ipeak_test==ipeak) then   ! if this particle needs to be checked for unbinding
                                        ! particle may be assigned to child/parent clump
            ! store the values for mergertrees
            call eparttot(ipeak, thispart, epart)
            particle_energy(ipart) = epart
            particle_energy_id(ipart) = thispart
          endif
          thispart=clmppart_next(thispart)
        enddo


        ! sort particles by lowest energy
        call quick_sort_real_int(particle_energy, particle_energy_id, nclmppart(ipeak))
        n = min(nclmppart(ipeak), nmost_bound)

        ! store bound and sorted particles
        do ipart = 1, n
          most_bound_energy(ipeak, ipart) = particle_energy(ipart)
          most_bound_pid(ipeak, ipart) = particle_energy_id(ipart)
        enddo

        ! count output  += 1 to estimate array sizes
        progenitorcount_written = progenitorcount_written + 1

        deallocate(particle_energy_id, particle_energy)

      endif ! final round and mergertrees
    endif   ! namegiver or not

  else
    hasatleastoneptcl(ipeak) = 0
  endif



end subroutine particle_unbinding
!############################################################################
!############################################################################
!############################################################################
!############################################################################
subroutine eparttot(ipeak, part_ind, epart)
  !-----------------------------------------------------------
  ! This function calculates the total energy of the particle
  ! with index part_ind assigned to clump ipeak.
  !-----------------------------------------------------------

  use pm_commons
  use clfind_commons
  implicit none

  integer, intent(in) :: ipeak, part_ind
  real(dp),intent(out):: epart
  real(dp) :: distance, kinetic_energy, minusphi
  real(dp),dimension(1:3) :: period
  integer :: i

  period=0d0
  if (periodical) then
    do i=1, ndim
      if (xp(part_ind,i)-peak_pos(ipeak,i) > 0.5*boxlen   ) period(i)=(-1.0)*boxlen
      if (xp(part_ind,i)-peak_pos(ipeak,i) < (-0.5*boxlen)) period(i)=boxlen
    enddo
  endif

#ifdef UNBINDINGCOM
  distance=(xp(part_ind,1)+period(1)-clmp_com_pb(ipeak,1))**2 + &
      (xp(part_ind,2)+period(2)-clmp_com_pb(ipeak,2))**2 + &
      (xp(part_ind,3)+period(3)-clmp_com_pb(ipeak,3))**2
#else
  distance=(xp(part_ind,1)+period(1)-peak_pos(ipeak,1))**2 + &
      (xp(part_ind,2)+period(2)-peak_pos(ipeak,2))**2 + &
      (xp(part_ind,3)+period(3)-peak_pos(ipeak,3))**2
#endif
  distance=sqrt(distance)

  kinetic_energy=0.5*((vp(part_ind,1)-clmp_vel_pb(ipeak,1))**2 + &
      (vp(part_ind,2)-clmp_vel_pb(ipeak,2))**2 + &
      (vp(part_ind,3)-clmp_vel_pb(ipeak,3))**2)

  call potential(ipeak, distance, minusphi)

  epart = kinetic_energy - minusphi
end subroutine eparttot
!########################################################
!########################################################
!########################################################
!########################################################
subroutine potential(ipeak, distance, pot)
  !------------------------------------------------------------------
  ! This function interpolates the potential of a particle for given distance
  ! It returns (-1)*phi (phi is expected to be <= 0 for gravity)
  !------------------------------------------------------------------

  use clfind_commons

  integer, intent(in) :: ipeak
  real(dp),intent(in) :: distance  ! is computed in function 'unbound', then passed
  real(dp),intent(out):: pot

  integer :: ibin, thisbin
  real(dp) :: a,b
  logical :: printedwarning

  ibin=1
  ! thisbin: the first cmp_distance which is greater than particle distance
  thisbin=0
  printedwarning = .false.

  do while (ibin<=nmassbins)
    if (distance<=cmp_distances(ipeak,ibin)) then
      thisbin=ibin
      exit
    endif
    ibin=ibin+1
  enddo

  if (thisbin==0) then
    ! for some reason, there can be minor differences in the distance calculation for intel compilers.
    ! (they do optimisations that are not value-safe if not specified otherwise)
    ! Don't crash, but report and assume it is outermost bin.
    if (.not. printedwarning) then
      write(*,'(A98,2(I9,x,A6,x),3E20.12)') &
        "(precision?) error in unbinding potential: distance > cmp_distances(ipeak,nmassbins) for ipeak", &
        ipeak, "on ID", myid, "diff:", cmp_distances(ipeak,nmassbins), distance, (distance-cmp_distances(ipeak, nmassbins))/distance
      printedwarning = .true.
    endif
    thisbin=nmassbins
  endif


  a=(phi_unb(thisbin)-phi_unb(thisbin-1))/(cmp_distances(ipeak,thisbin)-cmp_distances(ipeak,thisbin-1))
  b=phi_unb(thisbin-1)-a*cmp_distances(ipeak,thisbin-1)
  pot=(-1)*a*distance-b

end subroutine potential
!###############################################
!###############################################
!###############################################
subroutine compute_phi(ipeak)
  !-----------------------------------------------------------
  ! This subroutine computes the potential on each massbin
  ! It writes potential[r=ibin] into the array phi_unb[ibin]
  !-----------------------------------------------------------
  use clfind_commons
  use amr_commons!, only: dp
  integer, intent(in) :: ipeak
  real(dp) :: delta,add
  integer  :: i

  ! compute part of integral/sum for each bin
  phi_unb(nmassbins)=0.0
  do i=2,nmassbins
    delta=cmp_distances(ipeak,i)-cmp_distances(ipeak,i-1)
    phi_unb(i-1)=-0.5*GravConst*(cmp(ipeak,i)/cmp_distances(ipeak,i)**2+cmp(ipeak,i-1)/cmp_distances(ipeak,i-1)**2)*delta
  enddo
  delta=cmp_distances(ipeak,1)-cmp_distances(ipeak,0)
  phi_unb(0)=-0.5*GravConst*(cmp(ipeak,1)/cmp_distances(ipeak,1)**2)*delta

  ! sum bins up
  ! does not need to be done for i=nmassbins!
  add=-cmp(ipeak,nmassbins)/cmp_distances(ipeak,nmassbins)*GravConst !-G*M_tot/r_max
  do i=nmassbins-1,0,-1
    phi_unb(i)=phi_unb(i)+phi_unb(i+1) ! stops at phi_unb(1)
    phi_unb(i+1)=phi_unb(i+1)+add
  enddo
  phi_unb(0)=phi_unb(0)+add ! bypass division by 0, needed for interpolation.


end subroutine compute_phi
!######################################################################
!######################################################################
!######################################################################
subroutine dissolve_small_clumps(ilevel, for_halos, initial_cleanup)

  !----------------------------------------------------------------------
  ! Dissolve clumps with too small mass into parents/nothing.
  ! Clump is required to have at least mass_threshold number of
  ! its own particles.
  ! If mass_threshold <= 2, not removing clumps can be dangerous and
  ! lead to floating point exceptions and crashes further down in the
  ! code when computing clump properties.
  !----------------------------------------------------------------------

  use amr_commons
  use clfind_commons
  use mpi_mod

  implicit none

  integer, intent(in) :: ilevel
  logical, intent(in) :: for_halos ! whether to do it for halos or for subhalos
  logical, intent(in) :: initial_cleanup ! if it's the initial cleanup
  integer :: killed, appended
  integer :: ipeak, ipart, thispart, particle_local_id

#ifndef WITHOUTMPI
  integer, dimension(1:2) :: buf
  integer                 :: info
#endif



  !------------------------------------
  ! Kill or append too small clumps
  !------------------------------------

  if (.not. initial_cleanup) call get_exclusive_clump_mass(ilevel) ! subroutine further below

  killed = 0; appended = 0

  if (initial_cleanup) then
    ! remove all clumps that have no chance of having at least
    ! mass_threshold number of particles. Particle linked lists
    ! have been appended to parent clumps, so there should be no
    ! danger there. No need to care whether it's a namegiver or not.

    do ipeak=1, hfree-1
      if (lev_peak(ipeak) == ilevel) then

        ! if there are too few particles in there and clump is relevant (don't do unnecessary things for irrelevant stuff)
        ! for initial cleanup, look at inclusive clump mass, not exclusive
        if ( relevance(ipeak)    > relevance_threshold             .and. &
             clmp_mass_pb(ipeak) < (mass_threshold * partm_common) .and. &
             clmp_mass_pb(ipeak) > 0) then

          if (is_namegiver(ipeak)) then
            !--------------------------------
            ! if clump is namegiver, kill it
            !--------------------------------

            if(ipeak <= npeaks) killed = killed + 1 ! count only non-virtuals

            ! remove particles from clump
            thispart = clmppart_first(ipeak)
            do ipart = 1, nclmppart(ipeak)

              if (clmpidp(thispart) > 0) then
                ! no need to check that you're not removing particles of children:
                ! clmp_mass_pb at the first iteration includes all children particles as well
                ! so if halo has too little mass, children will as well
                ! otherwise, children's particles might be re-added to parent later
                clmpidp(thispart) = 0
              endif
              thispart = clmppart_next(thispart)

            enddo

            nclmppart(ipeak) = 0
            clmp_mass_pb(ipeak) = 0
            clmp_mass_exclusive(ipeak) = 0

          else
            !--------------------------------------------------
            ! if clump isn't namegiver, append it to parent
            !--------------------------------------------------
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

          endif ! namegiver or not
        endif ! if too small
      endif ! correct peak level
    enddo ! all peaks



  else ! not initial cleanup
    do ipeak=1, hfree-1
      if (lev_peak(ipeak) == ilevel) then

        ! if there are too few particles in there, but at least 1 (i.e. don't do it for noise)
        if (relevance(ipeak)           > relevance_threshold             .and. &
            clmp_mass_exclusive(ipeak) < (mass_threshold * partm_common) .and. &
            clmp_mass_exclusive(ipeak) > 0) then

          if (is_namegiver(ipeak) .and. for_halos) then
            !--------------------------------
            ! if clump is namegiver, kill it
            !--------------------------------

            if(ipeak <= npeaks) killed = killed + 1 ! count only non-virtuals

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

          endif ! namegiver or not
        endif ! if too small
      endif ! correct peak level
    enddo ! all peaks
  endif ! initial cleanup




  !---------------------------------------
  ! speak to me
  !---------------------------------------


  killed_tot = killed_tot + killed
  appended_tot = appended_tot + appended




  if (for_halos) then
    ! this is the last time this routine is called for this snapshot. Write some
    ! infos.

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
      write(*,'(A39,I12,A8)')  " Handling too small clumps: Dissolved ", killed_tot,   " halos;"
      write(*,'(A39,I12,A19)') "                               Merged ", appended_tot, " to their parents."
    endif

    ! reset values for next output step
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
      integer             :: ipeak, ipart, thispart!, i

      !----------------------------------------------------
      ! recompute clump properties after unbinding
      !----------------------------------------------------

      do ipeak=1, hfree-1 ! loop over all peaks

        ! reset values for virtual peaks to communicate multiple times
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
                  !   clmp_vel_exclusive(ipeak,i)=clmp_vel_exclusive(ipeak,i)+vp(thispart,i)*mp(thispart) ! get velocity sum
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
      call virtual_peak_dp(clmp_mass_exclusive,'sum')       ! collect
      call boundary_peak_dp(clmp_mass_exclusive)            ! scatter
      ! do i=1,3
      !   call virtual_peak_dp(clmp_vel_exclusive(1,i),'sum')  ! collect
      !   call boundary_peak_dp(clmp_vel_exclusive(1,i))       ! scatter
      ! enddo

    end subroutine get_exclusive_clump_mass

end subroutine dissolve_small_clumps






subroutine read_unbinding_params()
  use clfind_commons
  use mpi_mod
  implicit none

  namelist/unbinding_params/nmassbins,logbins,particlebased_clump_output &
       &, saddle_pot,iter_properties,conv_limit,repeat_max

  ! Read namelist file
  rewind(1)
  read(1,NML=unbinding_params,END=121)
  goto 122
121 if(myid==1)write(*,*)'You did not set up namelist &UNBDINGING_PARAMS in parameter file.'

122 rewind(1)

  if (unbind .and. .not. clumpfind) then
    if (myid==1) write(*,*) "You want particle unbinding, but didn't turn on clump finding."
    if (myid==1) write(*,*) "set clumpfind=.true. or unbind=.false. in your namelist."
    call clean_stop
  endif

  if (particlebased_clump_output .and..not. unbind) then
    ! if clumpfind = .false., we don't make it this far.
    if (myid==1) write(*,*) "You set particlebased_clump_output=.true., but not unbind=.true."
    if (myid==1) write(*,*) "I am setting unbind=.true."
    unbind=.true.
  endif
end subroutine read_unbinding_params






subroutine allocate_unbinding_arrays()
  use clfind_commons
  use pm_commons, only:npartmax
  implicit none

  !----------------------------------------------
  ! This subroutine allocates the necessary
  ! arrays and gives them initial values.
  !----------------------------------------------


  !-------------------
  ! Clump properties
  !-------------------
#ifdef UNBINDINGCOM
  allocate(clmp_com_pb(1:npeaks_max,1:3))
  clmp_com_pb=0d0
#endif
  allocate(clmp_vel_pb(1:npeaks_max,1:3))
  clmp_vel_pb=0d0
  allocate(clmp_mass_pb(1:npeaks_max))
  clmp_mass_pb=0d0
  allocate(cmp_distances(1:npeaks_max,0:nmassbins))
  cmp_distances=0d0
  allocate(cmp(1:npeaks_max,0:nmassbins))
  cmp=0d0
  ! careful with this! The first index of the second subscript
  ! of the cumulative mass aray (index 0) is there for reference
  ! for the enclosed mass interpolation.

  allocate(phi_unb(0:nmassbins)) ! array where to store the potential
  phi_unb=0d0

  if (saddle_pot) then
    allocate(closest_border(1:npeaks_max)) ! point of the closest border to CoM
    closest_border=3d0*boxlen**2
  endif

  allocate(to_iter(1:npeaks_max)) ! peak needs to be checked or not
  to_iter=.true.

  allocate(hasatleastoneptcl(1:npeaks_max))
  hasatleastoneptcl=1 ! initiate to yes

  allocate(is_namegiver(1:npeaks_max))
  is_namegiver=.false.

  !----------------------
  ! Particle linked list
  !----------------------
  allocate(clmpidp(1:npartmax))
  clmpidp=0

  allocate(clmppart_first(1:npeaks_max)) ! linked lists containing particles
  clmppart_first=0

  allocate(clmppart_last(1:npeaks_max))  ! linked lists containing particles
  clmppart_last=0

  allocate(clmppart_next(1:npartmax))    ! linked lists containing particles
  clmppart_next=0

  allocate(nclmppart(1:npeaks_max))      ! linked lists containing particles
  nclmppart=0

  allocate(contributes(1:npartmax))      ! particle contributes to clump properties or not
  contributes=.true.



  !------------------------
  ! Merger trees
  !------------------------

  if (make_mergertree) then
    allocate(most_bound_energy(1:npeaks_max, 1:nmost_bound))
    most_bound_energy = HUGE(0d0)

    allocate(most_bound_pid(1:npeaks_max, 1:nmost_bound))
    most_bound_pid = 0

    allocate(clmp_mass_exclusive(1:npeaks_max))
    clmp_mass_exclusive = 0d0

    ! allocate(clmp_vel_exclusive(1:npeaks_max, 1:3))
    ! clmp_vel_exclusive = 0
  endif


end subroutine allocate_unbinding_arrays
!########################################
!########################################
!########################################
subroutine deallocate_unbinding_arrays()
  use clfind_commons
  implicit none

#ifdef UNBINDINGCOM
  deallocate(clmp_com_pb)
#endif
  deallocate(clmp_vel_pb)
  deallocate(clmp_mass_pb)
  deallocate(cmp_distances)
  deallocate(cmp)

  deallocate(phi_unb)

  if(saddle_pot) deallocate(closest_border)

  deallocate(to_iter)

  deallocate(hasatleastoneptcl)
  deallocate(contributes)
  deallocate(is_namegiver)


  deallocate(clmpidp)
  deallocate(clmppart_last)
  deallocate(clmppart_first)
  deallocate(clmppart_next)
  deallocate(nclmppart)

  return

end subroutine deallocate_unbinding_arrays
!############################################
!############################################
!############################################
! endif: NDIM == 3
#endif
