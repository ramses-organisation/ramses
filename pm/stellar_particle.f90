subroutine make_stellar_from_sinks
  use pm_commons
  use amr_commons
  use sink_feedback_parameters
  use mpi_mod
  implicit none

  integer:: isink

  integer:: nbuf
  integer, parameter:: nbufmax = 1000
  integer, dimension(1:nbufmax):: buf_id
  real(dp):: mass_total
  integer:: iobj,nobj_new
  real(dp), dimension(1:nsink) :: dmfsink_sort

  if(.not. hydro) return
  if(ndim /= 3) return

  if(verbose) write(*,*) 'Entering make_stellar_from_sinks'

  nbuf = 0

  if(stellar_strategy=='local')then
     ! Check for each sink whether a stellar particle (or multiple) should be created
     do isink = 1, nsink
        do while(dmfsink(isink) .gt. stellar_msink_th)
           dmfsink(isink) = dmfsink(isink) - stellar_msink_th
           nbuf = nbuf + 1
           if(nbuf > nbufmax) then
              call create_stellar(nbufmax, nbufmax, buf_id)
              nbuf = 1
           end if
           ! save ID of sink
           buf_id(nbuf) = idsink(isink)
        end do
     end do
  else !global
     ! Determine how many new stellar objects must be created by summing all the mass
     ! within sinks and looking at how many objects have been created already.
     ! It then places them on the sinks which have the largest dmfsink, i.e. recently accreted gas

     !compare the total number of objects formed and the mass of the sinks
     mass_total = sum(dmfsink)

     !number of objects to be created
     nobj_new = mass_total / stellar_msink_th

     !order the sinks by recently accreted mass
     dmfsink_sort = dmfsink
     call quick_sort_dp(dmfsink_sort,idsink_sort,nsink)

     !loop over the sinks ranked by recently accreted mass
     !one object per sink is assumed
     !this assumes that the number of sinks is larger than the number of objects
     if(nobj_new .gt. nsink) then
        write(*,*) 'number of new objects is larger than the number of sinks ',nobj_new, nsink
        write(*,*) "use stellar_strategy='local'"
        stop
     endif

     do iobj = nsink - nobj_new + 1, nsink
        isink = idsink_sort(iobj)
        !note with this formulation dmfsink can be negative
        dmfsink(isink) = dmfsink(isink) - stellar_msink_th

        nbuf = nbuf + 1
        if(nbuf > nbufmax) then
           call create_stellar(nbufmax, nbufmax, buf_id)
           nbuf = 1
        end if

        buf_id(nbuf) = idsink(isink)
     end do

  endif

  call create_stellar(nbuf, nbufmax, buf_id)

  if (stellar_info)then
    call print_stellar_properties
  end if

end subroutine make_stellar_from_sinks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine create_stellar(ncreate, nbuf, id_new)
    use amr_commons, only: dp, myid, ncpu, ndim, t
    use sink_feedback_parameters
    use constants, only:M_sun
    use mpi_mod
    implicit none
    !------------------------------------------------------------------------
    ! Create new stellar objects
    !------------------------------------------------------------------------
    integer, intent(in):: ncreate, nbuf
    integer, dimension(1:nbuf), intent(in):: id_new
    integer:: ncreate_loc
    real(dp), dimension(1:ncreate):: mnew_loc, ltnew_loc
    real(dp), dimension(1:ncreate):: mnew, tnew, ltnew
#ifndef WITHOUTMPI
    integer, dimension(1:ncpu)::displ
    integer:: info, icpu, idim, isplit, nsplit
    integer, dimension(1:ncpu):: narr
#endif
    integer:: istellar

    if(ncreate == 0) return

    ! Check that there is enough space
    if(ncreate + nstellar > nstellarmax) then
        if(myid == 1) write(*, *) 'Not enough space for new stellar objects! Increase nstellarmax.'
        call clean_stop
    end if

#ifndef WITHOUTMPI
    ! Split work among processes (determining mass and lifetime)
    isplit = mod(ncreate, ncpu)
    nsplit = ncreate / ncpu
    narr(       1:isplit) = nsplit + 1
    narr(isplit+1:  ncpu) = nsplit
    displ(1) = 0
    do icpu = 2, ncpu
        displ(icpu) = displ(icpu - 1) + narr(icpu - 1)
    end do
    ncreate_loc = narr(myid)
#else
    ncreate_loc = ncreate
#endif

    ! Draw random masses from the IMF
    call sample_powerlaw(mnew_loc, imf_low, imf_high, imf_index, ncreate_loc)

    ! Compute lifetime
    ltnew_loc(1:ncreate_loc) = lt_t0 * exp(lt_a * (log(lt_m0 / mnew_loc))**lt_b)

    ! Communicate data
#ifndef WITHOUTMPI
    call MPI_ALLGATHERV(  mnew_loc, ncreate_loc, MPI_DOUBLE_PRECISION, mnew, narr, displ, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, info)
    call MPI_ALLGATHERV( ltnew_loc, ncreate_loc, MPI_DOUBLE_PRECISION, ltnew, narr, displ, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, info)
#else
    mnew = mnew_loc
    ltnew = ltnew_loc
#endif

    ! Add new objects to the arrays
    id_stellar(nstellar+1:nstellar+ncreate) = id_new(1:ncreate)
    ! Set birth time to current time
    tstellar(nstellar+1:nstellar+ncreate) = t

    ! Set stellar masses and lifetimes
    ! EITHER: use mnew
    ! OR: use mstellarini if this has non-zero values (sets first stellar objects to have specific mass)
    !TC: remove mstellarini? I can think of only a few cases where this is useful.
    !    Maybe not worth having this extra loop and checks
    do istellar = nstellar+1, nstellar+ncreate
       if(istellar .ge. nstellarini) then
          mstellar(istellar) = mnew(istellar-nstellar)
       else
          if (mstellarini(istellar).eq.0) then
            ! No initialised stellar mass? That's ok, use this one
             mstellar(istellar) = mnew(istellar-nstellar)
          else
            ! Current value? Leave it alone but overwrite the age with the correct one
             mstellar(istellar) = mstellarini(istellar)
             ltnew(istellar-nstellar) = lt_t0 * exp(lt_a * (log(lt_m0 / mstellar(istellar)))**lt_b)
         endif
       endif
    end do
    !mstellar(nstellar+1:nstellar+ncreate) = mnew
    ltstellar(nstellar+1:nstellar+ncreate) = ltnew

    if(myid == 1) then
        write(*, "('Created ', I5, ' stellar objects')") ncreate
    end if

    if(sn_direct) then
        ! explode immediately instead of after lifetime
        ltstellar(nstellar+1:nstellar+ncreate) = 0
    end if

    nstellar = nstellar + ncreate

end subroutine create_stellar
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine delete_stellar(flag_delete)
    use pm_commons
    use amr_commons
    use sink_feedback_parameters
    use mpi_mod
    implicit none
    !------------------------------------------------------------------------
    ! Delete flagged stellar objects
    !------------------------------------------------------------------------
    logical, dimension(1:nstellar), intent(in):: flag_delete !true if a particle should be deleted
    integer:: i, inew
    logical, dimension(1:nstellar):: flag_any

#ifndef WITHOUTMPI
    integer:: info

    ! Make sure every process deletes the same objects
    call MPI_ALLREDUCE(flag_delete, flag_any, nstellar, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, info)
#else
    flag_any = flag_delete
#endif

    inew = 1
    do i = 1, nstellar
        if(.not. flag_any(i)) then
            if(i > inew) then
                id_stellar(inew) = id_stellar(i)
                mstellar(inew) = mstellar(i)
                tstellar(inew) = tstellar(i)
                ltstellar(inew) = ltstellar(i)
            end if
            inew = inew + 1
        end if
    end do

    ! Update nstellar
    nstellar = inew - 1

end subroutine delete_stellar
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine sample_powerlaw(x, a, b, alpha, n)
    ! Sample from a power-law between a and b, with an index of alpha (for the PDF)
    use amr_commons  ,only:ncpu,myid
    use pm_commons   ,only:localseed,iseed
    use random
    implicit none
    real(8), dimension(1:n), intent(out):: x
    real(8), intent(in):: a, b, alpha
    integer, intent(in):: n
    integer, dimension(1:ncpu,1:IRandNumSize)::allseed
    real(8):: u, p, q
    integer:: i

    p = alpha + 1.0_8
    q = 1.0_8 / p

    ! If necessary, initialize random number generator
    if(localseed(1)==-1)then
        call rans(ncpu,iseed,allseed)
        localseed=allseed(myid,1:IRandNumSize)
    end if

    do i = 1, n
        call Ranf(localseed, u)
        ! u follows an uniform law between 0 and 1
        ! Scale it to b^p..a^p
        u = b**p + (a**p - b**p) * u
        x(i) = u**q
    end do

end subroutine sample_powerlaw
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine print_stellar_properties
    use amr_commons
    use sink_feedback_parameters
    use constants, only: M_sun, yr2sec
    implicit none
    integer::i,istellar
    real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m

    if(myid==1.and.nstellar>0.and.mod(nstep_coarse,ncontrol)==0) then
        ! Scaling factors
        call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
        scale_m=scale_d*scale_l**ndim

        ! sort by remaining lifetime
        !do istellar=1,nstellar
        !    time_remaining(istellar) = ltstellar(istellar) - (t-tstellar(istellar))
        !end do
        !call quick_sort_dp(time_remaining(1),idstellar_sort(1),nstellar)

        write(*,*)'Number of stellar objects = ',nstellar
        write(*, "('***********************************************')")
        write(*, "('   id    mass[Msol]     age[yr]    lifetime[yr]')")
        write(*, "('***********************************************')")
        do i=1,nstellar
            !istellar=idstellar_sort(i)
            write(*, "(I5,3(2X,1PE12.5))") id_stellar(i), &
                & mstellar(i)*scale_m/M_sun, (t-tstellar(i))*scale_t/yr2sec, ltstellar(i)*scale_t/yr2sec
        end do
        write(*,"('***********************************************')")
    end if

  end subroutine print_stellar_properties
