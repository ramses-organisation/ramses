subroutine make_stellar_from_sinks
  use pm_commons
  use amr_commons
  use sink_feedback_parameters
  use mpi_mod
  implicit none

  integer:: isink

  integer:: nbuf
  integer, parameter:: nbufmax = 1000
  real(dp), dimension(1:nbufmax, 1:ndim):: buf_x
  integer, dimension(1:nbufmax):: buf_id
  real(dp):: mass_total
  integer:: iobj,nobj_new
  real(dp), dimension(1:nsink) :: dmfsink_sort

  if(.not. hydro) return
  if(ndim /= 3) return

  if(verbose) write(*,*) 'Entering make_stellar_from_sinks'

  nbuf = 0

  if(stellar_strategy=='local')then
    ! Check for each sink whether a stellar particle should be created
    do isink = 1, nsink
        do while(dmfsink(isink) .gt. stellar_msink_th)
          dmfsink(isink) = dmfsink(isink) - stellar_msink_th
  
          nbuf = nbuf + 1
          if(nbuf > nbufmax) then
            call create_stellar(nbufmax, nbufmax, buf_x, buf_id, .true.)
            nbuf = 1
          end if
          ! save ID and position of sink
          buf_x(nbuf, 1:ndim) = xsink(isink, 1:ndim)
          buf_id(nbuf) = isink
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
       write(*,*) 'use make_stellar_from_sinks instead of make_stellar_from_sinks_glob or modify the code'
       stop
    endif
  
    do iobj = nsink - nobj_new + 1, nsink
       isink = idsink_sort(iobj)
       !note with this formulation dmfsink can be negative
       dmfsink(isink) = dmfsink(isink) - stellar_msink_th
  
       nbuf = nbuf + 1
       if(nbuf > nbufmax) then
          call create_stellar(nbufmax, nbufmax, buf_x, buf_id, .true.)
          nbuf = 1
       end if

       buf_x(nbuf, 1:ndim) = xsink(isink, 1:ndim)
       buf_id(nbuf) = isink
    end do

  endif

  call create_stellar(nbuf, nbufmax, buf_x, buf_id, .true.)

end subroutine make_stellar_from_sinks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine create_stellar(ncreate, nbuf, xnew, id_new, print_table)
    use amr_commons, only: dp, myid, ncpu, ndim, t
    use sink_feedback_parameters
    use constants, only:M_sun
    use mpi_mod
    implicit none
    !------------------------------------------------------------------------
    ! Create new stellar objects
    !------------------------------------------------------------------------
    integer, intent(in):: ncreate, nbuf
    real(dp), dimension(1:nbuf, 1:ndim), intent(in):: xnew
    integer, dimension(1:nbuf), intent(in):: id_new
    logical, intent(in):: print_table

    integer:: ncreate_loc
    real(dp), dimension(1:ncreate):: mnew_loc, tnew_loc, ltnew_loc
    real(dp), dimension(1:ncreate):: mnew, tnew, ltnew
    real(dp), dimension(1:ncreate, 1:ndim):: xnew_loc, xnew2
    integer, dimension(1:ncreate):: id_new_loc, id_new2
#ifndef WITHOUTMPI
    integer, dimension(1:ncpu)::displ
    integer:: info, icpu, idim, isplit, nsplit
    integer, dimension(1:ncpu):: narr
#endif
    integer:: istellar

    ! debugging
    write(*, *) 'ncreate', ncreate

    if(ncreate == 0) return
    
    ! Check that there is enough space
    if(ncreate + nstellar > nstellarmax) then
        if(myid == 1) write(*, *) 'Not enough space for new stellar objects! Increase nstellarmax.'
        call clean_stop
    end if

#ifndef WITHOUTMPI
    ! Split work among processes
    isplit = mod(ncreate, ncpu)
    nsplit = ncreate / ncpu
    narr(       1:isplit) = nsplit + 1
    narr(isplit+1:  ncpu) = nsplit
    displ(1) = 0
    do icpu = 2, ncpu
        displ(icpu) = displ(icpu - 1) + narr(icpu - 1)
    end do
    ncreate_loc = narr(myid)
    xnew_loc(1:ncreate_loc, 1:ndim) = xnew(displ(myid)+1:displ(myid)+ncreate_loc, 1:ndim)
    id_new_loc(1:ncreate_loc) = id_new(displ(myid)+1:displ(myid)+ncreate_loc)
#else
    ncreate_loc = ncreate
    xnew_loc(1:ncreate_loc, 1:ndim) = xnew(1:ncreate_loc, 1:ndim)
    id_new_loc(1:ncreate_loc) = id_new(1:ncreate_loc)
#endif

    ! Draw random masses from the IMF
    !if(imf_low==imf_high)then
    !    do istellar = 1, ncreate_loc
    !        mnew_loc(istellar) = imf_low !c.u.
    !    end do
    !else
        call sample_powerlaw(mnew_loc, imf_low, imf_high, imf_index, ncreate_loc)
    !endif

    ! Set birth time to current time
    tnew_loc(1:ncreate_loc) = t

    ! Compute lifetime
    ltnew_loc(1:ncreate_loc) = lt_t0 * exp(lt_a * (log(lt_m0 / mnew_loc))**lt_b)

    ! Communicate data
#ifndef WITHOUTMPI
    do idim = 1, ndim
        call MPI_ALLGATHERV(xnew_loc(:, idim), ncreate_loc, MPI_DOUBLE_PRECISION, xnew2(:, idim), narr, displ, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, info)
    end do
    call MPI_ALLGATHERV(id_new_loc, ncreate_loc, MPI_INTEGER, id_new2, narr, displ, MPI_INTEGER, MPI_COMM_WORLD, info)
    call MPI_ALLGATHERV(  mnew_loc, ncreate_loc, MPI_DOUBLE_PRECISION, mnew,    narr, displ, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, info)
    call MPI_ALLGATHERV(  tnew_loc, ncreate_loc, MPI_DOUBLE_PRECISION, tnew,    narr, displ, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, info)
    call MPI_ALLGATHERV( ltnew_loc, ncreate_loc, MPI_DOUBLE_PRECISION, ltnew,   narr, displ, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, info)
#else
    mnew = mnew_loc
    tnew = tnew_loc
    ltnew = ltnew_loc
#endif

    ! Add new objects to the arrays
    xstellar(nstellar+1:nstellar+ncreate, 1:ndim) = xnew2(1:ncreate, 1:ndim)
    id_stellar(nstellar+1:nstellar+ncreate) = id_new2(1:ncreate)

    ! Set stellar masses
    ! EITHER: use mnew
    ! OR: use mstellarini if this has non-zero values (sets first stellar objects to have specific mass)
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
    tstellar(nstellar+1:nstellar+ncreate) = tnew
    ltstellar(nstellar+1:nstellar+ncreate) = ltnew

    if(myid == 1) then
        write(*, "('Created ', I5, ' stellar objects:')") ncreate
        if(print_table) then
            write(*, "('===================================================================================================')")
            write(*, "('       x              y              z               Mass          Birth          LifeT       id   ')")
            write(*, "('===================================================================================================')")
            do istellar = nstellar + 1, nstellar + ncreate
                write(*, "(3F15.10, 2X, 3ES15.7,2X,i6)") xstellar(istellar, 1), xstellar(istellar, 2), xstellar(istellar, 3), &
                    & mstellar(istellar), tstellar(istellar), ltstellar(istellar), id_stellar(istellar)
            end do
        end if
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
                xstellar(inew, 1:ndim) = xstellar(i, 1:ndim)
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
    use pm_commons   ,only:localseed
    use pm_parameters,only:stellar_seed
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
        call rans(ncpu,stellar_seed,allseed)
        localseed=allseed(myid,1:IRandNumSize)
    end if

    do i = 1, n
        call Ranf(localseed, u)
        write(*,*) 'random number generated ', u
        !call random_number(u)
        ! u follows an uniform law between 0 and 1
        ! Scale it to b^p..a^p
        u = b**p + (a**p - b**p) * u
        ! Calculate x(i)
        x(i) = u**q
    end do

end subroutine sample_powerlaw
