!################################################################
!################################################################
!################################################################
!################################################################
module particle_defer
  implicit none
  ! NOTE: not amr_parameters' i8b, which is selected_int_kind(9) and therefore
  ! only 4 bytes here. The key packs two 32-bit fields, so it needs a real
  ! 8-byte integer.
  integer, parameter :: ikey = selected_int_kind(18)
  !--------------------------------------------------------------------------
  ! Deferred particle linked-list operations.
  !
  ! The routines that relocate particles (make_tree_fine, kill_tree_fine,
  ! virtual_tree_fine) traverse the grids in parallel but must splice the
  ! particles into shared linked lists. Doing that inside an !$omp critical is
  ! correct but not reproducible: mutual exclusion fixes *whether* threads
  ! collide, not the order in which they get the lock, so the resulting list
  ! order - and with it anything that walks the list, such as the MC tracer
  ! random draws - changes from run to run.
  !
  ! Instead each thread records its list operations here, tagged with the
  ! position the particle had in the traversal, and a single serial pass
  ! replays them in that order once the parallel region has ended. That
  ! reproduces the order a serial run would have produced, and removes the need
  ! for the named critical section altogether (which also avoids a known
  ! problem combining named critical regions with -ipo on the Intel compiler).
  !
  ! The traversal key is (icpu, jgrid): a grid is always processed in full by
  ! one thread, and within a grid the particles are visited in the order of the
  ! frozen headp_old/nextp_old list, so a *stable* sort on that key is enough
  ! to recover the full ordering without having to store the particle position
  ! within the grid as well.
  !--------------------------------------------------------------------------

  ! Operation codes
  integer, parameter :: DEFER_REMOVE_LIST = 1   ! unlink from grid linked list
  integer, parameter :: DEFER_ADD_LIST    = 2   ! append to grid linked list
  ! Note there is deliberately no DEFER_REMOVE_FREE: remove_free hands a slot
  ! back to its caller, which needs the index straight away to fill in the
  ! particle data, so it cannot be postponed. Deterministic slot allocation
  ! needs a different mechanism (reserving a block per rank up front).
  integer, parameter :: DEFER_ADD_FREE    = 3   ! give back to the free list

  type :: defer_buffer
     integer(ikey), allocatable :: key(:)
     integer,      allocatable :: op(:)
     integer,      allocatable :: part(:)
     integer,      allocatable :: grid(:)
     integer                   :: n = 0
  end type defer_buffer

  ! One buffer per thread. Indexed 1:nthr, allocated once on first use.
  type(defer_buffer), allocatable :: dbuf(:)

contains

  !----------------------------------------------------------------------
  ! Index of the calling thread, 1-based. 1 when built without OpenMP.
  !----------------------------------------------------------------------
  function defer_thread() result(ithr)
#ifdef _OPENMP
    use omp_lib
#endif
    implicit none
    integer :: ithr
#ifdef _OPENMP
    ithr = omp_get_thread_num() + 1
#else
    ithr = 1
#endif
  end function defer_thread

  !----------------------------------------------------------------------
  ! Encode the traversal position of a grid into a single sortable key.
  !----------------------------------------------------------------------
  pure function defer_key(icpu, jgrid) result(key)
    integer, intent(in) :: icpu, jgrid
    integer(ikey) :: key
    key = ishft(int(icpu, ikey), 32) + int(jgrid, ikey)
  end function defer_key

  !----------------------------------------------------------------------
  ! Reset every thread buffer, allocating them on first use.
  !----------------------------------------------------------------------
  subroutine defer_reset(nthr)
    integer, intent(in) :: nthr
    integer :: i
    if(.not.allocated(dbuf))then
       allocate(dbuf(1:nthr))
       do i=1,nthr
          allocate(dbuf(i)%key(1:1024), dbuf(i)%op(1:1024))
          allocate(dbuf(i)%part(1:1024), dbuf(i)%grid(1:1024))
       end do
    end if
    do i=1,nthr
       dbuf(i)%n = 0
    end do
  end subroutine defer_reset

  !----------------------------------------------------------------------
  ! Record one operation in the calling thread's buffer.
  !----------------------------------------------------------------------
  subroutine defer_push(ithr, key, op, ipart, igrid)
    integer,      intent(in) :: ithr, op, ipart, igrid
    integer(ikey), intent(in) :: key
    integer :: n
    n = dbuf(ithr)%n + 1
    if(n > size(dbuf(ithr)%key)) call defer_grow(ithr)
    dbuf(ithr)%key(n)  = key
    dbuf(ithr)%op(n)   = op
    dbuf(ithr)%part(n) = ipart
    dbuf(ithr)%grid(n) = igrid
    dbuf(ithr)%n = n
  end subroutine defer_push

  !----------------------------------------------------------------------
  ! Double the capacity of one thread buffer.
  !----------------------------------------------------------------------
  subroutine defer_grow(ithr)
    integer, intent(in) :: ithr
    integer(ikey), allocatable :: tkey(:)
    integer,      allocatable :: top(:), tpart(:), tgrid(:)
    integer :: old, new
    old = size(dbuf(ithr)%key)
    new = 2*old
    allocate(tkey(1:new), top(1:new), tpart(1:new), tgrid(1:new))
    tkey (1:old) = dbuf(ithr)%key (1:old)
    top  (1:old) = dbuf(ithr)%op  (1:old)
    tpart(1:old) = dbuf(ithr)%part(1:old)
    tgrid(1:old) = dbuf(ithr)%grid(1:old)
    call move_alloc(tkey , dbuf(ithr)%key )
    call move_alloc(top  , dbuf(ithr)%op  )
    call move_alloc(tpart, dbuf(ithr)%part)
    call move_alloc(tgrid, dbuf(ithr)%grid)
  end subroutine defer_grow

  !----------------------------------------------------------------------
  ! Replay every recorded operation in traversal order. Serial by design:
  ! this is the point at which the shared lists are actually modified.
  !
  ! Each thread's entries are already in ascending key order, because the
  ! !$omp do hands out chunks of increasing index and the icpu loop that
  ! encloses it is sequential within a thread. So the merge below only has to
  ! repeatedly pick the buffer with the smallest key, which is O(N*nthr) with
  ! no sorting and, being a stable merge, keeps the within-grid order intact.
  !----------------------------------------------------------------------
  subroutine defer_apply(nthr)
    integer, intent(in) :: nthr
    integer :: pos(1:nthr)
    integer :: i, ibest, total, done
    integer(ikey) :: best

    if(.not.allocated(dbuf)) return

    total = 0
    do i=1,nthr
       pos(i) = 1
       total = total + dbuf(i)%n
    end do

    do done=1,total
       ! Pick the thread whose next entry has the smallest key. Ties are
       ! broken by lowest thread index, but cannot actually occur: a grid is
       ! processed by exactly one thread.
       ibest = 0
       best = 0
       do i=1,nthr
          if(pos(i) > dbuf(i)%n) cycle
          if(ibest == 0 .or. dbuf(i)%key(pos(i)) < best)then
             ibest = i
             best = dbuf(i)%key(pos(i))
          end if
       end do
       if(ibest == 0) exit

       call defer_exec(dbuf(ibest)%op  (pos(ibest)), &
            &          dbuf(ibest)%part(pos(ibest)), &
            &          dbuf(ibest)%grid(pos(ibest)))
       pos(ibest) = pos(ibest) + 1
    end do

    do i=1,nthr
       dbuf(i)%n = 0
    end do
  end subroutine defer_apply

  !----------------------------------------------------------------------
  ! Apply a single recorded operation.
  !----------------------------------------------------------------------
  subroutine defer_exec(op, ipart, igrid)
    integer, intent(in) :: op, ipart, igrid
    select case(op)
    case(DEFER_REMOVE_LIST)
       call remove_list_one(ipart, igrid)
    case(DEFER_ADD_LIST)
       call add_list_one(ipart, igrid)
    case(DEFER_ADD_FREE)
       call add_free_one(ipart)
    end select
  end subroutine defer_exec

end module particle_defer
