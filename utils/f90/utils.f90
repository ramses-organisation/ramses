module utils
  use iso_fortran_env
  implicit none

#ifdef longint
  integer, parameter :: i8b = 8
#else
  integer, parameter :: i8b = 4
#endif

  interface hilbert3d
     subroutine hilbert3d_single(x,y,z,order,bit_length,npoint)
       integer, intent(in) :: bit_length, npoint
       integer, intent(in)::x,y,z
       real(kind=8),intent(out)::order
     end subroutine hilbert3d_single

     subroutine hilbert3d_vector(x,y,z,order,bit_length,npoint)
       integer, intent(in) :: bit_length, npoint
       integer, intent(in), dimension(1:1)::x,y,z
       real(kind=8),intent(out),dimension(1:npoint)::order
     end subroutine hilbert3d_vector
  end interface hilbert3d


  interface quick_sort
     subroutine quick_sort_int32(list, order, n)
       use iso_fortran_env
       integer, intent(in) :: n
       integer(int32), dimension (1:n), intent(inout)  :: list
       integer, dimension (1:n), intent(out)  :: order
     end subroutine quick_sort_int32

     subroutine quick_sort_int64(list, order, n)
       use iso_fortran_env
       integer, intent(in) :: n
       integer(int64), dimension (1:n), intent(inout)  :: list
       integer, dimension (1:n), intent(out)  :: order
     end subroutine quick_sort_int64

     subroutine quick_sort_real32(list, order, n)
       use iso_fortran_env
       integer, intent(in) :: n
       real(real32), dimension (1:n), intent(inout)  :: list
       integer, dimension (1:n), intent(out)  :: order
     end subroutine quick_sort_real32

     subroutine quick_sort_real64(list, order, n)
       use iso_fortran_env
       integer, intent(in) :: n
       real(real64), dimension (1:n), intent(inout)  :: list
       integer, dimension (1:n), intent(out)  :: order
     end subroutine quick_sort_real64

  end interface quick_sort

contains
  !=======================================================================
  subroutine title(n,nchar)
    integer::n
    character(5)::nchar

    character(1)::nchar1
    character(2)::nchar2
    character(3)::nchar3
    character(4)::nchar4
    character(5)::nchar5

    if(n.ge.10000)then
       write(nchar5,'(i5)') n
       nchar = nchar5
    elseif(n.ge.1000)then
       write(nchar4,'(i4)') n
       nchar = '0'//nchar4
    elseif(n.ge.100)then
       write(nchar3,'(i3)') n
       nchar = '00'//nchar3
    elseif(n.ge.10)then
       write(nchar2,'(i2)') n
       nchar = '000'//nchar2
    else
       write(nchar1,'(i1)') n
       nchar = '0000'//nchar1
    endif

  end subroutine title

  subroutine friedman(O_mat_0,O_vac_0,O_k_0,alpha,axp_min, &
       & axp_out,hexp_out,tau_out,t_out,ntable,age_tot)
    integer::ntable
    real(kind=8)::O_mat_0, O_vac_0, O_k_0
    real(kind=8)::alpha,axp_min,age_tot
    real(kind=8),dimension(0:ntable)::axp_out,hexp_out,tau_out,t_out
    ! ######################################################!
    ! This subroutine assumes that axp = 1 at z = 0 (today) !
    ! and that t and tau = 0 at z = 0 (today).              !
    ! axp is the expansion factor, hexp the Hubble constant !
    ! defined as hexp=1/axp*daxp/dtau, tau the conformal    !
    ! time, and t the look-back time, both in unit of 1/H0. !
    ! alpha is the required accuracy and axp_min is the     !
    ! starting expansion factor of the look-up table.       !
    ! ntable is the required size of the look-up table.     !
    ! ######################################################!
    real(kind=8)::axp_tau, axp_t
    real(kind=8)::axp_tau_pre, axp_t_pre
    real(kind=8)::dadtau, dadt
    real(kind=8)::dtau,dt
    real(kind=8)::tau,t
    integer::nstep,nout,nskip

    !  if( (O_mat_0+O_vac_0+O_k_0) .ne. 1.0D0 )then
    !     write(*,*)'Error: non-physical cosmological constants'
    !     write(*,*)'O_mat_0,O_vac_0,O_k_0=',O_mat_0,O_vac_0,O_k_0
    !     write(*,*)'The sum must be equal to 1.0, but '
    !     write(*,*)'O_mat_0+O_vac_0+O_k_0=',O_mat_0+O_vac_0+O_k_0
    !     stop
    !  end if

    axp_tau = 1.0D0
    axp_t = 1.0D0
    tau = 0.0D0
    t = 0.0D0
    nstep = 0

    do while ( (axp_tau .ge. axp_min) .or. (axp_t .ge. axp_min) )

       nstep = nstep + 1
       dtau = alpha * axp_tau / dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)
       axp_tau_pre = axp_tau - dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)*dtau/2.d0
       axp_tau = axp_tau - dadtau(axp_tau_pre,O_mat_0,O_vac_0,O_k_0)*dtau
       tau = tau - dtau

       dt = alpha * axp_t / dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
       axp_t_pre = axp_t - dadt(axp_t,O_mat_0,O_vac_0,O_k_0)*dt/2.d0
       axp_t = axp_t - dadt(axp_t_pre,O_mat_0,O_vac_0,O_k_0)*dt
       t = t - dt

    end do

    age_tot=-t
    write(*,666)-t
666 format(' Age of the Universe (in unit of 1/H0)=',1pe10.3)

    nskip=nstep/ntable

    axp_t = 1.d0
    t = 0.d0
    axp_tau = 1.d0
    tau = 0.d0
    nstep = 0
    nout=0
    t_out(nout)=t
    tau_out(nout)=tau
    axp_out(nout)=axp_tau
    hexp_out(nout)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau

    do while ( (axp_tau .ge. axp_min) .or. (axp_t .ge. axp_min) )

       nstep = nstep + 1
       dtau = alpha * axp_tau / dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)
       axp_tau_pre = axp_tau - dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)*dtau/2.d0
       axp_tau = axp_tau - dadtau(axp_tau_pre,O_mat_0,O_vac_0,O_k_0)*dtau
       tau = tau - dtau

       dt = alpha * axp_t / dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
       axp_t_pre = axp_t - dadt(axp_t,O_mat_0,O_vac_0,O_k_0)*dt/2.d0
       axp_t = axp_t - dadt(axp_t_pre,O_mat_0,O_vac_0,O_k_0)*dt
       t = t - dt

       if(mod(nstep,nskip)==0)then
          nout=nout+1
          t_out(nout)=t
          tau_out(nout)=tau
          axp_out(nout)=axp_tau
          hexp_out(nout)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau
       end if

    end do
    t_out(ntable)=t
    tau_out(ntable)=tau
    axp_out(ntable)=axp_tau
    hexp_out(ntable)=dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)/axp_tau

  end subroutine friedman

end module utils

subroutine hilbert3d_vector(x,y,z,order,bit_length,npoint)
  implicit none

  integer     ,INTENT(IN)                     ::bit_length,npoint
  integer     ,INTENT(IN) ,dimension(1:npoint)::x,y,z
  real(kind=8),INTENT(OUT),dimension(1:npoint)::order

  logical,dimension(0:3*bit_length-1)::i_bit_mask
  logical,dimension(0:1*bit_length-1)::x_bit_mask,y_bit_mask,z_bit_mask
  integer,dimension(0:7,0:1,0:11)::state_diagram
  integer::i,ip,cstate,nstate,b0,b1,b2,sdigit,hdigit

  if(bit_length>bit_size(bit_length))then
     write(*,*)'Maximum bit length=',bit_size(bit_length)
     write(*,*)'stop in hilbert3d'
     stop
  endif

  state_diagram = RESHAPE( (/   1, 2, 3, 2, 4, 5, 3, 5,&
       &   0, 1, 3, 2, 7, 6, 4, 5,&
       &   2, 6, 0, 7, 8, 8, 0, 7,&
       &   0, 7, 1, 6, 3, 4, 2, 5,&
       &   0, 9,10, 9, 1, 1,11,11,&
       &   0, 3, 7, 4, 1, 2, 6, 5,&
       &   6, 0, 6,11, 9, 0, 9, 8,&
       &   2, 3, 1, 0, 5, 4, 6, 7,&
       &  11,11, 0, 7, 5, 9, 0, 7,&
       &   4, 3, 5, 2, 7, 0, 6, 1,&
       &   4, 4, 8, 8, 0, 6,10, 6,&
       &   6, 5, 1, 2, 7, 4, 0, 3,&
       &   5, 7, 5, 3, 1, 1,11,11,&
       &   4, 7, 3, 0, 5, 6, 2, 1,&
       &   6, 1, 6,10, 9, 4, 9,10,&
       &   6, 7, 5, 4, 1, 0, 2, 3,&
       &  10, 3, 1, 1,10, 3, 5, 9,&
       &   2, 5, 3, 4, 1, 6, 0, 7,&
       &   4, 4, 8, 8, 2, 7, 2, 3,&
       &   2, 1, 5, 6, 3, 0, 4, 7,&
       &   7, 2,11, 2, 7, 5, 8, 5,&
       &   4, 5, 7, 6, 3, 2, 0, 1,&
       &  10, 3, 2, 6,10, 3, 4, 4,&
       &   6, 1, 7, 0, 5, 2, 4, 3 /), &
       & (/8 ,2, 12 /) )

  do ip=1,npoint

     ! convert to binary
     do i=0,bit_length-1
        x_bit_mask(i)=btest(x(ip),i)
        y_bit_mask(i)=btest(y(ip),i)
        z_bit_mask(i)=btest(z(ip),i)
     enddo

     ! interleave bits
     do i=0,bit_length-1
        i_bit_mask(3*i+2)=x_bit_mask(i)
        i_bit_mask(3*i+1)=y_bit_mask(i)
        i_bit_mask(3*i  )=z_bit_mask(i)
     end do

     ! build Hilbert ordering using state diagram
     cstate=0
     do i=bit_length-1,0,-1
        b2=0 ; if(i_bit_mask(3*i+2))b2=1
        b1=0 ; if(i_bit_mask(3*i+1))b1=1
        b0=0 ; if(i_bit_mask(3*i  ))b0=1
        sdigit=b2*4+b1*2+b0
        nstate=state_diagram(sdigit,0,cstate)
        hdigit=state_diagram(sdigit,1,cstate)
        i_bit_mask(3*i+2)=btest(hdigit,2)
        i_bit_mask(3*i+1)=btest(hdigit,1)
        i_bit_mask(3*i  )=btest(hdigit,0)
        cstate=nstate
     enddo

     ! save Hilbert key as double precision real
     order(ip)=0.
     do i=0,3*bit_length-1
        b0=0 ; if(i_bit_mask(i))b0=1
        order(ip)=order(ip)+dble(b0)*dble(2)**i
     end do

  end do

end subroutine hilbert3d_vector

subroutine hilbert3d_single(x,y,z,order,bit_length,npoint)
  implicit none

  integer     ,INTENT(in)::bit_length,npoint
  integer     ,INTENT(in)::x,y,z
  real(kind=8),intent(out)::order
  real(kind=8),dimension(1:1)::tmporder

  integer, dimension(1:1) :: xt, yt, zt

  xt(1) = x; yt(1) = y; zt(1) = z
  call hilbert3d_vector(xt, yt, zt, tmporder, bit_length, npoint)
  order = tmporder(1)

end subroutine hilbert3d_single

subroutine quick_sort_int32(list, order, n)
  use iso_fortran_env
  implicit none

  ! quick sort routine from:
  ! brainerd, w.s., goldberg, c.h. & adams, j.c. (1990) "programmer's guide to
  ! fortran 90", mcgraw-hill  isbn 0-07-000248-7, pages 149-150.
  ! modified by alan miller to include an associated integer array which gives
  ! the positions of the elements in the original order.
  integer, intent(in) :: n
  integer(int32), dimension (1:n), intent(inout)  :: list
  integer, dimension (1:n), intent(out)  :: order

  ! local variable
  integer :: i

  do i = 1, n
     order(i) = i
  end do

  call quick_sort_1(1, n)

contains

  recursive subroutine quick_sort_1(left_end, right_end)

    integer, intent(in) :: left_end, right_end

    !     local variables
    integer             :: i, j, itemp
    integer(int32)        :: reference, temp
    integer, parameter  :: max_simple_sort_size = 6

    if (right_end < left_end + max_simple_sort_size) then
       ! use interchange sort for small lists
       call interchange_sort(left_end, right_end)

    else
       ! use partition ("quick") sort
       reference = list((left_end + right_end)/2)
       i = left_end - 1; j = right_end + 1

       do
          ! scan list from left end until element >= reference is found
          do
             i = i + 1
             if (list(i) >= reference) exit
          end do
          ! scan list from right end until element <= reference is found
          do
             j = j - 1
             if (list(j) <= reference) exit
          end do


          if (i < j) then
             ! swap two out-of-order elements
             temp = list(i); list(i) = list(j); list(j) = temp
             itemp = order(i); order(i) = order(j); order(j) = itemp
          else if (i == j) then
             i = i + 1
             exit
          else
             exit
          end if
       end do

       if (left_end < j) call quick_sort_1(left_end, j)
       if (i < right_end) call quick_sort_1(i, right_end)
    end if

  end subroutine quick_sort_1

  subroutine interchange_sort(left_end, right_end)

    integer, intent(in) :: left_end, right_end

    !     local variables
    integer             :: i, j, itemp
    integer(int32)        :: temp

    do i = left_end, right_end - 1
       do j = i+1, right_end
          if (list(i) > list(j)) then
             temp = list(i); list(i) = list(j); list(j) = temp
             itemp = order(i); order(i) = order(j); order(j) = itemp
          end if
       end do
    end do

  end subroutine interchange_sort

end subroutine quick_sort_int32

subroutine quick_sort_int64(list, order, n)
  use iso_fortran_env
  implicit none

  ! quick sort routine from:
  ! brainerd, w.s., goldberg, c.h. & adams, j.c. (1990) "programmer's guide to
  ! fortran 90", mcgraw-hill  isbn 0-07-000248-7, pages 149-150.
  ! modified by alan miller to include an associated integer array which gives
  ! the positions of the elements in the original order.
  integer, intent(in) :: n
  integer(int64), dimension (1:n), intent(inout)  :: list
  integer, dimension (1:n), intent(out)  :: order

  ! local variable
  integer :: i

  do i = 1, n
     order(i) = i
  end do

  call quick_sort_1(1, n)

contains

  recursive subroutine quick_sort_1(left_end, right_end)

    integer, intent(in) :: left_end, right_end

    !     local variables
    integer             :: i, j, itemp
    integer(int64)        :: reference, temp
    integer, parameter  :: max_simple_sort_size = 6

    if (right_end < left_end + max_simple_sort_size) then
       ! use interchange sort for small lists
       call interchange_sort(left_end, right_end)

    else
       ! use partition ("quick") sort
       reference = list((left_end + right_end)/2)
       i = left_end - 1; j = right_end + 1

       do
          ! scan list from left end until element >= reference is found
          do
             i = i + 1
             if (list(i) >= reference) exit
          end do
          ! scan list from right end until element <= reference is found
          do
             j = j - 1
             if (list(j) <= reference) exit
          end do


          if (i < j) then
             ! swap two out-of-order elements
             temp = list(i); list(i) = list(j); list(j) = temp
             itemp = order(i); order(i) = order(j); order(j) = itemp
          else if (i == j) then
             i = i + 1
             exit
          else
             exit
          end if
       end do

       if (left_end < j) call quick_sort_1(left_end, j)
       if (i < right_end) call quick_sort_1(i, right_end)
    end if

  end subroutine quick_sort_1

  subroutine interchange_sort(left_end, right_end)

    integer, intent(in) :: left_end, right_end

    !     local variables
    integer             :: i, j, itemp
    integer(int64)      :: temp

    do i = left_end, right_end - 1
       do j = i+1, right_end
          if (list(i) > list(j)) then
             temp = list(i); list(i) = list(j); list(j) = temp
             itemp = order(i); order(i) = order(j); order(j) = itemp
          end if
       end do
    end do

  end subroutine interchange_sort

end subroutine quick_sort_int64

subroutine quick_sort_real32(list, order, n)
  use iso_fortran_env

  implicit none

  ! quick sort routine from:
  ! brainerd, w.s., goldberg, c.h. & adams, j.c. (1990) "programmer's guide to
  ! fortran 90", mcgraw-hill  isbn 0-07-000248-7, pages 149-150.
  ! modified by alan miller to include an associated integer array which gives
  ! the positions of the elements in the original order.
  integer, intent(in) :: n
  real(real32), dimension (1:n), intent(inout)  :: list
  integer, dimension (1:n), intent(out)  :: order

  ! local variable
  integer :: i

  do i = 1, n
     order(i) = i
  end do

  call quick_sort_1(1, n)

contains

  recursive subroutine quick_sort_1(left_end, right_end)

    integer, intent(in) :: left_end, right_end

    !     local variables
    integer             :: i, j, itemp
    real(real32)        :: reference, temp
    integer, parameter  :: max_simple_sort_size = 6

    if (right_end < left_end + max_simple_sort_size) then
       ! use interchange sort for small lists
       call interchange_sort(left_end, right_end)

    else
       ! use partition ("quick") sort
       reference = list((left_end + right_end)/2)
       i = left_end - 1; j = right_end + 1

       do
          ! scan list from left end until element >= reference is found
          do
             i = i + 1
             if (list(i) >= reference) exit
          end do
          ! scan list from right end until element <= reference is found
          do
             j = j - 1
             if (list(j) <= reference) exit
          end do


          if (i < j) then
             ! swap two out-of-order elements
             temp = list(i); list(i) = list(j); list(j) = temp
             itemp = order(i); order(i) = order(j); order(j) = itemp
          else if (i == j) then
             i = i + 1
             exit
          else
             exit
          end if
       end do

       if (left_end < j) call quick_sort_1(left_end, j)
       if (i < right_end) call quick_sort_1(i, right_end)
    end if

  end subroutine quick_sort_1

  subroutine interchange_sort(left_end, right_end)

    integer, intent(in) :: left_end, right_end

    !     local variables
    integer             :: i, j, itemp
    real(real32)        :: temp

    do i = left_end, right_end - 1
       do j = i+1, right_end
          if (list(i) > list(j)) then
             temp = list(i); list(i) = list(j); list(j) = temp
             itemp = order(i); order(i) = order(j); order(j) = itemp
          end if
       end do
    end do

  end subroutine interchange_sort

end subroutine quick_sort_real32

subroutine quick_sort_real64(list, order, n)
  use iso_fortran_env

  implicit none

  ! quick sort routine from:
  ! brainerd, w.s., goldberg, c.h. & adams, j.c. (1990) "programmer's guide to
  ! fortran 90", mcgraw-hill  isbn 0-07-000248-7, pages 149-150.
  ! modified by alan miller to include an associated integer array which gives
  ! the positions of the elements in the original order.
  integer, intent(in) :: n
  real(real64), dimension (1:n), intent(inout)  :: list
  integer, dimension (1:n), intent(out)  :: order

  ! local variable
  integer :: i

  do i = 1, n
     order(i) = i
  end do

  call quick_sort_1(1, n)

contains

  recursive subroutine quick_sort_1(left_end, right_end)

    integer, intent(in) :: left_end, right_end

    !     local variables
    integer             :: i, j, itemp
    real(real64)        :: reference, temp
    integer, parameter  :: max_simple_sort_size = 6

    if (right_end < left_end + max_simple_sort_size) then
       ! use interchange sort for small lists
       call interchange_sort(left_end, right_end)

    else
       ! use partition ("quick") sort
       reference = list((left_end + right_end)/2)
       i = left_end - 1; j = right_end + 1

       do
          ! scan list from left end until element >= reference is found
          do
             i = i + 1
             if (list(i) >= reference) exit
          end do
          ! scan list from right end until element <= reference is found
          do
             j = j - 1
             if (list(j) <= reference) exit
          end do


          if (i < j) then
             ! swap two out-of-order elements
             temp = list(i); list(i) = list(j); list(j) = temp
             itemp = order(i); order(i) = order(j); order(j) = itemp
          else if (i == j) then
             i = i + 1
             exit
          else
             exit
          end if
       end do

       if (left_end < j) call quick_sort_1(left_end, j)
       if (i < right_end) call quick_sort_1(i, right_end)
    end if

  end subroutine quick_sort_1

  subroutine interchange_sort(left_end, right_end)

    integer, intent(in) :: left_end, right_end

    !     local variables
    integer             :: i, j, itemp
    real(real64)        :: temp

    do i = left_end, right_end - 1
       do j = i+1, right_end
          if (list(i) > list(j)) then
             temp = list(i); list(i) = list(j); list(j) = temp
             itemp = order(i); order(i) = order(j); order(j) = itemp
          end if
       end do
    end do

  end subroutine interchange_sort

end subroutine quick_sort_real64

function dadtau(axp_tau,O_mat_0,O_vac_0,O_k_0)
  implicit none

  real(kind=8)::dadtau,axp_tau,O_mat_0,O_vac_0,O_k_0
  dadtau = axp_tau*axp_tau*axp_tau *  &
       &   ( O_mat_0 + &
       &     O_vac_0 * axp_tau*axp_tau*axp_tau + &
       &     O_k_0   * axp_tau )
  dadtau = sqrt(dadtau)
  return
end function dadtau

function dadt(axp_t,O_mat_0,O_vac_0,O_k_0)
  implicit none

  real(kind=8)::dadt,axp_t,O_mat_0,O_vac_0,O_k_0
  dadt   = (1.0D0/axp_t)* &
       &   ( O_mat_0 + &
       &     O_vac_0 * axp_t*axp_t*axp_t + &
       &     O_k_0   * axp_t )
  dadt = sqrt(dadt)
  return
end function dadt
