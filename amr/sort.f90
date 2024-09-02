SUBROUTINE quick_sort(list, order, n)

    ! Quick sort routine from:
    ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
    ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
    ! Modified by Alan Miller to include an associated integer array which gives
    ! the positions of the elements in the original order.

    use amr_parameters, ONLY: qdp
    IMPLICIT NONE
    integer :: n
    REAL(qdp), DIMENSION (1:n), INTENT(INOUT)  :: list
    INTEGER, DIMENSION (1:n), INTENT(OUT)  :: order

    ! Local variable
    integer :: i

    DO i = 1, n
        order(i) = i
    END DO

    CALL quick_sort_1(1, n)

CONTAINS

    RECURSIVE SUBROUTINE quick_sort_1(left_end, right_end)

        use amr_parameters, ONLY: qdp
        INTEGER, INTENT(IN) :: left_end, right_end

        ! Local variables
        integer             :: i, j, itemp
        real(qdp)           :: reference, temp
        INTEGER, PARAMETER  :: max_simple_sort_size = 6

        IF (right_end < left_end + max_simple_sort_size) THEN
            ! Use interchange sort for small lists
            CALL interchange_sort(left_end, right_end)

        ELSE
            ! Use partition ("quick") sort
            reference = list((left_end + right_end) / 2)
            i = left_end - 1; j = right_end + 1

            DO
                ! Scan list from left end until element >= reference is found
                DO
                    i = i + 1
                    IF (list(i) >= reference) EXIT
                END DO
                ! Scan list from right end until element <= reference is found
                DO
                    j = j - 1
                    IF (list(j) <= reference) EXIT
                END DO


                IF (i < j) THEN
                    ! Swap two out-of-order elements
                    temp = list(i); list(i) = list(j); list(j) = temp
                    itemp = order(i); order(i) = order(j); order(j) = itemp
                ELSE IF (i == j) THEN
                    i = i + 1
                    EXIT
                ELSE
                    EXIT
                END IF
            END DO

            IF (left_end < j) CALL quick_sort_1(left_end, j)
            IF (i < right_end) CALL quick_sort_1(i, right_end)
        END IF

    END SUBROUTINE quick_sort_1


    SUBROUTINE interchange_sort(left_end, right_end)

        use amr_parameters, ONLY: qdp
        INTEGER, INTENT(IN) :: left_end, right_end

        ! Local variables
        integer             :: i, j, itemp
        real(qdp)           :: temp

        DO i = left_end, right_end - 1
            DO j = i + 1, right_end
                IF (list(i) > list(j)) THEN
                    temp = list(i); list(i) = list(j); list(j) = temp
                    itemp = order(i); order(i) = order(j); order(j) = itemp
                END IF
            END DO
        END DO

    END SUBROUTINE interchange_sort

END SUBROUTINE quick_sort
!########################################################################
!########################################################################
!########################################################################
!########################################################################
!########################################################################
SUBROUTINE quick_sort_dp(list, order, n)
    use amr_parameters, ONLY: dp
    IMPLICIT NONE
    ! Quick sort routine from:
    ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
    ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
    ! Modified by Alan Miller to include an associated integer array which gives
    ! the positions of the elements in the original order.


    integer :: n
    REAL(dp), DIMENSION (1:n), INTENT(INOUT)  :: list
    INTEGER, DIMENSION (1:n), INTENT(OUT)  :: order

    ! Local variable
    integer :: i

    DO i = 1, n
        order(i) = i
    END DO

    CALL quick_sort_1_dp(1, n)

CONTAINS

    RECURSIVE SUBROUTINE quick_sort_1_dp(left_end, right_end)

        INTEGER, INTENT(IN) :: left_end, right_end

        ! Local variables
        integer             :: i, j, itemp
        real(kind=8)        :: reference, temp
        INTEGER, PARAMETER  :: max_simple_sort_size = 6

        IF (right_end < left_end + max_simple_sort_size) THEN
            ! Use interchange sort for small lists
            CALL interchange_sort_dp(left_end, right_end)

        ELSE
            ! Use partition ("quick") sort
            reference = list((left_end + right_end) / 2)
            i = left_end - 1; j = right_end + 1

            DO
                ! Scan list from left end until element >= reference is found
                DO
                    i = i + 1
                    IF (list(i) >= reference) EXIT
                END DO
                ! Scan list from right end until element <= reference is found
                DO
                    j = j - 1
                    IF (list(j) <= reference) EXIT
                END DO


                IF (i < j) THEN
                    ! Swap two out-of-order elements
                    temp = list(i); list(i) = list(j); list(j) = temp
                    itemp = order(i); order(i) = order(j); order(j) = itemp
                ELSE IF (i == j) THEN
                    i = i + 1
                    EXIT
                ELSE
                    EXIT
                END IF
            END DO
            IF (left_end < j) CALL quick_sort_1_dp(left_end, j)
            IF (i < right_end) CALL quick_sort_1_dp(i, right_end)
        END IF

    END SUBROUTINE quick_sort_1_dp


    SUBROUTINE interchange_sort_dp(left_end, right_end)

        INTEGER, INTENT(IN) :: left_end, right_end

        ! Local variables
        integer                :: i, j, itemp
        real(kind=8)           :: temp

        DO i = left_end, right_end - 1
            DO j = i + 1, right_end
                IF (list(i) > list(j)) THEN
                    temp = list(i); list(i) = list(j); list(j) = temp
                    itemp = order(i); order(i) = order(j); order(j) = itemp
                END IF
            END DO
        END DO

    END SUBROUTINE interchange_sort_dp

END SUBROUTINE quick_sort_dp
!########################################################################
!########################################################################
!########################################################################
!########################################################################
!########################################################################
SUBROUTINE quick_sort_real_int(list, order, n)

    ! ----------------------------------------------------------
    ! Sort array of reals (list), rearrange array of integers
    ! (order) in the same way.
    ! ----------------------------------------------------------

    use amr_parameters, ONLY: dp
    ! Quick sort routine from:
    ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
    ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
    ! Modified to sort the second given array by the same rules.

    IMPLICIT NONE
    integer :: n
    REAL(dp), DIMENSION (1:n), INTENT(INOUT)  :: list
    INTEGER, DIMENSION (1:n), INTENT(INOUT)   :: order


    CALL quick_sort_1_dp(1, n)

CONTAINS

    RECURSIVE SUBROUTINE quick_sort_1_dp(left_end, right_end)
        use amr_commons, only: dp

        INTEGER, INTENT(IN) :: left_end, right_end

        ! Local variables
        integer             :: i, j, itemp
        real(dp)            :: reference, temp
        INTEGER, PARAMETER  :: max_simple_sort_size = 6

        IF (right_end < left_end + max_simple_sort_size) THEN
            ! Use interchange sort for small lists
            CALL interchange_sort_dp(left_end, right_end)

        ELSE
            ! Use partition ("quick") sort
            reference = list((left_end + right_end) / 2)
            i = left_end - 1; j = right_end + 1

            DO
                ! Scan list from left end until element >= reference is found
                DO
                    i = i + 1
                    IF (list(i) >= reference) EXIT
                END DO
                ! Scan list from right end until element <= reference is found
                DO
                    j = j - 1
                    IF (list(j) <= reference) EXIT
                END DO


                IF (i < j) THEN
                    ! Swap two out-of-order elements
                    temp = list(i); list(i) = list(j); list(j) = temp
                    itemp = order(i); order(i) = order(j); order(j) = itemp
                ELSE IF (i == j) THEN
                    i = i + 1
                    EXIT
                ELSE
                    EXIT
                END IF
            END DO

            IF (left_end < j) CALL quick_sort_1_dp(left_end, j)
            IF (i < right_end) CALL quick_sort_1_dp(i, right_end)
        END IF

    END SUBROUTINE quick_sort_1_dp


    SUBROUTINE interchange_sort_dp(left_end, right_end)
        use amr_commons, only: dp

        INTEGER, INTENT(IN) :: left_end, right_end

        ! Local variables
        integer             :: i, j, itemp
        real(dp)            :: temp

        DO i = left_end, right_end - 1
            DO j = i + 1, right_end
                IF (list(i) > list(j)) THEN
                    temp = list(i); list(i) = list(j); list(j) = temp
                    itemp = order(i); order(i) = order(j); order(j) = itemp
                END IF
            END DO
        END DO

    END SUBROUTINE interchange_sort_dp

END SUBROUTINE quick_sort_real_int
!########################################################################
!########################################################################
!########################################################################
!########################################################################
!########################################################################
SUBROUTINE quick_sort_int_int(list, order, n)

    ! ------------------------------------------------------------
    ! Sort array of integers (list), rearrange array of integers
    ! (order) in the same way.
    ! ------------------------------------------------------------

    use amr_parameters, ONLY: i8b
    IMPLICIT NONE
    ! Quick sort routine from:
    ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
    ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
    ! Modified to sort the second given array by the same rules.


    integer :: n
    INTEGER(i8b), DIMENSION (1:n), INTENT(INOUT)  :: list
    INTEGER, DIMENSION (1:n), INTENT(INOUT)  :: order


    CALL quick_sort_1_int_int(1, n)

CONTAINS

    RECURSIVE SUBROUTINE quick_sort_1_int_int(left_end, right_end)

        INTEGER, INTENT(IN) :: left_end, right_end

        ! Local variables
        integer             :: i, j, itemp
        integer(i8b)        :: reference, temp
        INTEGER, PARAMETER  :: max_simple_sort_size = 6

        IF (right_end < left_end + max_simple_sort_size) THEN
            ! Use interchange sort for small lists
            CALL interchange_sort_int_int(left_end, right_end)

        ELSE
            ! Use partition ("quick") sort
            reference = list((left_end + right_end) / 2)
            i = left_end - 1; j = right_end + 1

            DO
                ! Scan list from left end until element >= reference is found
                DO
                    i = i + 1
                    IF (list(i) >= reference) EXIT
                END DO
                ! Scan list from right end until element <= reference is found
                DO
                    j = j - 1
                    IF (list(j) <= reference) EXIT
                END DO


                IF (i < j) THEN
                    ! Swap two out-of-order elements
                    temp = list(i); list(i) = list(j); list(j) = temp
                    itemp = order(i); order(i) = order(j); order(j) = itemp
                ELSE IF (i == j) THEN
                    i = i + 1
                    EXIT
                ELSE
                    EXIT
                END IF
            END DO
            IF (left_end < j) CALL quick_sort_1_int_int(left_end, j)
            IF (i < right_end) CALL quick_sort_1_int_int(i, right_end)
        END IF

    END SUBROUTINE quick_sort_1_int_int


    SUBROUTINE interchange_sort_int_int(left_end, right_end)

        INTEGER, INTENT(IN) :: left_end, right_end

        ! Local variables
        integer             :: i, j, itemp
        integer(i8b)        :: temp

        DO i = left_end, right_end - 1
            DO j = i + 1, right_end
                IF (list(i) > list(j)) THEN
                    temp = list(i); list(i) = list(j); list(j) = temp
                    itemp = order(i); order(i) = order(j); order(j) = itemp
                END IF
            END DO
        END DO

    END SUBROUTINE interchange_sort_int_int

END SUBROUTINE quick_sort_int_int
