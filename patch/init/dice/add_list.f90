!################################################################
!################################################################
!################################################################
!################################################################
subroutine add_list(ind_part, ind_grid, ok, np)
    use amr_commons
    use pm_commons
    implicit none
    integer, intent(in) :: np
    integer, dimension(1:nvector), intent(in) :: ind_part, ind_grid
    logical, dimension(1:nvector), intent(in) :: ok
    !
    ! Add particles to their new linked lists
    !
    integer :: j

    do j = 1, np
        if (ok(j)) then
            if (numbp(ind_grid(j)) > 0) then
                ! Add particle at the tail of its linked list
                nextp(tailp(ind_grid(j))) = ind_part(j)
                prevp(ind_part(j)) = tailp(ind_grid(j))
                nextp(ind_part(j)) = 0
                tailp(ind_grid(j)) = ind_part(j)
                numbp(ind_grid(j)) = numbp(ind_grid(j)) + 1
            else
                ! Initialise linked list
                headp(ind_grid(j)) = ind_part(j)
                tailp(ind_grid(j)) = ind_part(j)
                prevp(ind_part(j)) = 0
                nextp(ind_part(j)) = 0
                numbp(ind_grid(j)) = 1
            end if
        end if
    end do

end subroutine add_list
!################################################################
!################################################################
!################################################################
!################################################################
subroutine add_free(ind_part, np)
    use amr_commons
    use pm_commons
    use dice_commons
    implicit none
    integer, intent(in) :: np
    integer, dimension(1:nvector), intent(in) :: ind_part
    !
    ! Add particles to the free memory linked list
    ! and reset all particle variables
    !
    integer :: j, idim

    do idim = 1, ndim
        do j = 1, np
            xp(ind_part(j), idim) = 0
            vp(ind_part(j), idim) = 0
        end do
    end do
    do j = 1, np
        mp(ind_part(j)) = 0
        idp(ind_part(j)) = 0
        levelp(ind_part(j)) = 0
        typep(ind_part(j))%family = FAM_UNDEF
        typep(ind_part(j))%tag = 0
    end do
    if (star .or. sink) then
        do j = 1, np
            tp(ind_part(j)) = 0
        end do
        if (metal) then
            do j = 1, np
                zp(ind_part(j)) = 0
            end do
        end if
    end if
    ! DICE patch
    if (dice_init) then
        do j = 1, np
            up(ind_part(j)) = 0.0
        end do
    end if

    do j = 1, np
        if (numbp_free > 0) then
            ! Add particle at the tail of its linked list
            nextp(tailp_free) = ind_part(j)
            prevp(ind_part(j)) = tailp_free
            nextp(ind_part(j)) = 0
            tailp_free = ind_part(j)
            numbp_free = numbp_free + 1
        else
            ! Initialise linked list
            headp_free = ind_part(j)
            tailp_free = ind_part(j)
            prevp(ind_part(j)) = 0
            nextp(ind_part(j)) = 0
            numbp_free = 1
        end if
    end do
    npart = npartmax - numbp_free

end subroutine add_free
!################################################################
!################################################################
!################################################################
!################################################################
subroutine add_free_cond(ind_part, ok, np)
    use amr_commons
    use pm_commons
    use dice_commons
    implicit none
    integer :: np
    integer, dimension(1:nvector) :: ind_part
    logical, dimension(1:nvector) :: ok
    !
    ! Add particles to the free memory linked list
    ! and reset all particle variables
    !
    integer :: j, idim

    do idim = 1, ndim
        do j = 1, np
            if (ok(j)) then
                xp(ind_part(j), idim) = 0
                vp(ind_part(j), idim) = 0
            end if
        end do
    end do
    do j = 1, np
        if (ok(j)) then
            mp(ind_part(j)) = 0
            idp(ind_part(j)) = 0
            levelp(ind_part(j)) = 0
            typep(ind_part(j))%family = FAM_UNDEF
            typep(ind_part(j))%tag = 0
        end if
    end do
    if (star .or. sink) then
        do j = 1, np
            if (ok(j)) then
                tp(ind_part(j)) = 0
            end if
        end do
        if (metal) then
            do j = 1, np
                if (ok(j)) then
                    zp(ind_part(j)) = 0
                end if
            end do
        end if
    end if
    ! DICE patch
    if (dice_init) then
        do j = 1, np
            if (ok(j)) then
                up(ind_part(j)) = 0.0
            end if
        end do
    end if

    do j = 1, np
        if (ok(j)) then
            if (numbp_free > 0) then
                ! Add particle at the tail of its linked list
                nextp(tailp_free) = ind_part(j)
                prevp(ind_part(j)) = tailp_free
                nextp(ind_part(j)) = 0
                tailp_free = ind_part(j)
                numbp_free = numbp_free + 1
            else
                ! Initialise linked list
                headp_free = ind_part(j)
                tailp_free = ind_part(j)
                prevp(ind_part(j)) = 0
                nextp(ind_part(j)) = 0
                numbp_free = 1
            end if
        end if
    end do
    npart = npartmax - numbp_free

end subroutine add_free_cond
