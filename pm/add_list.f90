!################################################################
!################################################################
!################################################################
!################################################################
subroutine add_list(ind_part,ind_grid,ok,np)
  use amr_commons, only:nvector
  implicit none
  integer, intent(in)::np
  integer,dimension(1:nvector), intent(in)::ind_part,ind_grid
  logical,dimension(1:nvector), intent(in)::ok
  !----------------------------------------------------
  ! Add particles to their new linked lists
  !----------------------------------------------------
  integer::j

!$omp critical(omp_particle_list)
  do j=1,np
     if(ok(j)) call add_list_one(ind_part(j),ind_grid(j))
  end do
!$omp end critical(omp_particle_list)

end subroutine add_list
!################################################################
!################################################################
!################################################################
!################################################################
subroutine add_list_one(ipart,igrid)
  use pm_commons
  implicit none
  integer, intent(in)::ipart,igrid
  !----------------------------------------------------
  ! Append one particle to the tail of a grid linked list.
  ! OMP: Contains no synchronisation of its own: callers must either hold the
  ! particle-list critical section (add_list) or be running serially
  ! (apply_tree_moves).
  !----------------------------------------------------

  if (numbp(igrid) > 0) then
     ! Add particle at the tail of its linked list
     nextp(tailp(igrid)) = ipart
     prevp(ipart) = tailp(igrid)
     nextp(ipart) = 0
     tailp(igrid) = ipart
     numbp(igrid) = numbp(igrid) + 1
  else
     ! Initialise linked list
     headp(igrid) = ipart
     tailp(igrid) = ipart
     prevp(ipart) = 0
     nextp(ipart) = 0
     numbp(igrid) = 1
  end if

end subroutine add_list_one
!################################################################
!################################################################
!################################################################
!################################################################
subroutine add_free(ind_part,np)
  use amr_commons
  use pm_commons
  use dice_commons
  implicit none
  integer, intent(in)::np
  integer,dimension(1:nvector), intent(in)::ind_part
  !----------------------------------------------------
  ! Add particles to the free memory linked list
  ! and reset all particle variables
  !----------------------------------------------------
  integer::j,idim

  do idim=1,ndim
     do j=1,np
        xp(ind_part(j),idim)=0
        vp(ind_part(j),idim)=0
     end do
  end do
  do j=1,np
     mp(ind_part(j))=0
     idp(ind_part(j))=0
     levelp(ind_part(j))=0
     typep(ind_part(j))%family=FAM_UNDEF
     typep(ind_part(j))%tag=0
  end do
  if(star.or.sink)then
     do j=1,np
        tp(ind_part(j))=0
     end do
     if(metal)then
        do j=1,np
           zp(ind_part(j))=0
        end do
     end if
  end if
  if(dice_init) then
     do j=1,np
        up(ind_part(j))=0.0
     end do
  endif

!$omp critical(omp_particle_list)
  do j=1,np
     call add_free_one(ind_part(j))
  end do
!$omp end critical(omp_particle_list)

end subroutine add_free
!################################################################
!################################################################
!################################################################
!################################################################
subroutine add_free_one(ipart)
  use pm_commons
  implicit none
  integer, intent(in)::ipart
  !----------------------------------------------------
  ! Append one particle to the tail of the free memory linked list.
  ! Splice only: the caller is responsible for resetting the particle data.
  ! OMP: Contains no synchronisation of its own, callers must hold the
  ! particle-list critical section (add_free).
  !----------------------------------------------------

  if(numbp_free>0)then
     ! Add particle at the tail of its linked list
     nextp(tailp_free)=ipart
     prevp(ipart)=tailp_free
     nextp(ipart)=0
     tailp_free=ipart
     numbp_free=numbp_free+1
  else
     ! Initialise linked list
     headp_free=ipart
     tailp_free=ipart
     prevp(ipart)=0
     nextp(ipart)=0
     numbp_free=1
  end if
  npart=npartmax-numbp_free

end subroutine add_free_one
!################################################################
!################################################################
!################################################################
!################################################################
subroutine add_free_cond(ind_part,ok,np)
  use amr_commons
  use pm_commons
  use dice_commons
  implicit none
  integer::np
  integer,dimension(1:nvector)::ind_part
  logical,dimension(1:nvector)::ok
  !----------------------------------------------------
  ! Add particles to the free memory linked list
  ! and reset all particle variables
  !----------------------------------------------------
  integer::j,idim

  do idim=1,ndim
     do j=1,np
        if(ok(j))then
           xp(ind_part(j),idim)=0
           vp(ind_part(j),idim)=0
        endif
     end do
  end do
  do j=1,np
     if(ok(j))then
        mp(ind_part(j))=0
        idp(ind_part(j))=0
        levelp(ind_part(j))=0
        typep(ind_part(j))%family = FAM_UNDEF
        typep(ind_part(j))%tag = 0
     endif
  end do
  if(star.or.sink)then
     do j=1,np
        if(ok(j))then
           tp(ind_part(j))=0
        endif
     end do
     if(metal)then
        do j=1,np
           if(ok(j))then
              zp(ind_part(j))=0
           endif
        end do
     end if
  end if
  if(dice_init) then
     do j=1,np
        if(ok(j))then
           up(ind_part(j))=0.0
        endif
     end do
  endif

!$omp critical(omp_particle_list)
  do j=1,np
     if(ok(j)) call add_free_one(ind_part(j))
  end do
!$omp end critical(omp_particle_list)

end subroutine add_free_cond
