!################################################################
!################################################################
!################################################################
!################################################################
subroutine remove_list(ind_part,ind_grid,ok,np)
  use amr_commons, only:nvector
  implicit none
  integer, intent(in)::np
  integer,dimension(1:nvector), intent(in)::ind_part,ind_grid
  logical,dimension(1:nvector), intent(in)::ok
  !----------------------------------------------------
  ! Remove particles from their original linked lists
  !----------------------------------------------------
  integer::j

!$omp critical(omp_particle_list)
  do j=1,np
     if(ok(j)) call remove_list_one(ind_part(j),ind_grid(j))
  end do
!$omp end critical(omp_particle_list)

end subroutine remove_list
!################################################################
!################################################################
!################################################################
!################################################################
subroutine remove_list_one(ipart,igrid)
  use pm_commons
  implicit none
  integer, intent(in)::ipart,igrid
  !----------------------------------------------------
  ! Unlink one particle from its grid linked list.
  ! OMP: Contains no synchronisation of its own: callers must either hold the
  ! particle-list critical section (remove_list) or be running serially
  ! (apply_tree_moves).
  !----------------------------------------------------

  if(prevp(ipart) .ne. 0) then
     if( nextp(ipart) .ne. 0 )then
        nextp(prevp(ipart))=nextp(ipart)
        prevp(nextp(ipart))=prevp(ipart)
     else
        nextp(prevp(ipart))=0
        tailp(igrid)=prevp(ipart)
     end if
  else
     if(nextp(ipart) .ne. 0)then
        prevp(nextp(ipart))=0
        headp(igrid)=nextp(ipart)
     else
        headp(igrid)=0
        tailp(igrid)=0
     end if
  end if
  numbp(igrid)=numbp(igrid)-1

end subroutine remove_list_one
!################################################################
!################################################################
!################################################################
!################################################################
subroutine remove_free(ind_part,np)
  use amr_commons
  use pm_commons
  implicit none
  integer, intent(in) :: np
  integer, dimension(1:nvector), intent(out)::ind_part
  !-----------------------------------------------
  ! Get np particle from free memory linked list
  !-----------------------------------------------
  integer::j,ipart

!$omp critical(omp_particle_list)
  do j=1,np
     ipart=headp_free
     ind_part(j)=ipart
     numbp_free=numbp_free-1
     if(numbp_free<0)then
        write(*,*)'No more free memory'
        write(*,*)'in PE ',myid
        write(*,*)'Increase npartmax'
        call clean_stop
     end if
     headp_free=nextp(headp_free)
  end do
  npart=npartmax-numbp_free
!$omp end critical(omp_particle_list)

end subroutine remove_free
