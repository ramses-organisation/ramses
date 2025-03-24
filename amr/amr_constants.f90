module amr_constants
  implicit none
  integer, dimension(1:3,1:2,1:8),parameter :: iii=reshape( (/ &
                                              1,0,1,0,1,0,1,0, &
                                              0,2,0,2,0,2,0,2, &
                                              3,3,0,0,3,3,0,0, &
                                              0,0,4,4,0,0,4,4, &
                                              5,5,5,5,0,0,0,0, &
                                              0,0,0,0,6,6,6,6  &
                                       /), shape=(/3,2,8/), order=(/3,2,1/))
  integer, dimension(1:3,1:2,1:8),parameter :: jjj=reshape( (/ &
                                              2,1,4,3,6,5,8,7, &
                                              2,1,4,3,6,5,8,7, &
                                              3,4,1,2,7,8,5,6, &
                                              3,4,1,2,7,8,5,6, &
                                              5,6,7,8,1,2,3,4, &
                                              5,6,7,8,1,2,3,4  &
                                       /), shape=(/3,2,8/), order=(/3,2,1/))
end module amr_constants
