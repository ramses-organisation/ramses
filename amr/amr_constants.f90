module amr_constants
  use amr_parameters, only:twotondim,threetondim
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

#if NDIM==1
  integer,dimension(1:threetondim,1:twotondim),parameter::lll=reshape( (/ &
                                                                   2,1,1, &  !ind=1
                                                                   1,1,2  &  !ind=2
                                               /), shape=(/threetondim,twotondim/))

  integer,dimension(1:threetondim,1:twotondim),parameter::mmm=reshape( (/ &
                                                                   2,1,2, &  !ind=1
                                                                   1,2,1  &  !ind=2
                                               /), shape=(/threetondim,twotondim/))

#elif NDIM==2
  integer,dimension(1:threetondim,1:twotondim),parameter::lll=reshape( (/ &
                                                       4,3,3,2,1,1,2,1,1, &  !ind=1
                                                       3,3,4,1,1,2,1,1,2, &  !ind=2
                                                       2,1,1,2,1,1,4,3,3, &  !ind=3
                                                       1,1,2,1,1,2,3,3,4  &  !ind=4
                                               /), shape=(/threetondim,twotondim/))

  integer,dimension(1:threetondim,1:twotondim),parameter::mmm=reshape( (/ &
                                                       4,3,4,2,1,2,4,3,4, &  !ind=1
                                                       3,4,3,1,2,1,3,4,3, &  !ind=2
                                                       2,1,2,4,3,4,2,1,2, &  !ind=3
                                                       1,2,1,3,4,3,1,2,1  &  !ind=4
                                               /), shape=(/threetondim,twotondim/))

#elif NDIM==3
  integer,dimension(1:threetondim,1:twotondim),parameter::lll=reshape( (/ &
                   8,7,7,6,5,5,6,5,5,4,3,3,2,1,1,2,1,1,4,3,3,2,1,1,2,1,1, &  !ind=1
                   7,7,8,5,5,6,5,5,6,3,3,4,1,1,2,1,1,2,3,3,4,1,1,2,1,1,2, &  !ind=2
                   6,5,5,6,5,5,8,7,7,2,1,1,2,1,1,4,3,3,2,1,1,2,1,1,4,3,3, &  !ind=3
                   5,5,6,5,5,6,7,7,8,1,1,2,1,1,2,3,3,4,1,1,2,1,1,2,3,3,4, &  !ind=4
                   4,3,3,2,1,1,2,1,1,4,3,3,2,1,1,2,1,1,8,7,7,6,5,5,6,5,5, &  !ind=5
                   3,3,4,1,1,2,1,1,2,3,3,4,1,1,2,1,1,2,7,7,8,5,5,6,5,5,6, &  !ind=6
                   2,1,1,2,1,1,4,3,3,2,1,1,2,1,1,4,3,3,6,5,5,6,5,5,8,7,7, &  !ind=7
                   1,1,2,1,1,2,3,3,4,1,1,2,1,1,2,3,3,4,5,5,6,5,5,6,7,7,8  &  !ind=8
                                               /), shape=(/threetondim,twotondim/))

  integer,dimension(1:threetondim,1:twotondim),parameter::mmm=reshape( (/ &
                   8,7,8,6,5,6,8,7,8,4,3,4,2,1,2,4,3,4,8,7,8,6,5,6,8,7,8, &  !ind=1
                   7,8,7,5,6,5,7,8,7,3,4,3,1,2,1,3,4,3,7,8,7,5,6,5,7,8,7, &  !ind=2
                   6,5,6,8,7,8,6,5,6,2,1,2,4,3,4,2,1,2,6,5,6,8,7,8,6,5,6, &  !ind=3
                   5,6,5,7,8,7,5,6,5,1,2,1,3,4,3,1,2,1,5,6,5,7,8,7,5,6,5, &  !ind=4
                   4,3,4,2,1,2,4,3,4,8,7,8,6,5,6,8,7,8,4,3,4,2,1,2,4,3,4, &  !ind=5
                   3,4,3,1,2,1,3,4,3,7,8,7,5,6,5,7,8,7,3,4,3,1,2,1,3,4,3, &  !ind=6
                   2,1,2,4,3,4,2,1,2,6,5,6,8,7,8,6,5,6,2,1,2,4,3,4,2,1,2, &  !ind=7
                   1,2,1,3,4,3,1,2,1,5,6,5,7,8,7,5,6,5,1,2,1,3,4,3,1,2,1  &  !ind=8
                                               /), shape=(/threetondim,twotondim/))
#endif

  ! ---- Stencil constants ----
  ! in 3D:
  !            min      max    tot
  ! -----------------------------------------
  ! i,j,k1      0        2      3
  ! i,j,k2      0        1      2
  ! i,j,k3      1        2      2

  integer,parameter::i1min=0,i2min=0,i3min=1
  integer,parameter::j1min=0,j2min=0,j3min=1
  integer,parameter::k1min=0,k2min=0,k3min=1

  integer,parameter::i1max=2,i2max=1,i3max=2
#if NDIM==1
  integer,parameter::j1max=0,j2max=0,j3max=1
#else
  integer,parameter::j1max=2,j2max=1,j3max=2
#endif
#if NDIM==1 || NDIM==2
  integer,parameter::k1max=0,k2max=0,k3max=1
#else
  integer,parameter::k1max=2,k2max=1,k3max=2
#endif

end module amr_constants
