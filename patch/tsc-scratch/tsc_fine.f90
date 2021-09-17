! use amr_commons
! use pm_commons
! use poisson_commons
! use hydro_commons, ONLY: uold,unew,smallr,nvar,gamma
! implicit none
! integer::ng,np,ilevel,xtondim
! integer,dimension(1:nvector)::ind_grid
! integer,dimension(1:nvector)::ind_grid_part,ind_part
! !
! !
! !
! logical::error
! integer::i,j,ind,idim,nx_loc,isink,ivar_dust
! real(dp)::dx,scale,dx_loc,vol_loc
!
!
! ! Grid-based arrays
! real(dp),dimension(1:nvector,1:ndim),save::x0
! integer ,dimension(1:nvector),save::ind_cell
! integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
! integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
! ! Particle-based arrays
! logical ,dimension(1:nvector),save::ok
! real(dp),dimension(1:nvector,1:ndim),save::x,ff,new_vp
! real(dp),dimension(1:nvector,1:ndim),save::cl,cr,cc,wl,wr,wc
! integer ,dimension(1:nvector,1:ndim),save::igl,igr,igc,icl,icr,icc
! real(dp),dimension(1:nvector,1:threetondim),save::vol
! integer ,dimension(1:nvector,1:threetondim),save::igrid,icell,indp,kg
! integer::i,j,idim,ind,np,ng
if (ndim .ne. 3)then
   write(*,*)'TSC not supported for ndim neq 3'
   call clean_stop
end if

xtondim=threetondim

! Check for illegal moves
error=.false.
do idim=1,ndim
   do j=1,np
      if(x(j,idim)<1.0D0.or.x(j,idim)>5.0D0)error=.true.
   end do
end do
if(error)then
   write(*,*)'problem in tsc_fine'
   do idim=1,ndim
      do j=1,np
         if(x(j,idim)<1.0D0.or.x(j,idim)>5.0D0)then
            write(*,*)x(j,1:ndim)
         endif
      end do
   end do
   stop
end if

! TSC at level ilevel; a particle contributes
!     to three cells in each dimension
! cl: position of leftmost cell centre
! cc: position of central cell centre
! cr: position of rightmost cell centre
! wl: weighting function for leftmost cell
! wc: weighting function for central cell
! wr: weighting function for rightmost cell
do idim=1,ndim
   do j=1,np
      cl(j,idim)=dble(int(x(j,idim)))-0.5D0
      cc(j,idim)=dble(int(x(j,idim)))+0.5D0
      cr(j,idim)=dble(int(x(j,idim)))+1.5D0
      wl(j,idim)=0.50D0*(1.5D0-abs(x(j,idim)-cl(j,idim)))**2
      wc(j,idim)=0.75D0-          (x(j,idim)-cc(j,idim)) **2
      wr(j,idim)=0.50D0*(1.5D0-abs(x(j,idim)-cr(j,idim)))**2
   end do
end do

! Compute parent grids
do idim=1,ndim
   do j=1,np
      igl(j,idim)=(int(cl(j,idim)))/2
      igc(j,idim)=(int(cc(j,idim)))/2
      igr(j,idim)=(int(cr(j,idim)))/2
   end do
end do
! #if NDIM==3
do j=1,np
   kg(j,1 )=1+igl(j,1)+3*igl(j,2)+9*igl(j,3)
   kg(j,2 )=1+igc(j,1)+3*igl(j,2)+9*igl(j,3)
   kg(j,3 )=1+igr(j,1)+3*igl(j,2)+9*igl(j,3)
   kg(j,4 )=1+igl(j,1)+3*igc(j,2)+9*igl(j,3)
   kg(j,5 )=1+igc(j,1)+3*igc(j,2)+9*igl(j,3)
   kg(j,6 )=1+igr(j,1)+3*igc(j,2)+9*igl(j,3)
   kg(j,7 )=1+igl(j,1)+3*igr(j,2)+9*igl(j,3)
   kg(j,8 )=1+igc(j,1)+3*igr(j,2)+9*igl(j,3)
   kg(j,9 )=1+igr(j,1)+3*igr(j,2)+9*igl(j,3)
   kg(j,10)=1+igl(j,1)+3*igl(j,2)+9*igc(j,3)
   kg(j,11)=1+igc(j,1)+3*igl(j,2)+9*igc(j,3)
   kg(j,12)=1+igr(j,1)+3*igl(j,2)+9*igc(j,3)
   kg(j,13)=1+igl(j,1)+3*igc(j,2)+9*igc(j,3)
   kg(j,14)=1+igc(j,1)+3*igc(j,2)+9*igc(j,3)
   kg(j,15)=1+igr(j,1)+3*igc(j,2)+9*igc(j,3)
   kg(j,16)=1+igl(j,1)+3*igr(j,2)+9*igc(j,3)
   kg(j,17)=1+igc(j,1)+3*igr(j,2)+9*igc(j,3)
   kg(j,18)=1+igr(j,1)+3*igr(j,2)+9*igc(j,3)
   kg(j,19)=1+igl(j,1)+3*igl(j,2)+9*igr(j,3)
   kg(j,20)=1+igc(j,1)+3*igl(j,2)+9*igr(j,3)
   kg(j,21)=1+igr(j,1)+3*igl(j,2)+9*igr(j,3)
   kg(j,22)=1+igl(j,1)+3*igc(j,2)+9*igr(j,3)
   kg(j,23)=1+igc(j,1)+3*igc(j,2)+9*igr(j,3)
   kg(j,24)=1+igr(j,1)+3*igc(j,2)+9*igr(j,3)
   kg(j,25)=1+igl(j,1)+3*igr(j,2)+9*igr(j,3)
   kg(j,26)=1+igc(j,1)+3*igr(j,2)+9*igr(j,3)
   kg(j,27)=1+igr(j,1)+3*igr(j,2)+9*igr(j,3)
end do

do ind=1,threetondim
   do j=1,np
      igrid(j,ind)=son(nbors_father_cells(ind_grid_part(j),kg(j,ind)))
   end do
end do

! Check if particles are entirely in level ilevel
ok(1:np)=.true.
do ind=1,threetondim
   do j=1,np
      ok(j)=ok(j).and.igrid(j,ind)>0
   end do
end do

! If not, rescale position at level ilevel-1
do idim=1,ndim
   do j=1,np
      if(.not.ok(j))then
         x(j,idim)=x(j,idim)/2.0D0
      end if
   end do
end do
! If not, redo TSC at level ilevel-1
do idim=1,ndim
   do j=1,np
      if(.not.ok(j))then
        cl(j,idim)=dble(int(x(j,idim)))-0.5D0
        cc(j,idim)=dble(int(x(j,idim)))+0.5D0
        cr(j,idim)=dble(int(x(j,idim)))+1.5D0
        wl(j,idim)=0.50D0*(1.5D0-abs(x(j,idim)-cl(j,idim)))**2
        wc(j,idim)=0.75D0-          (x(j,idim)-cc(j,idim)) **2
        wr(j,idim)=0.50D0*(1.5D0-abs(x(j,idim)-cr(j,idim)))**2
      end if
   end do
end do

! Compute parent cell position
do idim=1,ndim
   do j=1,np
     if(ok(j))then
      icl(j,idim)=int(cl(j,idim))-2*igl(j,idim)
      icc(j,idim)=int(cc(j,idim))-2*igc(j,idim)
      icr(j,idim)=int(cr(j,idim))-2*igr(j,idim)
     else ! ERM: this else may or may not be correct? But I believe it is.
      icl(j,idim)=int(cl(j,idim))
      icc(j,idim)=int(cc(j,idim))
      icr(j,idim)=int(cr(j,idim))
     endif
   end do
end do

! #if NDIM==3
do j=1,np
  if(ok(j))then
   icell(j,1 )=1+icl(j,1)+2*icl(j,2)+4*icl(j,3)
   icell(j,2 )=1+icc(j,1)+2*icl(j,2)+4*icl(j,3)
   icell(j,3 )=1+icr(j,1)+2*icl(j,2)+4*icl(j,3)
   icell(j,4 )=1+icl(j,1)+2*icc(j,2)+4*icl(j,3)
   icell(j,5 )=1+icc(j,1)+2*icc(j,2)+4*icl(j,3)
   icell(j,6 )=1+icr(j,1)+2*icc(j,2)+4*icl(j,3)
   icell(j,7 )=1+icl(j,1)+2*icr(j,2)+4*icl(j,3)
   icell(j,8 )=1+icc(j,1)+2*icr(j,2)+4*icl(j,3)
   icell(j,9 )=1+icr(j,1)+2*icr(j,2)+4*icl(j,3)
   icell(j,10)=1+icl(j,1)+2*icl(j,2)+4*icc(j,3)
   icell(j,11)=1+icc(j,1)+2*icl(j,2)+4*icc(j,3)
   icell(j,12)=1+icr(j,1)+2*icl(j,2)+4*icc(j,3)
   icell(j,13)=1+icl(j,1)+2*icc(j,2)+4*icc(j,3)
   icell(j,14)=1+icc(j,1)+2*icc(j,2)+4*icc(j,3)
   icell(j,15)=1+icr(j,1)+2*icc(j,2)+4*icc(j,3)
   icell(j,16)=1+icl(j,1)+2*icr(j,2)+4*icc(j,3)
   icell(j,17)=1+icc(j,1)+2*icr(j,2)+4*icc(j,3)
   icell(j,18)=1+icr(j,1)+2*icr(j,2)+4*icc(j,3)
   icell(j,19)=1+icl(j,1)+2*icl(j,2)+4*icr(j,3)
   icell(j,20)=1+icc(j,1)+2*icl(j,2)+4*icr(j,3)
   icell(j,21)=1+icr(j,1)+2*icl(j,2)+4*icr(j,3)
   icell(j,22)=1+icl(j,1)+2*icc(j,2)+4*icr(j,3)
   icell(j,23)=1+icc(j,1)+2*icc(j,2)+4*icr(j,3)
   icell(j,24)=1+icr(j,1)+2*icc(j,2)+4*icr(j,3)
   icell(j,25)=1+icl(j,1)+2*icr(j,2)+4*icr(j,3)
   icell(j,26)=1+icc(j,1)+2*icr(j,2)+4*icr(j,3)
   icell(j,27)=1+icr(j,1)+2*icr(j,2)+4*icr(j,3)
 else
   icell(j,1 )=1+icl(j,1)+3*icl(j,2)+9*icl(j,3)
   icell(j,2 )=1+icc(j,1)+3*icl(j,2)+9*icl(j,3)
   icell(j,3 )=1+icr(j,1)+3*icl(j,2)+9*icl(j,3)
   icell(j,4 )=1+icl(j,1)+3*icc(j,2)+9*icl(j,3)
   icell(j,5 )=1+icc(j,1)+3*icc(j,2)+9*icl(j,3)
   icell(j,6 )=1+icr(j,1)+3*icc(j,2)+9*icl(j,3)
   icell(j,7 )=1+icl(j,1)+3*icr(j,2)+9*icl(j,3)
   icell(j,8 )=1+icc(j,1)+3*icr(j,2)+9*icl(j,3)
   icell(j,9 )=1+icr(j,1)+3*icr(j,2)+9*icl(j,3)
   icell(j,10)=1+icl(j,1)+3*icl(j,2)+9*icc(j,3)
   icell(j,11)=1+icc(j,1)+3*icl(j,2)+9*icc(j,3)
   icell(j,12)=1+icr(j,1)+3*icl(j,2)+9*icc(j,3)
   icell(j,13)=1+icl(j,1)+3*icc(j,2)+9*icc(j,3)
   icell(j,14)=1+icc(j,1)+3*icc(j,2)+9*icc(j,3)
   icell(j,15)=1+icr(j,1)+3*icc(j,2)+9*icc(j,3)
   icell(j,16)=1+icl(j,1)+3*icr(j,2)+9*icc(j,3)
   icell(j,17)=1+icc(j,1)+3*icr(j,2)+9*icc(j,3)
   icell(j,18)=1+icr(j,1)+3*icr(j,2)+9*icc(j,3)
   icell(j,19)=1+icl(j,1)+3*icl(j,2)+9*icr(j,3)
   icell(j,20)=1+icc(j,1)+3*icl(j,2)+9*icr(j,3)
   icell(j,21)=1+icr(j,1)+3*icl(j,2)+9*icr(j,3)
   icell(j,22)=1+icl(j,1)+3*icc(j,2)+9*icr(j,3)
   icell(j,23)=1+icc(j,1)+3*icc(j,2)+9*icr(j,3)
   icell(j,24)=1+icr(j,1)+3*icc(j,2)+9*icr(j,3)
   icell(j,25)=1+icl(j,1)+3*icr(j,2)+9*icr(j,3)
   icell(j,26)=1+icc(j,1)+3*icr(j,2)+9*icr(j,3)
   icell(j,27)=1+icr(j,1)+3*icr(j,2)+9*icr(j,3)
 endif
end do

! Compute parent cell adress
do ind=1,threetondim
   do j=1,np
     if(ok(j))then
      indp(j,ind)=ncoarse+(icell(j,ind)-1)*ngridmax+igrid(j,ind)
     else ! ERM: for AMR(?) there may be an issue with ind_grid_part(j) being used here.
       indp(j,ind)=nbors_father_cells(ind_grid_part(j),icell(j,ind))
     endif
   end do
end do

! Compute cloud volumes (NDIM==3)
do j=1,np
   vol(j,1 )=wl(j,1)*wl(j,2)*wl(j,3)
   vol(j,2 )=wc(j,1)*wl(j,2)*wl(j,3)
   vol(j,3 )=wr(j,1)*wl(j,2)*wl(j,3)
   vol(j,4 )=wl(j,1)*wc(j,2)*wl(j,3)
   vol(j,5 )=wc(j,1)*wc(j,2)*wl(j,3)
   vol(j,6 )=wr(j,1)*wc(j,2)*wl(j,3)
   vol(j,7 )=wl(j,1)*wr(j,2)*wl(j,3)
   vol(j,8 )=wc(j,1)*wr(j,2)*wl(j,3)
   vol(j,9 )=wr(j,1)*wr(j,2)*wl(j,3)
   vol(j,10)=wl(j,1)*wl(j,2)*wc(j,3)
   vol(j,11)=wc(j,1)*wl(j,2)*wc(j,3)
   vol(j,12)=wr(j,1)*wl(j,2)*wc(j,3)
   vol(j,13)=wl(j,1)*wc(j,2)*wc(j,3)
   vol(j,14)=wc(j,1)*wc(j,2)*wc(j,3)
   vol(j,15)=wr(j,1)*wc(j,2)*wc(j,3)
   vol(j,16)=wl(j,1)*wr(j,2)*wc(j,3)
   vol(j,17)=wc(j,1)*wr(j,2)*wc(j,3)
   vol(j,18)=wr(j,1)*wr(j,2)*wc(j,3)
   vol(j,19)=wl(j,1)*wl(j,2)*wr(j,3)
   vol(j,20)=wc(j,1)*wl(j,2)*wr(j,3)
   vol(j,21)=wr(j,1)*wl(j,2)*wr(j,3)
   vol(j,22)=wl(j,1)*wc(j,2)*wr(j,3)
   vol(j,23)=wc(j,1)*wc(j,2)*wr(j,3)
   vol(j,24)=wr(j,1)*wc(j,2)*wr(j,3)
   vol(j,25)=wl(j,1)*wr(j,2)*wr(j,3)
   vol(j,26)=wc(j,1)*wr(j,2)*wr(j,3)
   vol(j,27)=wr(j,1)*wr(j,2)*wr(j,3)
end do
