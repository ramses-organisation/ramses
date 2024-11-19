#if (EOS == 0)
  include 'eos_sf.f90'
#elif (EOS == 1)
  include 'eos_mg.f90'
#elif (EOS == 2)
  include 'eos_cc.f90'
#endif

!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmpdt(uu,grav,rr,dx,dt,ncell)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none
  integer::ncell
  real(dp)::dx,dt
  real(dp),dimension(1:nvector,1:nvar)::uu
  real(dp),dimension(1:nvector,1:ndim)::grav
  real(dp),dimension(1:nvector)::rr

  real(dp),dimension(1:nvector,1:npri),save::qq
  real(dp),dimension(1:nvector,1:nmat),save::ff,gg,kappa_matt
  real(dp),dimension(1:nvector),save::ekin,dtot,cc,st,pp,kappa_hatt
  real(dp)::dtcell,eps
  integer::k,idim,imat

  ! Convert to primitive variable

  ! Volume fraction and fluid density
  do imat = 1,nmat
     do k = 1,ncell
        ff(k,imat) = uu(k,imat+npri)
        gg(k,imat) = uu(k,imat+npri+nmat)
     end do
  end do

  ! Compute density
  do k = 1,ncell
     qq(k,1) = uu(k,1)
  end do

  ! Compute velocity and specific kinetic energy
  ekin(1:ncell)=0.0
  do idim = 1,ndim
     do k = 1,ncell
        qq(k,idim+1) = uu(k,idim+1)/uu(k,1)
        ekin(k) = ekin(k) + half*qq(k,idim+1)**2
     end do
  end do

  ! Compute total internal energy
  do k = 1,ncell
     qq(k,npri) = uu(k,npri) - uu(k,1)*ekin(k)
  end do

  ! Call eos routine
  call eos(ff,gg,qq,pp,cc,kappa_matt,kappa_hatt,ncell)

  ! Compute wave speed
  if(geom==3)then
     do k = 1,ncell
        eps = dx/two/rr(k)
        cc(k) = (abs(qq(k,2))+cc(k))*(one+eps)**2/(one+third*eps**2)
     end do
  else if(geom==2)then
     do k = 1,ncell
        eps = dx/two/rr(k)
        cc(k) = (abs(qq(k,2))+cc(k))*(one+eps)
     end do
  else
     do k = 1,ncell
        cc(k) = abs(qq(k,2))+cc(k)
     end do
  endif
  do idim = 2,ndim
     do k = 1,ncell
        cc(k) = cc(k) + abs(qq(k,idim+1))+cc(k)
     end do
  end do

  ! Compute gravity strength ratio
  do k = 1,ncell
     st(k) = zero
  end do
  do idim = 1,ndim
     do k = 1,ncell
        st(k) = st(k) + abs(grav(k,idim))
     end do
  end do
  do k = 1,ncell
     st(k) = st(k)*dx/cc(k)**2
     st(k) = MAX(st(k),0.0001_dp)
  end do

  ! Compute maximum time step for each authorized cell
  dt = courant_factor*dx/smallc
  do k = 1,ncell
     dtcell = dx/cc(k)*(sqrt(one+two*courant_factor*st(k))-one)/st(k)
     dt = min(dt,dtcell)
  end do

end subroutine cmpdt
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine hydro_refine(ug,um,ud,ok,current_dim,ncell)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none
  ! dummy arguments
  integer::ncell,current_dim
  real(dp),dimension(1:nvector,1:nvar)::ug,um,ud
  logical ,dimension(1:nvector)       ::ok

  integer::k,idim,imat
  real(dp),dimension(1:nvector,1:npri),save::qg,qm,qd
  real(dp),dimension(1:nvector,1:nmat),save::fg,fm,fd
  real(dp),dimension(1:nvector,1:nmat),save::gg,gm,gd
  real(dp),dimension(1:nvector,1:nmat),save::kappa_matg,kappa_matm,kappa_matd
  real(dp),dimension(1:nvector),save::eking,ekinm,ekind
  real(dp),dimension(1:nvector),save::pg,pm,pd
  real(dp),dimension(1:nvector),save::cg,cm,cd
  real(dp),dimension(1:nvector),save::kappa_hatg,kappa_hatm,kappa_hatd
  logical ,dimension(1:nvector),save::wg,wd,bg,bm,bd
  real(dp)::ffg,ffm,ffd,ddg,ddm,ddd
  real(dp)::ppg,ppm,ppd,vvg,vvm,vvd
  real(dp)::ccg,ccm,ccd,error

  ! Convert to primitive variables

  ! Volume fraction and fluid density
  do imat = 1,nmat
     do k = 1,ncell
        fg(k,imat) = ug(k,imat+npri)
        fm(k,imat) = um(k,imat+npri)
        fd(k,imat) = ud(k,imat+npri)
        gg(k,imat) = ug(k,imat+npri+nmat)
        gm(k,imat) = um(k,imat+npri+nmat)
        gd(k,imat) = ud(k,imat+npri+nmat)
     end do
  end do

  ! Detect embedded body
  if(static)then
     do k=1,ncell
        bg(k)=fg(k,1)>0.01
        bm(k)=fm(k,1)>0.01
        bd(k)=fd(k,1)>0.01
     end do
     do k=1,ncell
        if(bm(k))then
           fg(k,1:nmat)=fm(k,1:nmat)
           gg(k,1:nmat)=gm(k,1:nmat)
           ug(k,1:nvar)=um(k,1:nvar)
           fd(k,1:nmat)=fm(k,1:nmat)
           gd(k,1:nmat)=gm(k,1:nmat)
           ud(k,1:nvar)=um(k,1:nvar)
        else
           if(bg(k))then
              fg(k,1:nmat)=fm(k,1:nmat)
              gg(k,1:nmat)=gm(k,1:nmat)
              ug(k,1:nvar)=um(k,1:nvar)
              ug(k,current_dim+1)=-um(k,current_dim+1)
           endif
           if(bd(k))then
              fd(k,1:nmat)=fm(k,1:nmat)
              gd(k,1:nmat)=gm(k,1:nmat)
              ud(k,1:nvar)=um(k,1:nvar)
              ud(k,current_dim+1)=-um(k,current_dim+1)
           endif
        endif
     enddo
  endif

  ! Compute total density
  do k = 1,ncell
     qg(k,1) = ug(k,1)
     qm(k,1) = um(k,1)
     qd(k,1) = ud(k,1)
  end do

  ! Compute velocity and specific kinetic energy
  eking(1:ncell)=0.0
  ekinm(1:ncell)=0.0
  ekind(1:ncell)=0.0
  do idim = 1,ndim
     do k = 1,ncell
        qg(k,idim+1) = ug(k,idim+1)/ug(k,1)
        qm(k,idim+1) = um(k,idim+1)/um(k,1)
        qd(k,idim+1) = ud(k,idim+1)/ud(k,1)
        eking(k) = eking(k) + half*qg(k,idim+1)**2
        ekinm(k) = ekinm(k) + half*qm(k,idim+1)**2
        ekind(k) = ekind(k) + half*qd(k,idim+1)**2
     end do
  end do

  ! Compute total internal energy
  do k = 1,ncell
     qg(k,npri) = ug(k,npri) - qg(k,1)*eking(k)
     qm(k,npri) = um(k,npri) - qm(k,1)*ekinm(k)
     qd(k,npri) = ud(k,npri) - qd(k,1)*ekind(k)
  end do

  ! Call eos routine
  call eos(fg,gg,qg,pg,cg,kappa_matg,kappa_hatg,ncell)
  call eos(fm,gm,qm,pm,cm,kappa_matm,kappa_hatm,ncell)
  call eos(fd,gd,qd,pd,cd,kappa_matd,kappa_hatd,ncell)

  ! Compute errors
  if(err_grad_d >= 0.)then
     do k=1,ncell
        ddg=abs(qg(k,1)); ddm=abs(qm(k,1)); ddd=abs(qd(k,1))
        error=2.0d0*MAX( &
             & ABS((ddd-ddm)/(ddd+ddm+floor_d)) , &
             & ABS((ddm-ddg)/(ddm+ddg+floor_d)) )
        ok(k) = ok(k) .or. error > err_grad_d
     end do
  end if

  if(err_grad_f >= 0.)then
     do imat=1,nmat
        do k=1,ncell
           ffg=fg(k,imat); ffm=fm(k,imat); ffd=fd(k,imat)
           error=2.0d0*MAX( &
                & ABS((ffd-ffm)/(ffd+ffm+floor_f)) , &
                & ABS((ffm-ffg)/(ffm+ffg+floor_f)) )
           ok(k) = ok(k) .or. error > err_grad_f
        end do
     end do
  end if

  if(err_grad_p > -1.0)then
     do k=1,ncell
        ppg=pg(k); ppm=pm(k); ppd=pd(k)
        error=2.0d0*MAX( &
             & ABS((ppd-ppm)/(ppd+ppm+floor_p)), &
             & ABS((ppm-ppg)/(ppm+ppg+floor_p)) )
        ok(k) = ok(k) .or. error > err_grad_p
     end do
  end if

  if(err_grad_u >= 0.)then
     do idim = 1,ndim
        do k=1,ncell
           vvg=qg(k,idim+1); vvm=qm(k,idim+1); vvd=qd(k,idim+1)
           ccg=cg(k)       ; ccm=cm(k)       ; ccd=cd(k)
           error=2.0d0*MAX( &
                & ABS((vvd-vvm)/(ccd+ccm+ABS(vvd)+ABS(vvm)+floor_u)) , &
                & ABS((vvm-vvg)/(ccm+ccg+ABS(vvm)+ABS(vvg)+floor_u)) )
           ok(k) = ok(k) .or. error > err_grad_u
        end do
     end do
  end if

!!$  if(static)then
!!$     do k=1,ncell
!!$        if(wg(k).or.wd(k))then
!!$           ddg=abs(qg(k,1)); ddm=abs(qm(k,1)); ddd=abs(qd(k,1))
!!$           error=2.0d0*MAX( &
!!$                & ABS((ddd-ddm)/(ddd+ddm+floor_d)) , &
!!$                & ABS((ddm-ddg)/(ddm+ddg+floor_d)) )
!!$           write(*,*)wg(k),wd(k)
!!$           write(*,*)bg(k),bd(k)
!!$           write(*,*)ddg,ddm,ddd
!!$           write(*,*)error,ok(k)
!!$        endif
!!$     end do
!!$  endif

end subroutine hydro_refine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine riemann_acoustic(fl,fr,gl,gr,ql,qr,cl,cr,fgdnv,ggdnv,qgdnv,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  ! dummy arguments
  integer::ngrid
  real(dp),dimension(1:nvector,1:npri)::ql,qr,qgdnv
  real(dp),dimension(1:nvector,1:nmat)::fl,fr,fgdnv
  real(dp),dimension(1:nvector,1:nmat)::gl,gr,ggdnv
  real(dp),dimension(1:nvector,1:nmat)::kappa_matl,kappa_matr
  real(dp),dimension(1:nvector)::cl,cr
  real(dp),dimension(1:nvector)::kappa_hatl,kappa_hatr

  ! local variables
  integer::i,n,ir,ie,imat
  real(dp)::smallp,delp_p

  ! local arrays
  real(dp),dimension(1:nvector,1:npri),save::qo,qstar
  real(dp),dimension(1:nvector,1:nmat),save::fo,fstar
  real(dp),dimension(1:nvector,1:nmat),save::go,gstar
  real(dp),dimension(1:nvector,1:nmat),save::kappa_matstar
  real(dp),dimension(1:nvector),save::cstar,estar,pstar,ustar
  real(dp),dimension(1:nvector),save::kappa_hatstar
  real(dp),dimension(1:nvector),save::sgnm ,spin ,spout,ushock
  real(dp),dimension(1:nvector),save::frac,co,el,er
  logical ,dimension(1:nvector),save::wall
  real(dp)::wl,wr,ul,ur,pl,pr,dl,dr,dstar

  ! Sound speed
  call eosinv(fl,gl,ql,el,cl,kappa_matl,kappa_hatl,ngrid)
  call eosinv(fr,gr,qr,er,cr,kappa_matr,kappa_hatr,ngrid)

  ! Acoustic star state
  do i=1,ngrid
     dl = ql(i,1); dr = qr(i,1)
     ul = ql(i,2); ur = qr(i,2)
     pl = ql(i,npri); pr = qr(i,npri)
     wl = cl(i)*dl
     wr = cr(i)*dr
     pstar(i) = ( (wr*pl+wl*pr)+wl*wr*(ul-ur) ) / (wl+wr)
     ustar(i) = ( (wr*ur+wl*ul)+      (pl-pr) ) / (wl+wr)
  end do

  ! Left going or right going contact wave
  do i=1,ngrid
     sgnm(i) = sign(one,ustar(i))
  end do

  ! Left or right unperturbed state
  do i = 1,ngrid
     if(sgnm(i)==one)then
        fo(i,1:nmat) = fl(i,1:nmat)
        go(i,1:nmat) = gl(i,1:nmat)
        qo(i,1:npri) = ql(i,1:npri)
        co(i) = cl(i)
     else
        fo(i,1:nmat) = fr(i,1:nmat)
        go(i,1:nmat) = gr(i,1:nmat)
        qo(i,1:npri) = qr(i,1:npri)
        co(i) = cr(i)
     end if
  end do

  ! Star region density, internal energy and sound speed
  do i = 1,ngrid
     dstar = qo(i,1) + (pstar(i)-qo(i,npri))/co(i)**2
     dstar = max(dstar,smallr)
     fstar(i,1:nmat) = fo(i,1:nmat)
     gstar(i,1:nmat) = go(i,1:nmat)
     qstar(i,1) = dstar
     qstar(i,2) = ustar(i)
     qstar(i,npri) = pstar(i)
#if NDIM>1
     qstar(i,3) = qo(i,3)
#endif
#if NDIM>2
     qstar(i,4) = qo(i,4)
#endif
  end do
  call eosinv(fstar,gstar,qstar,estar,cstar,kappa_matstar,kappa_hatstar,ngrid)

  ! Head and tail speed of rarefaction
  do i=1,ngrid
     spout(i) = co   (i)-sgnm(i)*qo   (i,2)
     spin (i) = cstar(i)-sgnm(i)*qstar(i,2)
  end do

  ! Shock speed
  do i=1,ngrid
     ushock(i) = half*(spin(i)+spout(i))
     ushock(i) = max(ushock(i),-sgnm(i)*qstar(i,2))
  end do
  do i=1,ngrid
     if(pstar(i)>=qo(i,npri))then
        spout(i)=ushock(i)
        spin (i)=spout (i)
     end if
  end do

  ! Sample the solution at x/t=0
  do i=1,ngrid
     if(spout(i)<zero)then      ! Initial state
        qgdnv(i,1:npri) = qo(i,1:npri)
        fgdnv(i,1:nmat) = fo(i,1:nmat)
        ggdnv(i,1:nmat) = go(i,1:nmat)
     else if(spin(i)>=zero)then  ! Star region
        qgdnv(i,1:npri) = qstar(i,1:npri)
        fgdnv(i,1:nmat) = fstar(i,1:nmat)
        ggdnv(i,1:nmat) = gstar(i,1:nmat)
     else                        ! Rarefaction
        frac(i) = spout(i)/(spout(i)-spin(i))
        qgdnv(i,1:npri) = frac(i)*qstar(i,1:npri) + (one - frac(i))*qo(i,1:npri)
        fgdnv(i,1:nmat) = frac(i)*fstar(i,1:nmat) + (one - frac(i))*fo(i,1:nmat)
        ggdnv(i,1:nmat) = frac(i)*gstar(i,1:nmat) + (one - frac(i))*go(i,1:nmat)
     end if
  end do

end subroutine riemann_acoustic
!###########################################################
!###########################################################
!###########################################################
!###########################################################
