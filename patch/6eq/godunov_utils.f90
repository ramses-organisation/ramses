!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine eos(d,e,p,c,imat,inv,ncell)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none
  integer::imat
  integer::ncell
  logical::inv
  real(dp),dimension(1:nvector)::d,e,p,c
  ! Compute total internal energy and sound speed from total pressure
  ! On entry:
  !   d is the true density of each fluid
  !   imat is the identifier of the fluid species
  !   inv is a logical defining the case:
  !     inv=0:(d,e)-->(p,c) (e is given)
  !     inv=1:(d,p)-->(e,c) (p is given)
  ! On exit:
  !   c is the sound speed of each fluid
  !   e is the internal energy of each fluid
  !   p is the pressure of each fluid
  integer::k
  real(dp)::smallgamma,biggamma,p_0,rho_0,e_c,p_c,delpc,eta,q
  real(dp)::E_1,E_2,A_1,A_2,C_v,T_0,E_0,p_c_1,p_c_2
  do k=1,ncell

     if(eos_name =='stiffened gas')then

        ! Get Stiffened Gas EOS parameters
        smallgamma=eos_params(imat,1);q=eos_params(imat,2);p_0=eos_params(imat,3)
        ! Use the EOS to calculate the current pressure/internal energy in a given cell
        ! P  = (gamma - one) * (e - d*q) + gamma*P_inf
        if(inv .eqv. .false.)then
           p(k) = (smallgamma-1)*(e(k)-d(k)*q) - smallgamma*p_0
        else if(inv .eqv. .true.)then
           e(k) = (1/(smallgamma-1))*(p(k)+smallgamma*p_0) + d(k)*q
        end if

        ! Calculate the speed of sound
        c(k) = smallgamma*(p(k)+p_0)/d(k)
        c(k) = sqrt(max(c(k),smallc**2))

     else if(eos_name =='mie-grueneisen')then

        ! Get Mie-Grueneisen EOS parameters
        smallgamma=eos_params(imat,1);biggamma=eos_params(imat,2);p_0=eos_params(imat,3);rho_0=eos_params(imat,4)
        eta = d(k)/rho_0
        p_c = p_0 * eta**biggamma
        e_c = p_c / (biggamma-one)
        delpc = biggamma * p_c

        ! Use the EOS to calculate the current pressure/internal energy in a given cell
        ! P - P_c = (gamma - one) * (e - e_c) ; e = e_c + (P - P_c) / (gamma - one)
        if(inv .eqv. .false.)then
           p(k) = (smallgamma-1)*(e(k)-e_c) + p_c
        else if(inv .eqv. .true.)then
           e(k) = (1/(smallgamma-1))*(p(k)-p_c) + e_c
        end if

        ! Calculate the speed of sound
        ! c**2 = P_c' + gamma * (P - P_c) / rho
        c(k) = (delpc + smallgamma * (p(k)-p_c)) / d(k)
        c(k) = sqrt(max(c(k),smallc**2))

     else if(eos_name == 'cochran-chan')then

        ! Get Mie-Grueneisen EOS paramete
        smallgamma=eos_params(imat,1);rho_0=eos_params(imat,2)
        E_1=eos_params(imat,3);E_2=eos_params(imat,4)
        A_1=eos_params(imat,5);A_2=eos_params(imat,6)
        C_v=eos_params(imat,7);T_0=eos_params(imat,8)

        ! Define the Cochran-Chan constant term
        E_0 = A_1 / (E_1-one) - A_2 / (E_2-one) + rho_0 * C_v * T_0

        ! Update Mie-Gruneisen terms for each material
        eta   = d(k)/rho_0
        p_c_1 = A_1 * eta**E_1
        p_c_2 = A_2 * eta**E_2
        p_c   = p_c_1 - p_c_2
        e_c   = p_c_1 / (E_1-one) - p_c_2 / (E_2-one) - eta * E_0
        delpc = p_c_1 * E_1 - p_c_2 * E_2

        ! Use the EOS to calculate the current pressure/internal energy in a given cell
        ! P - P_c = (gamma - one) * (e - e_c) ; e = e_c + (P - P_c) / (gamma - one)
        if(inv .eqv. .false.)then
           p(k) = (smallgamma-1)*(e(k)-e_c) + p_c
        else if(inv .eqv. .true.)then
           e(k) = (1/(smallgamma-1))*(p(k)-p_c) + e_c
        end if

        ! Calculate the speed of sound of each fluid
        ! c**2 = P_c' + gamma * (P - P_c) / rho
        c(k) = (delpc + smallgamma * (p(k)-p_c) ) / d(k)
        c(k) = sqrt(max(c(k),smallc**2))

     end if
  end do
end subroutine eos
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine eos_s(d,e,s,imat,inv,ncell)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none
  integer::imat
  integer::ncell
  logical::inv
  real(dp),dimension(1:nvector)::d,e,s
  ! Compute entropy from internal energy
  ! On entry:
  !   d is the true density of each fluid
  !   imat is the identifier of the fluid species
  !   inv is a logical defining the case:
  !     inv=0:(d,e)-->(s) (e is given)
  !     inv=1:(d,s)-->(e) (s is given)
  ! On exit:
  !   e is the internal energy of each fluid
  !   s is the entropy of each fluid
  integer::k
  real(dp)::smallgamma,biggamma,p_0,rho_0,e_c,p_c,delpc,eta
  real(dp)::E_1,E_2,A_1,A_2,C_v,T_0,E_0,p_c_1,p_c_2

  do k=1,ncell

     if(eos_name=='mie-grueneisen')then

        ! Get Mie-Grueneisen EOS parameters
        smallgamma=eos_params(imat,1);biggamma=eos_params(imat,2);p_0=eos_params(imat,3);rho_0=eos_params(imat,4)
        eta = d(k)/rho_0
        p_c = p_0 * eta**biggamma
        e_c = p_c / (biggamma-one)

        ! Use the EOS to calculate the current entropy/internal energy in a given cell
        ! s = (e - e_c) / rho**gamma; e = e_c + s * rho**gamma
        if(inv .eqv. .false.)then
           s(k) = (e(k) - e_c) / d(k)**smallgamma
        else if(inv .eqv. .true.)then
           e(k) = e_c + s(k) * d(k)**smallgamma
        end if

     else if(eos_name == 'cochran-chan')then

        ! Cochran-Chan EOS written in terms of the Mie-Grueneisen EOS
        ! Get Mie-Grueneisen EOS parameters
        smallgamma=eos_params(imat,1);rho_0=eos_params(imat,2)
        E_1=eos_params(imat,3);E_2=eos_params(imat,4)
        A_1=eos_params(imat,5);A_2=eos_params(imat,6)
        C_v=eos_params(imat,7);T_0=eos_params(imat,8)

        ! Define the Cochran-Chan constant term
        E_0 = A_1 / (E_1-one) - A_2 / (E_2-one) + rho_0 * C_v * T_0

        ! Update Mie-Gruneisen terms for each material
        eta   = d(k)/rho_0
        p_c_1 = A_1 * eta**E_1
        p_c_2 = A_2 * eta**E_2
        p_c   = p_c_1 - p_c_2
        e_c   = p_c_1 / (E_1-one) - p_c_2 / (E_2-one) - eta * E_0

        ! Use the EOS to calculate the current entropy/internal energy in a given cell
        ! s = (e - e_c) / rho**gamma; e = e_c + s * rho**gamma
        if(inv .eqv. .false.)then
           s(k) = (e(k) - e_c) / d(k)**smallgamma
        else if(inv .eqv. .true.)then
           e(k) = e_c + s(k) * d(k)**smallgamma
        end if

     end if

  end do

end subroutine eos_s
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
  logical::inv
  real(dp)::dx,dt
  real(dp),dimension(1:nvector,1:nvar)::uu
  real(dp),dimension(1:nvector,1:ndim)::grav
  real(dp),dimension(1:nvector)::rr,dtot

  real(dp),dimension(1:nvector,1:npri),save::qq
  real(dp),dimension(1:nvector,1:nmat),save::ff,gg
  real(dp),dimension(1:nvector),save::gg_mat,ee_mat,pp_mat,cc_mat
  real(dp),dimension(1:nvector),save::ekin,cc,st
  real(dp)::dtcell,eps
  integer::k,idim,imat

  ! Convert to primitive variable

  ! Volume fraction and fluid density
  do imat = 1,nmat
     do k = 1,ncell
        ff(k,imat) = uu(k,imat)
        gg(k,imat) = uu(k,nmat+imat)/ff(k,imat)
     end do
  end do

  ! Compute total density
  dtot(1:ncell) = 0.0
  do imat=1,nmat
    do k = 1,ncell
       dtot(k)  = dtot(k) + uu(k,nmat+imat)
    end do
  end do

  ! Compute velocity and specific kinetic energy
  ekin(1:ncell)    = 0.0
  do idim = 1,ndim
     do k = 1,ncell
        qq(k,idim) = uu(k,2*nmat+idim)/dtot(k)
        ekin(k)    = ekin(k) + half*qq(k,idim)**2
     end do
  end do

  ! Compute partial internal energies
  do imat=1,nmat
    do k = 1,ncell
       qq(k,ndim+nmat+imat) = uu(k,2*nmat+ndim+imat)/ff(k,imat) - gg(k,imat)*ekin(k)
    end do
  end do

  ! Calculate the total speed of sound
  cc(1:ncell)=0
  inv=.false.
  do imat=1,nmat
    do k=1,ncell
      gg_mat(k) = gg(k,imat)
      ee_mat(k) = qq(k,ndim+nmat+imat)
    end do
    ! Call eos routine
    call eos(gg_mat,ee_mat,pp_mat,cc_mat,imat,inv,ncell)
    do k=1,ncell
      cc(k) = cc(k) + ff(k,imat)*gg(k,imat) * cc_mat(k)**2
    end do
  end do
  ! Convert rho c^2 to c
  cc(1:ncell)=sqrt(cc(1:ncell)/dtot(1:ncell))

  ! Compute wave speed
  do k = 1,ncell
     cc(k) = abs(qq(k,1))+cc(k)
  end do
  do idim = 2,ndim
     do k = 1,ncell
        cc(k) = cc(k) + abs(qq(k,idim))+cc(k)
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
  logical::inv
  real(dp),dimension(1:nvector,1:nvar)::ug,um,ud
  logical ,dimension(1:nvector)       ::ok

  integer::k,idim,imat
  real(dp),dimension(1:nvector,1:npri),save::qg,qm,qd
  real(dp),dimension(1:nvector,1:nmat),save::fg,fm,fd
  real(dp),dimension(1:nvector,1:nmat),save::gg,gm,gd
  real(dp),dimension(1:nvector),save::dtotg,dtotm,dtotd
  real(dp),dimension(1:nvector),save::gg_mat,gm_mat,gd_mat
  real(dp),dimension(1:nvector),save::eking,ekinm,ekind
  real(dp),dimension(1:nvector),save::eg_mat,em_mat,ed_mat
  real(dp),dimension(1:nvector),save::pg,pm,pd
  real(dp),dimension(1:nvector),save::pg_mat,pm_mat,pd_mat
  real(dp),dimension(1:nvector),save::cg,cm,cd
  real(dp),dimension(1:nvector),save::cg_mat,cm_mat,cd_mat
  logical ,dimension(1:nvector),save::wg,wd,bg,bm,bd
  real(dp)::ffg,ffm,ffd,ddg,ddm,ddd
  real(dp)::ppg,ppm,ppd,vvg,vvm,vvd
  real(dp)::ccg,ccm,ccd,error

  ! Convert to primitive variables

  ! Volume fraction and fluid density
  do imat = 1,nmat
     do k = 1,ncell
        fg(k,imat) = ug(k,imat)
        fm(k,imat) = um(k,imat)
        fd(k,imat) = ud(k,imat)
        gg(k,imat) = ug(k,imat+nmat)/fg(k,imat)
        gm(k,imat) = um(k,imat+nmat)/fm(k,imat)
        gd(k,imat) = ud(k,imat+nmat)/fd(k,imat)
     end do
  end do

  ! Compute total density
  dtotg(1:ncell)=0.0
  dtotm(1:ncell)=0.0
  dtotd(1:ncell)=0.0
  do imat=1,nmat
     do k = 1,ncell
        dtotg(k) = dtotg(k) + ug(k,nmat+imat)
        dtotm(k) = dtotm(k) + um(k,nmat+imat)
        dtotd(k) = dtotd(k) + ud(k,nmat+imat)
     end do
  end do

  ! Compute velocity and specific kinetic energy
  eking(1:ncell)=0.0
  ekinm(1:ncell)=0.0
  ekind(1:ncell)=0.0
  do idim = 1,ndim
     do k = 1,ncell
        qg(k,idim) = ug(k,2*nmat+idim)/dtotg(k)
        qm(k,idim) = um(k,2*nmat+idim)/dtotm(k)
        qd(k,idim) = ud(k,2*nmat+idim)/dtotd(k)
        eking(k) = eking(k) + half*qg(k,idim)**2
        ekinm(k) = ekinm(k) + half*qm(k,idim)**2
        ekind(k) = ekind(k) + half*qd(k,idim)**2
     end do
  end do

  ! Compute total internal energy
  do imat=1,nmat
     do k = 1,ncell
        qg(k,ndim+nmat+imat) = ug(k,2*nmat+ndim+imat)/ug(k,imat) - gg(k,imat)*eking(k)
        qm(k,ndim+nmat+imat) = um(k,2*nmat+ndim+imat)/um(k,imat) - gm(k,imat)*ekinm(k)
        qd(k,ndim+nmat+imat) = ud(k,2*nmat+ndim+imat)/ud(k,imat) - gd(k,imat)*ekind(k)
     end do
  end do

  ! Call eos routine to calculate the total pressure and the total speed of sound
  inv=.false.
  pg(1:ncell)=0
  pm(1:ncell)=0
  pd(1:ncell)=0
  cg(1:ncell)=0
  cm(1:ncell)=0
  cd(1:ncell)=0
  do imat=1,nmat
     do k=1,ncell
        gg_mat(k) = gg(k,imat)
        gm_mat(k) = gm(k,imat)
        gd_mat(k) = gd(k,imat)
        eg_mat(k) = qg(k,ndim+nmat+imat)
        em_mat(k) = qm(k,ndim+nmat+imat)
        ed_mat(k) = qd(k,ndim+nmat+imat)
     end do
     call eos(gg_mat,eg_mat,pg,cg_mat,imat,inv,ncell)
     call eos(gm_mat,em_mat,pm,cm_mat,imat,inv,ncell)
     call eos(gd_mat,ed_mat,pd,cd_mat,imat,inv,ncell)
     do k=1,ncell
        pg(k) = pg(k) + fg(k,imat) * pg_mat(k)
        pm(k) = pm(k) + fm(k,imat) * pg_mat(k)
        pd(k) = pd(k) + fd(k,imat) * pg_mat(k)
        cg(k) = cg(k) + fg(k,imat) * gg(k,imat) * cg_mat(k)**2
        cm(k) = cm(k) + fm(k,imat) * gm(k,imat) * cm_mat(k)**2
        cd(k) = cd(k) + fd(k,imat) * gd(k,imat) * cd_mat(k)**2
     end do
  end do
  ! Convert c^2 to c
  cg(1:ncell)=sqrt(cg(1:ncell)/dtotg(1:ncell))
  cm(1:ncell)=sqrt(cm(1:ncell)/dtotm(1:ncell))
  cd(1:ncell)=sqrt(cd(1:ncell)/dtotd(1:ncell))

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
           vvg=qg(k,idim); vvm=qm(k,idim); vvd=qd(k,idim)
           ccg=cg(k)     ; ccm=cm(k)     ; ccd=cd(k)
           error=2.0d0*MAX( &
                & ABS((vvd-vvm)/(ccd+ccm+ABS(vvd)+ABS(vvm)+floor_u)) , &
                & ABS((vvm-vvg)/(ccm+ccg+ABS(vvm)+ABS(vvg)+floor_u)) )
           ok(k) = ok(k) .or. error > err_grad_u
        end do
     end do
  end if

end subroutine hydro_refine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
! subroutine riemann_acoustic(fl,fr,gl,gr,ql,qr,cl,cr,fgdnv,ggdnv,qgdnv,ngrid)
!   use amr_parameters
!   use hydro_parameters
!   use const
!   implicit none

!   ! dummy arguments
!   integer::ngrid
!   logical::inv
!   real(dp),dimension(1:nvector,1:npri)::ql,qr,qgdnv
!   real(dp),dimension(1:nvector,1:nmat)::fl,fr,fgdnv
!   real(dp),dimension(1:nvector,1:nmat)::gl,gr,ggdnv
!   real(dp),dimension(1:nvector)::cl,cr
!   real(dp),dimension(1:nvector)::cl_mat,cr_mat

!   ! local variables
!   integer::i,n,ir,ie,imat
!   real(dp)::smallp,delp_p

!   ! local arrays
!   real(dp),dimension(1:nvector,1:npri),save::qo,qstar
!   real(dp),dimension(1:nvector,1:nmat),save::fo,fstar
!   real(dp),dimension(1:nvector,1:nmat),save::go,gstar
!   real(dp),dimension(1:nvector),save::cstar,estar,pstar,ustar
!   real(dp),dimension(1:nvector),save::cstar_mat,pstar_mat
!   real(dp),dimension(1:nvector),save::sgnm ,spin ,spout,ushock
!   real(dp),dimension(1:nvector),save::frac,co,el,er
!   logical ,dimension(1:nvector),save::wall
!   real(dp)::wl,wr,ul,ur,pl,pr,dl,dr,dstar
!   real(dp),dimension(1:nvector),save::::pl_mat,pr_mat ! Dummy variables

!   ! Sound speed
!   inv=.true.
!   cl(1:ngrid)=0
!   cr(1:ngrid)=0
!   do imat=1,nmat
!     call eos(gl(:,imat),el,ql(k,npri+imat-1),cl_mat,imat,inv,ngrid)
!     call eos(gr(:,imat),er,qr(k,npri+imat-1),cr_mat,imat,inv,ngrid)
!     do k=1,ngrid
!       cl(k) = cl(k) + (fl(k,imat)*gl(k,imat)/ql(k,1)) * cl_mat(k)
!       cr(k) = cr(k) + (fr(k,imat)*gr(k,imat)/qr(k,1)) * cr_mat(k)
!     end do
!   end do
!   ! Convert c^2 to c
!   cl(1:ngrid)=sqrt(cl(1:ngrid))
!   cr(1:ngrid)=sqrt(cr(1:ngrid))

!   ! Acoustic star state
!   do i=1,ngrid
!      dl = ql(i,1); dr = qr(i,1)
!      ul = ql(i,2); ur = qr(i,2)
!      pl = ql(i,npri); pr = qr(i,npri)
!      wl = cl(i)*dl
!      wr = cr(i)*dr
!      pstar(i) = ( (wr*pl+wl*pr)+wl*wr*(ul-ur) ) / (wl+wr)
!      ustar(i) = ( (wr*ur+wl*ul)+      (pl-pr) ) / (wl+wr)
!   end do

!   ! Left going or right going contact wave
!   do i=1,ngrid
!      sgnm(i) = sign(one,ustar(i))
!   end do

!   ! Left or right unperturbed state
!   do i = 1,ngrid
!      if(sgnm(i)==one)then
!         fo(i,1:nmat) = fl(i,1:nmat)
!         go(i,1:nmat) = gl(i,1:nmat)
!         qo(i,1:npri) = ql(i,1:npri)
!         co(i) = cl(i)
!      else
!         fo(i,1:nmat) = fr(i,1:nmat)
!         go(i,1:nmat) = gr(i,1:nmat)
!         qo(i,1:npri) = qr(i,1:npri)
!         co(i) = cr(i)
!      end if
!   end do

!   ! Star region density, internal energy and sound speed
!   do i = 1,ngrid
!      dstar = qo(i,1) + (pstar(i)-qo(i,npri))/co(i)**2
!      dstar = max(dstar,smallr)
!      fstar(i,1:nmat) = fo(i,1:nmat)
!      gstar(i,1:nmat) = go(i,1:nmat)
!      qstar(i,1) = dstar
!      qstar(i,2) = ustar(i)
!      qstar(i,npri) = pstar(i)
! #if NDIM>1
!      qstar(i,3) = qo(i,3)
! #endif
! #if NDIM>2
!      qstar(i,4) = qo(i,4)
! #endif
!   end do

!   ! Sound speed
!   inv=.true.
!   cl(1:ngrid)=0
!   do k=1,ngrid
!     do imat=1,nmat
!       call eos(gstar,qstar,estar,pstar_mat,cstar_mat,imat,inv,ngrid)
!       cstar(k) = cstar(k) + (fstar(k,imat)*gstar(k,imat)/qstar(k,1)) * cstar_mat(k)
!     end do
!   end do
!   ! Convert c^2 to c
!   cstar(1:ngrid)=sqrt(cstar(1:ngrid))


!   ! Head and tail speed of rarefaction
!   do i=1,ngrid
!      spout(i) = co   (i)-sgnm(i)*qo   (i,2)
!      spin (i) = cstar(i)-sgnm(i)*qstar(i,2)
!   end do

!   ! Shock speed
!   do i=1,ngrid
!      ushock(i) = half*(spin(i)+spout(i))
!      ushock(i) = max(ushock(i),-sgnm(i)*qstar(i,2))
!   end do
!   do i=1,ngrid
!      if(pstar(i)>=qo(i,npri))then
!         spout(i)=ushock(i)
!         spin (i)=spout (i)
!      end if
!   end do

!   ! Sample the solution at x/t=0
!   do i=1,ngrid
!      if(spout(i)<zero)then      ! Initial state
!         qgdnv(i,1:npri) = qo(i,1:npri)
!         fgdnv(i,1:nmat) = fo(i,1:nmat)
!         ggdnv(i,1:nmat) = go(i,1:nmat)
!      else if(spin(i)>=zero)then  ! Star region
!         qgdnv(i,1:npri) = qstar(i,1:npri)
!         fgdnv(i,1:nmat) = fstar(i,1:nmat)
!         ggdnv(i,1:nmat) = gstar(i,1:nmat)
!      else                        ! Rarefaction
!         frac(i) = spout(i)/(spout(i)-spin(i))
!         qgdnv(i,1:npri) = frac(i)*qstar(i,1:npri) + (one - frac(i))*qo(i,1:npri)
!         fgdnv(i,1:nmat) = frac(i)*fstar(i,1:nmat) + (one - frac(i))*fo(i,1:nmat)
!         ggdnv(i,1:nmat) = frac(i)*gstar(i,1:nmat) + (one - frac(i))*go(i,1:nmat)
!      end if
!   end do

! end subroutine riemann_acoustic
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine riemann_hllc(fl,fr,gl,gr,ql,qr,cl,cr,fgdnv,ugdnv,egdnv,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none
  ! HLLC Riemann solver (Toro)
  ! dummy arguments
  integer::ngrid
  real(dp),dimension(1:nvector,1:nmat)::fl,fr
  real(dp),dimension(1:nvector,1:nmat)::gl,gr
  real(dp),dimension(1:nvector,1:npri)::ql,qr
  real(dp),dimension(1:nvector,1:nmat)::cl,cr
  real(dp),dimension(1:nvector,1:nvar)::fgdnv
  real(dp),dimension(1:nvector,1:nmat)::egdnv
  real(dp),dimension(1:nvector)::ugdnv
  ! local variables
  REAL(dp)::SL,SR
  REAL(dp)::rl,ul,Ptotl
  real(dp),dimension(1:nmat)::Pl,ecinl,eintl,ekl
  real(dp),dimension(1:nmat)::Pr,ecinr,eintr,ekr
  REAL(dp)::rr,ur,Ptotr
  REAL(dp)::cfastl,rcl,rstarl,estarl
  REAL(dp)::cfastr,rcr,rstarr,estarr
  real(dp),dimension(1:nmat)::pkstarl,gkstarl,ekstarl
  real(dp),dimension(1:nmat)::pkstarr,gkstarr,ekstarr
  REAL(dp)::ustar,ptotstar
  REAL(dp)::ro,uo,ptoto,eo
  real(dp),dimension(1:nmat)::gko,fko,eko,pko
  INTEGER::ivar,i,imat
#if NVAR > NDIM + 3*NMAT
  integer::ipscal,npscal,n
#endif
  do i=1,ngrid
     ! Left variables
     ul    = ql(i,1)
     rl    = zero
     do imat=1,nmat
        rl          = rl + gl(i,imat)*fl(i,imat)
        Pl(imat)    = ql(i,ndim+imat)
        eintl(imat) = ql(i,ndim+nmat+imat)
        ecinl(imat) = half*gl(i,imat)*ul*ul
#if NDIM>1
        ecinl(imat) = ecinl(imat)+half*gl(i,imat)*ql(i,2)**2
#endif
#if NDIM>2
        ecinl(imat) = ecinl(imat)+half*gl(i,imat)*ql(i,3)**2
#endif
     end do
     ! We want to calculate the mixture speed of sound here
     cfastl = zero
     Ptotl  = zero
     do imat=1,nmat
        cfastl    = cfastl + fl(i,imat)*gl(i,imat)*cl(i,imat)**2  ! d.c^2
        Ptotl     = Ptotl  + fl(i,imat)*Pl(imat)
        ekl(imat) = eintl(imat) + ecinl(imat)
     end do
     ! Mixture speed of sound
     cfastl=sqrt(cfastl/rl)
     ! Right variables
     ur    = qr(i,1)
     rr    = zero
     do imat=1,nmat
        rr          = rr + gr(i,imat)*fr(i,imat)
        Pr(imat)    = qr(i,ndim+imat)
        eintr(imat) = qr(i,ndim+nmat+imat)
        ecinr(imat) = half*gr(i,imat)*ur*ur
#if NDIM>1
        ecinr(imat) = ecinr(imat)+half*gr(i,imat)*qr(i,2)**2
#endif
#if NDIM>2
        ecinr(imat) = ecinr(imat)+half*gr(i,imat)*qr(i,3)**2
#endif
     end do
     ! We want to calculate the mixture speed of sound here
     cfastr = zero
     Ptotr  = zero
     do imat=1,nmat
        cfastr    = cfastr + fr(i,imat)*gr(i,imat)*cr(i,imat)**2  ! d.c^2
        Ptotr     = Ptotr  + fr(i,imat)*Pr(imat)
        ekr(imat) = eintr(imat) + ecinr(imat)
     end do
     ! Mixture speed of sound
     cfastr=sqrt(cfastr/rr)
     ! Compute HLL wave speed
     SL=min(ul,ur)-max(cfastl,cfastr)
     SR=max(ul,ur)+max(cfastl,cfastr)
     ! Compute lagrangian sound speed
     rcl=rl*(ul-SL)
     rcr=rr*(SR-ur)
     ! Compute acoustic star state
     ustar    = ((Ptotl - Ptotr) + rcr*ur + rcl*ul)/(rcr+rcl)
     Ptotstar = (rcr*Ptotl + rcl*Ptotr + rcl*rcr*(ul-ur))/(rcr+rcl)
     ! Left star region variables
     rstarl = rl*(SL-ul)/(SL-ustar)
     do imat = 1,nmat
        gkstarl(imat) = gl(i,imat)*(ul-SL)/(ustar-SL)
        pkstarl(imat) = Pl(imat)-gl(i,imat)*(ul-SL)*(ustar-ul)
        ekstarl(imat) = ((SL-ul)*ekl(imat)-Pl(imat)*ul+pkstarl(imat)*ustar)/(SL-ustar)
     end do
     ! Right star region variables
     rstarr = rr*(SR-ur)/(SR-ustar)
     do imat=1,nmat
        gkstarr(imat) = gr(i,imat)*(ur-SR)/(ustar-SR)
        pkstarr(imat) = Pr(imat)-gr(i,imat)*(ur-SR)*(ustar-ur)
        ekstarr(imat) = ((SR-ur)*ekr(imat)-Pr(imat)*ur+pkstarr(imat)*ustar)/(SR-ustar)
     end do
     ! Sample the solution at x/t=0
     if(SL>0d0)then
        ro=rl
        uo=ul
        Ptoto=Ptotl
        do imat=1,nmat
           fko(imat) = fl(i,imat)
           gko(imat) = gl(i,imat)
           eko(imat) = ekl(imat)
           pko(imat) = Pl(imat)
        end do
     else if(ustar>0d0)then
        ro=rstarl
        uo=ustar
        Ptoto=Ptotstar
        do imat=1,nmat
           fko(imat) = fl(i,imat)
           gko(imat) = gkstarl(imat)
           eko(imat) = ekstarl(imat)
           pko(imat) = pkstarl(imat)
        end do
     else if(SR>0d0)then
        ro=rstarr
        uo=ustar
        Ptoto=Ptotstar
        do imat=1,nmat
           fko(imat) = fr(i,imat)
           gko(imat) = gkstarr(imat)
           eko(imat) = ekstarr(imat)
           pko(imat) = pkstarr(imat)
        end do
     else
        ro=rr
        uo=ur
        Ptoto=Ptotr
        do imat=1,nmat
           fko(imat) = fr(i,imat)
           gko(imat) = gr(i,imat)
           eko(imat) = ekr(imat)
           pko(imat) = Pr(imat)
        end do
     end if
     !=========================
     ! Compute the Godunov flux
     !=========================
     ! Volume fractions
     do imat=1,nmat
        fgdnv(i,imat)      = uo*fko(imat)
     end do
     ! Physical densities
     do imat=1,nmat
        fgdnv(i,nmat+imat) = uo*gko(imat)*fko(imat)
     end do
     ! Momentum flux
     ugdnv(i) = uo
     fgdnv(i,2*nmat+1)     = ro*uo*uo + Ptoto
     ! Transverse velocities
#if NDIM > 1
     if(ustar>0)then
        fgdnv(i,2*nmat+2)  = ro*uo*ql(i,2)
     else
        fgdnv(i,2*nmat+2)  = ro*uo*qr(i,2)
     endif
#endif
#if NDIM > 2
     if(ustar>0)then
        fgdnv(i,2*nmat+3)  = ro*uo*ql(i,3)
     else
        fgdnv(i,2*nmat+3)  = ro*uo*qr(i,3)
     endif
#endif
     ! Energy fluxes
     do imat=1,nmat
        egdnv(i,imat) = fko(imat)*pko(imat)*uo
        fgdnv(i,2*nmat+ndim+imat) = fko(imat)*(eko(imat)+pko(imat))*uo
     end do
     ! Passive scalars
#if NVAR > NDIM + 3*NMAT
     npscal = (nvar - ndim - 3*nmat) / nmat
     do imat = 1, nmat
        do ipscal = 1, npscal
           n = ndim + 3*nmat + npscal*(imat-1) + ipscal
           if(ustar>0)then
              fgdnv(i,n)  = uo*gko(imat)*fko(imat)*ql(i,n-nmat)
           else
              fgdnv(i,n)  = uo*gko(imat)*fko(imat)*qr(i,n-nmat)
           endif
        end do
     end do
#endif
  end do
end subroutine riemann_hllc
!###########################################################
!###########################################################
!###########################################################
!###########################################################
