! ---------------------------------------------------------------
!  UNSPLIT     Unsplit second order Godunov integrator for
!              real gases (arbitrary EOS) dynamics using
!              MUSCL-HANCOCK scheme
!              with various slope limiters.
!
!  inputs/outputs
!  uin         => (const)  input state
!  gravin      => (const)  input gravitational acceleration
!  iu1,iu2     => (const)  first and last index of input array,
!  ju1,ju2     => (const)  cell centered,
!  ku1,ku2     => (const)  including buffer cells.
!  flux       <=  (modify) return fluxes in the 3 coord directions
!  if1,if2     => (const)  first and last index of output array,
!  jf1,jf2     => (const)  edge centered,
!  kf1,kf2     => (const)  for active cells only.
!  dx,dy,dz    => (const)  (dx,dy,dz)
!  dt          => (const)  time step
!  ngrid       => (const)  number of sub-grids
!  ndim        => (const)  number of dimensions
! ----------------------------------------------------------------
subroutine unsplit(uin,gravin,flux,tmp,dx,dy,dz,dt,ngrid)
  use amr_parameters
  use const
  use hydro_parameters
  implicit none

  integer ::ngrid
  real(dp)::dx,dy,dz,dt

  ! Input states
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::uin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim)::gravin

  ! Output fluxes
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim)::flux
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:nmat+1,1:ndim)::tmp

  ! Primitive variables
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:npri),save::qin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nmat),save::fin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nmat),save::gin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nmat),save::cin

  ! Slopes
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:npri,1:ndim),save::dq
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nmat,1:ndim),save::df
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nmat,1:ndim),save::dg

  ! Left and right state arrays
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:npri,1:ndim),save::qm
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:npri,1:ndim),save::qp
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nmat,1:ndim),save::fm
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nmat,1:ndim),save::fp
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nmat,1:ndim),save::gm
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nmat,1:ndim),save::gp

  ! Intermediate fluxes
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),save::fx
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nmat+1),save::tx

  ! Local scalar variables
  integer::i,j,k,l,ivar,ir,ie,imat,iu
  integer::ilo,ihi,jlo,jhi,klo,khi
  real(dp),dimension(1:nvector)::ekin,flux_dtot
  real(dp)::fmil

  ilo=MIN(1,iu1+2); ihi=MAX(1,iu2-2)
  jlo=MIN(1,ju1+2); jhi=MAX(1,ju2-2)
  klo=MIN(1,ku1+2); khi=MAX(1,ku2-2)

  ! Translate to primitive variable (d,u,P) volume and mass fraction (f,g)
  call ctoprim(uin,qin,fin,gin,cin,gravin,dt,ngrid)

  ! Compute TVD slopes
  call qslope(qin,dq,dx,dt,ngrid)
  call fslope(fin,df,qin,dx,dt,ngrid)
  call fslope(gin,dg,qin,dx,dt,ngrid)

  ! Compute 3D traced-states in all three directions
#if NDIM==1
     call trace1d (qin,fin,gin,dq,df,dg,cin,qm,qp,fm,fp,gm,gp,dx      ,dt,ngrid)
#endif
#if NDIM==2
     call trace2d (qin,fin,gin,dq,df,dg,cin,qm,qp,fm,fp,gm,gp,dx,dy   ,dt,ngrid)
#endif
#if NDIM==3
     call trace3d (qin,fin,gin,dq,df,dg,cin,qm,qp,fm,fp,gm,gp,dx,dy,dz,dt,ngrid)
#endif

  ! Compute 1D flux in X direction
  call cmpgdnv(fm,gm,qm,iu1+1,iu2+1,ju1  ,ju2  ,ku1  ,ku2  , &
       &       fp,gp,qp,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , cin, &
       &          fx,tx,if1  ,if2  ,jlo  ,jhi  ,klo  ,khi  , 1,2,3,ngrid)
  ! Save flux in output array
  do i=if1,if2
  do j=jlo,jhi
  do k=klo,khi
     do ivar=1,nvar
        do l=1,ngrid
           flux(l,i,j,k,ivar,1)=fx(l,i,j,k,ivar)*dt/dx
        end do
     end do
     do ivar=1,nmat+1
        do l=1,ngrid
           tmp (l,i,j,k,ivar,1)=tx(l,i,j,k,ivar)*dt/dx
        end do
     end do
  end do
  end do
  end do

  ! Solve for 1D flux in Y direction
#if NDIM>1
  call cmpgdnv(fm,gm,qm,iu1  ,iu2  ,ju1+1,ju2+1,ku1  ,ku2  , &
       &       fp,gp,qp,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , cin, &
       &          fx,tx,ilo  ,ihi  ,jf1  ,jf2  ,klo  ,khi  , 2,1,3,ngrid)
  ! Save flux in output array
  do i=ilo,ihi
  do j=jf1,jf2
  do k=klo,khi
     do ivar=1,nvar
        do l=1,ngrid
           flux(l,i,j,k,ivar,2)=fx(l,i,j,k,ivar)*dt/dy
        end do
     end do
     do ivar=1,nmat+1
        do l=1,ngrid
           tmp (l,i,j,k,ivar,2)=tx(l,i,j,k,ivar)*dt/dy
        end do
     end do
  end do
  end do
  end do
#endif

  ! Solve for 1D flux in Z direction
#if NDIM>2
  call cmpgdnv(fm,gm,qm,iu1  ,iu2  ,ju1  ,ju2  ,ku1+1,ku2+1, &
       &       fp,gp,qp,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , cin, &
       &          fx,tx,ilo  ,ihi  ,jlo  ,jhi  ,kf1  ,kf2  , 3,1,2,ngrid)
  ! Save flux in output array
  do i=ilo,ihi
  do j=jlo,jhi
  do k=kf1,kf2
     do ivar=1,nvar
        do l=1,ngrid
           flux(l,i,j,k,ivar,3)=fx(l,i,j,k,ivar)*dt/dz
        end do
     end do
     do ivar=1,nmat+1
        do l=1,ngrid
           tmp (l,i,j,k,ivar,3)=tx(l,i,j,k,ivar)*dt/dz
        end do
     end do
  end do
  end do
  end do
#endif

end subroutine unsplit
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine ctoprim(uin,q,f,g,c,gravin,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none
  ! u: f_k, f_k.d_k, d.u, d.v, d.w, f_k.E_k           (Conservative variables)
  ! w: f_k,     d_k,   u,   v,   w, P_k, e_k, c_k     (Primitive variables + passive scalars)
  ! In this routine the primitive variables get broken down in to
  ! w = f_k + g_k + q + c_k ; f_k:Volume fractions, g_k:True densities, q:[u,v,w,P_k,e_k], c_k:Sound speeds
  integer ::ngrid
  logical ::inv
  real(dp)::dt
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::uin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim)::gravin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:npri)::q
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nmat)::f,g
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nmat)::c
  integer ::i, j, k, l, n, idim, imat
  real(dp),dimension(1:nvector),save::ekin
  real(dp),dimension(1:nvector),save::gg_mat,ee_mat,pp_mat,cc_mat,dtot
#if NVAR > NDIM + 3*NMAT
  integer::ipscal,npscal
#endif
  ! Convert to primitive variable
  do k = ku1, ku2
  do j = ju1, ju2
  do i = iu1, iu2

    ! Volume fraction and fluid density
    dtot(1:ngrid)=0
    do imat = 1,nmat
      do l = 1,ngrid
        ! Volume fractions
        f(l,i,j,k,imat) = uin(l,i,j,k,imat)
        ! True densities
        g(l,i,j,k,imat) = uin(l,i,j,k,nmat+imat)/f(l,i,j,k,imat)
        ! Total density
        dtot(l) = dtot(l) + uin(l,i,j,k,nmat+imat)
      end do
    end do

    ! Velocity components
    do idim=1,ndim
      do l=1,ngrid
        q(l,i,j,k,idim) = uin(l,i,j,k,2*nmat+idim)/dtot(l)
      end do
    end do

    ! Compute velocity and specific kinetic energy
    ekin(1:ngrid)=0.0
    do idim = 1,ndim
      do l = 1,ngrid
        ekin(l) = ekin(l) + half*q(l,i,j,k,idim)**2
      end do
    end do

    ! Compute individual internal energies
    do imat=1,nmat
      do l=1,ngrid
        q(l,i,j,k,ndim+nmat+imat) = uin(l,i,j,k,2*nmat+ndim+imat)/f(l,i,j,k,imat) - g(l,i,j,k,imat)*ekin(l)
      end do
    end do

    ! Calculate the total speed of sound and the total pressure from the EOS

    inv=.false.
    do imat=1,nmat
      do l=1,ngrid
      gg_mat(l) = g(l,i,j,k,imat)
      ee_mat(l) = q(l,i,j,k,ndim+nmat+imat)
      if(ee_mat(l)<0)then
         ! write(*,*)'ctoprim',imat,l,ee_mat(l)
      endif
      end do
      ! Call eos routine
      call eos(gg_mat,ee_mat,pp_mat,cc_mat,imat,inv,ngrid)
      do l=1,ngrid
        ! Individual pressures
        q(l,i,j,k,ndim+imat) = pp_mat(l)
        ! Indivual speed of sound
        c(l,i,j,k,imat) = cc_mat(l)
      end do
    end do

    ! Gravity predictor step (had to comment it out...)
    do idim = 1,ndim
!      q(1:ngrid,i,j,k,idim) = q(1:ngrid,i,j,k,idim) + gravin(1:ngrid,i,j,k,idim)*dt*half
    end do

  end do
  end do
  end do

#if NVAR > NDIM + 3*NMAT
  ! Passive scalar for each fluid
  npscal = (nvar - ndim - 3*nmat) / nmat
  do imat = 1, nmat
  do ipscal = 1, npscal
     n = ndim + 3*nmat + npscal*(imat-1) + ipscal

     do k = ku1, ku2
     do j = ju1, ju2
     do i = iu1, iu2
        do l = 1, ngrid
           q(l,i,j,k,n-nmat) = uin(l,i,j,k,n)/uin(l,i,j,k,nmat+imat)
           if(q(l,i,j,k,n-nmat)<0)then
              ! write(*,*)'ctoprim negative entropy',imat,ipscal,q(l,i,j,k,n-nmat)
           endif
        end do
     end do
     end do
     end do

  end do
  end do
#endif

end subroutine ctoprim
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine trace1d(qin,fin,gin,dq,df,dg,cin,qm,qp,fm,fp,gm,gp,dx,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer ::ngrid
  real(dp)::dx, dt

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:npri,1:ndim)::dq
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:npri,1:ndim)::qm
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:npri,1:ndim)::qp
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:npri)::qin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nmat,1:ndim)::df
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nmat,1:ndim)::fm
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nmat,1:ndim)::fp
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nmat)::fin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nmat,1:ndim)::dg
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nmat,1:ndim)::gm
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nmat,1:ndim)::gp
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nmat)::gin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nmat)::cin

  ! Local variables
  integer ::i, j, k, l, n, imat
  integer ::ilo,ihi,jlo,jhi,klo,khi
  integer ::iu
  real(dp)::dtdx, u ,dux, su0, curvilin
  real(dp)::r,dptotx
  real(dp),dimension(1:nmat)::f,g,p,e,csq
  real(dp),dimension(1:nmat)::dfx,dgx,dpx,dex
  real(dp),dimension(1:nmat)::sf0,sg0,sp0,se0
#if NVAR > NDIM + 3*NMAT
  real(dp)::a,dax,sa0
#endif

  dtdx = dt/dx

  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)
  iu=1

  do k = klo, khi
  do j = jlo, jhi
  do i = ilo, ihi

     do l = 1, ngrid
        r      = 0.0
        dptotx = 0.0
        ! Cell centered values
        do imat   = 1,nmat
        f(imat)   = fin(l,i,j,k,imat)
        g(imat)   = gin(l,i,j,k,imat)
        p(imat)   = qin(l,i,j,k,ndim+imat)
        e(imat)   = qin(l,i,j,k,ndim+nmat+imat)
        csq(imat) = cin(l,i,j,k,imat)**2
        r         = r + g(imat)*f(imat)
        end do
        u         = qin(l,i,j,k,iu)

        ! TVD slopes in X direction
        do imat     = 1,nmat
        dfx(imat)   = df(l,i,j,k,imat           ,1)
        dgx(imat)   = dg(l,i,j,k,imat           ,1)
        dpx(imat)   = dq(l,i,j,k,ndim+imat      ,1)
        dex(imat)   = dq(l,i,j,k,ndim+nmat+imat ,1)
        dptotx      = dptotx + dpx(imat)*f(imat) + dfx(imat)*p(imat)
        end do
        dux         = dq(l,i,j,k,iu             ,1)

        ! Source terms (including transverse derivatives and geometrical terms)
        do imat   = 1,nmat
        sf0(imat) = -u*dfx(imat)
        sg0(imat) = -u*dgx(imat) - (dux)*g(imat)
        sp0(imat) = -u*dpx(imat) - (dux)*g(imat)*csq(imat)
        se0(imat) = -u*dex(imat) - (dux)*(e(imat)+p(imat))
        end do
        su0       = -u*dux       - dptotx/r

        ! Right state
        do imat = 1,nmat
        fp(l,i,j,k,imat           ,1) = f(imat) -half*dfx(imat) +sf0(imat)*dtdx*half
        gp(l,i,j,k,imat           ,1) = g(imat) -half*dgx(imat) +sg0(imat)*dtdx*half
        qp(l,i,j,k,ndim+imat      ,1) = p(imat) -half*dpx(imat) +sp0(imat)*dtdx*half
        qp(l,i,j,k,ndim+nmat+imat ,1) = e(imat) -half*dex(imat) +se0(imat)*dtdx*half
        end do
        qp(l,i,j,k,iu             ,1) = u       -half*dux       +su0      *dtdx*half

        ! Left state
        do imat = 1,nmat
        fm(l,i,j,k,imat           ,1) = f(imat) +half*dfx(imat) +sf0(imat)*dtdx*half
        gm(l,i,j,k,imat           ,1) = g(imat) +half*dgx(imat) +sg0(imat)*dtdx*half
        qm(l,i,j,k,ndim+imat      ,1) = p(imat) +half*dpx(imat) +sp0(imat)*dtdx*half
        qm(l,i,j,k,ndim+nmat+imat ,1) = e(imat) +half*dex(imat) +se0(imat)*dtdx*half
        end do
        qm(l,i,j,k,iu             ,1) = u       +half*dux       +su0      *dtdx*half

     end do

  end do
  end do
  end do

#if NVAR > NDIM + 3*NMAT
  ! Passive scalars
  do n = ndim+2*nmat+1,nvar-nmat
     do k = klo, khi
     do j = jlo, jhi
     do i = ilo, ihi
        do l = 1, ngrid
           a   = qin(l,i,j,k,n)     ! Cell centered values
           u   = qin(l,i,j,k,iu)    ! x-velocity
           dax = dq(l,i,j,k,n,1)    ! TVD slopes
           sa0 = -u*dax             ! Source terms
           qp(l,i,j,k,n,1) = a - half*dax + sa0*dtdx*half   ! Right state
           qm(l,i,j,k,n,1) = a + half*dax + sa0*dtdx*half   ! Left state
        end do
     end do
     end do
     end do
  end do
#endif

end subroutine trace1d
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#if NDIM>1
subroutine trace2d(qin,fin,gin,dq,df,dg,cin,qm,qp,fm,fp,gm,gp,dx,dy,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer ::ngrid
  real(dp)::dx, dy, dt

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:npri,1:ndim)::dq
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:npri,1:ndim)::qm
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:npri,1:ndim)::qp
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:npri)::qin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nmat,1:ndim)::df
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nmat,1:ndim)::fm
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nmat,1:ndim)::fp
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nmat)::fin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nmat,1:ndim)::dg
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nmat,1:ndim)::gm
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nmat,1:ndim)::gp
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nmat)::gin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nmat)::cin

  ! Local variables
  integer ::i, j, k, l, n, imat
  integer ::ilo,ihi,jlo,jhi,klo,khi
  integer ::iu, iv
  real(dp)::dtdx, dtdy, u, v
  real(dp)::dux, dvx
  real(dp)::duy, dvy
  real(dp)::su0, sv0
  real(dp)::r,dptotx,dptoty
  real(dp),dimension(1:nmat)::f,g,p,e,csq
  real(dp),dimension(1:nmat)::sf0,sg0,sp0,se0
  real(dp),dimension(1:nmat)::dfx,dgx,dpx,dex
  real(dp),dimension(1:nmat)::dfy,dgy,dpy,dey
  real(dp)::smalle
#if NVAR > NDIM + 3*NMAT
  real(dp)::a,dax,day,sa0
#endif

  smalle = smallr**3
  dtdx = dt/dx
  dtdy = dt/dy

  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)
  iu=1; iv=ndim

  do k = klo, khi
  do j = jlo, jhi
  do i = ilo, ihi

     do l = 1, ngrid

        ! Cell centered values
        r = 0.0
        do imat = 1,nmat
           f(imat) = fin(l,i,j,k,imat)
           g(imat) = gin(l,i,j,k,imat)
           p(imat) = qin(l,i,j,k,ndim+imat)
           e(imat) = qin(l,i,j,k,ndim+nmat+imat)
           csq(imat) = cin(l,i,j,k,imat)**2
           r = r + g(imat)*f(imat)
        end do
        u = qin(l,i,j,k,iu)
        v = qin(l,i,j,k,iv)

        ! TVD slopes in the x-direction
        dptotx = 0.0
        do imat = 1,nmat
           dfx(imat) = df(l,i,j,k,imat,1)
           dgx(imat) = dg(l,i,j,k,imat,1)
           dpx(imat) = dq(l,i,j,k,ndim+imat,1)
           dex(imat) = dq(l,i,j,k,ndim+nmat+imat,1)
           dptotx = dptotx + dpx(imat)*f(imat) + dfx(imat)*p(imat)
        end do
        dux = dq(l,i,j,k,iu,1)
        dvx = dq(l,i,j,k,iv,1)

        ! TVD slopes in the y-direction
        dptoty = 0.0
        do imat     = 1,nmat
           dfy(imat)   = df(l,i,j,k,imat,2)
           dgy(imat)   = dg(l,i,j,k,imat,2)
           dpy(imat)   = dq(l,i,j,k,ndim+imat,2)
           dey(imat)   = dq(l,i,j,k,ndim+nmat+imat,2)
           dptoty      = dptoty + dpy(imat)*f(imat) + dfy(imat)*p(imat)
        end do
        duy         = dq(l,i,j,k,iu,2)
        dvy         = dq(l,i,j,k,iv,2)

        ! Source terms (with transverse derivatives)
        do imat   = 1,nmat
           sf0(imat) = -u*dfx(imat) -v*dfy(imat)
           sg0(imat) = -u*dgx(imat) -v*dgy(imat) -(dux+dvy)*g(imat)
           sp0(imat) = -u*dpx(imat) -v*dpy(imat) -(dux+dvy)*g(imat)*csq(imat)
           se0(imat) = -u*dex(imat) -v*dey(imat) -(dux+dvy)*(e(imat)+p(imat))
        end do
        su0 = -u*dux -v*duy -(dptotx)/max(r,smallr)
        sv0 = -u*dvx -v*dvy -(dptoty)/max(r,smallr)

        ! Right state at left interface
        do imat = 1,nmat

           fp(l,i,j,k,imat          ,1)  = f(imat) -half*dfx(imat) +sf0(imat)*dtdx*half
           gp(l,i,j,k,imat          ,1)  = g(imat) -half*dgx(imat) +sg0(imat)*dtdx*half
           qp(l,i,j,k,ndim+imat     ,1)  = p(imat) -half*dpx(imat) +sp0(imat)*dtdx*half
           qp(l,i,j,k,ndim+nmat+imat,1)  = e(imat) -half*dex(imat) +se0(imat)*dtdx*half

           if(fp(l,i,j,k,imat,1)<smallf)fp(l,i,j,k,imat,1)=f(imat)
           if(gp(l,i,j,k,imat,1)<smallr)gp(l,i,j,k,imat,1)=g(imat)
           if(qp(l,i,j,k,ndim+imat,1)<smalle)qp(l,i,j,k,ndim+imat,1)=p(imat)
           if(qp(l,i,j,k,ndim+nmat+imat,1)<smalle)qp(l,i,j,k,nmat+imat,1)=e(imat)

        end do
        qp(l,i,j,k,iu             ,1)  = u       -half*dux       +su0      *dtdx*half
        qp(l,i,j,k,iv             ,1)  = v       -half*dvx       +sv0      *dtdx*half

        ! Left state at right interface
        do imat = 1,nmat

           fm(l,i,j,k,imat          ,1)  = f(imat) +half*dfx(imat) +sf0(imat)*dtdx*half
           gm(l,i,j,k,imat          ,1)  = g(imat) +half*dgx(imat) +sg0(imat)*dtdx*half
           qm(l,i,j,k,ndim+imat     ,1)  = p(imat) +half*dpx(imat) +sp0(imat)*dtdx*half
           qm(l,i,j,k,ndim+nmat+imat,1)  = e(imat) +half*dex(imat) +se0(imat)*dtdx*half

           if(fm(l,i,j,k,imat,1)<smallf)fm(l,i,j,k,imat,1)=f(imat)
           if(gm(l,i,j,k,imat,1)<smallr)gm(l,i,j,k,imat,1)=g(imat)
           if(qm(l,i,j,k,ndim+imat,1)<smallr)qm(l,i,j,k,ndim+imat,1)=p(imat)
           if(qm(l,i,j,k,ndim+nmat+imat,1)<smallr)qm(l,i,j,k,ndim+nmat+imat,1)=e(imat)

        end do
        qm(l,i,j,k,iu             ,1)  = u       +half*dux       +su0      *dtdx*half
        qm(l,i,j,k,iv             ,1)  = v       +half*dvx       +sv0      *dtdx*half

        ! Top state at bottom interface
        do imat = 1,nmat

           fp(l,i,j,k,imat          ,2)  = f(imat) -half*dfy(imat) +sf0(imat)*dtdy*half
           gp(l,i,j,k,imat          ,2)  = g(imat) -half*dgy(imat) +sg0(imat)*dtdy*half
           qp(l,i,j,k,ndim+imat     ,2)  = p(imat) -half*dpy(imat) +sp0(imat)*dtdy*half
           qp(l,i,j,k,ndim+nmat+imat,2)  = e(imat) -half*dey(imat) +se0(imat)*dtdy*half

           if(fp(l,i,j,k,imat,2)<smallf)fp(l,i,j,k,imat,2)=f(imat)
           if(gp(l,i,j,k,imat,2)<smallr)gp(l,i,j,k,imat,2)=g(imat)
           if(qp(l,i,j,k,ndim+imat,2)<smalle)qp(l,i,j,k,ndim+imat,2)=p(imat)
           if(qp(l,i,j,k,ndim+nmat+imat,2)<smalle)qp(l,i,j,k,ndim+nmat+imat,2)=e(imat)

        end do
        qp(l,i,j,k,iu             ,2)  = u       -half*duy       +su0      *dtdy*half
        qp(l,i,j,k,iv             ,2)  = v       -half*dvy       +sv0      *dtdy*half

        ! Bottom state at top interface
        do imat = 1,nmat

           fm(l,i,j,k,imat          ,2)  = f(imat) +half*dfy(imat) +sf0(imat)*dtdy*half
           gm(l,i,j,k,imat          ,2)  = g(imat) +half*dgy(imat) +sg0(imat)*dtdy*half
           qm(l,i,j,k,ndim+imat     ,2)  = p(imat) +half*dpy(imat) +sp0(imat)*dtdy*half
           qm(l,i,j,k,ndim+nmat+imat,2)  = e(imat) +half*dey(imat) +se0(imat)*dtdy*half

           if(fm(l,i,j,k,imat,2)<smallf)fm(l,i,j,k,imat,2)=f(imat)
           if(gm(l,i,j,k,imat,2)<smallr)gm(l,i,j,k,imat,2)=g(imat)
           if(qm(l,i,j,k,ndim+imat,2)<smalle)qm(l,i,j,k,ndim+imat,2)=p(imat)
           if(qm(l,i,j,k,ndim+nmat+imat,2)<smalle)qm(l,i,j,k,ndim+nmat+imat,2)=e(imat)

        end do
        qm(l,i,j,k,iu             ,2)  = u       +half*duy       +su0      *dtdy*half
        qm(l,i,j,k,iv             ,2)  = v       +half*dvy       +sv0      *dtdy*half

     end do

  end do
  end do
  end do

#if NVAR > NDIM + 3*NMAT
  ! Passive scalars
  do n = ndim+2*nmat+1,nvar-nmat
     do k = klo, khi
     do j = jlo, jhi
     do i = ilo, ihi
        do l = 1, ngrid
           a   = qin(l,i,j,k,n)     ! Cell centered values
           u   = qin(l,i,j,k,iu)    ! x-velocity
           v   = qin(l,i,j,k,iv)    ! y-velocity
           dax = dq(l,i,j,k,n,1)    ! TVD slopes
           day = dq(l,i,j,k,n,2)
           sa0 = -u*dax-v*day       ! Source terms
           qp(l,i,j,k,n,1) = a - half*dax + sa0*dtdx*half   ! Right state
           qm(l,i,j,k,n,1) = a + half*dax + sa0*dtdx*half   ! Left state
           qp(l,i,j,k,n,2) = a - half*day + sa0*dtdy*half   ! Top state
           qm(l,i,j,k,n,2) = a + half*day + sa0*dtdy*half   ! Bottom state
           if(a<0)then
              write(*,*)'negative entropy',n,ndim+2*nmat,a
           endif
           if(qp(l,i,j,k,n,1)<0)qp(l,i,j,k,n,1)=a
           if(qm(l,i,j,k,n,1)<0)qm(l,i,j,k,n,1)=a
           if(qp(l,i,j,k,n,2)<0)qp(l,i,j,k,n,2)=a
           if(qm(l,i,j,k,n,2)<0)qm(l,i,j,k,n,2)=a
        end do
     end do
     end do
     end do
  end do
#endif

end subroutine trace2d
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#if NDIM>2
subroutine trace3d(qin,fin,gin,dq,df,dg,cin,qm,qp,fm,fp,ym,yp,dx,dy,dz,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer ::ngrid
  real(dp)::dx, dy, dz, dt

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:npri,1:ndim+1)::dq
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:npri,1:ndim+1)::qm
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:npri,1:ndim+1)::qp
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:npri)::qin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nmat,1:ndim)::df
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nmat,1:ndim)::fm
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nmat,1:ndim)::fp
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nmat)::fin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nmat,1:ndim)::dg
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nmat,1:ndim)::gm
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nmat,1:ndim)::gp
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nmat)::gin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nmat)::cin

  ! declare local variables
  integer ::i, j, k, l, n, imat
  integer ::ilo,ihi,jlo,jhi,klo,khi
  integer ::iu, iv, iw
  real(dp)::dtdx, dtdy, dtdz, u, v, w
  real(dp)::dux, dvx, dwx
  real(dp)::duy, dvy, dwy
  real(dp)::duz, dvz, dwz
  real(dp)::sr0, su0, sv0, sw0
  real(dp)::r,dptotx,dptoty,dptotz
  real(dp),dimension(1:nmat)::f,g,p,e,csq
  real(dp),dimension(1:nmat)::dfx,dgx,dpx,dex
  real(dp),dimension(1:nmat)::dfy,dgy,dpy,dey
  real(dp),dimension(1:nmat)::dfz,dgz,dpz,dez
  real(dp),dimension(1:nmat)::sf0,sg0,sp0, se0

  dtdx = dt/dx
  dtdy = dt/dy
  dtdz = dt/dz

  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)

  iu=1; iv=2; iw=ndim

  do k = klo, khi
  do j = jlo, jhi
  do i = ilo, ihi

     do l = 1, ngrid
        r      = 0.0
        dptotx = 0.0
        dptoty = 0.0
        dptotz = 0.0

        ! Cell centered values
        do imat   = 1,nmat
        f(imat)   = fin(l,i,j,k,imat)
        g(imat)   = gin(l,i,j,k,imat)
        p(imat)   = qin(l,i,j,k,ndim+imat)
        e(imat)   = qin(l,i,j,k,ndim+nmat+imat)
        csq(imat) = cin(l,i,j,k,imat)**2
        r         = r + g(imat)*f(imat)
        end do
        u         = qin(l,i,j,k,iu)
        v         = qin(l,i,j,k,iv)
        w         = qin(l,i,j,k,iw)

        ! TVD slopes in all directions
        do imat     = 1,nmat
        dfx(imat)   = df(l,i,j,k,imat           ,1)
        dgx(imat)   = dg(l,i,j,k,imat           ,1)
        dpx(imat)   = dq(l,i,j,k,ndim+imat      ,1)
        dex(imat)   = dq(l,i,j,k,ndim+nmat+imat ,1)
        dptotx      = dptotx + dpx(imat)*f(imat) + dfx(imat)*p(imat)
        end do
        dux         = dq(l,i,j,k,iu             ,1)
        dvx         = dq(l,i,j,k,iv             ,1)
        dwx         = dq(l,i,j,k,iw             ,1)
        do imat     = 1,nmat
        dfy(imat)   = df(l,i,j,k,imat           ,2)
        dgy(imat)   = dg(l,i,j,k,imat           ,2)
        dpy(imat)   = dq(l,i,j,k,ndim+imat      ,2)
        dey(imat)   = dq(l,i,j,k,ndim+nmat+imat ,2)
        dptoty      = dptoty + dpy(imat)*f(imat) + dfy(imat)*p(imat)
        end do
        duy         = dq(l,i,j,k,iu             ,2)
        dvy         = dq(l,i,j,k,iv             ,2)
        dwy         = dq(l,i,j,k,iw             ,2)
        do imat     = 1,nmat
        dfz(imat)   = df(l,i,j,k,imat           ,3)
        dgz(imat)   = dg(l,i,j,k,imat           ,3)
        dpz(imat)   = dq(l,i,j,k,ndim+imat      ,3)
        dez(imat)   = dq(l,i,j,k,ndim+nmat+imat ,3)
        dptotz      = dptotz + dpz(imat)*f(imat) + dfz(imat)*p(imat)
        end do
        duz         = dq(l,i,j,k,iu             ,3)
        dvz         = dq(l,i,j,k,iv             ,3)
        dwz         = dq(l,i,j,k,iw             ,3)

        ! Source terms (with transverse derivatives)
        do imat   = 1,nmat
        sf0(imat) = -u*dfx(imat) -v*dfy(imat) -w*dfz(imat)
        sg0(imat) = -u*dgx(imat) -v*dgy(imat) -w*dgz(imat) - (dux+dvy+dwz)*g(imat)
        sp0(imat) = -u*dpx(imat) -v*dpy(imat) -w*dpz(imat) - (dux+dvy+dwz)*g(imat)*csq(imat)
        se0(imat) = -u*dex(imat) -v*dey(imat) -w*dez(imat) - (dux+dvy+dwz)*(e(imat)+p(imat))
        end do
        su0       = -u*dux       -v*duy       -w*duz       - dptotx/max(r,smallr)
        sv0       = -u*dvx       -v*dvy       -w*dvz       - dptoty/max(r,smallr)
        sw0       = -u*dwx       -v*dwy       -w*dwz       - dptotz/max(r,smallr)

        ! Right state at left interface
        do imat = 1,nmat
        fp(l,i,j,k,imat           ,1) = f(imat) -half*dfx(imat) +sf0(imat)*dtdx*half
        gp(l,i,j,k,imat           ,1) = g(imat) -half*dgx(imat) +sgg(imat)*dtdx*half
        qp(l,i,j,k,ndim+imat      ,1) = p(imat) -half*dpx(imat) +sp0(imat)*dtdx*half
        qp(l,i,j,k,ndim+nmat+imat ,1) = e(imat) -half*dex(imat) +se0(imat)*dtdx*half
        end do
        qp(l,i,j,k,iu             ,1) = u       -half*dux       +su0      *dtdx*half
        qp(l,i,j,k,iv             ,1) = v       -half*dvx       +sv0      *dtdx*half
        qp(l,i,j,k,iw             ,1) = w       -half*dwx       +sw0      *dtdx*half

        ! Left state at right interface
        do imat = 1,nmat
        fm(l,i,j,k,imat           ,1) = f(imat) +half*dfx(imat) +sf0(imat)*dtdx*half
        gm(l,i,j,k,imat           ,1) = g(imat) +half*dgx(imat) +sg0(imat)*dtdx*half
        qm(l,i,j,k,ndim+imat      ,1) = p(imat) +half*dpx(imat) +sp0(imat)*dtdx*half
        qm(l,i,j,k,ndim+nmat+imat ,1) = e(imat) +half*dex(imat) +se0(imat)*dtdx*half
        end do
        qm(l,i,j,k,iu             ,1) = u       +half*dux       +su0      *dtdx*half
        qm(l,i,j,k,iv             ,1) = v       +half*dvx       +sv0      *dtdx*half
        qm(l,i,j,k,iw             ,1) = w       +half*dwx       +sw0      *dtdx*half

        ! Top state at bottom interface
        do imat = 1,nmat
        fp(l,i,j,k,imat           ,2) = f(imat) -half*dfy(imat) +sf0(imat)*dtdy*half
        gp(l,i,j,k,imat           ,2) = g(imat) -half*dgy(imat) +sg0(imat)*dtdy*half
        qp(l,i,j,k,ndim+imat      ,2) = p(imat) -half*dpy(imat) +sp0(imat)*dtdy*half
        qp(l,i,j,k,ndim+nmat+imat ,2) = e(imat) -half*dey(imat) +se0(imat)*dtdy*half
        end do
        qp(l,i,j,k,iu             ,2) = u       -half*duy       +su0      *dtdy*half
        qp(l,i,j,k,iv             ,2) = v       -half*dvy       +sv0      *dtdy*half
        qp(l,i,j,k,iw             ,2) = w       -half*dwy       +sw0      *dtdy*half

        ! Bottom state at top interface
        do imat = 1,nmat
        fm(l,i,j,k,imat           ,2) = f(imat) +half*dfy(imat) +sf0(imat)*dtdy*half
        gm(l,i,j,k,imat           ,2) = g(imat) +half*dgy(imat) +sg0(imat)*dtdy*half
        qm(l,i,j,k,ndim+imat      ,2) = p(imat) +half*dpy(imat) +sp0(imat)*dtdx*half
        qm(l,i,j,k,ndim+nmat+imat ,2) = e(imat) +half*dey(imat) +se0(imat)*dtdx*half
        end do
        qm(l,i,j,k,iu             ,2) = u       +half*duy       +su0      *dtdy*half
        qm(l,i,j,k,iv             ,2) = v       +half*dvy       +sv0      *dtdy*half
        qm(l,i,j,k,iw             ,2) = w       +half*dwy       +sw0      *dtdy*half

        ! Back state at front interface
        do imat = 1,nmat
        fp(l,i,j,k,imat           ,3) = f(imat) -half*dfz(imat) +sf0(imat)*dtdz*half
        gp(l,i,j,k,imat           ,3) = g(imat) -half*dgz(imat) +sg0(imat)*dtdz*half
        qp(l,i,j,k,ndim+imat      ,3) = p(imat) -half*dpz(imat) +sp0(imat)*dtdz*half
        qp(l,i,j,k,ndim+nmat+imat ,3) = e(imat) -half*dez(imat) +se0(imat)*dtdz*half
        end do
        qp(l,i,j,k,iu             ,3) = u       -half*duz       +su0      *dtdz*half
        qp(l,i,j,k,iv             ,3) = v       -half*dvz       +sv0      *dtdz*half
        qp(l,i,j,k,iw             ,3) = w       -half*dwz       +sw0      *dtdz*half

        ! Front state at back interface
        do imat = 1,nmat
        fm(l,i,j,k,imat           ,3) = f(imat) +half*dfz(imat) +sf0(imat)*dtdz*half
        gm(l,i,j,k,imat           ,3) = g(imat) +half*dgz(imat) +sg0(imat)*dtdz*half
        qm(l,i,j,k,ndim+imat      ,3) = p(imat) +half*dpz(imat) +sp0(imat)*dtdz*half
        qm(l,i,j,k,ndim+nmat+imat ,3) = e(imat) +half*dez(imat) +se0(imat)*dtdz*half
        end do
        qm(l,i,j,k,iu             ,3) = u       +half*duz       +su0      *dtdz*half
        qm(l,i,j,k,iv             ,3) = v       +half*dvz       +sv0      *dtdz*half
        qm(l,i,j,k,iw             ,3) = w       +half*dwz       +sw0      *dtdz*half

     end do
  end do
  end do
  end do

end subroutine trace3d
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmpgdnv(fm,gm,qm,im1,im2,jm1,jm2,km1,km2, &
     &             fp,gp,qp,ip1,ip2,jp1,jp2,kp1,kp2, c, &
     &              flx,tmp,ilo,ihi,jlo,jhi,klo,khi, &
     &                      ln,lt1,lt2,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer ::ngrid
  integer ::ln,lt1,lt2
  integer ::im1,im2,jm1,jm2,km1,km2
  integer ::ip1,ip2,jp1,jp2,kp1,kp2
  integer ::ilo,ihi,jlo,jhi,klo,khi
  real(dp),dimension(1:nvector,im1:im2,jm1:jm2,km1:km2,1:npri,1:ndim)::qm
  real(dp),dimension(1:nvector,ip1:ip2,jp1:jp2,kp1:kp2,1:npri,1:ndim)::qp
  real(dp),dimension(1:nvector,im1:im2,jm1:jm2,km1:km2,1:nmat,1:ndim)::fm
  real(dp),dimension(1:nvector,ip1:ip2,jp1:jp2,kp1:kp2,1:nmat,1:ndim)::fp
  real(dp),dimension(1:nvector,im1:im2,jm1:jm2,km1:km2,1:nmat,1:ndim)::gm
  real(dp),dimension(1:nvector,ip1:ip2,jp1:jp2,kp1:kp2,1:nmat,1:ndim)::gp
  real(dp),dimension(1:nvector,ip1:ip2,jp1:jp2,kp1:kp2,1:nvar)::flx
  real(dp),dimension(1:nvector,ip1:ip2,jp1:jp2,kp1:kp2,1:nmat+1)::tmp
  real(dp),dimension(1:nvector,ip1:ip2,jp1:jp2,kp1:kp2,1:nmat)::c

  ! local variables
  integer ::i, j, k, n, l, idim, jdim, imat, ivar
  real(dp),dimension(1:nvector,1:npri),save::qleft,qright
  real(dp),dimension(1:nvector,1:nmat),save::fleft,fright
  real(dp),dimension(1:nvector,1:nmat),save::gleft,gright
  real(dp),dimension(1:nvector,1:nvar),save::fgdnv
  real(dp),dimension(1:nvector),save::ugdnv
  real(dp),dimension(1:nvector,1:nmat),save::egdnv
  real(dp),dimension(1:nvector,1:nmat),save::cleft,cright

  idim= ln

  do k = klo, khi
  do j = jlo, jhi
  do i = ilo, ihi

     ! Normal velocity
     do l = 1, ngrid
        qleft (l,1) = qm(l,i,j,k,ln,idim)
        qright(l,1) = qp(l,i,j,k,ln,idim)
     end do

     ! Individual pressures
     do imat=1,nmat
        do l = 1, ngrid
           qleft (l,ndim+imat) = qm(l,i,j,k,ndim+imat,idim)
           qright(l,ndim+imat) = qp(l,i,j,k,ndim+imat,idim)
        end do
     end do

     ! Internal energies
     do imat=1,nmat
        do l=1,ngrid
           qleft (l,ndim+nmat+imat) = qm(l,i,j,k,ndim+nmat+imat,idim)
           qright(l,ndim+nmat+imat) = qp(l,i,j,k,ndim+nmat+imat,idim)
        end do
     end do
#if NVAR > NDIM + 3*NMAT
     ! Passive scalars
     do n=ndim+2*nmat+1, nvar-nmat
        do l=1,ngrid
           qleft (l,n) = qm(l,i,j,k,n,idim)
           qright(l,n) = qp(l,i,j,k,n,idim)
        end do
     end do
#endif
     ! Tangential velocity 1
#if NDIM>1
     do l = 1, ngrid
        qleft (l,2) = qm(l,i,j,k,lt1,idim)
        qright(l,2) = qp(l,i,j,k,lt1,idim)
     end do
#endif
     ! Tangential velocity 2
#if NDIM>2
     do l = 1, ngrid
        qleft (l,3) = qm(l,i,j,k,lt2,idim)
        qright(l,3) = qp(l,i,j,k,lt2,idim)
     end do
#endif
     ! Volume fractions
     do imat = 1,nmat
        do l = 1, ngrid
           fleft (l,imat) = fm(l,i,j,k,imat,idim)
           fright(l,imat) = fp(l,i,j,k,imat,idim)
        end do
     end do

     ! Physical densities
     do imat = 1,nmat
        do l = 1, ngrid
           gleft (l,imat) = gm(l,i,j,k,imat,idim)
           gright(l,imat) = gp(l,i,j,k,imat,idim)
        end do
     end do

     ! Sound speed
     if (ln==1)then
      do imat=1,nmat
        do l = 1,ngrid
           cleft (l,imat) = c(l,i-1,j,k,imat)
           cright(l,imat) = c(l,i  ,j,k,imat)
        end do
      end do
     else if(ln==2) then
      do imat=1,nmat
       do l = 1,ngrid
           cleft (l,imat) = c(l,i,j-1,k,imat)
           cright(l,imat) = c(l,i,j  ,k,imat)
        end do
      end do
     else
      do imat=1,nmat
        do l = 1,ngrid
           cleft (l,imat) = c(l,i,j,k-1,imat)
           cright(l,imat) = c(l,i,j,k  ,imat)
        end do
      end do
     end if

     ! Solve Riemann problem
     call riemann_hllc(fleft,fright,gleft,gright,qleft,qright,cleft,cright,fgdnv,ugdnv,egdnv,ngrid)

     ! Store fluxes
     do ivar = 1,nvar
        do l = 1, ngrid
           flx(l,i,j,k,ivar) = fgdnv(l,ivar)
        end do
     end do

     ! We need to reassign the positions of the velocity for the appropriate direction

     ! Normal velocity
     do l = 1, ngrid
        flx(l,i,j,k,2*nmat+ln)  = fgdnv(l,2*nmat+1)
     end do
     ! Tangential velocity 2
#if NDIM>1
     do l = 1, ngrid
        flx(l,i,j,k,2*nmat+lt1) = fgdnv(l,2*nmat+2)
     end do
#endif
     ! Tangential velocity 2
#if NDIM>2
     do l = 1, ngrid
        flx(l,i,j,k,2*nmat+lt2) = fgdnv(l,2*nmat+3)
     end do
#endif

     ! Store u for div(u)
     do l=1,ngrid
        tmp(l,i,j,k,1)=ugdnv(l)
     end do

     ! Store u for div(u)
     do imat=1,nmat
        do l=1,ngrid
           tmp(l,i,j,k,1+imat)=egdnv(l,imat)
        end do
     end do

  end do
  end do
  end do

end subroutine cmpgdnv
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine qslope(q,dq,dx,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer::ngrid
  real(dp)::dx,dt
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:npri)::q
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:npri,1:ndim)::dq

  ! local arrays
  integer::i, j, k, l, n
  real(dp)::dsgn, dlim, dcen, dlft, drgt, slop
  real(dp)::dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr
  real(dp)::dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl
  real(dp)::dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm
  real(dp)::dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr
  real(dp)::vmin,vmax,dfx,dfy,dfz,dff
  integer::ilo,ihi,jlo,jhi,klo,khi

  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)

  if(slope_type==0)then
     dq=zero
     return
  end if

#if NDIM==1
  do k = klo, khi
  do j = jlo, jhi
  do i = ilo, ihi
     do n = 1, npri
        if(slope_type==1.or.slope_type==2)then  ! minmod or average
           do l = 1, ngrid
              dlft = slope_type*(q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
              drgt = slope_type*(q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
              dcen = half*(dlft+drgt)/slope_type
              dsgn = sign(one, dcen)
              slop = min(abs(dlft),abs(drgt))
              dlim = slop
              if((dlft*drgt)<=zero)dlim=zero
              dq(l,i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
           end do
        else if(slope_type==3)then ! superbee
           do l = 1, ngrid
              dcen = q(l,i,j,k,1)*dt/dx
              dlft = two/(one+dcen)*(q(l,i,j,k,n)-q(l,i-1,j,k,n))
              drgt = two/(one-dcen)*(q(l,i+1,j,k,n)-q(l,i,j,k,n))
              dsgn = sign(one, dlft)
              slop = min(abs(dlft),abs(drgt))
              dlim = slop
              if((dlft*drgt)<=zero)dlim=zero
              dq(l,i,j,k,n,1) = dsgn*dlim
           end do
        else if(slope_type==4)then ! ultrabee
           do l = 1, ngrid
              dcen = q(l,i,j,k,1)*dt/dx
              if(dcen>=0)then
                 dlft = two/(zero+dcen+1d-10)*(q(l,i,j,k,n)-q(l,i-1,j,k,n))
                 drgt = two/(one -dcen      )*(q(l,i+1,j,k,n)-q(l,i,j,k,n))
              else
                 dlft = two/(one +dcen      )*(q(l,i,j,k,n)-q(l,i-1,j,k,n))
                 drgt = two/(zero-dcen+1d-10)*(q(l,i+1,j,k,n)-q(l,i,j,k,n))
              endif
              dsgn = sign(one, dlft)
              slop = min(abs(dlft),abs(drgt))
              dlim = slop
              if((dlft*drgt)<=zero)dlim=zero
              dq(l,i,j,k,n,1) = dsgn*dlim
           end do
        else
           write(*,*)'Unknown slope type'
           stop
        end if
     end do
  end do
  end do
  end do
#endif

#if NDIM==2
  do k = klo, khi
  do j = jlo, jhi
  do i = ilo, ihi
     if(slope_type==1.or.slope_type==2)then  ! minmod or average
        do n = 1, npri
           ! slopes in first coordinate direction
           do l = 1, ngrid
              dlft = slope_type*(q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
              drgt = slope_type*(q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
              dcen = half*(dlft+drgt)/slope_type
              dsgn = sign(one, dcen)
              slop = min(abs(dlft),abs(drgt))
              dlim = slop
              if((dlft*drgt)<=zero)dlim=zero
              dq(l,i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
           end do
           ! slopes in second coordinate direction
           do l = 1, ngrid
              dlft = slope_type*(q(l,i,j  ,k,n) - q(l,i,j-1,k,n))
              drgt = slope_type*(q(l,i,j+1,k,n) - q(l,i,j  ,k,n))
              dcen = half*(dlft+drgt)/slope_type
              dsgn = sign(one,dcen)
              slop = min(abs(dlft),abs(drgt))
              dlim = slop
              if((dlft*drgt)<=zero)dlim=zero
              dq(l,i,j,k,n,2) = dsgn*min(dlim,abs(dcen))
           end do
        end do
     else if(slope_type==3)then ! positivity preserving 2d unsplit slope
        do n = 1, npri
           do l = 1, ngrid
              dfll = q(l,i-1,j-1,k,n)-q(l,i,j,k,n)
              dflm = q(l,i-1,j  ,k,n)-q(l,i,j,k,n)
              dflr = q(l,i-1,j+1,k,n)-q(l,i,j,k,n)
              dfml = q(l,i  ,j-1,k,n)-q(l,i,j,k,n)
              dfmm = q(l,i  ,j  ,k,n)-q(l,i,j,k,n)
              dfmr = q(l,i  ,j+1,k,n)-q(l,i,j,k,n)
              dfrl = q(l,i+1,j-1,k,n)-q(l,i,j,k,n)
              dfrm = q(l,i+1,j  ,k,n)-q(l,i,j,k,n)
              dfrr = q(l,i+1,j+1,k,n)-q(l,i,j,k,n)
              vmin = min(dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr)
              vmax = max(dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr)
              dfx  = half*(q(l,i+1,j,k,n)-q(l,i-1,j,k,n))
              dfy  = half*(q(l,i,j+1,k,n)-q(l,i,j-1,k,n))
              dff  = half*(abs(dfx)+abs(dfy))
              if(dff>zero)then
                 slop = min(one,min(abs(vmin),abs(vmax))/dff)
              else
                 slop = one
              endif
              dlim = slop
              dq(l,i,j,k,n,1) = dlim*dfx
              dq(l,i,j,k,n,2) = dlim*dfy
           end do
        end do
     else
        write(*,*)'Unknown slope type'
        stop
     endif
  end do
  end do
  end do
#endif

#if NDIM==3
  do k = klo, khi
  do j = jlo, jhi
  do i = ilo, ihi
     if(slope_type==1.or.slope_type==2)then  ! minmod or average
        do n = 1, npri
           ! slopes in first coordinate direction
           do l = 1, ngrid
              dlft = slope_type*(q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
              drgt = slope_type*(q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
              dcen = half*(dlft+drgt)/slope_type
              dsgn = sign(one, dcen)
              slop = min(abs(dlft),abs(drgt))
              dlim = slop
              if((dlft*drgt)<=zero)dlim=zero
              dq(l,i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
           end do
           ! slopes in second coordinate direction
           do l = 1, ngrid
              dlft = slope_type*(q(l,i,j  ,k,n) - q(l,i,j-1,k,n))
              drgt = slope_type*(q(l,i,j+1,k,n) - q(l,i,j  ,k,n))
              dcen = half*(dlft+drgt)/slope_type
              dsgn = sign(one,dcen)
              slop = min(abs(dlft),abs(drgt))
              dlim = slop
              if((dlft*drgt)<=zero)dlim=zero
              dq(l,i,j,k,n,2) = dsgn*min(dlim,abs(dcen))
           end do
           ! slopes in third coordinate direction
           do l = 1, ngrid
              dlft = slope_type*(q(l,i,j,k  ,n) - q(l,i,j,k-1,n))
              drgt = slope_type*(q(l,i,j,k+1,n) - q(l,i,j,k  ,n))
              dcen = half*(dlft+drgt)/slope_type
              dsgn = sign(one,dcen)
              slop = min(abs(dlft),abs(drgt))
              dlim = slop
              if((dlft*drgt)<=zero)dlim=zero
              dq(l,i,j,k,n,3) = dsgn*min(dlim,abs(dcen))
           end do
        end do
     else if(slope_type==3)then ! positivity preserving 3d unsplit slope
        do n = 1, npri
           do l = 1, ngrid
              dflll = q(l,i-1,j-1,k-1,n)-q(l,i,j,k,n)
              dflml = q(l,i-1,j  ,k-1,n)-q(l,i,j,k,n)
              dflrl = q(l,i-1,j+1,k-1,n)-q(l,i,j,k,n)
              dfmll = q(l,i  ,j-1,k-1,n)-q(l,i,j,k,n)
              dfmml = q(l,i  ,j  ,k-1,n)-q(l,i,j,k,n)
              dfmrl = q(l,i  ,j+1,k-1,n)-q(l,i,j,k,n)
              dfrll = q(l,i+1,j-1,k-1,n)-q(l,i,j,k,n)
              dfrml = q(l,i+1,j  ,k-1,n)-q(l,i,j,k,n)
              dfrrl = q(l,i+1,j+1,k-1,n)-q(l,i,j,k,n)
              dfllm = q(l,i-1,j-1,k  ,n)-q(l,i,j,k,n)
              dflmm = q(l,i-1,j  ,k  ,n)-q(l,i,j,k,n)
              dflrm = q(l,i-1,j+1,k  ,n)-q(l,i,j,k,n)
              dfmlm = q(l,i  ,j-1,k  ,n)-q(l,i,j,k,n)
              dfmmm = q(l,i  ,j  ,k  ,n)-q(l,i,j,k,n)
              dfmrm = q(l,i  ,j+1,k  ,n)-q(l,i,j,k,n)
              dfrlm = q(l,i+1,j-1,k  ,n)-q(l,i,j,k,n)
              dfrmm = q(l,i+1,j  ,k  ,n)-q(l,i,j,k,n)
              dfrrm = q(l,i+1,j+1,k  ,n)-q(l,i,j,k,n)
              dfllr = q(l,i-1,j-1,k+1,n)-q(l,i,j,k,n)
              dflmr = q(l,i-1,j  ,k+1,n)-q(l,i,j,k,n)
              dflrr = q(l,i-1,j+1,k+1,n)-q(l,i,j,k,n)
              dfmlr = q(l,i  ,j-1,k+1,n)-q(l,i,j,k,n)
              dfmmr = q(l,i  ,j  ,k+1,n)-q(l,i,j,k,n)
              dfmrr = q(l,i  ,j+1,k+1,n)-q(l,i,j,k,n)
              dfrlr = q(l,i+1,j-1,k+1,n)-q(l,i,j,k,n)
              dfrmr = q(l,i+1,j  ,k+1,n)-q(l,i,j,k,n)
              dfrrr = q(l,i+1,j+1,k+1,n)-q(l,i,j,k,n)
              vmin = min(dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl, &
                   &     dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm, &
                   &     dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr)
              vmax = max(dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl, &
                   &     dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm, &
                   &     dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr)
              dfx  = half*(q(l,i+1,j,k,n)-q(l,i-1,j,k,n))
              dfy  = half*(q(l,i,j+1,k,n)-q(l,i,j-1,k,n))
              dfz  = half*(q(l,i,j,k+1,n)-q(l,i,j,k-1,n))
              dff  = half*(abs(dfx)+abs(dfy)+abs(dfz))
              if(dff>zero)then
                 slop = min(one,min(abs(vmin),abs(vmax))/dff)
              else
                 slop = one
              endif
              dlim = slop
              dq(l,i,j,k,n,1) = dlim*dfx
              dq(l,i,j,k,n,2) = dlim*dfy
              dq(l,i,j,k,n,3) = dlim*dfz
           end do
        end do
     else
        write(*,*)'Unknown slope type'
        stop
     endif
  end do
  end do
  end do
#endif

end subroutine qslope
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine fslope(f,df,q,dx,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer::ngrid
  real(dp)::dx,dt
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:npri)::q
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nmat)::f
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nmat,1:ndim)::df

  ! local arrays
  integer::i, j, k, l, n
  real(dp)::dsgn, dcen, dlft, drgt, slop
  real(dp)::dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr
  real(dp)::dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl
  real(dp)::dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm
  real(dp)::dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr
  real(dp)::vmin,vmax,dfx,dfy,dfz,dff
  integer::ilo,ihi,jlo,jhi,klo,khi
  real(dp),dimension(1:nvector),save::dlim

  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)

  if(slope_type==0)then
     df=zero
     return
  end if

#if NDIM==1
  do k = klo, khi
  do j = jlo, jhi
  do i = ilo, ihi
     if(slope_type==1.or.slope_type==2)then  ! minmod or average
        do n = 1, nmat
           do l = 1, ngrid
              dlim(l)=one
           end do
        end do
        do n = 1, nmat
           do l = 1, ngrid
              dlft = slope_type*(f(l,i,j,k,n)-f(l,i-1,j,k,n))
              drgt = slope_type*(f(l,i+1,j,k,n)-f(l,i,j,k,n))
              dcen = half*(dlft+drgt)/slope_type
              dsgn = sign(one, dcen)
              slop = min(abs(dlft),abs(drgt))
              if((dlft*drgt)<=zero)then
                 slop = zero
              else
                 slop = min(slop/abs(dcen),one)
              endif
              dlim(l) = min(slop,dlim(l))
           end do
        end do
        do n = 1, nmat
           do l = 1, ngrid
              df(l,i,j,k,n,1) = dlim(l)*half*(f(l,i+1,j,k,n)-f(l,i-1,j,k,n))
           end do
        end do
     else if(slope_type==3)then ! superbee
        do n = 1, nmat
           do l = 1, ngrid
              dcen = q(l,i,j,k,1)*dt/dx
              dlft = two/(one+dcen)*(f(l,i,j,k,n)-f(l,i-1,j,k,n))
              drgt = two/(one-dcen)*(f(l,i+1,j,k,n)-f(l,i,j,k,n))
              dsgn = sign(one, dlft)
              slop = min(abs(dlft),abs(drgt))
              dlim(l) = slop
              if((dlft*drgt)<=zero)dlim(l)=zero
              df(l,i,j,k,n,1) = dsgn*dlim(l)
           end do
        end do
     else if(slope_type==4)then ! ultrabee
        do n = 1, nmat
           do l = 1, ngrid
              dcen = q(l,i,j,k,1)*dt/dx
              if(dcen>=0)then
                 dlft = two/(zero+dcen+1d-10)*(f(l,i,j,k,n)-f(l,i-1,j,k,n))
                 drgt = two/(one -dcen      )*(f(l,i+1,j,k,n)-f(l,i,j,k,n))
              else
                 dlft = two/(one +dcen      )*(f(l,i,j,k,n)-f(l,i-1,j,k,n))
                 drgt = two/(zero-dcen+1d-10)*(f(l,i+1,j,k,n)-f(l,i,j,k,n))
              endif
              dsgn = sign(one, dlft)
              slop = min(abs(dlft),abs(drgt))
              dlim(l) = slop
              if((dlft*drgt)<=zero)dlim(l)=zero
              df(l,i,j,k,n,1) = dsgn*dlim(l)
           end do
        end do
     else
        write(*,*)'Unknown slope type'
        stop
     end if
  end do
  end do
  end do
#endif

#if NDIM==2
  do k = klo, khi
  do j = jlo, jhi
  do i = ilo, ihi
     if(slope_type==1.or.slope_type==2)then  ! minmod or average
        ! slopes in first coordinate direction
        do l = 1, ngrid
           dlim(l)=one
        enddo
        do n = 1, nmat
           do l = 1, ngrid
              dlft = slope_type*(f(l,i,j,k,n)-f(l,i-1,j,k,n))
              drgt = slope_type*(f(l,i+1,j,k,n)-f(l,i,j,k,n))
              dcen = half*(dlft+drgt)/slope_type
              dsgn = sign(one, dcen)
              slop = min(abs(dlft),abs(drgt))
              if((dlft*drgt)<=zero)then
                 slop=zero
              else
                 slop=min(slop/abs(dcen),one)
              endif
              dlim(l)=min(slop,dlim(l))
           end do
        end do
        do n = 1, nmat
           do l = 1, ngrid
              df(l,i,j,k,n,1) = dlim(l)*half*(f(l,i+1,j,k,n)-f(l,i-1,j,k,n))
           end do
        end do
        ! slopes in second coordinate direction
        do l = 1, ngrid
           dlim(l)=one
        enddo
        do n = 1, nmat
           do l = 1, ngrid
              dlft = slope_type*(f(l,i,j,k,n)-f(l,i,j-1,k,n))
              drgt = slope_type*(f(l,i,j+1,k,n)-f(l,i,j,k,n))
              dcen = half*(dlft+drgt)/slope_type
              dsgn = sign(one, dcen)
              slop = min(abs(dlft),abs(drgt))
              if((dlft*drgt)<=zero)then
                 slop=zero
              else
                 slop=min(slop/abs(dcen),one)
              endif
              dlim(l)=min(slop,dlim(l))
           end do
        end do
        do n = 1, nmat
           do l = 1, ngrid
              df(l,i,j,k,n,2) = dlim(l)*half*(f(l,i,j+1,k,n)-f(l,i,j-1,k,n))
           end do
        end do
     else if(slope_type==3)then ! positivity preserving 2d unsplit slope
        do l = 1, ngrid
           dlim(l) = one
        enddo
        do n = 1, nmat
           do l = 1, ngrid
              dfll = f(l,i-1,j-1,k,n)-f(l,i,j,k,n)
              dflm = f(l,i-1,j  ,k,n)-f(l,i,j,k,n)
              dflr = f(l,i-1,j+1,k,n)-f(l,i,j,k,n)
              dfml = f(l,i  ,j-1,k,n)-f(l,i,j,k,n)
              dfmm = f(l,i  ,j  ,k,n)-f(l,i,j,k,n)
              dfmr = f(l,i  ,j+1,k,n)-f(l,i,j,k,n)
              dfrl = f(l,i+1,j-1,k,n)-f(l,i,j,k,n)
              dfrm = f(l,i+1,j  ,k,n)-f(l,i,j,k,n)
              dfrr = f(l,i+1,j+1,k,n)-f(l,i,j,k,n)
              vmin = min(dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr)
              vmax = max(dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr)
              dfx  = half*(f(l,i+1,j,k,n)-f(l,i-1,j,k,n))
              dfy  = half*(f(l,i,j+1,k,n)-f(l,i,j-1,k,n))
              dff  = half*(abs(dfx)+abs(dfy))
              if(dff>zero)then
                 slop = min(one,min(abs(vmin),abs(vmax))/dff)
              else
                 slop = one
              endif
              dlim(l) = min(slop,dlim(l))
           end do
        end do
        ! Computed 2nd order limited slopes
        do n = 1, nmat
           do l = 1, ngrid
              df(l,i,j,k,n,1) = dlim(l)*half*(f(l,i+1,j,k,n)-f(l,i-1,j,k,n))
              df(l,i,j,k,n,2) = dlim(l)*half*(f(l,i,j+1,k,n)-f(l,i,j-1,k,n))
           end do
        end do
     else
        write(*,*)'Unkown slope type'
        stop
     endif
  end do
  end do
  end do
#endif

#if NDIM==3
  do k = klo, khi
  do j = jlo, jhi
  do i = ilo, ihi
     if(slope_type==1.or.slope_type==2)then  ! minmod or average
        ! slopes in first coordinate direction
        do l = 1, ngrid
           dlim(l)=one
        enddo
        do n = 1, nmat
           do l = 1, ngrid
              dlft = slope_type*(f(l,i,j,k,n)-f(l,i-1,j,k,n))
              drgt = slope_type*(f(l,i+1,j,k,n)-f(l,i,j,k,n))
              dcen = half*(dlft+drgt)/slope_type
              dsgn = sign(one, dcen)
              slop = min(abs(dlft),abs(drgt))
              if((dlft*drgt)<=zero)then
                 slop=zero
              else
                 slop=min(slop/abs(dcen),one)
              endif
              dlim(l)=min(slop,dlim(l))
           end do
        end do
        do n = 1, nmat
           do l = 1, ngrid
              df(l,i,j,k,n,1) = dlim(l)*half*(f(l,i+1,j,k,n)-f(l,i-1,j,k,n))
           end do
        end do
        ! slopes in second coordinate direction
        do l = 1, ngrid
           dlim(l)=one
        enddo
        do n = 1, nmat
           do l = 1, ngrid
              dlft = slope_type*(f(l,i,j,k,n)-f(l,i,j-1,k,n))
              drgt = slope_type*(f(l,i,j+1,k,n)-f(l,i,j,k,n))
              dcen = half*(dlft+drgt)/slope_type
              dsgn = sign(one, dcen)
              slop = min(abs(dlft),abs(drgt))
              if((dlft*drgt)<=zero)then
                 slop=zero
              else
                 slop=min(slop/abs(dcen),one)
              endif
              dlim(l)=min(slop,dlim(l))
           end do
        end do
        do n = 1, nmat
           do l = 1, ngrid
              df(l,i,j,k,n,2) = dlim(l)*half*(f(l,i,j+1,k,n)-f(l,i,j-1,k,n))
           end do
        end do
        ! slopes in third coordinate direction
        do l = 1, ngrid
           dlim(l)=one
        enddo
        do n = 1, nmat
           do l = 1, ngrid
              dlft = slope_type*(f(l,i,j,k,n)-f(l,i,j,k-1,n))
              drgt = slope_type*(f(l,i,j,k+1,n)-f(l,i,j,k,n))
              dcen = half*(dlft+drgt)/slope_type
              dsgn = sign(one, dcen)
              slop = min(abs(dlft),abs(drgt))
              if((dlft*drgt)<=zero)then
                 slop=zero
              else
                 slop=min(slop/abs(dcen),one)
              endif
              dlim(l)=min(slop,dlim(l))
           end do
        end do
        do n = 1, nmat
           do l = 1, ngrid
              df(l,i,j,k,n,3) = dlim(l)*half*(f(l,i,j,k+1,n)-f(l,i,j,k-1,n))
           end do
        end do
     else if(slope_type==3)then ! positivity preserving 3d unsplit slope
        do l = 1, ngrid
           dlim(l)=one
        enddo
        do n = 1, nmat
           do l = 1, ngrid
              dflll = f(l,i-1,j-1,k-1,n)-f(l,i,j,k,n)
              dflml = f(l,i-1,j  ,k-1,n)-f(l,i,j,k,n)
              dflrl = f(l,i-1,j+1,k-1,n)-f(l,i,j,k,n)
              dfmll = f(l,i  ,j-1,k-1,n)-f(l,i,j,k,n)
              dfmml = f(l,i  ,j  ,k-1,n)-f(l,i,j,k,n)
              dfmrl = f(l,i  ,j+1,k-1,n)-f(l,i,j,k,n)
              dfrll = f(l,i+1,j-1,k-1,n)-f(l,i,j,k,n)
              dfrml = f(l,i+1,j  ,k-1,n)-f(l,i,j,k,n)
              dfrrl = f(l,i+1,j+1,k-1,n)-f(l,i,j,k,n)
              dfllm = f(l,i-1,j-1,k  ,n)-f(l,i,j,k,n)
              dflmm = f(l,i-1,j  ,k  ,n)-f(l,i,j,k,n)
              dflrm = f(l,i-1,j+1,k  ,n)-f(l,i,j,k,n)
              dfmlm = f(l,i  ,j-1,k  ,n)-f(l,i,j,k,n)
              dfmmm = f(l,i  ,j  ,k  ,n)-f(l,i,j,k,n)
              dfmrm = f(l,i  ,j+1,k  ,n)-f(l,i,j,k,n)
              dfrlm = f(l,i+1,j-1,k  ,n)-f(l,i,j,k,n)
              dfrmm = f(l,i+1,j  ,k  ,n)-f(l,i,j,k,n)
              dfrrm = f(l,i+1,j+1,k  ,n)-f(l,i,j,k,n)
              dfllr = f(l,i-1,j-1,k+1,n)-f(l,i,j,k,n)
              dflmr = f(l,i-1,j  ,k+1,n)-f(l,i,j,k,n)
              dflrr = f(l,i-1,j+1,k+1,n)-f(l,i,j,k,n)
              dfmlr = f(l,i  ,j-1,k+1,n)-f(l,i,j,k,n)
              dfmmr = f(l,i  ,j  ,k+1,n)-f(l,i,j,k,n)
              dfmrr = f(l,i  ,j+1,k+1,n)-f(l,i,j,k,n)
              dfrlr = f(l,i+1,j-1,k+1,n)-f(l,i,j,k,n)
              dfrmr = f(l,i+1,j  ,k+1,n)-f(l,i,j,k,n)
              dfrrr = f(l,i+1,j+1,k+1,n)-f(l,i,j,k,n)
              vmin = min(dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl, &
                   &     dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm, &
                   &     dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr)
              vmax = max(dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl, &
                   &     dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm, &
                   &     dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr)
              dfx  = half*(f(l,i+1,j,k,n)-f(l,i-1,j,k,n))
              dfy  = half*(f(l,i,j+1,k,n)-f(l,i,j-1,k,n))
              dfz  = half*(f(l,i,j,k+1,n)-f(l,i,j,k-1,n))
              dff  = half*(abs(dfx)+abs(dfy)+abs(dfz))
              if(dff>zero)then
                 slop = min(one,min(abs(vmin),abs(vmax))/dff)
              else
                 slop = one
              endif
              dlim(l) = min(slop,dlim(l))
           end do
        end do
        do n = 1, nmat
           do l = 1, ngrid
              df(l,i,j,k,n,1) = dlim(l)*half*(f(l,i+1,j,k,n)-f(l,i-1,j,k,n))
              df(l,i,j,k,n,2) = dlim(l)*half*(f(l,i,j+1,k,n)-f(l,i,j-1,k,n))
              df(l,i,j,k,n,3) = dlim(l)*half*(f(l,i,j,k+1,n)-f(l,i,j,k-1,n))
           end do
        end do
     else
        write(*,*)'Unkown slope type'
        stop
     endif
  end do
  end do
  end do
#endif

end subroutine fslope
!###########################################################
!###########################################################
!###########################################################
!###########################################################
