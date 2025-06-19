! ---------------------------------------------------------------
!  UNSPLIT     Unsplit second order Godunov integrator for
!              polytropic gas dynamics using either
!              MUSCL-HANCOCK scheme or Collela's PLMDE scheme
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
subroutine unsplit(uin,gravin,pin,flux,tmp,dx,dy,dz,dt,ngrid)
  use amr_parameters
  use const
  use hydro_parameters
  implicit none

  integer ::ngrid
  real(dp)::dx,dy,dz,dt
  real(dp)::dtdx

  ! Input states
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::pin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::uin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim)::gravin

  ! Output fluxes
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim)::flux
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:2   ,1:ndim)::tmp

  ! Primitive variables
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),save::qin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2       ),save::cin

  ! Slopes
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim),save::dq

  ! Left and right state arrays
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim),save::qm
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim),save::qp

  ! Velocity divergence
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2)::divu

  ! Local scalar variables
  integer::i,j,k,l,ivar
  integer::ilo,ihi,jlo,jhi,klo,khi

  ilo=MIN(1,iu1+2); ihi=MAX(1,iu2-2)
  jlo=MIN(1,ju1+2); jhi=MAX(1,ju2-2)
  klo=MIN(1,ku1+2); khi=MAX(1,ku2-2)
  dtdx = dt/dx

  ! Translate to primitive variables, compute sound speeds
  call ctoprim(uin,qin,cin,gravin,dt,ngrid)

  ! Compute TVD slopes
  call uslope(qin,dq,dx,dt,ngrid)

  ! Compute 3D traced-states in all three directions
  if(scheme=='muscl')then
#if NDIM==1
     call trace(qin,dq,qm,qp,dx      ,dt,ngrid)
#endif
#if NDIM==2
    call trace2d(qin,dq,qm,qp,dx,dy   ,dt,ngrid)
#endif
!#if NDIM==3
!     call trace3d(qin,dq,qm,qp,dx,dy,dz,dt,ngrid)
!#endif
  endif
  if(scheme=='plmde')then
#if NDIM==1
     call tracex  (qin,dq,cin,qm,qp,dx      ,dt,ngrid)
#endif
#if NDIM==2
     call tracexy (qin,dq,cin,qm,qp,dx,dy   ,dt,ngrid)
#endif
#if NDIM==3
     call tracexyz(qin,dq,cin,qm,qp,dx,dy,dz,dt,ngrid)
#endif
  endif

  ! Solve for 1D flux in X direction
  call cmpflxm(qm,iu1+1,iu2+1,ju1  ,ju2  ,ku1  ,ku2  , &
       &       qp,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
       &          if1  ,if2  ,jlo  ,jhi  ,klo  ,khi  , &
       &       2,3,4,flux,tmp,1,dtdx,ngrid)

  ! Solve for 1D flux in Y direction
#if NDIM>1
  call cmpflxm(qm,iu1  ,iu2  ,ju1+1,ju2+1,ku1  ,ku2  , &
       &       qp,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
       &          ilo  ,ihi  ,jf1  ,jf2  ,klo  ,khi  , &
       &       3,2,4,flux,tmp,2,dtdx,ngrid)
#endif

  ! Solve for 1D flux in Z direction
#if NDIM>2
  call cmpflxm(qm,iu1  ,iu2  ,ju1  ,ju2  ,ku1+1,ku2+1, &
       &       qp,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
       &          ilo  ,ihi  ,jlo  ,jhi  ,kf1  ,kf2  , &
       &       4,2,3,flux,tmp,3,dtdx,ngrid)
#endif

  if(difmag>0.0)then
    call cmpdivu(qin,divu,dx,dy,dz,ngrid)
    call consup(uin,flux,divu,dt,ngrid)
  endif

end subroutine unsplit
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine trace(q,dq,qm,qp,dx,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer ::ngrid
  real(dp)::dx, dt

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qm
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qp

  ! Local variables
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::oneoverr
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::S0
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::dvel_sum
  integer ::i, j, k, l, ivar, idim, ivel
  integer ::ilo,ihi,jlo,jhi,klo,khi
  integer,parameter ::ir=1,ip=neul
  real(dp)::dtdx
  real(dp)::r, u, p,vel
  real(dp)::dvar, source
  real(dp)::sr0, sp0
#if NENER>0
  integer::irad
  real(dp)::e
#endif

  dtdx = dt/dx
  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)

  ! Initialize qd and qm to q for all variables
  ! and apply TVD slopes
  do idim=1,ndim
  do ivar=1,nvar
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              do l = 1, ngrid
                 dvar = half*dq(l,i,j,k,ivar,idim)
                 qp(l,i,j,k,ivar,idim) = q(l,i,j,k,ivar) - dvar
                 qm(l,i,j,k,ivar,idim) = q(l,i,j,k,ivar) + dvar
               end do
           end do
        end do
     end do
  end do
  end do

  ! Calculate first term of source (sum over directions)
  S0 = 0d0
  do idim=1,ndim
  do ivar=1,nvar
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              do l = 1, ngrid
                 S0(l,i,j,k,ivar) = S0(l,i,j,k,ivar) - q(l,i,j,k,1+idim)*dq(l,i,j,k,ivar,idim)
               end do
           end do
        end do
     end do
  end do
  end do

  ! precalc 1/rho
  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid
              r   =  q(l,i,j,k,ir)
              oneoverr(l,i,j,k) = 1d0/r
            end do
         end do
      end do
   end do

  ! add second term of source for velocities
  do ivel=1,ndim
  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid
              source = - (dq(l,i,j,k,ip,ivel))*oneoverr(l,i,j,k)
#if NENER>0
              do irad=1,nener
                 source = source - dq(l,i,j,k,ip+irad,ivel)*oneoverr(l,i,j,k)
              end do
#endif

               S0(l,i,j,k,1+ivel) = S0(l,i,j,k,1+ivel) + source
            end do
         end do
      end do
   end do
  end do

  ! Apply first term of source for all variables
  do idim=1,ndim
  do ivar=1,nvar
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              do l = 1, ngrid
                 source = S0(l,i,j,k,ivar) * dtdx*half
                 qp(l,i,j,k,ivar,idim) = qp(l,i,j,k,ivar,idim) + source
                 qm(l,i,j,k,ivar,idim) = qm(l,i,j,k,ivar,idim) + source
               end do
           end do
        end do
     end do
  end do
  end do



  ! Transverse derivatives for velocities
  do ivel=1,ndim
  do idim=1,ndim
  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid
              source = -dq(l,i,j,k,ip,idim)
#if NENER>0
              ! correct source vel for nener
              do irad=1,nener
                 dvar = dq(l,i,j,k,ip+irad,idim)
                 source = source - dvar
              end do
#endif
              ! apply source
              source = source*oneoverr(l,i,j,k) * dtdx*half
              qp(l,i,j,k,1+ivel,idim) = qp(l,i,j,k,1+ivel,idim) + source
              qm(l,i,j,k,1+ivel,idim) = qm(l,i,j,k,1+ivel,idim) + source

            end do
         end do
      end do
   end do
   end do
   end do

   ! calc transverse term for rho and pressure
   dvel_sum = 0
  do idim=1,ndim
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              do l = 1, ngrid
                 dvel_sum(l,i,j,k) = dvel_sum(l,i,j,k) + dq(l,i,j,k,1+idim,idim)
               end do
           end do
        end do
     end do
  end do

  ! Transverse derivatives for density and pressure
  do idim=1,ndim
  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid

              ! Cell centered values for density and pressure
              r   =  q(l,i,j,k,ir)
              p   =  q(l,i,j,k,ip)

              ! Source transverse derivatives
              ! density:  - (dux+dvy+dwz)/r
              ! pressure: - (dux+dvy+dwz)*gamma*p
              sr0 = - ( dvel_sum(l,i,j,k))*r       * dtdx*half
              sp0 = - ( dvel_sum(l,i,j,k))*gamma*p * dtdx*half

              ! Right state
              qp(l,i,j,k,ir,idim) = qp(l,i,j,k,ir,idim) + sr0
              qp(l,i,j,k,ip,idim) = qp(l,i,j,k,ip,idim) + sp0

              ! Left state
              qm(l,i,j,k,ir,idim) = qm(l,i,j,k,ir,idim) + sr0
              qm(l,i,j,k,ip,idim) = qm(l,i,j,k,ip,idim) + sp0
           end do
        end do
     end do
  end do
  end do

  ! Transverse derivatives for nener
#if NENER>0
  do idim=1,ndim
  do irad=1,nener
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              do l = 1, ngrid
                 e = q(l,i,j,k,ip+irad)
                 dvar = dq(l,i,j,k,1+idim,idim)
                 source = - (dvar)*gamma_rad(irad)*e * dtdx*half
                 qp(l,i,j,k,ip+irad,idim) = qp(l,i,j,k,ip+irad,idim) + source
                 qm(l,i,j,k,ip+irad,idim) = qm(l,i,j,k,ip+irad,idim) + source
              end do
           end do
        end do
     end do
  end do
  end do
#endif

  ! Limit for density
  do idim=1,ndim
  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid
              r   =  q(l,i,j,k,ir)
              if(qp(l,i,j,k,ir,idim)<smallr)qp(l,i,j,k,ir,idim)=r
              if(qm(l,i,j,k,ir,idim)<smallr)qm(l,i,j,k,ir,idim)=r
           end do
        end do
     end do
  end do
  end do

end subroutine trace
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#if NDIM>1
subroutine trace2d(q,dq,qm,qp,dx,dy,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer ::ngrid
  real(dp)::dx, dy, dt

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qm
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::qp

  ! declare local variables
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::oneoverr
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::S0
  integer ::i, j, k, l, ivar, idim, ivel
  integer ::ilo,ihi,jlo,jhi,klo,khi
  integer ::ir, iu, iv, ip
  real(dp)::dtdx
  real(dp)::r, u, v, p, vel, dvar, source
  real(dp)::drx, dux, dvx, dpx
  real(dp)::dry, duy, dvy, dpy
  real(dp)::sr0, su0, sv0, sp0
  real(dp)::sr0_bis, su0_bis, sv0_bis, sp0_bis
#if NENER>0
  integer ::irad
  real(dp),dimension(1:nener)::e, dex, dey, se0
#endif
#if NVAR > NHYDRO + NENER
  integer ::n
  real(dp)::a, dax, day, sa0
#endif

  dtdx = dt/dx
  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)
  ir=1; iu=2; iv=3; ip=4

  ! Initialize qd and qm to q for all variables
  ! and apply TVD slopes
  do idim=1,ndim
  do ivar=1,nvar
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              do l = 1, ngrid
                 dvar = half*dq(l,i,j,k,ivar,idim)
                 qp(l,i,j,k,ivar,idim) = q(l,i,j,k,ivar) - dvar! + source
                 qm(l,i,j,k,ivar,idim) = q(l,i,j,k,ivar) + dvar! + source
               end do
           end do
        end do
     end do
  end do
  end do

  ! Calculate first term of source (sum over directions)
  S0 = 0d0
  do idim=1,ndim
  do ivar=1,nvar
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              do l = 1, ngrid
                 S0(l,i,j,k,ivar) = S0(l,i,j,k,ivar) - q(l,i,j,k,1+idim)*dq(l,i,j,k,ivar,idim)
               end do
           end do
        end do
     end do
  end do
  end do

  ! precalc 1/rho
  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid
              r   =  q(l,i,j,k,ir)
              oneoverr(l,i,j,k) = 1d0/r
            end do
         end do
      end do
   end do

  ! add second term of source for velocities
  do ivel=1,ndim
  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid
              source = - (dq(l,i,j,k,ip,ivel))*oneoverr(l,i,j,k)
#if NENER>0
              do irad=1,nener
                 source = source - dq(l,i,j,k,ip+irad,ivel)*oneoverr(l,i,j,k)
              end do
#endif

               S0(l,i,j,k,1+ivel) = S0(l,i,j,k,1+ivel) + source
            end do
         end do
      end do
   end do
  end do

  ! Apply first term of source for all variables
  do idim=1,ndim
  do ivar=1,nvar
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              do l = 1, ngrid
                 source = S0(l,i,j,k,ivar) * dtdx*half
                 qp(l,i,j,k,ivar,idim) = qp(l,i,j,k,ivar,idim) + source
                 qm(l,i,j,k,ivar,idim) = qm(l,i,j,k,ivar,idim) + source
               end do
           end do
        end do
     end do
  end do
  end do

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid

              ! TVD slopes in all directions
              dux = dq(l,i,j,k,iu,1)
              dvy = dq(l,i,j,k,iv,2)
              r   =  q(l,i,j,k,ir)
              p   =  q(l,i,j,k,ip)
              sr0 = - (dux+dvy)*r
              sp0 = - (dux+dvy)*gamma*p

              do idim=1,ndim

              ! Right state at left interface
              qp(l,i,j,k,ir,idim) = qp(l,i,j,k,ir,idim) + sr0*dtdx*half
              qp(l,i,j,k,ip,idim) = qp(l,i,j,k,ip,idim) + sp0*dtdx*half
              if(qp(l,i,j,k,ir,idim)<smallr)qp(l,i,j,k,ir,idim)=r

              ! Left state at right interface
              qm(l,i,j,k,ir,idim) = qm(l,i,j,k,ir,idim) + sr0*dtdx*half
              qm(l,i,j,k,ip,idim) = qm(l,i,j,k,ip,idim) + sp0*dtdx*half
              if(qm(l,i,j,k,ir,idim)<smallr)qm(l,i,j,k,ir,idim)=r
              end do

#if NENER>0
              do irad=1,nener
                 e(irad) = q(l,i,j,k,ip+irad)
                 se0(irad) = - (dux+dvy)*gamma_rad(irad)*e(irad)
                 qp(l,i,j,k,ip+irad,1) = e(irad) + se0(irad)*dtdx*half
                 qm(l,i,j,k,ip+irad,1) = e(irad) + se0(irad)*dtdx*half
                 qp(l,i,j,k,ip+irad,2) = e(irad) + se0(irad)*dtdx*half
                 qm(l,i,j,k,ip+irad,2) = e(irad)  + se0(irad)*dtdx*half
              end do
#endif

           end do
        end do
     end do
  end do

end subroutine trace2d
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#if NDIM>2
subroutine trace3d(q,dq,qm,qp,dx,dy,dz,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer,intent(in)::ngrid
  real(dp),intent(in)::dx, dy, dz, dt

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar),intent(in)::q
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim),intent(in)::dq
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim),intent(out)::qm
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim),intent(out)::qp

  ! declare local variables
  integer ::i, j, k, l, ivar, idim
  integer ::ilo,ihi,jlo,jhi,klo,khi
  integer,parameter::ir=1,iu=2,iv=3,iw=4,ip=5
  real(dp)::dtdx
  real(dp)::r, u, v, w, p, oneoverr
  real(dp)::dvar, dvarx, dvary, dvarz, dvel_sum
  real(dp)::base_source,source_r,source_u,source_v,source_w,source_p
#if NENER>0
  integer ::irad
  real(dp)::e,source_e
#endif

  dtdx = dt/dx
  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)

  ! Initialize qd and qm to q for all variables
  do idim=1,ndim
  do ivar=1,nvar
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              do l = 1, ngrid
                 qp(l,i,j,k,ivar,idim) = q(l,i,j,k,ivar)
                 qm(l,i,j,k,ivar,idim) = q(l,i,j,k,ivar)
               end do
           end do
        end do
     end do
  end do
  end do

  ! Apply TVD slopes for all variables
  do idim=1,ndim
  do ivar=1,nvar
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              do l = 1, ngrid
                 ! Half TVD slope
                 dvar = half*dq(l,i,j,k,ivar,idim)
                 qp(l,i,j,k,ivar,idim) = qp(l,i,j,k,ivar,idim) - dvar
                 qm(l,i,j,k,ivar,idim) = qm(l,i,j,k,ivar,idim) + dvar
               end do
            end do
         end do
      end do
   end do
   end do

  ! Apply universal source for all variables
  do ivar=1,nvar
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              do l = 1, ngrid
                 u   =  q(l,i,j,k,iu)
                 v   =  q(l,i,j,k,iv)
                 w   =  q(l,i,j,k,iw)
                 dvarx = dq(l,i,j,k,ivar,1)
                 dvary = dq(l,i,j,k,ivar,2)
                 dvarz = dq(l,i,j,k,ivar,3)
                 base_source = -u*dvarx-v*dvary-w*dvarz
                 base_source = base_source * dtdx * half
                 qp(l,i,j,k,ivar,1) = qp(l,i,j,k,ivar,1) + base_source
                 qm(l,i,j,k,ivar,1) = qm(l,i,j,k,ivar,1) + base_source
                 qp(l,i,j,k,ivar,2) = qp(l,i,j,k,ivar,2) + base_source
                 qm(l,i,j,k,ivar,2) = qm(l,i,j,k,ivar,2) + base_source
                 qp(l,i,j,k,ivar,3) = qp(l,i,j,k,ivar,3) + base_source
                 qm(l,i,j,k,ivar,3) = qm(l,i,j,k,ivar,3) + base_source
               end do
            end do
         end do
      end do
   end do

  ! Transverse derivatives for velocities
   do k = klo, khi
      do j = jlo, jhi
         do i = ilo, ihi
            do l = 1, ngrid
               r   =  q(l,i,j,k,ir)
               oneoverr = 1d0/r

               source_u = dq(l,i,j,k,ip,1)
               source_v = dq(l,i,j,k,ip,2)
               source_w = dq(l,i,j,k,ip,3)

#if NENER>0
              ! correct source u,v,w for nener
              do irad=1,nener
                 dvarx = dq(l,i,j,k,ip+irad,1)
                 dvary = dq(l,i,j,k,ip+irad,2)
                 dvarz = dq(l,i,j,k,ip+irad,3)
                 source_u = source_u + dvarx
                 source_v = source_v + dvary
                 source_w = source_w + dvarz
              end do
#endif
               source_u = source_u * oneoverr * dtdx * half
               source_v = source_v * oneoverr * dtdx * half
               source_w = source_w * oneoverr * dtdx * half

              ! apply source
              ! Right state at left interface
              qp(l,i,j,k,iu,1) = qp(l,i,j,k,iu,1) - source_u
              qp(l,i,j,k,iv,1) = qp(l,i,j,k,iv,1) - source_v
              qp(l,i,j,k,iw,1) = qp(l,i,j,k,iw,1) - source_w

              ! Left state at right interface
              qm(l,i,j,k,iu,1) = qm(l,i,j,k,iu,1) - source_u
              qm(l,i,j,k,iv,1) = qm(l,i,j,k,iv,1) - source_v
              qm(l,i,j,k,iw,1) = qm(l,i,j,k,iw,1) - source_w
              
              ! Top state at bottom interface
              qp(l,i,j,k,iu,2) = qp(l,i,j,k,iu,2) - source_u
              qp(l,i,j,k,iv,2) = qp(l,i,j,k,iv,2) - source_v
              qp(l,i,j,k,iw,2) = qp(l,i,j,k,iw,2) - source_w

              ! Bottom state at top interface
              qm(l,i,j,k,iu,2) = qm(l,i,j,k,iu,2) - source_u
              qm(l,i,j,k,iv,2) = qm(l,i,j,k,iv,2) - source_v
              qm(l,i,j,k,iw,2) = qm(l,i,j,k,iw,2) - source_w

              ! Back state at front interface
              qp(l,i,j,k,iu,3) = qp(l,i,j,k,iu,3) - source_u
              qp(l,i,j,k,iv,3) = qp(l,i,j,k,iv,3) - source_v
              qp(l,i,j,k,iw,3) = qp(l,i,j,k,iw,3) - source_w

              ! Front state at back interface
              qm(l,i,j,k,iu,3) = qm(l,i,j,k,iu,3) - source_u
              qm(l,i,j,k,iv,3) = qm(l,i,j,k,iv,3) - source_v
              qm(l,i,j,k,iw,3) = qm(l,i,j,k,iw,3) - source_w

            end do
         end do
      end do
   end do
   
  ! Transverse derivatives for density and pressures
   do k = klo, khi
      do j = jlo, jhi
         do i = ilo, ihi
            do l = 1, ngrid
               r   =  q(l,i,j,k,ir)
               p   =  q(l,i,j,k,ip)

               dvarx = dq(l,i,j,k,iu,1)
               dvary = dq(l,i,j,k,iv,2)
               dvarz = dq(l,i,j,k,iw,3)
               dvel_sum = (dvarx+dvary+dvarz)

               source_r = dvel_sum*r       * dtdx*half
               source_p = dvel_sum*gamma*p * dtdx*half

              ! apply source
              ! Right state at left interface
               qp(l,i,j,k,ir,1) = qp(l,i,j,k,ir,1) - source_r
               qp(l,i,j,k,ip,1) = qp(l,i,j,k,ip,1) - source_p
 
               ! Left state at right interface
               qm(l,i,j,k,ir,1) = qm(l,i,j,k,ir,1) - source_r
               qm(l,i,j,k,ip,1) = qm(l,i,j,k,ip,1) - source_p
               
               ! Top state at bottom interface
               qp(l,i,j,k,ir,2) = qp(l,i,j,k,ir,2) - source_r
               qp(l,i,j,k,ip,2) = qp(l,i,j,k,ip,2) - source_p
 
               ! Bottom state at top interface
               qm(l,i,j,k,ir,2) = qm(l,i,j,k,ir,2) - source_r
               qm(l,i,j,k,ip,2) = qm(l,i,j,k,ip,2) - source_p
 
               ! Back state at front interface
               qp(l,i,j,k,ir,3) = qp(l,i,j,k,ir,3) - source_r
               qp(l,i,j,k,ip,3) = qp(l,i,j,k,ip,3) - source_p
 
               ! Front state at back interface
               qm(l,i,j,k,ir,3) = qm(l,i,j,k,ir,3) - source_r
               qm(l,i,j,k,ip,3) = qm(l,i,j,k,ip,3) - source_p

               ! Limit rho
               if(qp(l,i,j,k,ir,1)<smallr)qp(l,i,j,k,ir,1)=r
               if(qm(l,i,j,k,ir,1)<smallr)qm(l,i,j,k,ir,1)=r
               if(qp(l,i,j,k,ir,2)<smallr)qp(l,i,j,k,ir,2)=r
               if(qm(l,i,j,k,ir,2)<smallr)qm(l,i,j,k,ir,2)=r
               if(qp(l,i,j,k,ir,3)<smallr)qp(l,i,j,k,ir,3)=r
               if(qm(l,i,j,k,ir,3)<smallr)qm(l,i,j,k,ir,3)=r

#if NENER>0
              ! nener source
              do irad=1,nener
                 e = q(l,i,j,k,ip+irad)
                 source_e = dvel_sum*gamma_rad(irad)*e * dtdx*half
                 qp(l,i,j,k,ip+irad,1) = qp(l,i,j,k,ip+irad,1) - source_e
                 qm(l,i,j,k,ip+irad,1) = qm(l,i,j,k,ip+irad,1) - source_e
                 qp(l,i,j,k,ip+irad,2) = qp(l,i,j,k,ip+irad,2) - source_e
                 qm(l,i,j,k,ip+irad,2) = qm(l,i,j,k,ip+irad,2) - source_e
                 qp(l,i,j,k,ip+irad,3) = qp(l,i,j,k,ip+irad,3) - source_e
                 qm(l,i,j,k,ip+irad,3) = qm(l,i,j,k,ip+irad,3) - source_e
              end do
#endif

            end do
         end do
      end do
   end do

end subroutine trace3d

!###########################################################

!subroutine update_qd_qm_dvar(qp,qm,ivar,dvarx,dvary,dvarz)
!   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim),intent(inout)::qp
!   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim),intent(inout)::qm
!   integer,intent(in)::ivar
!   real(dp),intent(in)::dvarx,dvary,dvarz
   ! Right state at left interface
!   qp(l,i,j,k,ivar,1) = qp(l,i,j,k,ivar,1) - dvarx
   ! Left state at right interface
!   qm(l,i,j,k,ivar,1) = qp(l,i,j,k,ivar,1) + dvarx
   ! Top state at bottom interface
!   qp(l,i,j,k,ivar,2) = qp(l,i,j,k,ivar,2) - dvary
   ! Bottom state at top interface
!   qm(l,i,j,k,ivar,2) = qp(l,i,j,k,ivar,2) + dvary
   ! Back state at front interface
!   qp(l,i,j,k,ivar,3) = qp(l,i,j,k,ivar,3) - dvarz
   ! Front state at back interface
!   qm(l,i,j,k,ivar,3) = qp(l,i,j,k,ivar,3) + dvarz
!end subroutine update_qd_qm_dvar

!###########################################################

!subroutine update_qd_qm_source(qp,qm,ivar,source)
!   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim),intent(inout)::qp
!   real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim),intent(inout)::qm
!   integer,intent(in)::ivar
!   real(dp),intent(in)::source
   ! Right state at left interface
!   qp(l,i,j,k,ivar,1) = qp(l,i,j,k,ivar,1) + source
   ! Left state at right interface
!   qm(l,i,j,k,ivar,1) = qp(l,i,j,k,ivar,1) + source

   ! Top state at bottom interface
!   qp(l,i,j,k,ivar,2) = qp(l,i,j,k,ivar,2) + source
   ! Bottom state at top interface
!   qm(l,i,j,k,ivar,2) = qp(l,i,j,k,ivar,2) + source

   ! Back state at front interface
!   qp(l,i,j,k,ivar,3) = qp(l,i,j,k,ivar,3) + source
   ! Front state at back interface
!   qm(l,i,j,k,ivar,3) = qp(l,i,j,k,ivar,3) + source
!end subroutine update_qd_qm_dvar

#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmpflxm(qm,im1,im2,jm1,jm2,km1,km2, &
     &             qp,ip1,ip2,jp1,jp2,kp1,kp2, &
     &                ilo,ihi,jlo,jhi,klo,khi, ln,lt1,lt2, &
     &             flux,tmp,idim,dtdx,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  real(dp)::dtdx
  integer ::idim,ngrid
  integer ::ln,lt1,lt2
  integer ::im1,im2,jm1,jm2,km1,km2
  integer ::ip1,ip2,jp1,jp2,kp1,kp2
  integer ::ilo,ihi,jlo,jhi,klo,khi
  real(dp),dimension(1:nvector,im1:im2,jm1:jm2,km1:km2,1:nvar,1:ndim)::qm
  real(dp),dimension(1:nvector,ip1:ip2,jp1:jp2,kp1:kp2,1:nvar,1:ndim)::qp
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim)::flux
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:2,   1:ndim)::tmp

  ! local variables
  integer ::i, j, k, l, xdim
  real(dp)::entho
  real(dp),dimension(1:nvector,1:nvar),save::qleft,qright
  real(dp),dimension(1:nvector,1:nvar+1),save::fgdnv
#if NVAR > NHYDRO
  integer ::n
#endif

  entho=one/(gamma-one)
  xdim=ln-1

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi

           ! Mass density
           do l = 1, ngrid
              qleft (l,1) = qm(l,i,j,k,1,xdim)
              qright(l,1) = qp(l,i,j,k,1,xdim)
           end do

           ! Normal velocity
           do l = 1, ngrid
              qleft (l,2) = qm(l,i,j,k,ln,xdim)
              qright(l,2) = qp(l,i,j,k,ln,xdim)
           end do

           ! Pressure
           do l = 1, ngrid
              qleft (l,3) = qm(l,i,j,k,neul,xdim)
              qright(l,3) = qp(l,i,j,k,neul,xdim)
           end do

           ! Tangential velocity 1
#if NDIM>1
           do l = 1, ngrid
              qleft (l,4) = qm(l,i,j,k,lt1,xdim)
              qright(l,4) = qp(l,i,j,k,lt1,xdim)
           end do
#endif
           ! Tangential velocity 2
#if NDIM>2
           do l = 1, ngrid
              qleft (l,5) = qm(l,i,j,k,lt2,xdim)
              qright(l,5) = qp(l,i,j,k,lt2,xdim)
           end do
#endif
#if NVAR > NHYDRO
           ! Other advected quantities
           do n = nhydro+1, nvar
              do l = 1, ngrid
                 qleft (l,n) = qm(l,i,j,k,n,xdim)
                 qright(l,n) = qp(l,i,j,k,n,xdim)
              end do
           end do
#endif
           ! Solve Riemann problem
           if(riemann.eq.'acoustic')then
              call riemann_acoustic(qleft,qright,fgdnv,ngrid)
           else if (riemann.eq.'exact')then
              call riemann_approx  (qleft,qright,fgdnv,ngrid)
           else if (riemann.eq.'llf')then
              call riemann_llf     (qleft,qright,fgdnv,ngrid)
           else if (riemann.eq.'hllc')then
              call riemann_hllc    (qleft,qright,fgdnv,ngrid)
           else if (riemann.eq.'hll')then
              call riemann_hll     (qleft,qright,fgdnv,ngrid)
           else
              write(*,*)'unknown Riemann solver'
              stop
           end if

           ! Compute fluxes

           ! Mass density
           do l = 1, ngrid
              flux(l,i,j,k,1,idim) = fgdnv(l,1) * dtdx
           end do

           ! Normal momentum
           do l = 1, ngrid
              flux(l,i,j,k,ln,idim) = fgdnv(l,2) * dtdx
           end do

           ! Transverse momentum 1
#if NDIM>1
           do l = 1, ngrid
              flux(l,i,j,k,lt1,idim) = fgdnv(l,4) * dtdx
           end do
#endif
           ! Transverse momentum 2
#if NDIM>2
           do l = 1, ngrid
              flux(l,i,j,k,lt2,idim) = fgdnv(l,5) * dtdx
           end do
#endif
           ! Total energy
           do l = 1, ngrid
              flux(l,i,j,k,neul,idim) = fgdnv(l,3) * dtdx
           end do

#if NVAR > NHYDRO
           ! Other advected quantities
           do n = nhydro+1, nvar
              do l = 1, ngrid
                 flux(l,i,j,k,n,idim) = fgdnv(l,n) * dtdx
              end do
           end do
#endif
           ! Normal velocity
           do l = 1, ngrid
              tmp(l,i,j,k,1,idim) = half*(qleft(l,2)+qright(l,2)) * dtdx
           end do
           ! Internal energy flux
           do l = 1,ngrid
              tmp(l,i,j,k,2,idim) = fgdnv(l,nvar+1) * dtdx
           end do

        end do
     end do
  end do

end subroutine cmpflxm
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine ctoprim(uin,q,c,gravin,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer ::ngrid
  real(dp)::dt
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::uin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim)::gravin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::c

  integer ::i, j, k, l
  real(dp)::eint, smalle, dtxhalf, oneoverrho
  real(dp)::eken, erad
#if NVAR > NHYDRO + NENER
  integer ::n
#endif
#if NENER>0
  integer ::irad
#endif

  smalle = smallc**2/gamma/(gamma-one)
  dtxhalf = dt*half

  ! Convert to primitive variable
  do k = ku1, ku2
     do j = ju1, ju2
        do i = iu1, iu2
           do l = 1, ngrid

              ! Compute density
              q(l,i,j,k,1) = max(uin(l,i,j,k,1),smallr)

              ! Compute velocities
              oneoverrho = one/q(l,i,j,k,1)
              q(l,i,j,k,2) = uin(l,i,j,k,2)*oneoverrho
#if NDIM>1
              q(l,i,j,k,3) = uin(l,i,j,k,3)*oneoverrho
#endif
#if NDIM>2
              q(l,i,j,k,4) = uin(l,i,j,k,4)*oneoverrho
#endif

              ! Compute specific kinetic energy
              eken = half*q(l,i,j,k,2)*q(l,i,j,k,2)
#if NDIM>1
              eken = eken + half*q(l,i,j,k,3)*q(l,i,j,k,3)
#endif
#if NDIM>2
              eken = eken + half*q(l,i,j,k,4)*q(l,i,j,k,4)
#endif
              ! Compute non-thermal pressure
              erad = zero
#if NENER>0
              do irad = 1,nener
                 q(l,i,j,k,nhydro+irad) = (gamma_rad(irad)-one)*uin(l,i,j,k,nhydro+irad)
                 erad = erad+uin(l,i,j,k,nhydro+irad)*oneoverrho
              enddo
#endif
              ! Compute thermal pressure
              eint = MAX(uin(l,i,j,k,neul)*oneoverrho-eken-erad,smalle)
              q(l,i,j,k,neul) = (gamma-one)*q(l,i,j,k,1)*eint

              ! Compute sound speed
              c(l,i,j,k)=gamma*q(l,i,j,k,neul)
#if NENER>0
              do irad=1,nener
                 c(l,i,j,k)=c(l,i,j,k)+gamma_rad(irad)*q(l,i,j,k,nhydro+irad)
              enddo
#endif
              c(l,i,j,k)=sqrt(c(l,i,j,k)*oneoverrho)

              ! Gravity predictor step
              q(l,i,j,k,2) = q(l,i,j,k,2) + gravin(l,i,j,k,1)*dtxhalf
#if NDIM>1
              q(l,i,j,k,3) = q(l,i,j,k,3) + gravin(l,i,j,k,2)*dtxhalf
#endif
#if NDIM>2
              q(l,i,j,k,4) = q(l,i,j,k,4) + gravin(l,i,j,k,3)*dtxhalf
#endif

           end do
        end do
     end do
  end do

#if NVAR > NHYDRO + NENER
  ! Passive scalar
  do n = ndim+nener+3, nvar
     do k = ku1, ku2
        do j = ju1, ju2
           do i = iu1, iu2
              do l = 1, ngrid
                 oneoverrho = one/q(l,i,j,k,1)
                 q(l,i,j,k,n) = uin(l,i,j,k,n)*oneoverrho
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
subroutine uslope(q,dq,dx,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer::ngrid
  real(dp)::dx,dt
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar,1:ndim)::dq

  ! local arrays
  integer::i, j, k, l, n
  real(dp)::dsgn, dlim, dcen, dlft, drgt, slop
#if NDIM==2
  real(dp)::dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr
#endif
#if NDIM==3
  real(dp)::dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl
  real(dp)::dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm
  real(dp)::dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr
  real(dp)::dfz
#endif
#if NDIM>1
  real(dp)::vmin,vmax,dfx,dfy,dff
#endif
  integer::ilo,ihi,jlo,jhi,klo,khi

  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)

  if(slope_type==0)then
     dq=zero
     return
  end if

#if NDIM==1
  do n = 1, nvar
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              if(slope_type==1.or.slope_type==2.or.slope_type==3)then  ! minmod or average
                 do l = 1, ngrid
                    dlft = MIN(slope_type,2)*(q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                    drgt = MIN(slope_type,2)*(q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                    dcen = half*(dlft+drgt)/MIN(slope_type,2)
                    dsgn = sign(one, dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
                 end do
              else if(slope_type==4)then ! superbee
                 do l = 1, ngrid
                    dcen = q(l,i,j,k,2)*dt/dx
                    dlft = two/(one+dcen)*(q(l,i,j,k,n)-q(l,i-1,j,k,n))
                    drgt = two/(one-dcen)*(q(l,i+1,j,k,n)-q(l,i,j,k,n))
                    dcen = half*(q(l,i+1,j,k,n)-q(l,i-1,j,k,n))
                    dsgn = sign(one, dlft)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,1) = dsgn*dlim !min(dlim,abs(dcen))
                 end do
              else if(slope_type==5)then ! ultrabee
                 if(n==1)then
                    do l = 1, ngrid
                       dcen = q(l,i,j,k,2)*dt/dx
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
                       dcen = half*(q(l,i+1,j,k,n)-q(l,i-1,j,k,n))
                       if((dlft*drgt)<=zero)dlim=zero
                       dq(l,i,j,k,n,1) = dsgn*dlim !min(dlim,abs(dcen))
                    end do
                 else
                    do l = 1, ngrid
                       dq(l,i,j,k,n,1) = 0
                    end do
                 end if
              else if(slope_type==6)then ! unstable
                 if(n==1)then
                    do l = 1, ngrid
                       dlft = (q(l,i,j,k,n)-q(l,i-1,j,k,n))
                       drgt = (q(l,i+1,j,k,n)-q(l,i,j,k,n))
                       slop = 0.5d0*(dlft+drgt)
                       dlim = slop
                       dq(l,i,j,k,n,1) = dlim
                    end do
                 else
                    do l = 1, ngrid
                       dq(l,i,j,k,n,1) = 0
                    end do
                 end if
              else if(slope_type==7)then ! van Leer
                 do l = 1, ngrid
                    dlft = (q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                    drgt = (q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                    if((dlft*drgt)<=zero) then
                       dq(l,i,j,k,n,1)=zero
                    else
                       dq(l,i,j,k,n,1)=(2*dlft*drgt/(dlft+drgt))
                    end if
                 end do
              else if(slope_type==8)then ! generalized moncen/minmod parameterisation (van Leer 1979)
                 do l = 1, ngrid
                    dlft = (q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                    drgt = (q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                    dcen = half*(dlft+drgt)
                    dsgn = sign(one, dcen)
                    slop = min(slope_theta*abs(dlft),slope_theta*abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
                 end do
              else
                 write(*,*)'Unknown slope type',dx,dt
                 stop
              end if
           end do
        end do
     end do
  end do
#endif

#if NDIM==2
  if(slope_type==1.or.slope_type==2)then  ! minmod or average
     do n = 1, nvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
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
           end do
        end do
     end do
  else if(slope_type==3)then ! positivity preserving 2d unsplit slope
     do n = 1, nvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
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
           end do
        end do
     end do
  else if(slope_type==7)then ! van Leer
     do n = 1, nvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 ! slopes in first coordinate direction
                 do l = 1, ngrid
                    dlft = (q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                    drgt = (q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                    if((dlft*drgt)<=zero) then
                       dq(l,i,j,k,n,1)=zero
                    else
                       dq(l,i,j,k,n,1)=(2*dlft*drgt/(dlft+drgt))
                       end if
                 end do
                 ! slopes in second coordinate direction
                 do l = 1, ngrid
                    dlft = (q(l,i,j  ,k,n) - q(l,i,j-1,k,n))
                    drgt = (q(l,i,j+1,k,n) - q(l,i,j  ,k,n))
                    if((dlft*drgt)<=zero) then
                       dq(l,i,j,k,n,2)=zero
                    else
                       dq(l,i,j,k,n,2)=(2*dlft*drgt/(dlft+drgt))
                    end if
                 end do
              end do
           end do
        end do
     end do
  else if(slope_type==8)then ! generalized moncen/minmod parameterisation (van Leer 1979)
     do n = 1, nvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 ! slopes in first coordinate direction
                 do l = 1, ngrid
                    dlft = (q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                    drgt = (q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                    dcen = half*(dlft+drgt)
                    dsgn = sign(one, dcen)
                    slop = min(slope_theta*abs(dlft),slope_theta*abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
                 end do
                 ! slopes in second coordinate direction
                 do l = 1, ngrid
                    dlft = (q(l,i,j  ,k,n) - q(l,i,j-1,k,n))
                    drgt = (q(l,i,j+1,k,n) - q(l,i,j  ,k,n))
                    dcen = half*(dlft+drgt)
                    dsgn = sign(one,dcen)
                    slop = min(slope_theta*abs(dlft),slope_theta*abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,2) = dsgn*min(dlim,abs(dcen))
                 end do
              end do
           end do
        end do
     end do
  else
     write(*,*)'Unknown slope type',dx,dt
     stop
  endif
#endif

#if NDIM==3
  if(slope_type==1)then  ! minmod
     do n = 1, nvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 ! slopes in first coordinate direction
                 do l = 1, ngrid
                    dlft = q(l,i  ,j,k,n) - q(l,i-1,j,k,n)
                    drgt = q(l,i+1,j,k,n) - q(l,i  ,j,k,n)
                    if((dlft*drgt)<=zero) then
                       dq(l,i,j,k,n,1) = zero
                    else if(dlft>0) then
                       dq(l,i,j,k,n,1) = min(dlft,drgt)
                    else
                       dq(l,i,j,k,n,1) = max(dlft,drgt)
                    end if
                 end do
                 ! slopes in second coordinate direction
                 do l = 1, ngrid
                    dlft = q(l,i,j  ,k,n) - q(l,i,j-1,k,n)
                    drgt = q(l,i,j+1,k,n) - q(l,i,j  ,k,n)
                    if((dlft*drgt)<=zero) then
                       dq(l,i,j,k,n,2) = zero
                    else if(dlft>0) then
                       dq(l,i,j,k,n,2) = min(dlft,drgt)
                    else
                       dq(l,i,j,k,n,2) = max(dlft,drgt)
                    end if
                 end do
                 ! slopes in third coordinate direction
                 do l = 1, ngrid
                    dlft = q(l,i,j,k  ,n) - q(l,i,j,k-1,n)
                    drgt = q(l,i,j,k+1,n) - q(l,i,j,k  ,n)
                    if((dlft*drgt)<=zero) then
                       dq(l,i,j,k,n,3) = zero
                    else if(dlft>0) then
                       dq(l,i,j,k,n,3) = min(dlft,drgt)
                    else
                       dq(l,i,j,k,n,3) = max(dlft,drgt)
                    end if
                 end do
              end do
           end do
        end do
     end do
  else if(slope_type==2)then ! moncen
     do n = 1, nvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
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
           end do
        end do
     end do
  else if(slope_type==3)then ! positivity preserving 3d unsplit slope
     do n = 1, nvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
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
           end do
        end do
     end do
  else if(slope_type==7)then ! van Leer
     do n = 1, nvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 ! slopes in first coordinate direction
                 do l = 1, ngrid
                    dlft = (q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                    drgt = (q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                    if((dlft*drgt)<=zero) then
                       dq(l,i,j,k,n,1)=zero
                    else
                       dq(l,i,j,k,n,1)=(2*dlft*drgt/(dlft+drgt))
                    end if
                 end do
                 ! slopes in second coordinate direction
                 do l = 1, ngrid
                    dlft = (q(l,i,j  ,k,n) - q(l,i,j-1,k,n))
                    drgt = (q(l,i,j+1,k,n) - q(l,i,j  ,k,n))
                    if((dlft*drgt)<=zero) then
                       dq(l,i,j,k,n,2)=zero
                    else
                       dq(l,i,j,k,n,2)=(2*dlft*drgt/(dlft+drgt))
                    end if
                 end do
                 ! slopes in third coordinate direction
                 do l = 1, ngrid
                    dlft = (q(l,i,j,k  ,n) - q(l,i,j,k-1,n))
                    drgt = (q(l,i,j,k+1,n) - q(l,i,j,k  ,n))
                    if((dlft*drgt)<=zero) then
                       dq(l,i,j,k,n,3)=zero
                    else
                       dq(l,i,j,k,n,3)=(2*dlft*drgt/(dlft+drgt))
                    end if
                 end do
              end do
           end do
        end do
     end do
  else if(slope_type==8)then ! generalized moncen/minmod parameterisation (van Leer 1979)
     do n = 1, nvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 ! slopes in first coordinate direction
                 do l = 1, ngrid
                    dlft = (q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                    drgt = (q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                    dcen = half*(dlft+drgt)
                    dsgn = sign(one, dcen)
                    slop = min(slope_theta*abs(dlft),slope_theta*abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
                 end do
                 ! slopes in second coordinate direction
                 do l = 1, ngrid
                    dlft = (q(l,i,j  ,k,n) - q(l,i,j-1,k,n))
                    drgt = (q(l,i,j+1,k,n) - q(l,i,j  ,k,n))
                    dcen = half*(dlft+drgt)
                    dsgn = sign(one,dcen)
                    slop = min(slope_theta*abs(dlft),slope_theta*abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,2) = dsgn*min(dlim,abs(dcen))
                 end do
                 ! slopes in third coordinate direction
                 do l = 1, ngrid
                    dlft = (q(l,i,j,k  ,n) - q(l,i,j,k-1,n))
                    drgt = (q(l,i,j,k+1,n) - q(l,i,j,k  ,n))
                    dcen = half*(dlft+drgt)
                    dsgn = sign(one,dcen)
                    slop = min(slope_theta*abs(dlft),slope_theta*abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,3) = dsgn*min(dlim,abs(dcen))
                 end do
              end do
           end do
        end do
     end do
  else
     write(*,*)'Unknown slope type',dx,dt
     stop
  endif
#endif

end subroutine uslope
