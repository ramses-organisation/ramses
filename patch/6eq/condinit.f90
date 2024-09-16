!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(x,u,dx,nn)
  use amr_parameters
  use hydro_parameters
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1:nmat): f_k, U(i,nmat:2*nmat): f_k.d_k,
  ! U(i,2*nmat:2*nmat+ndim): d.u,d.v,d.w, U(i,2*nmat+ndim:3*+ndim):f_k.E_k
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1:ndim): u,v,w, Q(i,ndim:nmat+ndim): P_k,
  ! Q(i,nmat+ndim:2*nmat+ndim): e_k,
  ! If nvar >= ndim+3, remaining variables are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
  integer::ivar,imat,idim,i,k,id,iu,iv,iw,ip1,ip2,io
  real(dp),dimension(1:nvector,1:npri),save::q             ! Primitive variables
  real(dp),dimension(1:nvector,1:nmat),save::f,g           ! Volume fraction and densities
  real(dp),dimension(1:nvector),save::ekin,dtot,eint,p,cs,ftot
  real(dp),dimension(1:nvector),save::g_mat,eint_mat,p_mat,cs_mat
  logical,save::read_flag=.false.
  logical::inv
  integer,parameter::nrows=10000,ncols=3          ! CSV file parameters
  real(dp),dimension(1:nrows, 1:ncols),save::xx   ! Lane-Emden solutions (r, d, P)
  integer,save::nmax
  real(dp),save::rmax,dmax,pmax
  real(dp)::v0,rr,etot_mat
#if NVAR > NDIM + 3*NMAT
  integer::ipscal,npscal
  real(dp),dimension(1:nvector),save::s_mat
#endif

#if NDIM==1
  ! Call built-in initial condition generator
  call region_condinit(x,q,f,g,dx,nn)
#endif

#if NDIM==2
  iu=1; iv=2

  if (.not. read_flag) then
     ! Read Lane-Emden solutions into an array
     xx=0d0
     open (unit=10,file="patch/planet/le_0.5.csv",action="read",status="old")
     nmax=0
     do
        nmax=nmax+1
        read (10,*,iostat=io) xx(nmax,:)
        if(io.ne.0)exit
     end do
     read_flag = .true.
     ! Maximum radius for the planet
     rmax=xx(nmax-1,1)
     ! Minimum density for the planet
     dmax=xx(nmax-1,2)
     write(*,*)'Lane Emden file read'
     write(*,*)'nmax=',nmax,' rmax=',rmax,' dmin=',dmax
  end if

  ! Loop over all the cells
  do k=1,nn
     q(k,iu) = 0.0d0
     q(k,iv) = 0.0d0
     g(k,1)  = 1d-4
     g(k,2)  = 1d-4
     do imat=1,nmat
        q(k,ndim+imat) = g(k,1)**3 + 0.005*g(k,1)**3
     end do
     f(k,1)  = 1e-8
     f(k,2)  = 1.0

     rr=((x(k,1)-boxlen/2.0)**2+(x(k,2)-boxlen/2.0)**2)**(1.0/2)
     if(rr<rmax)then
        i=(rr/rmax)*(nmax-1)
        do imat=1,nmat
           g(k,imat)      = xx(i,2) + (rr-xx(i,1))*((xx(i+1,2)-xx(i,2))/(xx(i+1,1)-xx(i,1)))
           q(k,ndim+imat) = g(k,imat)**3 + 0.005*g(k,1)**3
        end do
        q(k,iu) = 0.0
        q(k,iv) = 0.0
        f(k,1)  = 1.0
        f(k,2)  = 1e-8
     endif

     ! rr=((x(k,1)-2.3*boxlen/4.0)**2+(x(k,2)-boxlen/2.0)**2)**(1.0/2)
     ! if(rr<rmax)then
     !    i=(rr/rmax)*(nmax-1)
     !    do imat=1,nmat
     !       g(k,imat)      = xx(i,2) + (rr-xx(i,1))*((xx(i+1,2)-xx(i,2))/(xx(i+1,1)-xx(i,1)))
     !       q(k,ndim+imat) = g(k,imat)**3 + 0.005*g(k,1)**3
     !    end do
     !    q(k,iu) = 0.0
     !    q(k,iv) = 0.0
     !    f(k,1)  = 1.0
     !    f(k,2)  = 1e-8
     ! endif
  end do
#endif

  ! normalize volume fraction
  ftot(1:nn)=0.0
  do imat=1,nmat
     do k=1,nn
        ftot(k)=ftot(k)+f(k,imat)
     end do
  end do
  do imat=1,nmat
     do k=1,nn
        f(k,imat)=f(k,imat)/ftot(k)
     end do
  end do

  ! compute kinetic energy
  ekin(1:nn)=0.0
  do k=1,nn
    do idim=1,ndim
      ekin(k) = ekin(k) + 0.5*q(k,idim)**2
    end do
  end do

  ! pressures --> partial total energies
  inv=.true.
  dtot(1:nn)=0.0
  do imat=1,nmat
    do k=1,nn
      g_mat(k) = g(k,imat)
      p_mat(k) = q(k,ndim+imat)
    end do
    ! call inverse eos routine (g,p) -> (e,c)
    call eos(g_mat,eint_mat,p_mat,cs_mat,imat,inv,nn)
    do k=1,nn
      etot_mat = eint_mat(k) + g(k,imat)* ekin(k) ! E_k
      u(k,2*nmat+ndim+imat) = etot_mat * f(k,imat) ! E_k * f_k
      dtot(k) = dtot(k) + f(k,imat) * g(k,imat)    ! total density
    end do
#if NVAR > NDIM + 3*NMAT
    call eos_s(g_mat,eint_mat,s_mat,imat,.false.,nn)
    do k=1,nn
       q(k,2*nmat+ndim+imat) = s_mat(k)
    end do
#endif
  end do

  ! volume fractions --> volume fractions
  u(1:nn,1:nmat)      = f(1:nn,1:nmat)
  ! physical densities -> partial masses (m_k = rho_k * f_k)
  do imat=1,nmat
    u(1:nn,nmat+imat) = f(1:nn,imat)*g(1:nn,imat)
  end do
  ! velocity -> momentum
  u(1:nn,2*nmat+1)    = dtot(1:nn)*q(1:nn,1)
#if NDIM>1
  u(1:nn,2*nmat+2)    = dtot(1:nn)*q(1:nn,2)
#endif
#if NDIM>2
  u(1:nn,2*nmat+ndim) = dtot(1:nn)*q(1:nn,ndim)
#endif
  ! passive scalars for each fluid
#if NVAR > NDIM + 3*NMAT
  npscal = (nvar - ndim - 3*nmat) / nmat
  do imat = 1, nmat
     do ipscal = 1, npscal
        ivar = ndim + 3*nmat + npscal*(imat-1) + ipscal
        u(1:nn,ivar) = f(1:nn,imat)*g(1:nn,imat)*q(1:nn,ivar-nmat)
     end do
  end do
#endif

end subroutine condinit
