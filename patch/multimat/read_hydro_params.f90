subroutine read_hydro_params(nml_ok)
  use amr_commons
  use hydro_commons
  use mpi_mod
  implicit none
  logical::nml_ok
  !--------------------------------------------------
  ! Local variables
  !--------------------------------------------------
  integer::i,idim,imat
  integer ,dimension(1:MAXBOUND)::bound_type
  real(dp)::scale,dtot,ftot
  real(dp),dimension(1:nvector,1:npri)::q
  real(dp),dimension(1:nvector,1:nmat)::f,g,kappa_mat
  real(dp),dimension(1:nvector)::ek_bound,eint,cs,kappa_hat

  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------
  namelist/init_params/initfile,nregion,region_type &
       & ,x_center,y_center,z_center &
       & ,length_x,length_y,length_z,exp_region &
       & ,d1_region,d2_region,d3_region,d4_region &
       & ,f1_region,f2_region,f3_region,f4_region &
       & ,u_region,v_region,w_region,p_region
  namelist/hydro_params/courant_factor,smallr,smallc,smallf,slope_type,scheme,eos_name
  namelist/refine_params/x_refine,y_refine,z_refine,r_refine,a_refine,b_refine,exp_refine &
       & ,m_refine,mass_sph,err_grad_f,err_grad_d,err_grad_p,err_grad_u &
       & ,floor_f,floor_d,floor_u,floor_p &
       & ,interpol_var,interpol_type
  namelist/boundary_params/nboundary,bound_type &
       & ,ibound_min,ibound_max,jbound_min,jbound_max,kbound_min,kbound_max &
       & ,d1_bound,d2_bound,d3_bound,d4_bound &
       & ,f1_bound,f2_bound,f3_bound,f4_bound &
       & ,u_bound,v_bound,w_bound,p_bound
  namelist/material_params/eos_type,eos_file,eos_params

  ! Read namelist file
  rewind(1)
  read(1,NML=init_params)
  rewind(1)
  read(1,NML=material_params)
  rewind(1)
  if(nlevelmax>levelmin)read(1,NML=refine_params)
  rewind(1)
  if(hydro)read(1,NML=hydro_params)
  rewind(1)
  read(1,NML=boundary_params,END=103)
  simple_boundary=.true.
  goto 104
103 simple_boundary=.false.
104 if(nboundary>MAXBOUND)then
    write(*,*) 'Error: nboundary>MAXBOUND'
    call clean_stop
  end if
  rewind(1)
  !-------------------------------------------------
  ! This section deals with hydro boundary conditions
  !-------------------------------------------------
  if(simple_boundary.and.nboundary==0)then
     simple_boundary=.false.
  endif

  if (simple_boundary)then

     ! Compute new coarse grid boundaries
     do i=1,nboundary
        if(ibound_min(i)*ibound_max(i)==1.and.ndim>0)then
           nx=nx+1
           if(ibound_min(i)==-1)then
              icoarse_min=icoarse_min+1
              icoarse_max=icoarse_max+1
           end if
        end if
     end do
     do i=1,nboundary
        if(jbound_min(i)*jbound_max(i)==1.and.ndim>1)then
           ny=ny+1
           if(jbound_min(i)==-1)then
              jcoarse_min=jcoarse_min+1
              jcoarse_max=jcoarse_max+1
           end if
        end if
     end do
     do i=1,nboundary
        if(kbound_min(i)*kbound_max(i)==1.and.ndim>2)then
           nz=nz+1
           if(kbound_min(i)==-1)then
              kcoarse_min=kcoarse_min+1
              kcoarse_max=kcoarse_max+1
           end if
        end if
     end do

     ! Compute boundary geometry
     do i=1,nboundary
        if(ibound_min(i)*ibound_max(i)==1.and.ndim>0)then
           if(ibound_min(i)==-1)then
              ibound_min(i)=icoarse_min+ibound_min(i)
              ibound_max(i)=icoarse_min+ibound_max(i)
              if(bound_type(i)==1)boundary_type(i)=1
              if(bound_type(i)==2)boundary_type(i)=11
              if(bound_type(i)==3)boundary_type(i)=21
           else
              ibound_min(i)=icoarse_max+ibound_min(i)
              ibound_max(i)=icoarse_max+ibound_max(i)
              if(bound_type(i)==1)boundary_type(i)=2
              if(bound_type(i)==2)boundary_type(i)=12
              if(bound_type(i)==3)boundary_type(i)=22
           end if
           if(ndim>1)jbound_min(i)=jcoarse_min+jbound_min(i)
           if(ndim>1)jbound_max(i)=jcoarse_max+jbound_max(i)
           if(ndim>2)kbound_min(i)=kcoarse_min+kbound_min(i)
           if(ndim>2)kbound_max(i)=kcoarse_max+kbound_max(i)
        else if(jbound_min(i)*jbound_max(i)==1.and.ndim>1)then
           ibound_min(i)=icoarse_min+ibound_min(i)
           ibound_max(i)=icoarse_max+ibound_max(i)
           if(jbound_min(i)==-1)then
              jbound_min(i)=jcoarse_min+jbound_min(i)
              jbound_max(i)=jcoarse_min+jbound_max(i)
              if(bound_type(i)==1)boundary_type(i)=3
              if(bound_type(i)==2)boundary_type(i)=13
              if(bound_type(i)==3)boundary_type(i)=23
           else
              jbound_min(i)=jcoarse_max+jbound_min(i)
              jbound_max(i)=jcoarse_max+jbound_max(i)
              if(bound_type(i)==1)boundary_type(i)=4
              if(bound_type(i)==2)boundary_type(i)=14
              if(bound_type(i)==3)boundary_type(i)=24
           end if
           if(ndim>2)kbound_min(i)=kcoarse_min+kbound_min(i)
           if(ndim>2)kbound_max(i)=kcoarse_max+kbound_max(i)
        else if(kbound_min(i)*kbound_max(i)==1.and.ndim>2)then
           ibound_min(i)=icoarse_min+ibound_min(i)
           ibound_max(i)=icoarse_max+ibound_max(i)
           jbound_min(i)=jcoarse_min+jbound_min(i)
           jbound_max(i)=jcoarse_max+jbound_max(i)
           if(kbound_min(i)==-1)then
              kbound_min(i)=kcoarse_min+kbound_min(i)
              kbound_max(i)=kcoarse_min+kbound_max(i)
              if(bound_type(i)==1)boundary_type(i)=5
              if(bound_type(i)==2)boundary_type(i)=15
              if(bound_type(i)==3)boundary_type(i)=25
           else
              kbound_min(i)=kcoarse_max+kbound_min(i)
              kbound_max(i)=kcoarse_max+kbound_max(i)
              if(bound_type(i)==1)boundary_type(i)=6
              if(bound_type(i)==2)boundary_type(i)=16
              if(bound_type(i)==3)boundary_type(i)=26
           end if
        end if
     end do
     do i=1,nboundary
        ! Check for errors
        if( (ibound_min(i)<0.or.ibound_max(i)>(nx-1)) .and. (ndim>0) )then
           if(myid==1)write(*,*)'Error in the namelist'
           if(myid==1)write(*,*)'Check boundary conditions along X direction',i
           nml_ok=.false.
        end if
        if( (jbound_min(i)<0.or.jbound_max(i)>(ny-1)) .and. (ndim>1) )then
           if(myid==1)write(*,*)'Error in the namelist'
           if(myid==1)write(*,*)'Check boundary conditions along Y direction',i
           nml_ok=.false.
        end if
        if( (kbound_min(i)<0.or.kbound_max(i)>(nz-1)) .and. (ndim>2) )then
           if(myid==1)write(*,*)'Error in the namelist'
           if(myid==1)write(*,*)'Check boundary conditions along Z direction',i
           nml_ok=.false.
        end if
     end do
  end if

  !--------------------------------------------------
  ! Compute boundary conservative variables
  !--------------------------------------------------
  do i=1,nboundary
     imat=1
     f(i,imat)=max(f1_bound(i),smallf)
     g(i,imat)=max(d1_bound(i),smallr)
     if(nmat>1)then
        imat=imat+1
        f(i,imat)=max(f2_bound(i),smallf)
        g(i,imat)=max(d2_bound(i),smallr)
     end if
     if(nmat>2)then
        imat=imat+1
        f(i,imat)=max(f3_bound(i),smallf)
        g(i,imat)=max(d3_bound(i),smallr)
     end if
     if(nmat>3)then
        imat=imat+1
        f(i,imat)=max(f4_bound(i),smallf)
        g(i,imat)=max(d4_bound(i),smallr)
     end if
     ! Normalize volume fraction
     ftot=1d-15
     do imat=1,nmat
        ftot=ftot+f(i,imat)
     end do
     do imat=1,nmat
        f(i,imat)=f(i,imat)/ftot
     end do
     ! Compute total density
     q(i,1)=0.0
     do imat=1,nmat
        q(i,1)=q(i,1)+f(i,imat)*g(i,imat)
     end do
     ! Compute velocity
     q(i,2)=u_bound(i)
#if NDIM>1
     q(i,3)=v_bound(i)
#endif
#if NDIM>2
     q(i,4)=w_bound(i)
#endif
     q(i,npri)=p_bound(i)
  end do

  if(nboundary>0)call eosinv(f,g,q,eint,cs,kappa_mat,kappa_hat,nboundary)

  ! density -> density
  do i=1,nboundary
     boundary_var(i,1)=q(i,1)
  end do
  ! velocity -> momentum
  ek_bound(1:nboundary)=0.0
  do idim=1,ndim
     do i=1,nboundary
        boundary_var(i,idim+1)=q(i,1)*q(i,idim+1)
        ek_bound(i)=ek_bound(i)+0.5*q(i,idim+1)**2
     end do
  end do
  ! total energy
  do i=1,nboundary
     boundary_var(i,npri)=eint(i)+q(i,1)*ek_bound(i)
  end do
  ! volume fraction
  do imat=1,nmat
     do i=1,nboundary
        boundary_var(i,imat+npri)=f(i,imat)
     end do
  end do
  ! fluid density
  do imat=1,nmat
     do i=1,nboundary
        boundary_var(i,imat+npri+nmat)=g(i,imat)
     end do
  end do

end subroutine read_hydro_params
