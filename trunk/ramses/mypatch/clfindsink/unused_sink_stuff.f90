subroutine create_just_one_sink(pos,vel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use cooling_module, ONLY: XH=>X, rhoc, mH
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  real(dp),dimension(1:3)::pos,vel
  !----------------------------------------------------------------------                                                      
  ! Description: This subroutine creates one sink paricle at a given place                                                     
  ! It only updates the sink array and doesn't create true RAMSES particles.                                                   
  ! Andreas Bleuler, August 2011                                                                                               
  !----------------------------------------------------------------------                                                      
  ! local constants                                                                                                            

  integer ::ncache,nnew,ivar,ngrid,icpu,index_sink,index_sink_tot,icloud
  integer ::igrid,ix,iy,iz,ind,i,j,n,iskip,isink,inew,nx_loc
  integer ::ii,jj,kk,ind_cloud,ncloud
  integer ::info
  logical ::ok_free
#ifdef SOLVERhydro
  integer ::imetal=ndim+3
#endif
#ifdef SOLVERmhd
  integer ::imetal=9
#endif

  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m
  real(dp)::d,x,y,z,u,v,w,e,temp,zg,factG
  real(dp)::dxx,dyy,dzz,drr
  real(dp)::rsink_max2,rmax,rmax2
  real(dp)::velc,uc,vc,wc,l_jeans,d_jeans,d_thres,d_sink
  real(dp)::birth_epoch,xx,yy,zz,rr
  real(dp),dimension(1:3)::skip_loc
  real(dp)::dx,dx_loc,scale,vol_loc,dx_min,vol_min
  real(dp)::bx1,bx2,by1,by2,bz1,bz2


  if(.not. hydro)return
  if(ndim.ne.3)return



  ! Set new sink variables to zero                                                                                             
  msink_new=0d0; tsink_new=0d0; delta_mass_new=0d0; xsink_new=0d0; vsink_new=0d0; oksink_new=0d0; idsink_new=0




  if(myid .eq. 1)then

     print*,pos(1),pos(2),pos(3)


     if(verbose)write(*,*)' creating a sink from a clump '


     ! Conversion factor from user units to cgs units                                                                          
     call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
     scale_m=scale_d*scale_l**3

     ! Minimum radius to create a new sink from any other                                                                      
     rsink_max=10.0 ! in kpc                                                                                                   
     rsink_max=rsink_max*3.08d21/scale_l
     rsink_max2=rsink_max**2

     ! Maximum value for the initial sink mass                                                                                 
     msink_max=1d5 ! in Msol                                                                                                   
     msink_max=msink_max*2d33/scale_m

     ! Gravitational constant                                                                                                  
     factG=1
     if(cosmo)factG=3d0/8d0/3.1415926*omega_m*aexp

     ! Density threshold for sink particle creation                                                                            
     d_sink=n_sink/scale_nH

     ! Mesh spacing in that level                                                                                              
     dx=0.5D0**nlevelmax
     nx_loc=(icoarse_max-icoarse_min+1)
     skip_loc=(/0.0d0,0.0d0,0.0d0/)
     if(ndim>0)skip_loc(1)=dble(icoarse_min)
     if(ndim>1)skip_loc(2)=dble(jcoarse_min)
     if(ndim>2)skip_loc(3)=dble(kcoarse_min)
     scale=boxlen/dble(nx_loc)
     dx_loc=dx*scale
     vol_loc=dx_loc**ndim
     dx_min=(0.5D0**nlevelmax)*scale
     vol_min=dx_min**ndim

     rmax=dble(ir_cloud)*dx_min/aexp ! Linking length in physical units                                                        
     rmax2=rmax*rmax

     ! Birth epoch                                                                                                             
     birth_epoch=t


     !------------------------------                                                                                           
     ! Create new sink particles                                                                                               
     !------------------------------                                                                                           
     ! Starting identity number                                                                                                
     index_sink=nsink+1
     index_sink_tot=nindsink+1

     ! Mass of the new sink                   
     msink_new(index_sink)=vol_loc*d_sink/100.
     delta_mass_new(index_sink)=0d0

     ! Global index of the new sink 
     oksink_new(index_sink)=1d0
     idsink_new(index_sink)=index_sink_tot

     ! Store properties of the new sink
     tsink_new(index_sink)=birth_epoch
     xsink_new(index_sink,1)=pos(1)
     xsink_new(index_sink,2)=pos(2)
     xsink_new(index_sink,3)=pos(3)
     vsink_new(index_sink,1)=vel(1)
     vsink_new(index_sink,2)=vel(2)
     vsink_new(index_sink,3)=vel(3)
endif

nsink=nsink+1
nindsink=nindsink+1


#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(oksink_new,oksink_all,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(idsink_new,idsink_all,nsinkmax,MPI_INTEGER         ,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(msink_new ,msink_all ,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(tsink_new ,tsink_all ,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(xsink_new ,xsink_all ,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(vsink_new ,vsink_all ,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(delta_mass_new,delta_mass_all,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)

#else
  oksink_all=oksink_new
  idsink_all=idsink_new
  msink_all=msink_new
  tsink_all=tsink_new
  xsink_all=xsink_new
  vsink_all=vsink_new
  delta_mass_all=delta_mass_new
#endif
  do isink=1,nsink
     if(oksink_all(isink)==1)then
        idsink(isink)=idsink_all(isink)
        msink(isink)=msink_all(isink)
        tsink(isink)=tsink_all(isink)
        xsink(isink,1:ndim)=xsink_all(isink,1:ndim)
        vsink(isink,1:ndim)=vsink_all(isink,1:ndim)
        delta_mass(isink)=delta_mass_all(isink)
        if (myid==1)then
           print*,isink
           print*,idsink(isink)
           print*,msink(isink)
           print*,tsink(isink)
           print*,xsink(isink,1),xsink(isink,2),xsink(isink,3)
           print*,vsink(isink,1),vsink(isink,2),vsink(isink,3)
        endif
     endif
  end do



end subroutine create_just_one_sink
!################################################################
!################################################################
!################################################################
!################################################################ 
