subroutine init_sink
  use amr_commons
  use pm_commons
  use clfind_commons
  use mpi_mod
  use amr_parameters, only:levelmin
  implicit none
#ifndef WITHOUTMPI
  integer,parameter::tag=1112,tag2=1113
  integer::dummy_io,info2
#endif
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  integer::isink, nsinkold
  logical::eof,ic_sink=.false.
  character(LEN=80)::filename
  character(LEN=80)::fileloc
  character(LEN=5)::nchar,ncharcpu

  integer::sid,slevel
  real(dp)::sm1,sx1,sx2,sx3,sv1,sv2,sv3,sl1,sl2,sl3
  real(dp)::stform,sacc_rate,sacc_mass,srho_gas,sc2_gas,seps_sink,svg1,svg2,svg3,sm2
  character::co
  character(LEN=200)::comment_line

  ! Allocate all sink related quantities...
  allocate(idsink(1:nsinkmax))
  idsink=0 ! Important: need to set idsink to zero
  allocate(msink(1:nsinkmax))
  allocate(msmbh(1:nsinkmax))
  allocate(xsink(1:nsinkmax,1:ndim))
  xsink=boxlen/2.
  allocate(vsink(1:nsinkmax,1:ndim))
  allocate(lsink(1:nsinkmax,1:ndim))
  allocate(delta_mass(1:nsinkmax))

  allocate(tsink(1:nsinkmax))
  allocate(vsold(1:nsinkmax,1:ndim,levelmin:nlevelmax))
  allocate(vsnew(1:nsinkmax,1:ndim,levelmin:nlevelmax))
  allocate(fsink_partial(1:nsinkmax,1:ndim,levelmin:nlevelmax))
  allocate(fsink(1:nsinkmax,1:ndim))

  allocate(msum_overlap(1:nsinkmax))
  msum_overlap=0.0
  allocate(rho_sink_tff(levelmin:nlevelmax))

  ! Temporary sink variables
  allocate(wden(1:nsinkmax))
  allocate(wmom(1:nsinkmax,1:ndim))
  allocate(weth(1:nsinkmax))
  allocate(wvol(1:nsinkmax))
  allocate(wdiv(1:nsinkmax))
  allocate(wden_new(1:nsinkmax))
  allocate(wmom_new(1:nsinkmax,1:ndim))
  allocate(weth_new(1:nsinkmax))
  allocate(wvol_new(1:nsinkmax))
  allocate(wdiv_new(1:nsinkmax))
  allocate(msink_new(1:nsinkmax))
  allocate(msmbh_new(1:nsinkmax))
  allocate(msmbh_all(1:nsinkmax))
  allocate(msink_all(1:nsinkmax))
  allocate(tsink_new(1:nsinkmax))
  allocate(tsink_all(1:nsinkmax))
  allocate(idsink_new(1:nsinkmax))
  allocate(idsink_all(1:nsinkmax))
  allocate(idsink_old(1:nsinkmax))
  allocate(vsink_new(1:nsinkmax,1:ndim))
  allocate(vsink_all(1:nsinkmax,1:ndim))
  allocate(fsink_new(1:nsinkmax,1:ndim))
  allocate(fsink_all(1:nsinkmax,1:ndim))
  allocate(lsink_new(1:nsinkmax,1:ndim))
  allocate(lsink_all(1:nsinkmax,1:ndim))
  allocate(xsink_new(1:nsinkmax,1:ndim))
  allocate(xsink_all(1:nsinkmax,1:ndim))
  allocate(sink_jump(1:nsinkmax,1:ndim,levelmin:nlevelmax))
  sink_jump=0.d0
  allocate(dMsink_overdt(1:nsinkmax))
  allocate(dMBHoverdt(1:nsinkmax))
  allocate(dMsmbh_overdt(1:nsinkmax))
  allocate(dMBHoverdt_smbh(1:nsinkmax))
  allocate(eps_sink(1:nsinkmax))
  eps_sink=0.d0
  allocate(volume_gas(1:nsinkmax))
  volume_gas=0.d0
  allocate(vel_gas(1:nsinkmax,1:ndim))
  vel_gas=0.d0
  allocate(rho_gas(1:nsinkmax))
  rho_gas=0.d0
  allocate(c2sink(1:nsinkmax))
  allocate(weighted_density(1:nsinkmax,1:nlevelmax))
  weighted_density = 0.d0
  allocate(weighted_volume(1:nsinkmax,1:nlevelmax))
  weighted_volume = 0.d0
  allocate(weighted_ethermal(1:nsinkmax,1:nlevelmax))
  weighted_ethermal = 0.d0
  allocate(weighted_momentum(1:nsinkmax,1:nlevelmax,1:ndim))
  weighted_momentum = 0.d0
  allocate(weighted_divergence(1:nsinkmax,1:nlevelmax))
  weighted_divergence = 0.d0
  allocate(oksink_new(1:nsinkmax))
  allocate(oksink_all(1:nsinkmax))
  allocate(idsink_sort(1:nsinkmax))
  allocate(xmsink(1:nsinkmax))
  allocate(delta_mass_new(1:nsinkmax),delta_mass_all(1:nsinkmax))
  allocate(ok_blast_agn(1:nsinkmax),ok_blast_agn_all(1:nsinkmax))
  allocate(direct_force_sink(1:nsinkmax))
  allocate(new_born(1:nsinkmax),new_born_all(1:nsinkmax),new_born_new(1:nsinkmax))

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Loading sinks from the restart
  if(nrestart>0)then

     call title(nrestart,nchar)

     if(IOGROUPSIZEREP>0)then
        call title(((myid-1)/IOGROUPSIZEREP)+1,ncharcpu)
        fileloc='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/sink_'//TRIM(nchar)//'.csv'
     else
        fileloc='output_'//TRIM(nchar)//'/sink_'//TRIM(nchar)//'.csv'
     endif

     ! Wait for the token
#ifndef WITHOUTMPI
     if(IOGROUPSIZE>0) then
        if (mod(myid-1,IOGROUPSIZE)/=0) then
           call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
                & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
        end if
     endif
#endif

     nsink=0
     open(10,file=fileloc,form='formatted')
     eof=.false.
     ! scrolling over the comment lines
     read(10,'(A200)')comment_line
     read(10,'(A200)')comment_line
     do
        read(10,'(I10,20(A1,ES20.10),A1,I10)',end=104)sid,co, sm1,co,&
                           sx1,co,sx2,co,sx3,co, &
                           sv1,co,sv2,co,sv3,co, &
                           sl1,co,sl2,co,sl3,co, &
                           stform,co, sacc_rate,co, &
                           sacc_mass,co, &
                           srho_gas,co, sc2_gas,co, seps_sink,co, &
                           svg1,co,svg2,co,svg3,co, &
                           sm2,co,slevel
        nsink=nsink+1
        nindsink=nindsink+1
        idsink(nsink)=sid
        msink(nsink)=sm1
        xsink(nsink,1)=sx1
        xsink(nsink,2)=sx2
        xsink(nsink,3)=sx3
        vsink(nsink,1)=sv1
        vsink(nsink,2)=sv2
        vsink(nsink,3)=sv3
        lsink(nsink,1)=sl1
        lsink(nsink,2)=sl2
        lsink(nsink,3)=sl3
        tsink(nsink)=stform
        dMBHoverdt(nsink)=sacc_rate
        delta_mass(nsink)=sacc_mass
        rho_gas(nsink)=srho_gas
        c2sink(nsink)=sc2_gas
        eps_sink(nsink)=seps_sink
        vel_gas(nsink,1)=svg1
        vel_gas(nsink,2)=svg2
        vel_gas(nsink,3)=svg3
        new_born(nsink)=.false. ! this is a restart
        msmbh(nsink)=sm2
     end do
     sinkint_level=slevel
104  continue
     close(10)

     ! Send the token
#ifndef WITHOUTMPI
     if(IOGROUPSIZE>0) then
        if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
           dummy_io=1
           call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
                & MPI_COMM_WORLD,info2)
        end if
     endif
#endif

  end if

  ! Loading sinks from the ICs (ic_sink or ic_sink_restart)
  if (nrestart>0)then
     nsinkold=nsink
     if(TRIM(initfile(levelmin)).NE.' ')then
        filename=TRIM(initfile(levelmin))//'/ic_sink_restart'
     else
        filename='ic_sink_restart'
     end if
     INQUIRE(FILE=filename, EXIST=ic_sink)
     if (myid==1)write(*,*)'Looking for file ic_sink_restart: ',filename
     if (.not. ic_sink)then
        filename='ic_sink_restart'
        INQUIRE(FILE=filename, EXIST=ic_sink)
     end if
  else
     nsink=0
     nindsink=0
     nsinkold=0
     if(TRIM(initfile(levelmin)).NE.' ')then
        filename=TRIM(initfile(levelmin))//'/ic_sink'
     else
        filename='ic_sink'
     end if
     INQUIRE(FILE=filename, EXIST=ic_sink)
     if (myid==1)write(*,*)'Looking for file ic_sink: ',filename
     if (.not. ic_sink)then
        filename='ic_sink'
        INQUIRE(FILE=filename, EXIST=ic_sink)
     end if
  end if

  if (ic_sink)then

     ! Wait for the token
#ifndef WITHOUTMPI
     if(IOGROUPSIZE>0) then
        if (mod(myid-1,IOGROUPSIZE)/=0) then
           call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag2,&
                & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
        end if
     endif
#endif

     open(10,file=filename,form='formatted')
     eof=.false.
     do
        read(10,*,end=103)sm1,sx1,sx2,sx3,sv1,sv2,sv3,sl1,sl2,sl3,sm2
        nsink=nsink+1
        nindsink=nindsink+1
        idsink(nsink)=nindsink
        msink(nsink)=sm1
        xsink(nsink,1)=sx1+boxlen/2.0
        xsink(nsink,2)=sx2+boxlen/2.0
        xsink(nsink,3)=sx3+boxlen/2.0
        vsink(nsink,1)=sv1
        vsink(nsink,2)=sv2
        vsink(nsink,3)=sv3
        lsink(nsink,1)=sl1
        lsink(nsink,2)=sl2
        lsink(nsink,3)=sl3
        tsink(nsink)=t
        new_born(nsink)=.false.
        msmbh(nsink)=sm2
     end do
     sinkint_level=levelmin
103  continue
     close(10)

     ! Send the token
#ifndef WITHOUTMPI
     if(IOGROUPSIZE>0) then
        if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
           dummy_io=1
           call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag2, &
                & MPI_COMM_WORLD,info2)
        end if
     endif
#endif

  end if

  ! Compute number of cloud particles within sink sphere
  call compute_ncloud_sink

  ! Output sink properties to screen
  if (myid==1.and.nsink-nsinkold>0)then
     write(*,*)'sinks read from file '//filename
     write(*,'("   Id           M             x             y             z            vx            vy            vz            lx            ly            lz       ")')
     write(*,'("======================================================================================================================================================")')
     do isink=nsinkold+1,nsink
        write(*,'(I8,2X,10(2X,E12.5))')idsink(isink),msink(isink),xsink(isink,1:ndim),&
             vsink(isink,1:ndim),lsink(isink,1:ndim)
     end do
  end if

  ! Set direct force boolean
  if(mass_sink_direct_force .ge. 0.0)then
     do isink=1,nsink
        direct_force_sink(isink)=(msink(isink) .ge. mass_sink_direct_force*2d33/(scale_d*scale_l**3))
     end do
  else
     do isink=1,nsink
        direct_force_sink(isink)=.False.
     end do
  endif

end subroutine init_sink
!################################################################
!################################################################
!################################################################
!################################################################
subroutine compute_ncloud_sink
  use amr_commons, only:dp
  use pm_commons, only:ir_cloud,ir_cloud_massive,ncloud_sink,ncloud_sink_massive
  real(dp)::xx,yy,zz,rr
  integer::ii,jj,kk
  ! Compute number of cloud particles
  ncloud_sink=0
  ncloud_sink_massive=0
  do kk=-2*ir_cloud,2*ir_cloud
     zz=dble(kk)/2.0
     do jj=-2*ir_cloud,2*ir_cloud
        yy=dble(jj)/2.0
        do ii=-2*ir_cloud,2*ir_cloud
           xx=dble(ii)/2.0
           rr=sqrt(xx*xx+yy*yy+zz*zz)
           if(rr<=dble(ir_cloud))ncloud_sink=ncloud_sink+1
           if(rr<=dble(ir_cloud_massive))ncloud_sink_massive=ncloud_sink_massive+1
        end do
     end do
  end do
end subroutine compute_ncloud_sink
