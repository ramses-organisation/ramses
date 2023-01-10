subroutine output_sink_csv(filename)
  use amr_commons
  use pm_commons
  implicit none
  character(LEN=80)::filename,fileloc

  integer::isink

  if(verbose)write(*,*)'Entering output_sink_csv'

  fileloc=TRIM(filename)
  open(unit=123,file=TRIM(fileloc),form='formatted',status='replace', recl=500)
  !======================
  ! Write sink properties
  !======================
  write(123,'(" # id,msink,x,y,z,vx,vy,vz,lx,ly,lz,tform,acc_rate,del_mass,rho_gas,cs**2,etherm,vx_gas,vy_gas,vz_gas,mbh,dmfsink,level ")')
  write(123,'(" # 1,m,l,l,l,l t**-1,l t**-1,l t**-1,m l**2 t**-1,m l**2 t**-1,m l**2 t**-1,t,m t**-1,m,m l**-3,l**2 t**-2,m l**2 t**-2,l t**-1,l t**-1,l t**-1,m,m,1")')
  do isink=1,nsink
     write(123,'(I10,21(A1,ES21.10),A1,I10)')idsink(isink),',',msink(isink),&
          ',',xsink(isink,1),',',xsink(isink,2),',',xsink(isink,3),&
          ',',vsink(isink,1),',',vsink(isink,2),',',vsink(isink,3),&
          ',',lsink(isink,1),',',lsink(isink,2),',',lsink(isink,3),&
          ',',tsink(isink),',',dMBHoverdt(isink),&
          ',',delta_mass(isink),&
          ',',rho_gas(isink),',',c2sink(isink),',',eps_sink(isink),&
          ',',vel_gas(isink,1),',',vel_gas(isink,2),',',vel_gas(isink,3),&
          ',',msmbh(isink),',',dmfsink(isink),',',sinkint_level
  end do

  close(123)

end subroutine output_sink_csv
