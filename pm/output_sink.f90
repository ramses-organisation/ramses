subroutine output_sink_csv(filename)
  use amr_commons
  use pm_commons
  implicit none
  character(LEN=80)::filename,fileloc

  integer::isink

  if(verbose)write(*,*)'Entering output_sink_csv'

  fileloc=TRIM(filename)
  ! recl must hold the whole record: I10 + 21*(A1+ES25.16) + A1 + I10 = 567.
  ! ES25.16 is 17 significant digits, which is what an IEEE double needs to
  ! restart exactly.
  ! The explicit E3 matters. Without it Fortran drops the 'E' when the exponent
  ! needs three digits, writing e.g. 1.2345678901234567-123. Fortran reads that back
  ! correctly, but numpy cannot, so visu_ramses turns any such value into NaN
  ! when the test suite reads the file. E3 always emits the E.
  open(unit=123,file=TRIM(fileloc),form='formatted',status='replace', recl=1000)
  !======================
  ! Write sink properties
  !======================
  write(123,'(" # id,msink,x,y,z,vx,vy,vz,lx,ly,lz,tform,acc_rate,del_mass,rho_gas,cs**2,etherm,vx_gas,vy_gas,vz_gas,mbh,dmfsink,level ")')
  write(123,'(" # 1,m,l,l,l,l t**-1,l t**-1,l t**-1,m l**2 t**-1,m l**2 t**-1,m l**2 t**-1,t,m t**-1,m,m l**-3,l**2 t**-2,m l**2 t**-2,l t**-1,l t**-1,l t**-1,m,m,1")')
  do isink=1,nsink
     write(123,'(I10,21(A1,ES25.16E3),A1,I10)')idsink(isink),',',msink(isink),&
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
