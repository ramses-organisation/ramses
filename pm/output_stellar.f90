subroutine output_stellar_csv(filename)
  use amr_commons
  use pm_commons
  use sink_feedback_parameters
  implicit none
  character(LEN=80)::filename,fileloc
  integer::istellar

  if(verbose)write(*,*)'Entering output_stellar_csv'

  fileloc=TRIM(filename)
  open(unit=123,file=TRIM(fileloc),form='formatted',status='replace', recl=500)
  !======================
  ! Write stellar properties
  !======================
  write(123,'(" # id,mstellar,tform,tlife ")')
  write(123,'(" # 1,m,t,t")')
  do istellar=1,nstellar
     write(123,'(I10,3(A1,ES21.10))')id_stellar(istellar),',',mstellar(istellar),&
          ',',tstellar(istellar),',',ltstellar(istellar)
  end do

  close(123)

end subroutine output_stellar_csv
