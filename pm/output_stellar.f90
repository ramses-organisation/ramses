subroutine output_stellar_csv(filename)
  use amr_commons
  use pm_commons

  use sink_feedback_module

  implicit none
  character(LEN=80)::filename,fileloc

  integer::ilun,icpu,istellar

  if(verbose)write(*,*)'Entering output_stellar_csv'

  ilun=2*ncpu+myid+10

  fileloc=TRIM(filename)
  open(unit=ilun,file=TRIM(fileloc),form='formatted',status='replace', recl=500)
  !======================
  ! Write stellar properties
  !======================
  do istellar=1,nstellar

     write(ilun,'(I10,6(A1,ES20.10))')id_stellar(istellar),',',mstellar(istellar),&
          ',',xstellar(istellar,1),',',xstellar(istellar,2),',',xstellar(istellar,3),&
          ',',tstellar(istellar),',',ltstellar(istellar)
  end do

  close(ilun)

end subroutine output_stellar_csv
