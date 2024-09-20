#if USE_TURB==1
subroutine read_turb_fields
   use turb_commons
   implicit none

   integer            :: ilun, nturbtemp_tmp
   character(LEN=256) :: filedir
   logical            :: ok
   character(LEN=5)   :: nchar
   character(len=16)  :: precision_str_tmp
   real(dp)           :: aux
   character(len=17)  :: turb_last_time_z ! turb_last_time in hex
   character(len=17)  :: turb_next_time_z ! turb_next_time in hex

   ilun = myid+10
   call title(nrestart,nchar)
   filedir = 'output_'//TRIM(nchar)//'/'
   turb_file_header = trim(filedir)//'turb_fields.dat'
   inquire(file=turb_file_header, exist=ok)
   if (.not. ok)then
      write(*,*)'Restart failed:'
      write(*,*)'File '//TRIM(turb_file_header)//' not found'
      call clean_stop
   end if

   ! Read header file of important information
   open(ilun,file=turb_file_header,status="old",form="formatted")
   read(ilun,*) nturbtemp_tmp ! Not used, SEREN compatibility
   read(ilun,*) precision_str_tmp
   read(ilun,*) turb_file_last
   read(ilun,*) aux, turb_last_time_z
   read(ilun,*) turb_file_next
   read(ilun,*) aux, turb_next_time_z
   read(ilun,*); read(ilun,*); read(ilun,*); read(ilun,*); read(ilun,*)
   read(ilun,*) kiss64_state
   close(ilun)
   ! We don't care about some of the stuff in the file
   ! as it should be set in the params file

   read(turb_last_time_z, '(Z16)') turb_last_time
   read(turb_next_time_z, '(Z16)') turb_next_time

   ! Read turbulent fields
   open(ilun,file=trim(filedir)//trim(turb_file_last),status="unknown",&
        &access="stream",form="unformatted")
   read(ilun) turb_last
   close(ilun)

   open(ilun,file=trim(filedir)//trim(turb_file_next),status="unknown",&
        access="stream",form="unformatted")
   read(ilun) turb_next
   close(ilun)

end subroutine read_turb_fields
#endif
