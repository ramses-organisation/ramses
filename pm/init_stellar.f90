subroutine init_stellar
    use amr_commons
    use pm_commons
    use sink_feedback_parameters
    use mpi_mod
    implicit none
#ifndef WITHOUTMPI
    integer, parameter:: tag=1112
    integer:: dummy_io, info2
#endif
    integer:: ilun
    logical::eof=.false.
    character(len=80):: fileloc
    character(len=5):: nchar, ncharcpu
    integer:: idim
    integer::sid
    real(dp)::sm,sx1,sx2,sx3,stform,stlife
    character::co
    character(LEN=200)::comment_line

    if(.not. stellar) return

    ! Allocate all stellar object related quantities
    allocate(xstellar(1:nstellarmax, 1:ndim))
    allocate(mstellar(1:nstellarmax))
    allocate(tstellar(1:nstellarmax))
    allocate(ltstellar(1:nstellarmax))
    allocate(id_stellar(1:nstellarmax))
    
    ! Load stellar particles from the restart
    if(nrestart > 0) then
        ilun = 4*ncpu + myid + 11
        call title(nrestart, nchar)

        if(IOGROUPSIZEREP > 0) then
            call title(((myid - 1) / IOGROUPSIZEREP) + 1, ncharcpu)
            fileloc='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/stellar_'//TRIM(nchar)//'.csv'
        else
            fileloc='output_'//TRIM(nchar)//'/stellar_'//TRIM(nchar)//'.csv'
        end if

        ! Wait for the token                                                                                                                                                                    
#ifndef WITHOUTMPI
        if(IOGROUPSIZE > 0) then
            if(mod(myid - 1, IOGROUPSIZE) /= 0) then
                call MPI_RECV(dummy_io, 1, MPI_INTEGER, myid-1-1, tag, &
                    & MPI_COMM_WORLD, MPI_STATUS_IGNORE, info2)
            end if
        end if
#endif

        nstellar=0
        open(ilun, file=fileloc, form='formatted')
        eof=.false.
        ! scrolling over the comment lines
        read(ilun,'(A200)')comment_line
        read(ilun,'(A200)')comment_line
        do
            read(ilun,'(I10,6(A1,ES21.10))',end=104)sid,co,sm,co,&
                                sx1,co,sx2,co,sx3,co, &
                                stform,co,stlife
            nstellar=nstellar+1
            id_stellar(nstellar)=sid
            mstellar(nstellar)=sm
            xstellar(nstellar,1)=sx1
            xstellar(nstellar,2)=sx2
            xstellar(nstellar,3)=sx3
            tstellar(nstellar)=stform
            ltstellar(nstellar)=stlife
        end do
   104  continue
        close(ilun)

        ! Send the token                                                                                                                                                                        
#ifndef WITHOUTMPI
        if(IOGROUPSIZE > 0) then
            if(mod(myid, IOGROUPSIZE) /=0 .and. (myid < ncpu)) then
                dummy_io = 1
                call MPI_SEND(dummy_io, 1, MPI_INTEGER, myid-1+1, tag, &
                    & MPI_COMM_WORLD, info2)
            end if
        end if
#endif
    end if

    ! TODO: add possibility to manually add stellar particles on restart as for sinks?
    !       replaces fixes sources in namelist

end subroutine init_stellar
