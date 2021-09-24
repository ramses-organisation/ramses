subroutine init_stellar
    use amr_commons
    use pm_commons
    use sink_feedback_module
    use mpi_mod
    implicit none

    integer:: ilun
    character(len=80):: fileloc
    character(len=5):: nchar, ncharcpu
    integer:: dummy_io, info2
    real(dp), allocatable, dimension(:):: xdp
    integer, allocatable, dimension(:):: xin
    integer, parameter:: tag = 1112
    integer:: nstellar_var, nstellar_var_tmp
    integer:: idim

    nstellar_var = ndim + 3 ! positions, mass, birth and life times

    ! Allocate all stellar object related quantities
    allocate(xstellar(1:nstellarmax, 1:ndim))
    allocate(mstellar(1:nstellarmax))
    allocate(tstellar(1:nstellarmax))
    allocate(ltstellar(1:nstellarmax))
    allocate(id_stellar(1:nstellarmax))
    
    ! Read restart variables from output files
    if(nrestart > 0) then
        ilun = 4*ncpu + myid + 11
        call title(nrestart, nchar)

        if(IOGROUPSIZEREP > 0) then
            call title(((myid - 1) / IOGROUPSIZEREP) + 1, ncharcpu)
            fileloc='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/stellar_'//TRIM(nchar)//'.out'
        else
            fileloc='output_'//TRIM(nchar)//'/stellar_'//TRIM(nchar)//'.out'
        end if

        call title(myid, nchar)
        fileloc = TRIM(fileloc) // TRIM(nchar)

        ! Wait for the token                                                                                                                                                                    
#ifndef WITHOUTMPI
        if(IOGROUPSIZE > 0) then
            if(mod(myid - 1, IOGROUPSIZE) /= 0) then
                call MPI_RECV(dummy_io, 1, MPI_INTEGER, myid-1-1, tag, &
                    & MPI_COMM_WORLD, MPI_STATUS_IGNORE, info2)
            end if
        end if
#endif

        ! to do: read csv instead
        open(unit=ilun, file=fileloc, form='unformatted')
        rewind(ilun)
        read(ilun) nstellar_var_tmp
        ! TODO: check that nstellar_var_tmp == nstellar_var
        read(ilun) nstellar

!        read(ilun) nstellar_tot

        if(nstellar > 0) then
            allocate(xdp(1:nstellar))
            allocate(xin(1:nstellar))

            ! Read stellar object position
            do idim = 1, ndim
                read(ilun) xdp
                xstellar(1:nstellar, idim) = xdp
            end do

            ! Read stellar object mass
            read(ilun) xdp
            mstellar(1:nstellar) = xdp

            ! Read stellar object birth time
            read(ilun) xdp
            tstellar(1:nstellar) = xdp

            ! Read stellar object life time
            read(ilun) xdp
            ltstellar(1:nstellar) = xdp

            ! Read stellar object sink particle id
            read(ilun) xin
            id_stellar(1:nstellar) = xin
        end if

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

    ! T.C. This should be in log file or in snapshots as a list until the time of snapshot
    ! Create file for HII region feedback logging
    if(myid == 1 .and. nrestart == 0) then
        open(104, file='hii.txt', form='formatted', status='unknown', position='append')
        write(104,*) 't ', 'x ', 'y ', 'z ', 'st_mass ', 'p_inj ', 'p_exp ', 'e_exp '
        close(104)
    end if

end subroutine init_stellar
