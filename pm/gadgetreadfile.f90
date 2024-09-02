MODULE gadgetreadfilemod
    !
    ! Routines to read L-Gadget2 particle data files in Fortran
    ! Adapted from http://astro.dur.ac.uk/~jch/password_pages/index.html
    !
    ! CHARACTER(LEN=*) :: basefile - the base snapshot name. It may get a file id added
    !
    !
    !
    !
    ! Data type corresponding to gadget file header
    TYPE gadgetheadertype
    INTEGER(4), DIMENSION(6) :: npart
    REAL(8), DIMENSION(6) :: mass
    real(8) :: time
    real(8) :: redshift
    integer(4) :: flag_sfr
    integer(4) :: flag_feedback
    INTEGER(4), DIMENSION(6) :: nparttotal
    integer(4) :: flag_cooling
    integer(4) :: numfiles
    real(8) :: boxsize
    real(8) :: omega0
    real(8) :: omegalambda
    real(8) :: hubbleparam
    integer(4) :: flag_stellarage
    integer(4) :: flag_metals
    INTEGER(4), DIMENSION(6)  :: totalhighword
    integer(4) :: flag_entropy_instead_u
    integer(4) :: flag_doubleprecision
    integer(4) :: flag_ic_info
    real(4) :: lpt_scalingfactor
    CHARACTER, DIMENSION(48) :: unused
    END TYPE gadgetheadertype

CONTAINS

    ! ---------------------------------------------------------------------------

    SUBROUTINE gadgetreadheader(basename, ifile, header, ok)
        !
        ! Read and return the gadget file header for the specified file
        !

#ifndef WITHOUTMPI
        use amr_commons, only:myid, IOGROUPSIZE, ncpu
#endif
        use mpi_mod
        implicit none
        ! Input parameters
        CHARACTER(LEN=*), INTENT(IN) :: basename
        INTEGER, INTENT(IN) :: ifile
        ! Header to return
        TYPE (gadgetheadertype), INTENT(OUT) :: header
        logical, INTENT(OUT) :: ok
        ! Internal variables
        character(len=256) :: filename
        character(len=6) :: fileno
        integer, parameter :: tag=1104
#ifndef WITHOUTMPI
        integer :: dummy_io, info2
#endif
        filename = TRIM(basename)
        INQUIRE(file = filename, exist = ok)
        if (.not. ok) then
            ! Generate the number to go on the end of the filename
            IF (ifile < 10) THEN
                WRITE(fileno,'(".",1i1.1)') ifile
            ELSE IF (ifile < 100) THEN
                WRITE(fileno,'(".",1i2.2)') ifile
            ELSE IF (ifile < 1000) THEN
                WRITE(fileno,'(".",1i3.3)') ifile
            ELSE IF (ifile < 10000) THEN
                WRITE(fileno,'(".",1i4.4)') ifile
            ELSE
                WRITE(fileno,'(".",1i5.5)') ifile
            END IF
            filename = TRIM(basename) // fileno
            INQUIRE(file = filename, exist = ok)
            if (.not. ok) then
                write(*, *) 'No file '// basename //' or '// filename
                RETURN
            end if
        end if

        ! Wait for the token
#ifndef WITHOUTMPI
        if (IOGROUPSIZE > 0) then
            if (mod(myid - 1, IOGROUPSIZE) /= 0) then
                call MPI_RECV(dummy_io, 1, MPI_INTEGER, myid - 1 - 1, tag,&
                    & MPI_COMM_WORLD, MPI_STATUS_IGNORE, info2)
            end if
        end if
#endif

        OPEN(unit=1, file=filename, status='old', action='read', form='unformatted')
        ! Byte swapping doesn't work if you just do READ(1)header
        READ(1) header%npart, header%mass, header%time, header%redshift, &
            header%flag_sfr, header%flag_feedback, header%nparttotal, &
            header%flag_cooling, header%numfiles, header%boxsize, &
            header%omega0, header%omegalambda, header%hubbleparam, &
            header%flag_stellarage, header%flag_metals, header%totalhighword, &
            header%flag_entropy_instead_u, header%flag_doubleprecision, &
            header%flag_ic_info, header%lpt_scalingfactor
        CLOSE(1)

        ! Send the token
#ifndef WITHOUTMPI
        if (IOGROUPSIZE > 0) then
            if (mod(myid, IOGROUPSIZE) /= 0 .and. (myid < ncpu)) then
                dummy_io = 1
                call MPI_SEND(dummy_io, 1, MPI_INTEGER, myid - 1 + 1, tag, &
                    & MPI_COMM_WORLD, info2)
            end if
        end if
#endif


    END SUBROUTINE gadgetreadheader

    ! ---------------------------------------------------------------------------

    SUBROUTINE gadgetreadfile(basename, ifile, header, pos, vel, id)
        !
        ! Read and return all data from the specified file. Output arrays must
        ! already be allocated. Use readheader to get particle numbers to do this.
        !
#ifndef WITHOUTMPI
        use amr_commons, only:myid, IOGROUPSIZE, ncpu
#endif
        use mpi_mod

        IMPLICIT NONE
        ! Input parameters
        CHARACTER(LEN=*), INTENT(IN) :: basename
        INTEGER, INTENT(IN) :: ifile
        ! Header and hash table to return
        TYPE (gadgetheadertype) :: header
        ! Particle data
        REAL, DIMENSION(3,*) :: pos, vel
#ifndef LONGINT
        INTEGER(4), DIMENSION(*) :: id
#else
        INTEGER(8), DIMENSION(*) :: id
#endif
        ! Internal variables
        character(len=256) :: filename
        character(len=6) :: fileno
        integer :: np
        logical :: ok
        integer, parameter :: tag=1105
#ifndef WITHOUTMPI
        integer :: dummy_io, info2
#endif
        ! Generate the number to go on the end of the filename
        IF (ifile < 10) THEN
            WRITE(fileno,'(".",1i1.1)') ifile
        ELSE IF (ifile < 100) THEN
            WRITE(fileno,'(".",1i2.2)') ifile
        ELSE IF (ifile < 1000) THEN
            WRITE(fileno,'(".",1i3.3)') ifile
        ELSE IF (ifile < 10000) THEN
            WRITE(fileno,'(".",1i4.4)') ifile
        ELSE
            WRITE(fileno,'(".",1i5.5)') ifile
        END IF

        filename = TRIM(basename) // fileno

        INQUIRE(file = filename, exist = ok)

        if (.not. ok) then
            write(*, *) 'No file '// filename
            RETURN
        end if

        ! Wait for the token (this token might be moved to init_part for best performance)
#ifndef WITHOUTMPI
        if (IOGROUPSIZE > 0) then
            if (mod(myid - 1, IOGROUPSIZE) /= 0) then
                call MPI_RECV(dummy_io, 1, MPI_INTEGER, myid - 1 - 1, tag,&
                    & MPI_COMM_WORLD, MPI_STATUS_IGNORE, info2)
            end if
        end if
#endif


        OPEN(unit=1, file=filename, status='old', action='read', form='unformatted')
        ! Byte swapping doesn't appear to work if you just do READ(1)header
        READ(1) header%npart, header%mass, header%time, header%redshift, &
            header%flag_sfr, header%flag_feedback, header%nparttotal, &
            header%flag_cooling, header%numfiles, header%boxsize, &
            header%omega0, header%omegalambda, header%hubbleparam, &
            header%flag_stellarage, header%flag_metals, header%totalhighword, &
            header%flag_entropy_instead_u, header%flag_doubleprecision, &
            header%flag_ic_info, header%lpt_scalingfactor
        np = header%npart(2)
        READ(1) pos(1:3, 1:np)
        READ(1) vel(1:3, 1:np)
        READ(1) id(1:np)
        CLOSE(1)

        ! Send the token
#ifndef WITHOUTMPI
        if (IOGROUPSIZE > 0) then
            if (mod(myid, IOGROUPSIZE) /= 0 .and. (myid < ncpu)) then
                dummy_io = 1
                call MPI_SEND(dummy_io, 1, MPI_INTEGER, myid - 1 + 1, tag, &
                    & MPI_COMM_WORLD, info2)
            end if
        end if
#endif

    END SUBROUTINE gadgetreadfile

    ! ---------------------------------------------------------------------------

    SUBROUTINE gadgetwritefile(basename, ifile, header, pos, vel, id)
        !
        ! Read and return all data from the specified file. Output arrays must
        ! already be allocated. Use readheader to get particle numbers to do this.
        !
#ifndef WITHOUTMPI
        use amr_commons, only:myid, IOGROUPSIZE, ncpu
#endif
        use mpi_mod

        IMPLICIT NONE
        ! Input parameters
        CHARACTER(LEN=*), INTENT(IN) :: basename
        INTEGER, INTENT(IN) :: ifile
        ! Header and hash table to return
        TYPE (gadgetheadertype) :: header
        ! Particle data
        REAL, DIMENSION(3,*) :: pos, vel
#ifndef LONGINT
        INTEGER(4), DIMENSION(*) :: id
#else
        INTEGER(8), DIMENSION(*) :: id
#endif
        ! Internal variables
        character(len=256) :: filename
        character(len=6) :: fileno
        integer :: np
        integer, parameter :: tag=1106
#ifndef WITHOUTMPI
        integer :: dummy_io, info2
#endif
        ! Generate the number to go on the end of the filename
        IF (ifile < 10) THEN
            WRITE(fileno,'(".",1i1.1)') ifile
        ELSE IF (ifile < 100) THEN
            WRITE(fileno,'(".",1i2.2)') ifile
        ELSE IF (ifile < 1000) THEN
            WRITE(fileno,'(".",1i3.3)') ifile
        ELSE IF (ifile < 10000) THEN
            WRITE(fileno,'(".",1i4.4)') ifile
        ELSE
            WRITE(fileno,'(".",1i5.5)') ifile
        END IF

        filename = TRIM(basename) // fileno

        ! Wait for the token
#ifndef WITHOUTMPI
        if (IOGROUPSIZE > 0) then
            if (mod(myid - 1, IOGROUPSIZE) /= 0) then
                call MPI_RECV(dummy_io, 1, MPI_INTEGER, myid - 1 - 1, tag,&
                    & MPI_COMM_WORLD, MPI_STATUS_IGNORE, info2)
            end if
        end if
#endif

        OPEN(unit=1, file=filename, status='unknown', action='write', form='unformatted')
        WRITE(1) header%npart, header%mass, header%time, header%redshift, &
            header%flag_sfr, header%flag_feedback, header%nparttotal, &
            header%flag_cooling, header%numfiles, header%boxsize, &
            header%omega0, header%omegalambda, header%hubbleparam, &
            header%flag_stellarage, header%flag_metals, header%totalhighword, &
            header%flag_entropy_instead_u, header%flag_doubleprecision, &
            header%flag_ic_info, header%lpt_scalingfactor, header%unused
        np = header%npart(2)
        WRITE(1) pos(1:3, 1:np)
        WRITE(1) vel(1:3, 1:np)
        WRITE(1) id(1:np)

        CLOSE(1)

        ! Send the token
#ifndef WITHOUTMPI
        if (IOGROUPSIZE > 0) then
            if (mod(myid, IOGROUPSIZE) /= 0 .and. (myid < ncpu)) then
                dummy_io = 1
                call MPI_SEND(dummy_io, 1, MPI_INTEGER, myid - 1 + 1, tag, &
                    & MPI_COMM_WORLD, info2)
            end if
        end if
#endif

    END SUBROUTINE gadgetwritefile
END MODULE gadgetreadfilemod
