module buildinfo
  implicit none

  ! Variables for executable identification
  ! Set at compile time (via -D flags)
  character(len=300),parameter::builddate    = BUILDDATE
  character(len=300),parameter::buildcommand = BUILDCOMMAND
  character(len=300),parameter::patchdir     = PATCH
  character(len=300),parameter::gitrepo      = GITREPO
  character(len=300),parameter::gitbranch    = GITBRANCH
  character(len=300),parameter::githash      = GITHASH

contains

  subroutine write_gitinfo(unit)
    integer,intent(in)::unit

    write(unit,'(" compile date    = ",A)')TRIM(builddate)
    write(unit,'(" compile command = ",A)')TRIM(buildcommand)
    write(unit,'(" patch dir       = ",A)')TRIM(patchdir)
    write(unit,'(" remote repo     = ",A)')TRIM(gitrepo)
    write(unit,'(" local branch    = ",A)')TRIM(gitbranch)
    write(unit,'(" last commit     = ",A)')TRIM(githash)

  end subroutine write_gitinfo

end module buildinfo
