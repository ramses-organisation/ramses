subroutine write_gitinfo
  use amr_commons, ONLY:builddate,patchdir,gitrepo,gitbranch,githash

  builddate = BUILDDATE
  patchdir  = PATCH
  gitrepo   = GITREPO
  gitbranch = GITBRANCH
  githash   = GITHASH

  write(*,*)' '
  write(*,'(" compile date = ",A)')TRIM(builddate)
  write(*,'(" patch dir    = ",A)')TRIM(patchdir)
  write(*,'(" remote repo  = ",A)')TRIM(gitrepo)
  write(*,'(" local branch = ",A)')TRIM(gitbranch)
  write(*,'(" last commit  = ",A)')TRIM(githash)
  write(*,*)' '

end subroutine write_gitinfo
