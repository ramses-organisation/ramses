!NAME:  I/O PATCH FOR RAMSES
!
!PURPOSE: Optimization of I/O for run on large number of MPI tasks (>10 000).
!
!METHOD: MPI tasks are grouped to avoid concommitent writing/reading. Only one mpi task per group can write/read at a time. 
!
!CATEGORY: Patch
!
!CALLING: Adjust IO parameters in amr_parameters (later on I will put this parameters in the namelist)
!
!INPUT: Input can optionnally be read in subfolder named group_xxxxx
!      
!OPTIONAL: Can write or synchronize when IO are done. For large run it's better to create folders beforehand.
!          
!OUTPUT: Input could be read in subfolder named group_xxxxx
!       
!COMMON: 
!
!EXAMPLE:
!  integer::IOGROUPSIZE=64         ! Main snapshot (default 0)
!  integer::IOGROUPSIZECONE=8      ! Lightcone (default 0)
!  integer::IOGROUPSIZEREP=64      ! Subfolder size (default 0)
!  logical::withoutmkdir=.true.    !If true mkdir should be done before the run (default false)
!  logical::print_when_io=.true.   !If true print when IO (default false)
!  logical::synchro_when_io=.true. !If true synchronize when IO (default false)
!  
!
!HISTORY: P.Wautelet 2007-2010
!         Update and current version Y. Rasera 2010-2014
!
!REFERENCE: The MPI ticket was used in Horizon and DEUS grand challenge run. 
!           Teyssier et al,2009, Rasera et al, 2010, Alimi et al, 2012 
! 
!WARNING: For mhd run you should adapt the Makefile with init_flow_fine_mhd.f90,
!         init_hydro_mhd.f90, output_hydro_mhd.f90
!
!
!remarks to be removed: 
!-I changed the makefile this is not the default one
