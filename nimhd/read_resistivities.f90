subroutine read_resistivities
   use nimhd_parameters
   use nimhd_commons
   use constants, ONLY:mH

   integer::i,j,k
   real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
   !------------------------------------------
   ! Read resistivity tables for non-ideal MHD
   !------------------------------------------
   ! this should be refactored
   ! add parameter to choose how to compute the resistivity
   ! 1) fixed resistivity
   ! 2) analytical model resitivity(rho,T), Shu? Default option? Can be easily patched by the user
   ! 3) a table (like from Marchand 2016)
   ! make separate routine that computes the resistivity based on the chosen option
   ! move use_x3d to patch. We don't want to include a large table in ramses
   !TODO: What about if TEST?
   if(ntestDADM==1) return

   if(use_x3d==1)then
      open(42,file='marchand2016_table.dat',status='old')
      read(42,*) nchimie, tchimie, xichimie, nvarchimie
      read(42,*)
      read(42,*)
      allocate(resistivite_chimie_x(-2:nvarchimie+4,nchimie,tchimie,xichimie))
      do k=1,xichimie
      do i=1,tchimie
         do j=1,nchimie
            read(42,*)resistivite_chimie_x(-2:nvarchimie+4,j,i,k)
         end do
         read(42,*)
      end do
      end do
      close(42)

      rho_threshold=max(rho_threshold,resistivite_chimie_x(-2,1,1,1)*(mu_gas*mH)/scale_d) ! input in part/cc, output in code units
      nminchimie=(resistivite_chimie_x(-2,1,1,1))
      dnchimie=(log10(resistivite_chimie_x(-2,nchimie,1,1))-log10(resistivite_chimie_x(-2,1,1,1)))/&
               &(nchimie-1)
      tminchimie=(resistivite_chimie_x(-1,1,1,1))
      dtchimie=(log10(resistivite_chimie_x(-1,1,tchimie,1))-log10(resistivite_chimie_x(-1,1,1,1)))/&
               &(tchimie-1)
      ximinchimie=(resistivite_chimie_x(0,1,1,1))
      dxichimie=(log10(resistivite_chimie_x(0,1,1,xichimie))-log10(resistivite_chimie_x(0,1,1,1)))/&
               &(xichimie-1)
      call rq_3d
      call nimhd_4dtable
   else
      print*, 'must choose an input for abundances or resistivities'
      stop
   endif

end subroutine read_resistivities