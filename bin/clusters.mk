#############################################################################
# This part of the Makefile contains machine specific compiler options
#############################################################################

# Default "generic" cluster: no special tuning
ifeq ($(MACHINE),generic)
   MPIF90 = mpif90
endif

####### EuroHPC clusters #######

#
ifeq ($(MACHINE),discoverer)
   MPIF90 = mpif90
endif

#
ifeq ($(MACHINE),karolina)
   MPIF90 = mpif90
endif

#
ifeq ($(MACHINE),leonardo)
   MPIF90 = mpiifort
endif

#
ifeq ($(MACHINE),lumi)
   MPIF90 = ftn
   # mismatch flag needed to circumvent MPI_ALLREDUCE error
   FFLAGS_BASE += -fallow-argument-mismatch
endif

#
ifeq ($(MACHINE),marenostrum)
   ifeq ($(COMPILER),INTEL)
      MPIF90 = mpif90
   else ifeq ($(COMPILER),INTELX)
      MPIF90 = mpiifx
   endif
endif

#
ifeq ($(MACHINE),meluxina)
   ifeq ($(COMPILER),GNU)
      MPIF90 = mpif90
      FFLAGS_OPT += -march=native
   else ifeq ($(COMPILER),INTELX)
      MPIF90 = mpiifx
   endif
   # old intel ifort compiler not available on meluxina
endif

#
ifeq ($(MACHINE),vega)
   MPIF90 = mpif90
endif
