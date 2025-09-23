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
   MPIF90 = mpif90
endif

#
ifeq ($(MACHINE),meluxina)
   MPIF90 = mpif90
   ifeq ($(COMPILER),GNU)
      FFLAGS_OPT += -march=native
   endif
endif

#
ifeq ($(MACHINE),vega)
   MPIF90 = mpif90
endif
