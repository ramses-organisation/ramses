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
   ifeq ($(COMPILER),GNU)
      FFLAGS_OPT += -march=native
      # link-time optimization accross files
      LDFLAGS += -flto -fwhole-program
      FFLAGS_OPT += -flto -fwhole-program
   endif
endif

#
ifeq ($(MACHINE),karolina)
   MPIF90 = mpif90
   ifeq ($(COMPILER),GNU)
      FFLAGS_OPT += -march=native
      # link-time optimization accross files
      LDFLAGS += -flto -fwhole-program
      FFLAGS_OPT += -flto -fwhole-program
   endif
endif

#
ifeq ($(MACHINE),leonardo)
   MPIF90 = mpiifort
   ifeq ($(COMPILER),INTEL)
      # link-time optimization and inlining accross files
      LDFLAGS += -ipo -inline-forceinline
      FFLAGS_OPT += -ipo -inline-forceinline
   endif
endif

#
ifeq ($(MACHINE),lumi)
   MPIF90 = ftn
   # mismatch flag needed to circumvent MPI_ALLREDUCE error
   FFLAGS_BASE += -fallow-argument-mismatch
   # link-time optimization accross files
   LDFLAGS += -flto -fwhole-program
   FFLAGS_OPT += -flto -fwhole-program
endif

#
ifeq ($(MACHINE),marenostrum)
   MPIF90 = mpif90
   ifeq ($(COMPILER),INTEL)
      # link-time optimization and inlining accross files
      LDFLAGS += -ipo -inline-forceinline
      FFLAGS_OPT += -ipo -inline-forceinline
   endif
endif

#
ifeq ($(MACHINE),meluxina)
   MPIF90 = mpif90
   ifeq ($(COMPILER),GNU)
      FFLAGS_OPT += -march=native
      # link-time optimization accross files
      LDFLAGS += -flto -fwhole-program
      FFLAGS_OPT += -flto -fwhole-program
   endif
endif

#
ifeq ($(MACHINE),vega)
   MPIF90 = mpif90
   ifeq ($(COMPILER),GNU)
      FFLAGS_OPT += -march=native
      # link-time optimization accross files
      LDFLAGS += -flto -fwhole-program
      FFLAGS_OPT += -flto -fwhole-program
   endif
endif
