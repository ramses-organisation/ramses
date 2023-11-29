#############################################################################
# If you have problems with this makefile, contact Romain.Teyssier@gmail.com
#############################################################################
# Compilation time parameters

# Do we want a debug build? 1=Yes, 0=No
DEBUG = 0
# Compiler flavor: GNU or INTEL
COMPILER = GNU
# Size of vector cache
NVECTOR = 32
# Number of dimensions
NDIM = 1
# Float precision size
NPRE = 8
# hydro/mhd/rhd solver
SOLVER = hydro
# Patch
PATCH =
# Use RT? 1=Yes, 0=No
RT = 0
# Use turbulence? 1=Yes, 0=No (requires fftw3)
USE_TURB = 0
# Use MPI? 1=Yes, 0=No
MPI = 0
MPIF90 = mpif90
# Root name of executable
EXEC = ramses
# Number of additional energies
NENER = 0
# Use Grackle cooling? 1=Yes, 0=No
GRACKLE = 0
# Number of metal species
NMETALS = 0
# Use ATON? 1=Yes, 0=No
ATON = 0
# Number of ions for RT
# use 1 for HI+HII, +1 for added H2, +2 for added HeI, HeII, HeIII
NIONS = 0
# Number of photon groups for RT
NGROUPS = 0
# Number of passive scalars
NPSCAL = 0
# Light MPI communicator structure (for > 10k MPI processes)
LIGHT_MPI_COMM = 0

# Compute NVAR
NVAR = 2+$(NDIM)
ifeq ($(SOLVER),mhd)
   NVAR = 8
endif
NVAR := $(NVAR)+$(NENER)+$(NPSCAL)+$(NMETALS)
ifeq ($(RT),1)
   NVAR := $(NVAR)+$(NIONS)
endif

# Set to one to use 'include "mpif.h"' instead of more recent "use mpi"
OLD_MPI_SUPPORT = 0
#############################################################################
GITBRANCH = $(shell git rev-parse --abbrev-ref HEAD)
GITHASH = $(shell git log --pretty=format:'%H' -n 1)
GITREMOTE = $(shell git config --get branch.$(GITBRANCH).remote)
GITREPO = $(shell git config --get remote.$(GITREMOTE).url)
BUILDDATE = $(shell date +"%D-%T")
DEFINES = -DNVECTOR=$(NVECTOR) -DNDIM=$(NDIM) -DNPRE=$(NPRE) -DNENER=$(NENER) \
          -DNVAR=$(NVAR) -DSOLVER$(SOLVER)
ifeq ($(ATON),1)
   DEFINES += -DATON
endif
ifeq ($(GRACKLE),1)
   DEFINES += -Dgrackle
endif
ifeq ($(OLD_MPI_SUPPORT),1)
   DEFINES += -DMPI_OLD
endif
ifeq ($(RT),1)
   DEFINES += -DRT -DNIONS=$(NIONS) -DNGROUPS=$(NGROUPS)
endif
ifeq ($(USE_TURB),1)
   DEFINES += -DUSE_TURB
endif
ifeq ($(LIGHT_MPI_COMM), 1)
   DEFINES += -DLIGHT_MPI_COMM
endif
#############################################################################
# Fortran compiler options and directives

# GNU compiler (gfortran)
ifeq ($(COMPILER),GNU)
   FFLAGS = -cpp $(DEFINES)
   ifeq ($(MPI),1)
      F90 = $(MPIF90)
      FC = gfortran
   else
      F90 = gfortran
      FC = gfortran
      FFLAGS += -DWITHOUTMPI
   endif
   F90 += -ffree-line-length-none -fimplicit-none
   FC += -ffree-line-length-none -fimplicit-none
   ifeq ($(DEBUG),1)
      F90 += -g -O0 -fbacktrace -fbounds-check -Wuninitialized -Wall
      FFLAGS += -ffpe-trap=zero,underflow,overflow,invalid -finit-real=nan
   else
      F90 += -O3
   endif
endif

# Intel compiler
ifeq ($(COMPILER),INTEL)
   FFLAGS = -cpp $(DEFINES)
   ifeq ($(MPI),1)
      F90 = $(MPIF90)
      FC = ifort
      FFLAGS += -DNOSYSTEM
   else
      F90 = ifort
      FC = ifort
      FFLAGS += -DWITHOUTMPI
   endif
   F90 += -fp-model source
   FC += -fp-model source
   ifeq ($(DEBUG),1)
      F90 += -warn all -O0 -g -traceback -check bounds
      FFLAGS += -fpe0 -ftrapuv -init=zero -init=snan -init=arrays
   else
      F90 += -O3
   endif
endif

#############################################################################
MOD = mod
#############################################################################
# MPI librairies
LIBMPI =
#LIBMPI = -lfmpi -lmpi -lelan

# --- CUDA libraries, for Titane ---
LIBCUDA = -L/opt/cuda/lib  -lm -lcuda -lcudart

ifeq ($(GRACKLE),1)
   # Add include and library install path for grackle and hdf5 here
   LIBS_GRACKLE = -L$(HOME)/local/lib -lgrackle -lhdf5 -lz -lgfortran -ldl
   LIBS_OBJ     = -I$(HOME)/local/include -DCONFIG_BFLOAT_8 -DH5_USE_16_API -fPIC
endif
ifeq ($(USE_TURB),1)
   LIBS_TURB = -L$(HOME)/local/lib -lfftw3
   LIBS_OBJ_TURB = -I$(HOME)/local/include -lfftw3
endif
LIBS = $(LIBMPI) $(LIBS_GRACKLE) $(LIBS_TURB)
#############################################################################
# Sources directories are searched in this exact order
VPATH = $(PATCH):../$(SOLVER):../aton:
ifeq ($(RT),1)
   VPATH += ../rt:
endif
ifeq ($(USE_TURB),1)
   VPATH += ../turb:
endif
VPATH += ../hydro:../pm:../poisson:../amr:../io
#############################################################################
# All objects
MODOBJ = mpi_mod.o amr_parameters.o amr_commons.o random.o pm_parameters.o sink_feedback_parameters.o \
         pm_commons.o poisson_parameters.o dump_utils.o constants.o file_module.o
ifeq ($(GRACKLE),1)
   MODOBJ += grackle_parameters.o
endif
MODOBJ += poisson_commons.o hydro_parameters.o hydro_commons.o \
          cooling_module.o bisection.o sparse_mat.o clfind_commons.o \
          gadgetreadfile.o write_makefile.o write_patch.o write_gitinfo.o \
          sink_sn_feedback.o
ifeq ($(RT),1)
   MODOBJ += rt_parameters.o rt_hydro_commons.o coolrates_module.o \
             rt_spectra.o rt_cooling_module.o rt_flux_module.o
endif
ifeq ($(USE_TURB),1)
   MODOBJ += turb_parameters.o turb_commons.o
endif

AMROBJ = read_params.o init_amr.o init_time.o init_refine.o tracer_utils.o adaptive_loop.o \
         amr_step.o update_time.o output_amr.o flag_utils.o \
         physical_boundaries.o virtual_boundaries.o refine_utils.o \
         nbors_utils.o hilbert.o load_balance.o title.o sort.o cooling_fine.o \
         eos.o units.o light_cone.o movie.o memory.o end.o
# Particle-Mesh objects
PMOBJ = init_part.o output_part.o rho_fine.o synchro_fine.o move_fine.o \
        newdt_fine.o particle_tree.o add_list.o remove_list.o star_formation.o \
        sink_particle.o feedback.o clump_finder.o clump_merger.o output_clump.o \
        flag_formation_sites.o init_sink.o output_sink.o \
        unbinding.o merger_tree.o move_tracer.o init_tracer.o \
        read_sink_feedback_params.o sink_rt_feedback.o \
        stellar_particle.o init_stellar.o output_stellar.o
# Poisson solver objects
POISSONOBJ = init_poisson.o phi_fine_cg.o interpol_phi.o force_fine.o \
             multigrid_coarse.o multigrid_fine_commons.o multigrid_fine_fine.o \
             multigrid_fine_coarse.o gravana.o boundary_potential.o rho_ana.o \
             output_poisson.o
# Hydro objects
HYDROOBJ = init_hydro.o init_flow_fine.o write_screen.o output_hydro.o \
           courant_fine.o godunov_fine.o uplmde.o umuscl.o interpol_hydro.o \
           godunov_utils.o condinit.o hydro_flag.o hydro_boundary.o boundana.o \
           read_hydro_params.o synchro_hydro_fine.o cooling_module_ism.o
# RT objects
RTOBJ = rt_init_hydro.o rt_init_xion.o rt_init.o rt_init_flow_fine.o \
        rt_output_hydro.o rt_godunov_fine.o rt_interpol_hydro.o \
        rt_godunov_utils.o rt_condinit.o rt_hydro_flag.o rt_hydro_boundary.o \
        rt_boundana.o rt_units.o
# Turbulence objects
TURBOBJ = turb_force_utils.o read_turb_params.o init_turb.o mpi_share_turb_fields.o \
          turb_next_field.o turb_check_time.o write_turb_fields.o read_turb_fields.o \
          add_turb_forcing.o
# Patch objects
sinclude $(PATCH)/Makefile

# All objects
AMRLIB = $(AMROBJ) $(HYDROOBJ) $(PMOBJ) $(POISSONOBJ)
ifeq ($(RT),1)
   AMRLIB += $(RTOBJ)
endif
ifeq ($(USE_TURB),1)
   AMRLIB += $(TURBOBJ)
endif
# ATON objects
ATON_MODOBJ = timing.o radiation_commons.o rad_step.o
ATON_OBJ = observe.o init_radiation.o rad_init.o rad_boundary.o rad_stars.o \
           rad_backup.o ../aton/atonlib/libaton.a
#############################################################################
ramses:	$(MODOBJ) $(AMRLIB) ramses.o
	$(F90) $(MODOBJ) $(AMRLIB) ramses.o -o $(EXEC)$(NDIM)d $(LIBS)
	rm write_makefile.f90
	rm write_patch.f90
ramses_aton: $(MODOBJ) $(ATON_MODOBJ) $(AMRLIB) $(ATON_OBJ) ramses.o
	$(F90) $(MODOBJ) $(ATON_MODOBJ) $(AMRLIB) $(ATON_OBJ) ramses.o \
		-o $(EXEC)$(NDIM)d $(LIBS) $(LIBCUDA)
	rm write_makefile.f90
	rm write_patch.f90
#############################################################################
write_gitinfo.o: FORCE
	$(FC) -O0 -cpp -DPATCH=\'$(PATCH)\' -DGITBRANCH=\'$(GITBRANCH)\' \
		-DGITHASH=\'"$(GITHASH)"\' -DGITREPO=\'$(GITREPO)\' \
		-DBUILDDATE=\'"$(BUILDDATE)"\' -c ../amr/write_gitinfo.f90 -o $@
write_makefile.o: FORCE
	../utils/scripts/cr_write_makefile.sh $(MAKEFILE_LIST)
	$(FC) -O0 -c write_makefile.f90 -o $@
write_patch.o: FORCE
	../utils/scripts/cr_write_patch.sh $(PATCH)
	$(FC) -O0 -c write_patch.f90 -o $@
%.o:%.F
	$(F90) $(FFLAGS) -c $^ -o $@ $(LIBS_OBJ) $(LIBS_OBJ_TURB)
%.o:%.f90
	$(F90) $(FFLAGS) -c $^ -o $@ $(LIBS_OBJ) $(LIBS_OBJ_TURB)
FORCE:
#############################################################################
clean:
	rm -f *.o *.$(MOD) *.i
#############################################################################
