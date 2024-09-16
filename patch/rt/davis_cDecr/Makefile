#############################################################################
# If you have problems with this makefile, contact Romain.Teyssier@cea.fr
#############################################################################
# Compilation time parameters
NVECTOR = 64
NDIM = 2
NPRE = 8
NVAR = 9
NENER=1
SOLVER = hydro
PATCH = ../patch/rt/rt_pressure ../patch/rt/davis  ../patch/rt/davis_cDecr
EXEC = ramses
ATON_FLAGS = #-DATON  # Uncomment to enable ATON.
NGROUPS = 1
DEFINES = -DNVECTOR=$(NVECTOR) -DNDIM=$(NDIM) -DNPRE=$(NPRE) -DNVAR=$(NVAR) -DSOLVER$(SOLVER) -DRT -DNGROUPS=$(NGROUPS) -DNENER=$(NENER) $(ATON_FLAGS)
#############################################################################
# Fortran compiler options and directives
# Fortran compiler options and directives
F90 = mpif90
FFLAGS = -O3 -fpp $(DEFINES) -DNOSYSTEM

#############################################################################
MOD = mod
#############################################################################
# MPI librairies
LIBMPI =
#LIBMPI = -lfmpi -lmpi -lelan

# --- CUDA libraries, for Titane ---
LIBCUDA = -L/opt/cuda/lib  -lm -lcuda -lcudart

LIBS = $(LIBMPI)
#############################################################################
# Sources directories are searched in this exact order
VPATH = $(PATCH):../$(SOLVER):../aton:../rt:../hydro:../pm:../poisson:../amr
#############################################################################
# All objects
MODOBJ = amr_parameters.o amr_commons.o random.o pm_parameters.o pm_commons.o poisson_parameters.o poisson_commons.o hydro_parameters.o hydro_commons.o cooling_module.o bisection.o sparse_mat.o clfind_commons.o gadgetreadfile.o rt_parameters.o rt_hydro_commons.o coolrates_module.o rt_spectra.o rt_cooling_module.o rt_flux_module.o
AMROBJ = read_params.o init_amr.o init_time.o init_refine.o adaptive_loop.o amr_step.o update_time.o output_amr.o flag_utils.o physical_boundaries.o virtual_boundaries.o refine_utils.o nbors_utils.o hilbert.o load_balance.o title.o sort.o cooling_fine.o units.o light_cone.o movie.o
# Particle-Mesh objects
PMOBJ = init_part.o output_part.o rho_fine.o synchro_fine.o move_fine.o newdt_fine.o particle_tree.o add_list.o remove_list.o star_formation.o sink_particle.o feedback.o clump_finder.o clump_merger.o flag_formation_sites.o init_sink.o output_sink.o
# Poisson solver objects
POISSONOBJ = init_poisson.o phi_fine_cg.o interpol_phi.o force_fine.o multigrid_coarse.o multigrid_fine_commons.o multigrid_fine_fine.o multigrid_fine_coarse.o gravana.o boundary_potential.o rho_ana.o output_poisson.o
# Hydro objects
HYDROOBJ = init_hydro.o init_flow_fine.o write_screen.o output_hydro.o courant_fine.o godunov_fine.o uplmde.o umuscl.o interpol_hydro.o godunov_utils.o condinit.o hydro_flag.o hydro_boundary.o boundana.o read_hydro_params.o synchro_hydro_fine.o
# RT objects
RTOBJ = rt_init_hydro.o rt_init_xion.o rt_init.o rt_init_flow_fine.o rt_output_hydro.o rt_godunov_fine.o rt_interpol_hydro.o rt_godunov_utils.o rt_condinit.o rt_hydro_flag.o rt_hydro_boundary.o rt_boundana.o rt_read_hydro_params.o rt_units.o
#AGN objects
AGNOBJ = #stellar_wind.o init_stellar.o
# All objects
AMRLIB = $(AMROBJ) $(HYDROOBJ) $(PMOBJ) $(POISSONOBJ) $(AGNOBJ) $(RTOBJ)
# ATON objects
ATON_MODOBJ = timing.o radiation_commons.o rad_step.o
ATON_OBJ = observe.o init_radiation.o rad_init.o rad_boundary.o rad_stars.o rad_backup.o ../aton/atonlib/libaton.a
#############################################################################
ramses:	$(MODOBJ) $(AMRLIB) ramses.o
	$(F90) $(MODOBJ) $(AMRLIB) ramses.o -o $(EXEC)$(NDIM)d $(LIBS)
ramses_aton: $(MODOBJ) $(ATON_MODOBJ) $(AMRLIB) $(ATON_OBJ) ramses.o
	$(F90) $(MODOBJ) $(ATON_MODOBJ) $(AMRLIB) $(ATON_OBJ) ramses.o -o $(EXEC)$(NDIM)d $(LIBS) $(LIBCUDA)
#############################################################################
%.o:%.f90
	$(F90) $(FFLAGS) -c $^ -o $@
#############################################################################
clean :
	rm *.o *.$(MOD)
#############################################################################
