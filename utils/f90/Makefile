# Makefile for RAMSES utils
BINDIR=.

# Ifort configuration
F90=ifort
FLAGS=-ffree-line-length-none

# # Gfortran configuration
# F90=gfortran
# FLAGS=-ffree-line-length-none -std=f2003 -Ofast

MODOBJ=random.o utils.o io_ramses.o

# Make all targets
all: amr2prof amr2cylprof ramses2tipsy amr2map part2map part2prof part2cube
amr2prof: $(MODOBJ) amr2prof.o
	$(F90) $(FLAGS) $^ -o $(BINDIR)/$@
amr2cylprof: $(MODOBJ) amr2cylprof.o
	$(F90) $(FLAGS) $^ -o $(BINDIR)/$@
ramses2tipsy: $(MODOBJ) ramses2tipsy.o
	$(F90) $(FLAGS) $^ -o $(BINDIR)/$@
amr2map: $(MODOBJ) amr2map.o
	$(F90) $(FLAGS) $^ -o $(BINDIR)/$@
part2map: $(MODOBJ) utils.o part2map.o
	$(F90) $(FLAGS) $^ -o $(BINDIR)/$@
part2cube: $(MODOBJ) utils.o part2cube.o
	$(F90) $(FLAGS) $^ -o $(BINDIR)/$@
part2prof: $(MODOBJ) utils.o part2prof.o
	$(F90) $(FLAGS) $^ -o $(BINDIR)/$@

# Make a specific object file
%.o: %.f90
	$(F90) $(FLAGS) -c $^
clean:
	rm *.o *.mod
