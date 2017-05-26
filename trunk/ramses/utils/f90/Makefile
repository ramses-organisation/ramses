# Makefile for RAMSES utils
F90=ifort 
#gfortran
BINDIR=.
FLAGS=-ffree-line-length-none

MODOBJ=random.o io_ramses.o

# Make all targets
all: amr2prof amr2cylprof ramses2tipsy amr2map part2map part2prof
amr2prof: $(MODOBJ) amr2prof.o
	$(F90) $(MODOBJ) $(FLAGS) amr2prof.o -o $(BINDIR)/amr2prof
amr2cylprof: $(MODOBJ) amr2cylprof.o
	$(F90) $(MODOBJ) $(FLAGS) amr2cylprof.o -o $(BINDIR)/amr2cylprof
ramses2tipsy: $(MODOBJ) ramses2tipsy.o
	$(F90) $(MODOBJ) $(FLAGS) ramses2tipsy.o -o $(BINDIR)/ramses2tipsy
amr2map: $(MODOBJ) amr2map.o
	$(F90) $(FLAGS) amr2map.o -o $(BINDIR)/amr2map
part2map: $(MODOBJ) part2map.o
	$(F90) $(FLAGS) part2map.o -o $(BINDIR)/part2map
part2prof: $(MODOBJ) part2prof.o
	$(F90) $(FLAGS) part2prof.o -o $(BINDIR)/part2prof

# Make a specific object file
%.o: %.f90
	$(F90) $(FLAGS) -c $^
clean:
	rm *.o *.mod
