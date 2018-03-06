# Makefile for RAMSES utils
BINDIR=.

# Ifort configuration
F90=ifort
CFLAGS=-Ofast -cpp -traceback

# # Gfortran configuration
#F90=gfortran
#CFLAGS=-ffree-line-length-none -std=f2008 -Ofast -Wall -x f95-cpp-input

LFLAGS=

MODOBJ=random.o utils.o io_ramses.o

# Make all targets
all: amr2prof amr2cylprof ramses2tipsy amr2map part2map part2prof \
     part2cube part2list partcenter part2birth part2sfr

amr2prof: $(MODOBJ) amr2prof.o
	$(F90) $(LFLAGS) $^ -o $(BINDIR)/$@
amr2cylprof: $(MODOBJ) amr2cylprof.o
	$(F90) $(LFLAGS) $^ -o $(BINDIR)/$@
amr2map: $(MODOBJ) amr2map.o
	$(F90) $(LFLAGS) $^ -o $(BINDIR)/$@

ramses2tipsy: $(MODOBJ) ramses2tipsy.o
	$(F90) $(LFLAGS) $^ -o $(BINDIR)/$@

part2map: $(MODOBJ) utils.o part2map.o
	$(F90) $(LFLAGS) $^ -o $(BINDIR)/$@
part2cube: $(MODOBJ) utils.o part2cube.o
	$(F90) $(LFLAGS) $^ -o $(BINDIR)/$@
part2prof: $(MODOBJ) utils.o part2prof.o
	$(F90) $(LFLAGS) $^ -o $(BINDIR)/$@
part2list: $(MODOBJ) utils.o part2list.f90
	$(F90) $(LFLAGS) $^ -o $(BINDIR)/$@

partcenter: $(MODOBJ) utils.o partcenter.f90
	$(F90) $(LFLAGS) $^ -o $(BINDIR)/$@
part2birth: $(MODOBJ) utils.o part2birth.f90
	$(F90) $(LFLAGS) $^ -o $(BINDIR)/$@
part2sfr: $(MODOBJ) part2sfr.f90
	$(F90) $(LFLAGS) $^ -o $(BINDIR)/$@

# Make a specific object file
%.o: %.f90
	$(F90) $(CFLAGS) -c $^ -o $@
clean:
	rm *.o *.mod
