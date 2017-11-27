# Makefile for RAMSES utils
BINDIR=.

# Ifort configuration
F90=ifort
CFLAGS=-Ofast -cpp

# Gfortran configuration
# F90=gfortran
# CFLAGS=-ffree-line-length-none -std=f2008 -Ofast -Wall -x f95-cpp-input -J..

LFLAGS=

MODOBJ=

# Make all targets
all: get_music_refmask geticmask

get_music_refmask: $(MODOBJ) ../utils.o get_music_refmask.o
	$(F90) $(LFLAGS) $^ -o $(BINDIR)/$@

geticmask: $(MODOBJ) ../utils.o geticmask.o
	$(F90) $(LFLAGS) $^ -o $(BINDIR)/$@

# Make a specific object file
%.o: %.f90
	$(F90) $(CFLAGS) -c $^ -o $@
clean:
	rm *.o *.mod
