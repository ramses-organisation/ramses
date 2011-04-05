# Makefile for RAMSES utils
F90=ifort
BINDIR=.

# Make all targets
all: $(patsubst %.f90, %, $(wildcard *.f90))

# Make a specific f90 file
%: %.f90
	$(F90) $< -o $(BINDIR)/$@
