# Mayall:
C_LIBS = -lm -lcuda  -L/usr/local/cuda/lib64 -L/home/itp/tstranex/cuda/sdk/C/lib -lcudart -I/usr/local/cuda/include
# C_LIBS = -lm -lcuda  -L/opt/cuda/lib -lcudart
CUDA_LIBS =

F_OBJS=
C_OBJS=
CC_OBJS= cpu_numerical.o cpu_memory.o validation.o aton_fortran.o gpu_init.o aton_cpp.o
CUDA_OBJS= gpu_memory.o gpu_numerical.o
OBJECTS = $(C_OBJS) $(CUDA_OBJS) $(F_OBJS) $(CC_OBJS)
VPATH=internal

DEFINES = -DS_100000  # -DTEST7_RAYTRACING_TRICKS
PUBLIC_INCLUDES=$(shell pwd)


NVCC= nvcc -arch=sm_13 #--device-emulation -g
C = mpicc
CC = mpic++ -O2
F90 = mpif90

EXECUTABLE = 

libaton.a: $(OBJECTS)
	ar rcs libaton.a $(OBJECTS)

aton_cpu_test: $(OBJECTS) internal/aton_cpu_test.o
	$(F90) -o aton_cpu_test internal/aton_cpu_test.o $(OBJECTS) $(C_LIBS)

stromgren_test: $(OBJECTS) examples/stromgren_test.o
	$(CC) -o stromgren_test stromgren_test.o $(OBJECTS) $(C_LIBS)

%.o:%.f90
	$(F90) -c $*.f90
.c.o:
	$(C) $(DEFINES) -I. -c $<
.cc.o:
	$(CC) $(DEFINES) $(C_LIBS) -I$(PUBLIC_INCLUDES) -c $^
%.o:%.cu
	$(NVCC) $(DEFINES) $(C_LIBS) $(CUDA_LIBS) -I. -c $^

.PHONY : all
all:$(C_OBJS) $(CUDA_OBJS) $(F_OBJS)

clean:
	rm -f *.o *.cudafe1.* *.cudafe2.* *.hash *.ptx *fatbin.c *.cubin *.cpp* $(EXECUTABLE) libaton.a stromgren_test *~ *.mod
