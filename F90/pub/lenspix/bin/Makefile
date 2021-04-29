#now using v11+ intel ifort

#Intel. Note heap-arrays needs v10+ of compiler, and avoids Seg Faults for large arrays
F90C     = ifort
#F90C     = mpif90

#FFLAGS = -openmp -openmp-report2 -O3 -ip -fpp -error-limit 5 -heap-arrays -check all -ftrapuv -traceback -g
FFLAGS = -qopenmp -O3 -ip -fpp -error-limit 5 -heap-arrays -ftrapuv
#DMPIPIX and DMPI is not included in version 11.0

Lib = ../
MyLib = ../../../
HEALPIX = ${MyLib}/pub/Healpix/
cfitsio = ${MyLib}/pub/cfitsio/

F90FLAGS = $(FFLAGS) -I$(MyLib)/mod -L$(MyLib)/lib -I$(Lib)/mod -L$(Lib)/lib -I$(HEALPIX)/include -L$(cfitsio) -L$(healpix)/lib -lutils -llenspix -lhealpix -lcfitsio

OBJFILES= main.o

default: exe
all: exe

.f.o:
	f77 $(F90FLAGS) -c $<

%.o: %.f90
	$(F90C) $(F90FLAGS) -c $*.f90

%.o: %.F90
	$(F90C) $(F90FLAGS) -c $*.F90

exe: $(OBJFILES) 	
	$(F90C) -o exe $(OBJFILES) $(F90FLAGS)

clean:
	rm -f *.o* *.e* *.mod *.d *.pc *.obj core* *.il

