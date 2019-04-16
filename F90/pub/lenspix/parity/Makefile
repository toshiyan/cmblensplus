
exename=exe

## Complier and options ##
FC = ifort
#FLAG = -O3
FLAG = -openmp -openmp-report2 -O3 -ip -fpp -error-limit 5 -heap-arrays -ftrapuv
#DBAG = -check all -std -gen_interfaces -fpe0 -ftrapuv -traceback
FLAGS = $(FLAG) $(DBAG)

libroot = ../../../../../Lib/
MyLib = ${libroot}/mylib
Healpix = ${libroot}/Healpix
cfitsio = ${libroot}/cfitsio/
Lenspix = ${libroot}/lenspix/

MODD = -I$(MyLib)/mod -I$(Lenspix)/mod -I$(Healpix)/include
LIBD = -L$(MyLib)/lib -L$(Lenspix)/lib -L$(cfitsio) -L$(Healpix)/lib
LINK = -llenspix -lhealpix -lcfitsio
OPTION = $(MODD) $(LIBD) $(LINK)

OBJ = main.o

$(exename): $(OBJ)
	$(FC) $(FLAGS) $(OBJ) $(OPTION) -o $@

%.o: %.f90
	$(FC) $(FLAGS) $(OPTION) -c $*.f90

clean:
	rm -f *.o* *.e* *.mod *.d *.pc *.obj core* *.il

