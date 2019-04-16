# * Makefile

## libfile ##
target=libanalytic.a

## Complier and options ##
FC = ifort
FLAG = -O3 -fpp -fPIC
#DBAG = -check all -std -gen_interfaces -fpe0 -ftrapuv -traceback
FLAGS = $(FLAG) $(DBAG)

## Directories ##
libroot = ../../F90/
libloc = ../
#LAPACK = ${libroot}/pub/LAPACK95/

MODD = -I$(libroot)/mod
LIBD = -L$(libroot)/lib
LINK = -lutils 
#MODD = -I$(LAPACK)/mod -I$(libroot)/mod
#LIBD = -L$(LAPACK)/lib -L$(libroot)/lib
#LINK = -llinalg -lutils -llapack95 -llapack -lrefblas

OPTION = $(MODD) $(LIBD) $(LINK)

OBJ = bb_flat.o al_lens.o delensing.o al_rot.o al_tau.o al_src.o al_flat.o

$(target): $(OBJ)
	ar cr $@ $^

%.o: %.f90
	$(FC) $(FLAGS) $(OPTION) -c $*.f90

install:
	mv *.mod ${libloc}/mod; mv $(target) ${libloc}/lib/

clean:
	rm -f *.o*

