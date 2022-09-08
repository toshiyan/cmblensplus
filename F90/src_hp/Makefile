#////////////////////////////////////////////////////#
# * Makefile
#////////////////////////////////////////////////////#

## libfile ##
target = libhp.a

## Complier and options ##
FC = ifort
FLAG = -O3 -fpp -fPIC -qopenmp
#DBAG = -check all -std -gen_interfaces -fpe0 -ftrapuv -traceback
FLAGS = $(FLAG) $(DBAG)

## Directories ##
flib = ../
healpix = ${flib}/pub/Healpix
cfitsio = ${flib}/pub/cfitsio
LAPACK  = ${flib}/pub/LAPACK95/

MODD= -I${healpix}/include -I$(LAPACK)/mod -I${flib}/mod
LIBD= -L${healpix}/lib -L${cfitsio} -L$(LAPACK)/lib -L${flib}/lib
#LINK= -lutils -lhealpix -lsharp -lcfitsio -liomp5 -lpthread -llapack95 -llapack -lrefblas
LINK= -lutils -lhealpix -lcfitsio -liomp5 -lpthread -llapack95 -llapack -lrefblas

OPTION = $(MODD) $(LIBD) $(LINK)

SRC=$(wildcard *.f90)
#OBJ=$(patsubst %.f90, %.o, $(SRC))
OBJ = hp_utils.o hp_udgrade.o hp_spht.o hp_bsp.o hp_cgd.o

$(target): $(OBJ)
	ar cr $@ $^

%.o: %.f90
	$(FC) $(FLAGS) $(OPTION) -c $*.f90

install:
	mv *.mod ${flib}/mod; mv $(target) ${flib}/lib/

clean:
	rm -f *.o *.mod *.a

remove:
	rm -f $(patsubst %.f90, ${flib}/mod/%.mod, $(SRC))


