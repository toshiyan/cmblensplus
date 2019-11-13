#////////////////////////////////////////////////////#
# * Makefile
#////////////////////////////////////////////////////#

## libfile ##
target = libhp.a

## Complier and options ##
FC = ifort
FLAG = -O3 -fpp -fPIC
#DBAG = -check all -std -gen_interfaces -fpe0 -ftrapuv -traceback
FLAGS = $(FLAG) $(DBAG)

## Directories ##
flib = ../
healpix = ${flib}/pub/Healpix
cfitsio = ${flib}/pub/cfitsio

MODD= -I${healpix}/include -I${flib}/mod
LIBD= -L${healpix}/lib -L${cfitsio} -L${flib}/lib
LINK= -lutils -lhealpix -lcfitsio -liomp5 -lpthread

OPTION = $(MODD) $(LIBD) $(LINK)

OBJ = spht.o hp_bsp.o

$(target): $(OBJ)
	ar cr $@ $^

%.o: %.f90
	$(FC) $(FLAGS) $(OPTION) -c $*.f90

install:
	mv *.mod ${flib}/mod; mv $(target) ${flib}/lib/

clean:
	rm -f *.o* *.mod *.a

