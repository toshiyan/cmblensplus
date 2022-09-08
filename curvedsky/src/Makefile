#////////////////////////////////////////////////////#
# * Makefile
#////////////////////////////////////////////////////#

## libfile ##
target=libcurvedsky.a

## Complier and options ##
FC = ifort
FLAG = -qopenmp -fpp -fPIC
FLAGS = $(FLAG) $(DBAG)

## Directories ##
flibloc = ../
flibutl = ../../F90/
healpix = ${flibutl}/pub/Healpix
cfitsio = ${flibutl}/pub/cfitsio
lenspix = ${flibutl}/pub/lenspix
LAPACK  = ${flib}/pub/LAPACK95/

MODD= -I${lenspix}/mod -I${healpix}/include -I$(LAPACK)/mod -I${flibutl}/mod
LIBD= -L${lenspix}/lib -L${healpix}/lib -L${cfitsio} -L$(LAPACK)/lib -L${flibutl}/lib
LINK= -llenspix -lutils -lhealpix -lcfitsio -liomp5 -lpthread -llapack95 -llapack -lrefblas

OPTION = $(MODD) $(LIBD) $(LINK) -Wl,-R${healpix}/lib

OBJ = utils.o bispec.o rec_lens.o delens.o rec_tau.o rec_src.o

$(target): $(OBJ)
	ar cr $@ $^

%.o: %.f90
	$(FC) $(FLAGS) $(OPTION) -c $*.f90

install:
	mv *.mod ${flibloc}/mod; mv $(target) ${flibloc}/lib/

clean:
	rm -f *.o*

