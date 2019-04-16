#////////////////////////////////////////////////////#
# * Makefile
#////////////////////////////////////////////////////#

## libfile ##
target=libmatrix.a

# complier and options
FC = ifort
FLAG = -O3 -fpp -fPIC
#DBAG = -check all -std -gen_interfaces -fpe0 -ftrapuv -traceback
FLAGS = $(FLAG) $(DBAG)

# directories
flibutils = ../
flibloc = ../
LAPACK = ${flibutils}/pub/LAPACK95/

# compile links
MODD = -I$(flibutils)/mod -I$(LAPACK)/mod
LIBD = -L$(flibutils)/lib -L$(LAPACK)/lib
LINK = -lutils -llapack95 -llapack -lrefblas
OPTION = $(MODD) $(LIBD) $(LINK)

# files to be compiled
OBJ = blas.o lapack.o lapack95.o

$(target): $(OBJ)
	ar crs $@ $^

%.o: %.f90
	$(FC) $(FLAGS) $(OPTION) -c $*.f90

install:
	mv *.mod ${flibloc}/mod; mv $(target) ${flibloc}/lib/

clean:
	rm -f *.o* *.mod *.a

# file dependencies

