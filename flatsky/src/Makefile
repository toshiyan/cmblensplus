#////////////////////////////////////////////////////#
# * Makefile
#////////////////////////////////////////////////////#

## libfile ##
target = libflat.a

## Complier and options ##
FC = ifort
FLAG = -O3 -fpp -fPIC
#DBAG = -check all -std -gen_interfaces -fpe0 -ftrapuv -traceback
FLAGS = $(FLAG) $(DBAG)

## Directories ##
floc = ../
flib = ../../F90/
FFTW = ${flib}/pub/FFTW/
MODD = -I$(FFTW)/include -I$(flib)/mod
LIBD = -L$(FFTW)/lib -L$(flib)/lib
LINK = -lfttw -lutils -lfftw3

OPTION = $(MODD) $(LIBD) $(LINK)

OBJ = utils.o ffttools.o rec_lens.o norm_lens.o rec_rot.o rec_tau.o delens.o inp.o remap.o bispec.o

$(target): $(OBJ)
	ar cr $@ $^

%.o: %.f90
	$(FC) $(FLAGS) $(OPTION) -c $*.f90

install:
	mv *.mod ${floc}/mod; mv $(target) ${floc}/lib/

clean:
	rm -f *.o* *.mod *.a


