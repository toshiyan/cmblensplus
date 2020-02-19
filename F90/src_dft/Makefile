#////////////////////////////////////////////////////#
# * Makefile
#////////////////////////////////////////////////////#

## libfile ##
target = libdft.a

## Complier and options ##
FC = ifort
FLAG = -O3 -fpp -fPIC
#DBAG = -check all -std -gen_interfaces -fpe0 -ftrapuv -traceback
FLAGS = $(FLAG) $(DBAG)

## Directories ##
flib = ../
FFTW = ${flib}/pub/FFTW/
MODD = -I$(FFTW)/include -I$(flib)/mod
LIBD = -L$(FFTW)/lib -L$(flib)/lib
LINK = -lutils -lfftw3

OPTION = $(MODD) $(LIBD) $(LINK)

OBJ = fftw.o #fft_remap.o

$(target): $(OBJ)
	ar cr $@ $^

%.o: %.f90
	$(FC) $(FLAGS) $(OPTION) -c $*.f90

install:
	mv *.mod ${flib}/mod; mv $(target) ${flib}/lib/

clean:
	rm -f *.o *.mod *.a

SRC=$(wildcard *.f90)
remove:
	rm -f $(patsubst %.f90, ${flib}/mod/%.mod, $(SRC))

