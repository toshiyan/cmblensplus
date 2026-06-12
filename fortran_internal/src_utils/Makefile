#////////////////////////////////////////////////////#
# * Makefile
#////////////////////////////////////////////////////#

## libfile ##
target = libutils.a

# compiler and options
FC ?= ifort

# Select Fortran preprocessor flag depending on compiler.
# - ifort / ifx: use -fpp
# - gfortran: use -cpp
FC_BASENAME := $(notdir $(FC))

ifeq ($(findstring gfortran,$(FC_BASENAME)),gfortran)
  FPPFLAG = -cpp
else ifeq ($(findstring ifx,$(FC_BASENAME)),ifx)
  FPPFLAG = -fpp
else ifeq ($(findstring ifort,$(FC_BASENAME)),ifort)
  FPPFLAG = -fpp
else
  # Default to gfortran-style preprocessing.
  FPPFLAG = -cpp
endif

FLAG = -O3 $(FPPFLAG) -fPIC
#DBAG = -check all -std -gen_interfaces -fpe0 -ftrapuv -traceback
FLAGS = $(FLAG) $(DBAG)


# directories
flib = ../

# compile links
MODD = -I$(flib)/mod
OPTION = $(MODD)

#//// files to be compiled ////#
# independent modules
OBJ = constants.o wigner_f90.o funcs.o random.o io.o general.o cosmofunc.o pstool.o bstool.o grid2d.o alkernel.o delensing.o alflat.o bbflat.o utilsgal.o utilsgw.o

$(target): $(OBJ)
	ar crs $@ $^

%.o: %.f90
	$(FC) $(FLAGS) $(OPTION) -c $*.f90

install:
	cp *.mod ${flib}/mod; cp $(target) ${flib}/lib/

clean:
	rm -f *.o *.mod *.a

SRC=$(wildcard *.f90)
remove:
	rm -f $(patsubst %.f90, ${flib}/mod/%.mod, $(SRC)) ${flib}/lib/${target}


