#### Makefile ####
## library ##
target=liblenspix.a

# complier and options
FC = ifort
FLAG = -qopenmp -ip -fpp -fPIC
#DBAG = -check all -std -gen_interfaces -fpe0 -ftrapuv -traceback
FLAGS = $(FLAG) $(DBAG)

# directories
libroot = ../../
healpix = ${libroot}/Healpix/
cfitsio = ${libroot}/cfitsio/

# compile links
MODD = -I$(healpix)/include
LIBD = -L$(cfitsio)/ -L$(healpix)/lib
LINK = -lhealpix -lcfitsio
OPTION = $(MODD) $(LIBD) $(LINK)

# files to be compiled
#OBJ = toms760.o inifile.o utils.o spin_alm_tools.o remap.o #HealpixObj.o HealpixVis.o
OBJ = utils.o toms760.o spin_alm_tools.o remap.o

$(target): $(OBJ)
	ar crs $@ $^

%.o: %.f90
	$(FC) $(FLAGS) $(OPTION) -c $*.f90

install:
	mv *.mod ../mod; mv $(target) ../lib/

clean:
	rm -f *.o* *.mod *.a

# file dependencies
spin_alm_tools.o:  utils.o toms760.o

