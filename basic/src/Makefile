#////////////////////////////////////////////////////#
# * Makefile
#////////////////////////////////////////////////////#

## libfile ##
target = libbasic.a

# complier and options
FC = ifort
FLAG = -O3 -fpp -fPIC
FLAGS = $(FLAG) $(DBAG)

# directories
flibloc = ../
flibutl = ../../utils/

MODD = -I${flibutl}/mod
LIBD = -L${flibutl}/lib
LINK = -lutils
OPTION = $(MODD) $(LIBD) $(LINK)

#//// files to be compiled ////#
OBJ = cosmofunc.o pstool.o bstool.o utilsgal.o utilsGW.o

$(target): $(OBJ)
	ar crs $@ $^

%.o: %.f90
	$(FC) $(FLAGS) $(OPTION) -c $*.f90

install:
	mv *.mod ${flibloc}/mod; mv $(target) ${flibloc}/lib/

clean:
	rm -f *.o* *.mod *.a


