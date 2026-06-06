#////////////////////////////////////////////////////#
# * Makefile for libdft.a
#////////////////////////////////////////////////////#

target = libdft.a

# compiler and options
FC ?= ifort

# Select Fortran preprocessor flag depending on compiler.
# - ifort / ifx: use -fpp
# - gfortran: use -cpp
FC_BASENAME := $(notdir $(FC))

ifeq ($(findstring gfortran,$(FC_BASENAME)),gfortran)
  FPPFLAG ?= -cpp
else ifeq ($(findstring ifx,$(FC_BASENAME)),ifx)
  FPPFLAG ?= -fpp
else ifeq ($(findstring ifort,$(FC_BASENAME)),ifort)
  FPPFLAG ?= -fpp
else
  FPPFLAG ?= -cpp
endif

OPTFLAGS ?= -O3
PICFLAG ?= -fPIC
INCLUDES ?= -I../mod

# Meson may override this:
#   make FC=... FFLAGS="..."
FFLAGS ?= $(OPTFLAGS) $(FPPFLAG) $(PICFLAG) $(INCLUDES) $(DBAG)


# Optional extra include/library flags.
LDFLAGS ?=
LIBS ?=

# Directories.
flib = ../
MODDIR = $(flib)/mod
LIBDIR = $(flib)/lib

OBJ = fftw.o fft_utils.o

.PHONY: all install clean remove

all: $(target)

$(target): $(OBJ)
	ar crs $@ $^

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

install: $(target)
	cp *.mod $(MODDIR)/
	cp $(target) $(LIBDIR)/

clean:
	rm -f *.o *.mod *.a

SRC = $(wildcard *.f90)

remove:
	rm -f $(patsubst %.f90,$(MODDIR)/%.mod,$(SRC))
	rm -f $(LIBDIR)/$(target)
