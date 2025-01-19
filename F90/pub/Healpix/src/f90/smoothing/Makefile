# -*- Makefile -*-

fac     = smoothing
inc 	= smo_sub_inc.f90
VPATH	= $(BINDIR)
exec	= $(BINDIR)/$(fac)
#src	= $(fac).f90
obj 	= $(fac).o

MKDIR   = mkdir -p
OD :=$(BUILDDIR)/$(fac)
SD :=$(shell pwd)
OBJ:=$(obj:%=$(OD)/%)
LIB:=$(LIBDIR)/libhealpix$(LIBSUFFIX)


%_mkdir:
	@if [ ! -d $* ]; then $(MKDIR) $* ; fi 

$(OD)/%.o: $(SD)/%.f90 $(inc) $(OD)_mkdir
	$(FC) $(FFLAGS) -c $< -o $@ $(MODDIR)$(OD)

default:  $(exec) 

$(exec): $(OBJ) $(LIB)
	$(FC) $(FFLAGS) -o $(exec) $(OBJ) $(LDFLAGS)
