# -*- Makefile -*-

fac     = synfast
inc     = syn_sub_inc.f90
VPATH	= $(BINDIR)
exec	= $(BINDIR)/$(fac)
#src	= $(fac).f90
obj 	= $(fac).o

MKDIR   = mkdir -p
OD :=$(BUILDDIR)/$(fac)
SD :=$(shell pwd)
OBJ:=$(obj:%=$(OD)/%)
LIB:=$(LIBDIR)/libhealpix$(LIBSUFFIX)

# executable creation
.DEFAULT_GOAL := all

all: $(OD)_mkdir $(exec)

$(exec): $(OBJ) $(LIB)
	$(FC) $(FFLAGS) -o $(exec) $(OBJ) $(LDFLAGS)

# create build sub-directory
%_mkdir:
	@if [ ! -d $* ]; then $(MKDIR) $* ; fi 

# objects creation
$(OD)/%.o: $(SD)/%.f90 $(inc)
	$(FC) $(FFLAGS) -c $< -o $@ $(MODDIR)$(OD)


